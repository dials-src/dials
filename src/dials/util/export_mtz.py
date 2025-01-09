from __future__ import annotations

import logging
import time
import warnings
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass, field
from math import isclose

import gemmi
import numpy as np
import pandas as pd

from cctbx import uctbx
from dxtbx import flumpy
from iotbx import mtz
from libtbx import env
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from scitbx import matrix
from scitbx.math import r3_rotation_axis_and_angle_from_matrix

import dials.util.ext
from dials.algorithms.scaling.scaling_library import (
    MergedHalfDatasets,
    determine_best_unit_cell,
)
from dials.array_family import flex
from dials.util.batch_handling import (
    assign_batches_to_reflections,
    calculate_batch_offsets,
    get_image_ranges,
)
from dials.util.filter_reflections import filter_reflection_table
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.reindex import reindex_experiments, reindex_reflections
from dials.util.version import dials_version

logger = logging.getLogger(__name__)


class MTZWriterBase:
    """Helper for adding metadata, crystals and datasets to an mtz file object."""

    def __init__(self, space_group, unit_cell=None):
        """If a unit cell is provided, will be used as default unless specified
        for each crystal."""
        warnings.warn(
            "MTZWriterBase classes (MergedMTZWriter and MADMergedMTZWriter) are deprecated. Use MergedMTZCreator instead.\n",
            DeprecationWarning,
            stacklevel=2,
        )
        mtz_file = mtz.object()
        mtz_file.set_title(f"From {env.dispatcher_name}")
        date_str = time.strftime("%Y-%m-%d at %H:%M:%S %Z")
        if time.strftime("%Z") != "GMT":
            date_str += time.strftime("  (%Y-%m-%d at %H:%M:%S %Z)", time.gmtime())
        mtz_file.add_history(f"From {dials_version()}, run on {date_str}")
        mtz_file.set_space_group_info(space_group.info())
        self.mtz_file = mtz_file
        if unit_cell:
            self.unit_cell = unit_cell
        self.current_crystal = None
        self.current_dataset = None
        self.n_crystals = 0
        self.n_datasets = 0

    def add_crystal(self, crystal_name=None, project_name=None, unit_cell=None):
        """Add a crystal to the mtz file object."""
        if not unit_cell:
            if not self.unit_cell:
                raise ValueError("Unit cell must be provided.")
            else:
                unit_cell = self.unit_cell
        if not crystal_name:
            crystal_name = f"crystal_{self.n_crystals + 1}"
        if not project_name:
            project_name = "DIALS"
        self.current_crystal = self.mtz_file.add_crystal(
            crystal_name, project_name, unit_cell.parameters()
        )
        self.n_crystals += 1

    def add_empty_dataset(self, wavelength, name=None):
        """Add an empty dataset object to the mtz file."""
        if not name:
            name = "FROMDIALS"
        self.current_dataset = self.current_crystal.add_dataset(name, wavelength)
        self.n_datasets += 1


class MergedMTZWriter(MTZWriterBase):
    """Mtz writer for merged data."""

    def add_dataset(
        self,
        merged_array=None,
        anom_array=None,
        amplitudes=None,
        anom_amplitudes=None,
        dano=None,
        multiplicities=None,
        anom_multiplicities=None,
        suffix=None,
        half_datasets: MergedHalfDatasets | None = None,
        r_free_array=None,
    ):
        """Add merged data to the most recent dataset.

        Args:
            merged_array: A merged miller array of IMEAN intensities
            wavelength: The wavelength of the dataset
            anom_array (Optional): An anomalous merged miller array
            amplitudes (Optional): A merged miller array of amplitudes
            anom_amplitudes (Optional): An anomalous merged array of amplitudes
            suffix (Optional[str]): Column name suffix to use for this dataset.
        """
        if not suffix:
            suffix = ""
        if merged_array:
            self.current_dataset.add_miller_array(merged_array, "IMEAN" + suffix)
        if multiplicities:
            self.current_dataset.add_miller_array(multiplicities, "N" + suffix)
        if amplitudes:
            self.current_dataset.add_miller_array(amplitudes, "F" + suffix)
        if anom_array:
            self.current_dataset.add_miller_array(anom_array, "I" + suffix)
        if anom_multiplicities:
            self.current_dataset.add_miller_array(anom_multiplicities, "N" + suffix)
        if anom_amplitudes:
            self.current_dataset.add_miller_array(anom_amplitudes, "F" + suffix)
        if dano:
            self.current_dataset.add_miller_array(
                dano, "DANO" + suffix, column_types="DQ"
            )
        if half_datasets:
            self.current_dataset.add_miller_array(
                half_datasets.data1, "IHALF1" + suffix, column_types="JQ"
            )
            self.current_dataset.add_miller_array(
                half_datasets.data2, "IHALF2" + suffix, column_types="JQ"
            )
            self.current_dataset.add_miller_array(
                half_datasets.multiplicity1,
                "NHALF1" + suffix,
            )
            self.current_dataset.add_miller_array(
                half_datasets.multiplicity2,
                "NHALF2" + suffix,
            )
        if r_free_array:
            self.current_dataset.add_miller_array(
                r_free_array, column_root_label="FreeR_flag", column_types="I"
            )


class MADMergedMTZWriter(MergedMTZWriter):
    """Mtz writer for multi-wavelength merged data."""

    def add_dataset(
        self,
        merged_array=None,
        anom_array=None,
        amplitudes=None,
        anom_amplitudes=None,
        dano=None,
        multiplicities=None,
        anom_multiplicities=None,
        suffix=None,
        half_datasets: MergedHalfDatasets | None = None,
        r_free_array=None,
    ):
        if not suffix:
            suffix = f"_WAVE{self.n_datasets}"
        super().add_dataset(
            merged_array,
            anom_array,
            amplitudes,
            anom_amplitudes,
            dano,
            multiplicities,
            anom_multiplicities,
            suffix,
            half_datasets,
            r_free_array=r_free_array,
        )


def add_batch_list(
    mtz,
    image_range,
    experiment,
    wavelength,
    dataset_id,
    batch_offset,
    force_static_model,
):
    """Add batch metadata to the gemmi mtz object."""

    # Recalculate useful numbers and references here
    n_batches = image_range[1] - image_range[0] + 1
    phi_start = flex.float(n_batches, 0)
    phi_range = flex.float(n_batches, 0)
    umat_array = flex.float(flex.grid(n_batches, 9))
    cell_array = flex.float(flex.grid(n_batches, 6))

    # Reciprocal lattice vectors in the lab frame at zero scan angle
    if experiment.goniometer:
        S = matrix.sqr(experiment.goniometer.get_setting_rotation())
        F = matrix.sqr(experiment.goniometer.get_fixed_rotation())
        UBlab = S * F * matrix.sqr(experiment.crystal.get_A())

        axis = matrix.col(experiment.goniometer.get_rotation_axis())
        axis_datum = matrix.col(experiment.goniometer.get_rotation_axis_datum())

    else:
        UBlab = matrix.sqr(experiment.crystal.get_A())

    i0 = image_range[0]
    for i in range(n_batches):
        if experiment.scan and experiment.scan.get_oscillation()[1] != 0.0:
            phi_start[i], phi_range[i] = experiment.scan.get_image_oscillation(i + i0)

        # Unit cell and UB matrix for the centre of the image for scan-varying model
        if (
            not force_static_model
            and experiment.crystal.num_scan_points > 0
            and experiment.goniometer
        ):
            # Get the index of the image in the sequence e.g. first => 0, second => 1
            image_index = i + i0 - experiment.scan.get_image_range()[0]

            # Find the U matrix at the frame centre by calculating the linear transform
            # that goes from the start of the frame to the end, and then applying half of
            # that to the start value
            U0 = matrix.sqr(experiment.crystal.get_U_at_scan_point(image_index))
            U1 = matrix.sqr(experiment.crystal.get_U_at_scan_point(image_index + 1))
            M = U1 * U0.inverse()
            (
                angle_M,
                axis_M,
            ) = M.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
                deg=False
            )
            M_half = axis_M.axis_and_angle_as_r3_rotation_matrix(angle_M / 2, deg=False)
            Ucentre = M_half * U0

            # Find the B matrix at the frame centre by interpolation
            B0 = matrix.sqr(experiment.crystal.get_B_at_scan_point(image_index))
            B1 = matrix.sqr(experiment.crystal.get_B_at_scan_point(image_index + 1))
            Bcentre = (B0 + B1) / 2

            # Unit cell at the frame centre
            unit_cell = uctbx.unit_cell(
                orthogonalization_matrix=Bcentre.transpose().inverse()
            )

            # Get full lab frame UB then unwind to zero scan angle
            phi_centre = phi_start[i] + phi_range[i] / 2
            R = matrix.sqr(
                axis_datum.axis_and_angle_as_r3_rotation_matrix(phi_centre, deg=False)
            )
            Rlab_inv = matrix.sqr(
                axis.axis_and_angle_as_r3_rotation_matrix(-phi_centre, deg=False)
            )
            _UBlab = Rlab_inv * S * R * F * Ucentre * Bcentre

        else:
            unit_cell = experiment.crystal.get_unit_cell()
            _UBlab = UBlab

        # We assume a single-axis goniometer as it is not clear that multi-
        # axis goniometry was ever fully supported in MTZ format. Orientation
        # will be taken from the laboratory frame for this image.
        U = matrix.sqr(dials.util.ext.ub_to_mosflm_u(_UBlab, unit_cell))

        # FIXME need to get what was refined and what was constrained from the
        # crystal model - see https://github.com/dials/dials/issues/355
        _unit_cell_params = unit_cell.parameters()
        for j in range(6):
            cell_array[i, j] = _unit_cell_params[j]
        # Transpose to put in column-major order for MTZ export
        U_t_elements = U.transpose().elems
        for j in range(9):
            umat_array[i, j] = U_t_elements[j]

    # We ignore panels beyond the first one, at the moment
    panel = experiment.detector[0]
    panel_size = panel.get_image_size()
    panel_distance = panel.get_directed_distance()

    if experiment.goniometer:
        axis = flex.float(experiment.goniometer.get_rotation_axis())
    else:
        axis = flex.float((0.0, 0.0, 0.0))

    source = flex.float(experiment.beam.get_sample_to_source_direction())

    # get the mosaic spread though today it may not actually be set - should
    # this be in the BATCH headers?
    try:
        mosaic = experiment.crystal.get_mosaicity()
    except AttributeError:
        mosaic = 0.0

    max_batch_number = 0
    if mtz.batches:
        max_batch_number = mtz.batches[-1].number

    batch_offset += image_range[0] - 1
    if max_batch_number > batch_offset:
        batch_offset = max_batch_number

    batch = gemmi.Mtz.Batch()

    # Setting fields that are the same for all batches
    batch.dataset_id = dataset_id
    batch.wavelength = wavelength
    batch.ints[12] = 1  # ncryst
    batch.ints[14] = 2  # ldtype 3D
    batch.ints[15] = 1  # jsaxs - goniostat scan axis number
    batch.ints[17] = 1  # ngonax - number of goniostat axes
    batch.ints[19] = 1  # ndet

    batch.floats[21] = mosaic  # crydat[0]
    for j in range(3):
        batch.floats[38 + j] = axis[j]  # scanax
    batch.floats[43] = 1.0  # bscale (batch scale)
    for j in range(3):
        batch.floats[59 + j] = axis[j]  # e1
    batch.floats[80 + flex.min_index(source)] = -1.0  # idealised source vector
    for j in range(3):
        batch.floats[83 + j] = source[j]  # source including tilts
    batch.floats[111] = panel_distance  # dx
    batch.floats[114] = panel_size[0]  # NX
    batch.floats[116] = panel_size[1]  # NY
    batch.axes = ["AXIS"]  # gonlab[0]

    # Setting fields that differ
    for i_batch in range(n_batches):
        batch.number = batch_offset + i_batch + 1
        batch.title = f"Batch {batch.number}"
        for j in range(6):
            batch.floats[j] = cell_array[i_batch, j]  # cell
        for j in range(9):
            batch.floats[6 + j] = umat_array[i_batch, j]  # Umat
        batch.floats[36] = phi_start[i_batch]  # phistt
        batch.floats[37] = phi_start[i_batch] + phi_range[i_batch]  # phiend
        batch.floats[47] = phi_range[i_batch]  # phirange

        # Append this batch
        mtz.batches.append(batch)

    return


def write_columns(mtz, reflection_table):
    """Write the column definitions AND data to the current dataset."""

    nref = len(reflection_table["miller_index"])
    assert nref
    xdet, ydet, _ = (
        flex.double(x) for x in reflection_table["xyzobs.px.value"].parts()
    )

    type_table = {
        "H": "H",
        "K": "H",
        "L": "H",
        "I": "J",
        "SIGI": "Q",
        "IPR": "J",
        "SIGIPR": "Q",
        "BG": "R",
        "SIGBG": "R",
        "XDET": "R",
        "YDET": "R",
        "BATCH": "B",
        "BGPKRATIOS": "R",
        "WIDTH": "R",
        "MPART": "I",
        "M_ISYM": "Y",
        "FLAG": "I",
        "LP": "R",
        "FRACTIONCALC": "R",
        "ROT": "R",
        "QE": "R",
    }

    mtz_data = pd.DataFrame(
        flumpy.to_numpy(reflection_table["miller_index"]).astype("float32"),
        columns=["H", "K", "L"],
    )
    mtz_data.insert(3, "M/ISYM", np.zeros(nref, dtype="float32"))

    # H, K, L are in the base dataset, but we have to add M/ISYM
    mtz.add_column("M/ISYM", type_table["M_ISYM"])
    mtz.add_column("BATCH", type_table["BATCH"])
    mtz_data.insert(
        4, "BATCH", flumpy.to_numpy(reflection_table["batch"]).astype("float32")
    )

    # if intensity values used in scaling exist, then just export these as I, SIGI
    if "intensity.scale.value" in reflection_table:
        I_scaling = reflection_table["intensity.scale.value"]
        V_scaling = reflection_table["intensity.scale.variance"]
        assert V_scaling.all_gt(0)  # Trap negative variances
        mtz.add_column("I", type_table["I"])
        mtz_data.insert(
            len(mtz_data.columns), "I", flumpy.to_numpy(I_scaling).astype("float32")
        )
        mtz.add_column("SIGI", type_table["SIGI"])
        mtz_data.insert(
            len(mtz_data.columns),
            "SIGI",
            flumpy.to_numpy(flex.sqrt(V_scaling)).astype("float32"),
        )
        mtz.add_column("SCALEUSED", "R")
        mtz_data.insert(
            len(mtz_data.columns),
            "SCALEUSED",
            flumpy.to_numpy(reflection_table["inverse_scale_factor"]).astype("float32"),
        )
        mtz.add_column("SIGSCALEUSED", "R")
        mtz_data.insert(
            len(mtz_data.columns),
            "SIGSCALEUSED",
            flumpy.to_numpy(
                flex.sqrt(reflection_table["inverse_scale_factor_variance"])
            ).astype("float32"),
        )
    else:
        if "intensity.prf.value" in reflection_table:
            if "intensity.sum.value" in reflection_table:
                col_names = ("IPR", "SIGIPR")
            else:
                col_names = ("I", "SIGI")
            I_profile = reflection_table["intensity.prf.value"]
            V_profile = reflection_table["intensity.prf.variance"]
            assert V_profile.all_gt(0)  # Trap negative variances
            mtz.add_column(col_names[0], type_table["I"])
            mtz_data.insert(
                len(mtz_data.columns),
                col_names[0],
                flumpy.to_numpy(I_profile.as_float()).astype("float32"),
            )
            mtz.add_column(col_names[1], type_table["SIGI"])
            mtz_data.insert(
                len(mtz_data.columns),
                col_names[1],
                flumpy.to_numpy(flex.sqrt(V_profile)).astype("float32"),
            )

        if "intensity.sum.value" in reflection_table:
            I_sum = reflection_table["intensity.sum.value"]
            V_sum = reflection_table["intensity.sum.variance"]
            assert V_sum.all_gt(0)  # Trap negative variances
            mtz.add_column("I", type_table["I"])
            mtz_data.insert(
                len(mtz_data.columns), "I", flumpy.to_numpy(I_sum).astype("float32")
            )
            mtz.add_column("SIGI", type_table["SIGI"])
            mtz_data.insert(
                len(mtz_data.columns),
                "SIGI",
                flumpy.to_numpy(flex.sqrt(V_sum)).astype("float32"),
            )

    if (
        "background.sum.value" in reflection_table
        and "background.sum.variance" in reflection_table
    ):
        bg = reflection_table["background.sum.value"]
        varbg = reflection_table["background.sum.variance"]
        assert (varbg >= 0).count(False) == 0
        sigbg = flex.sqrt(varbg)
        mtz.add_column("BG", type_table["BG"])
        mtz_data.insert(
            len(mtz_data.columns), "BG", flumpy.to_numpy(bg).astype("float32")
        )
        mtz.add_column("SIGBG", type_table["SIGBG"])
        mtz_data.insert(
            len(mtz_data.columns), "SIGBG", flumpy.to_numpy(sigbg).astype("float32")
        )

    mtz.add_column("FRACTIONCALC", type_table["FRACTIONCALC"])
    mtz_data.insert(
        len(mtz_data.columns),
        "FRACTIONCALC",
        flumpy.to_numpy(reflection_table["fractioncalc"]).astype("float32"),
    )

    mtz.add_column("XDET", type_table["XDET"])
    mtz_data.insert(
        len(mtz_data.columns), "XDET", flumpy.to_numpy(xdet).astype("float32")
    )
    mtz.add_column("YDET", type_table["YDET"])
    mtz_data.insert(
        len(mtz_data.columns), "YDET", flumpy.to_numpy(ydet).astype("float32")
    )
    mtz.add_column("ROT", type_table["ROT"])
    mtz_data.insert(
        len(mtz_data.columns),
        "ROT",
        flumpy.to_numpy(reflection_table["ROT"]).astype("float32"),
    )
    if "lp" in reflection_table:
        mtz.add_column("LP", type_table["LP"])
        mtz_data.insert(
            len(mtz_data.columns),
            "LP",
            flumpy.to_numpy(reflection_table["lp"]).astype("float32"),
        )
    if "qe" in reflection_table:
        mtz.add_column("QE", type_table["QE"])
        mtz_data.insert(
            len(mtz_data.columns),
            "QE",
            flumpy.to_numpy(reflection_table["qe"]).astype("float32"),
        )
    elif "dqe" in reflection_table:
        mtz.add_column("QE", type_table["QE"])
        mtz_data.insert(
            len(mtz_data.columns),
            "QE",
            flumpy.to_numpy(reflection_table["dqe"]).astype("float32"),
        )
    else:
        mtz.add_column("QE", type_table["QE"])
        mtz_data.insert(len(mtz_data.columns), "QE", np.ones(nref).astype("float32"))

    mtz.switch_to_original_hkl()
    mtz.set_data(mtz_data.to_numpy())


def export_mtz(
    reflection_table,
    experiment_list,
    intensity_choice,
    filename,
    best_unit_cell=None,
    partiality_threshold=0.4,
    combine_partials=True,
    min_isigi=-5,
    filter_ice_rings=False,
    d_min=None,
    force_static_model=False,
    crystal_name=None,
    project_name=None,
    wavelength_tolerance=1e-4,
):
    """Export data from reflection_table corresponding to experiment_list to an
    MTZ file hklout."""

    # First get the experiment identifier information out of the data
    expids_in_table = reflection_table.experiment_identifiers()
    if not list(expids_in_table.keys()):
        reflection_tables = parse_multiple_datasets([reflection_table])
        experiment_list, refl_list = assign_unique_identifiers(
            experiment_list, reflection_tables
        )
        reflection_table = flex.reflection_table()
        for reflections in refl_list:
            reflection_table.extend(reflections)
        expids_in_table = reflection_table.experiment_identifiers()
    reflection_table.assert_experiment_identifiers_are_consistent(experiment_list)
    expids_in_list = list(experiment_list.identifiers())

    # Convert geometry to the Cambridge frame
    experiment_list = convert_to_cambridge(experiment_list)

    # Validate multi-experiment assumptions
    if len(experiment_list) > 1:
        # All experiments should match crystals, or else we need multiple crystals/datasets
        if not all(
            x.crystal == experiment_list[0].crystal for x in experiment_list[1:]
        ):
            logger.warning(
                "Experiment crystals differ. Using first experiment crystal for file-level data."
            )

        # At least, all experiments must have the same space group
        if len({x.crystal.get_space_group().make_tidy() for x in experiment_list}) != 1:
            raise ValueError("Experiments do not have a unique space group")

    # Reindex to a tabulated setting if necessary
    sg = experiment_list[0].crystal.get_space_group()
    if sg.match_tabulated_settings().number() == 0:
        logger.warning(
            "The data will be reindexed to a tabulated setting of the space group"
        )
        cb_op = sg.info().change_of_basis_op_to_reference_setting()
        experiment_list = reindex_experiments(experiment_list, cb_op)
        reflection_table = reindex_reflections(
            [
                reflection_table,
            ],
            cb_op,
        )

    # Convert experiment_list to a real python list or else identity assumptions
    # fail like:
    #   assert experiment_list[0] is experiment_list[0]
    # And assumptions about added attributes break
    experiment_list = list(experiment_list)

    wavelengths = match_wavelengths(experiment_list, wavelength_tolerance)
    for w in wavelengths.values():
        w.calculate_weighted_mean([reflection_table])

    if len(wavelengths) > 1:
        identifiers_list = [e.identifier for e in experiment_list]
        logger.info(
            "Multiple wavelengths found: \n%s",
            "\n".join(
                "  Wavelength: {:.5f}, experiment numbers: {} ".format(
                    v.weighted_mean,
                    ",".join(
                        map(str, [identifiers_list.index(i) for i in v.identifiers])
                    ),
                )
                for v in wavelengths.values()
            ),
        )

    # also only work correctly with one panel (for the moment)
    if any(len(experiment.detector) != 1 for experiment in experiment_list):
        logger.warning("Ignoring multiple panels in output MTZ")

    if best_unit_cell is None:
        best_unit_cell = determine_best_unit_cell(experiment_list)
    reflection_table["d"] = best_unit_cell.d(reflection_table["miller_index"])

    # Clean up the data with the passed in options
    reflection_table = filter_reflection_table(
        reflection_table,
        intensity_choice=intensity_choice,
        partiality_threshold=partiality_threshold,
        combine_partials=combine_partials,
        min_isigi=min_isigi,
        filter_ice_rings=filter_ice_rings,
        d_min=d_min,
    )

    # get batch offsets and image ranges - even for scanless experiments
    batch_offsets = [
        expt.scan.get_batch_offset()
        for expt in experiment_list
        if expt.scan is not None
    ]
    unique_offsets = set(batch_offsets)
    if len(unique_offsets) != len(batch_offsets):
        logger.debug("Calculating new batches")
        batch_offsets = calculate_batch_offsets(experiment_list)
        batch_starts = [
            e.scan.get_image_range()[0] if e.scan else 0 for e in experiment_list
        ]
        effective_offsets = [o + s for o, s in zip(batch_offsets, batch_starts)]
        unique_offsets = set(effective_offsets)
    else:
        logger.debug("Keeping existing batches")
    image_ranges = get_image_ranges(experiment_list)
    if len(unique_offsets) != len(batch_offsets):
        raise ValueError(
            "Duplicate batch offsets detected: %s"
            % ", ".join(
                str(item) for item, count in Counter(batch_offsets).items() if count > 1
            )
        )

    # Create the mtz file
    mtz = gemmi.Mtz(with_base=True)
    mtz.title = f"From {env.dispatcher_name}"
    date_str = time.strftime("%Y-%m-%d at %H:%M:%S %Z")
    if time.strftime("%Z") != "GMT":
        date_str += time.strftime("  (%Y-%m-%d at %H:%M:%S %Z)", time.gmtime())
    mtz.history += [
        f"From {dials_version()}, run on {date_str}",
    ]

    # Create the right gemmi spacegroup from the crystal's cctbx space_group
    # via a Hall symbol
    hall = experiment_list[0].crystal.get_space_group().type().hall_symbol()
    ops = gemmi.symops_from_hall(hall)
    mtz.spacegroup = gemmi.find_spacegroup_by_ops(ops)

    # FIXME TODO for more than one experiment into an MTZ file:
    #
    # - add an epoch (or recover an epoch) from the scan and add this as an extra
    #   column to the MTZ file for scaling, so we know that the two lattices were
    #   integrated at the same time
    # ✓ decide a sensible BATCH increment to apply to the BATCH value between
    #   experiments and add this

    for id_ in expids_in_table.keys():
        # Grab our subset of the data
        loc = expids_in_list.index(
            expids_in_table[id_]
        )  # get strid and use to find loc in list
        experiment = experiment_list[loc]
        identifier = experiment.identifier

        if len(wavelengths) > 1:
            for i, wl in enumerate(wavelengths.values()):
                if identifier in wl.identifiers:
                    wavelength = wl.weighted_mean
                    dataset_id = i + 1
                    break
        else:
            wavelength = list(wavelengths.values())[0].weighted_mean
            dataset_id = 1

        reflections = reflection_table.select(reflection_table["id"] == id_)
        batch_offset = batch_offsets[loc]
        image_range = image_ranges[loc]
        reflections = assign_batches_to_reflections([reflections], [batch_offset])[0]
        experiment.data = dict(reflections)

        s0n = matrix.col(experiment.beam.get_s0()).normalize().elems
        logger.debug("Beam vector: {:.4f} {:.4f} {:.4f}".format(*s0n))

        add_batch_list(
            mtz,
            image_range,
            experiment,
            wavelength,
            dataset_id,
            batch_offset=batch_offset,
            force_static_model=force_static_model,
        )

        # Create the batch offset array. This gives us an experiment (id)-dependent
        # batch offset to calculate the correct batch from image number.
        experiment.data["batch_offset"] = flex.int(
            len(experiment.data["id"]), batch_offset
        )

        # Calculate whether we have a ROT value for this experiment, and set the column
        _, _, z = experiment.data["xyzcal.px"].parts()
        if experiment.scan:
            experiment.data["ROT"] = experiment.scan.get_angle_from_array_index(z)
        else:
            experiment.data["ROT"] = z

    mtz.set_cell_for_all(gemmi.UnitCell(*best_unit_cell.parameters()))

    # For multi-wave unmerged mtz, we add an empty dataset for each wavelength,
    # but only write the data into the final dataset (for unmerged the batches
    # link the unmerged data to the individual wavelengths).
    for wavelength in wavelengths.values():
        ds = mtz.add_dataset("FROMDIALS")
        ds.crystal_name = crystal_name
        ds.project_name = project_name
        ds.wavelength = wavelength.weighted_mean

    # Combine all of the experiment data columns before writing
    combined_data = {k: v.deep_copy() for k, v in experiment_list[0].data.items()}
    for experiment in experiment_list[1:]:
        for k, v in experiment.data.items():
            combined_data[k].extend(v)
    # ALL columns must be the same length
    assert len({len(v) for v in combined_data.values()}) == 1, "Column length mismatch"
    assert len(combined_data["id"]) == len(
        reflection_table["id"]
    ), "Lost rows in split/combine"

    # Write all the data and columns to the mtz file
    write_columns(mtz, combined_data)

    # Switch to ASU indices and sort file in standard order
    mtz.switch_to_asu_hkl()
    mtz.sort(5)

    logger.info(
        "Saving %s integrated reflections to %s", len(combined_data["id"]), filename
    )
    mtz.write_to_file(filename)
    log_summary(mtz)

    return mtz


def log_summary(mtz):
    """Log a summary of an MTZ object, based on the output of `gemmi mtz --dump`"""

    logger.info("Title: " + mtz.title)
    logger.info(f"Total Number of Datasets = {len(mtz.datasets)}\n")
    for ds in mtz.datasets:
        logger.info(
            f"Dataset {ds.id:4d}   {ds.project_name} > {ds.crystal_name} > {ds.dataset_name}:"
        )
        logger.info(
            "        cell  {:7g} {:7g} {:7g}  {:6g} {:6g} {:6g}".format(
                *ds.cell.parameters
            )
        )
        logger.info(f"  wavelength  {ds.wavelength:g}")
    logger.info(f"\nNumber of Columns = {len(mtz.columns)}")
    logger.info(
        f"Number of Reflections = {mtz.nreflections}",
    )
    logger.info(f"Number of Batches = {len(mtz.batches)}")
    logger.info(f"Missing values marked as: {mtz.valm}")
    logger.info(
        "Global Cell (obsolete): {:7g} {:7g} {:7g}  {:6g} {:6g} {:6g}".format(
            *mtz.cell.parameters
        )
    )
    mtz.update_reso()
    logger.info(
        f"Resolution: {mtz.resolution_high():.2f} - {mtz.resolution_low():.2f} A"
    )
    logger.info("Sort Order: {:d} {:d} {:d} {:d} {:d}".format(*mtz.sort_order))
    logger.info(f"Space Group: {mtz.spacegroup.hm}")
    logger.info(f"Space Group Number: {mtz.spacegroup.ccp4}")
    logger.info("Header info:")
    logger.info("Column    Type  Dataset    Min        Max")
    for col in mtz.columns:
        # col.min_value and col.max_value are not set, so we have to calculate them here
        logger.info(
            f"{col.label:<12s} {col.type} {col.dataset_id:2d} {np.nanmin(col.array):12.6g} {np.nanmax(col.array):10.6g}"
        )
    logger.info(f"History ({len(mtz.history)} lines):")
    for line in mtz.history:
        logger.info(line)


@dataclass
class WavelengthGroup:
    min_wl: float
    max_possible_wl: float
    identifiers: list[str] = field(default_factory=list)
    exp_nos: list[int] = field(default_factory=list)
    wavelengths: list[float] = field(default_factory=list)
    weighted_mean: float = 0

    def add_experiment(self, identifier: str, loc_in_list: int, wl: float) -> None:
        self.identifiers.append(identifier)
        self.exp_nos.append(loc_in_list)
        self.wavelengths.append(wl)

    def calculate_weighted_mean(
        self, reflection_tables: list[flex.reflection_table]
    ) -> None:
        n, nw = (0, 0)
        for i, w in zip(self.identifiers, self.wavelengths):
            for table in reflection_tables:
                refls = table.select_on_experiment_identifiers([i])
                n_this = refls.select(
                    refls.get_flags(refls.flags.integrated, all=False)
                ).size()
                if n_this:
                    n += n_this
                    nw += n_this * w
                    break
        if n:
            self.weighted_mean = nw / n


def match_wavelengths(experiments, absolute_tolerance=1e-4):
    wavelengths = {}
    for i, x in enumerate(experiments):
        w = x.beam.get_wavelength()
        matches = [isclose(w, k, abs_tol=absolute_tolerance) for k in wavelengths]
        if not any(matches):
            wavelengths[w] = WavelengthGroup(w, w + absolute_tolerance)
            wavelengths[w].add_experiment(x.identifier, i, w)
        else:
            match_w = list(wavelengths.keys())[matches.index(True)]
            wavelengths[match_w].add_experiment(x.identifier, i, w)
    return wavelengths


def convert_to_cambridge(experiments):
    """Rotate the geometry of an experiment list to match the Cambridge frame,
    in which X is along the idealized X-ray beam and Z is along the primary
    rotation axis"""

    # First handle the potential for shared experiment models - we don't
    # want to apply multiple transformations to shared models, simplest way is
    # to make a copy for each experiment.
    n_expt = len(experiments)
    if len(experiments.crystals()) < n_expt:
        for expt in experiments:
            expt.crystal = deepcopy(expt.crystal)
    if len(experiments.beams()) < n_expt:
        for expt in experiments:
            expt.beam = deepcopy(expt.beam)
    if len(experiments.detectors()) < n_expt:
        for expt in experiments:
            expt.detector = deepcopy(expt.detector)
    if any(experiments.goniometers()) and len(experiments.goniometers()) < n_expt:
        for expt in experiments:
            expt.goniometer = deepcopy(expt.goniometer)

    for expt in experiments:
        if expt.goniometer:
            primary_axis = matrix.col(expt.goniometer.get_rotation_axis_datum())
        else:
            primary_axis = matrix.col((1.0, 0.0, 0.0))
        us0 = expt.beam.get_unit_s0()

        R = align_reference_frame(primary_axis, (0, 0, 1), us0, (1, 0, 0))
        axis_angle = r3_rotation_axis_and_angle_from_matrix(R)
        axis = matrix.col(axis_angle.axis)
        angle = axis_angle.angle()
        logger.debug(
            f"Rotating experiment{'s' if len(experiments) else ''} about axis {axis.elems} by {np.degrees(angle):.2f}°"
        )
        expt.detector.rotate_around_origin(axis, angle, deg=False)
        expt.beam.rotate_around_origin(axis, angle, deg=False)

        # For the goniometer, each component needs transformation
        if expt.goniometer:
            F = matrix.sqr(expt.goniometer.get_fixed_rotation())
            expt.goniometer.set_fixed_rotation(R * F * R.transpose())
            expt.goniometer.set_rotation_axis_datum(R * primary_axis)
            S = matrix.sqr(expt.goniometer.get_setting_rotation())
            expt.goniometer.set_setting_rotation(R * S * R.transpose())

        if expt.crystal is not None:
            expt.crystal = rotate_crystal(expt.crystal, R, axis, angle)

    return experiments


def rotate_crystal(crystal, Rmat, axis, angle):
    Amats = []
    if crystal.num_scan_points > 0:
        scan_pts = list(range(crystal.num_scan_points))
        Amats = [
            Rmat
            * matrix.sqr(crystal.get_U_at_scan_point(t))
            * matrix.sqr(crystal.get_B_at_scan_point(t))
            for t in scan_pts
        ]

    crystal.rotate_around_origin(axis, angle, deg=False)
    if Amats:
        crystal.set_A_at_scan_points(Amats)

    return crystal
