from __future__ import annotations

import json
import logging
import math
import os
import pathlib
import sys
from dataclasses import dataclass, field
from multiprocessing import Pool
from typing import Any

import numpy as np

from dxtbx.model import Crystal, Experiment, ExperimentList
from libtbx import Auto, phil

from dials.algorithms.indexing import DialsIndexError
from dials.algorithms.indexing.indexer import Indexer
from dials.algorithms.indexing.max_cell import find_max_cell
from dials.array_family import flex
from dials.util.combine_experiments import CombineWithReference

RAD2DEG = 180 / math.pi

logger = logging.getLogger("dials")


@dataclass
class InputToIndex:
    reflection_table: flex.reflection_table = None
    experiment: Experiment = None
    parameters: Any = None
    image_identifier: str = ""
    image_no: int = 0
    method_list: list[str] = field(default_factory=list)
    known_crystal_models: list[Crystal] | None = None
    imageset_no: int = 0


@dataclass
class IndexingResult:
    image: str
    image_no: int
    reflection_table: flex.reflection_table = None
    experiments: ExperimentList = None
    n_strong: int = 0
    n_indexed: list[int] = field(default_factory=list)
    identifiers: list[str] = field(default_factory=list)
    rmsd_x: list[float] = field(default_factory=list)
    rmsd_y: list[float] = field(default_factory=list)
    rmsd_dpsi: list[float] = field(default_factory=list)
    imageset_no: int = 0
    unindexed_experiment: Experiment = None


loggers_to_disable = [
    "dials.algorithms.refinement.reflection_processor",
    "dials.algorithms.refinement.refiner",
    "dials.algorithms.refinement.reflection_manager",
    "dials.algorithms.refinement.outlier_detection.outlier_base",
    "dials.algorithms.indexing.stills_indexer",
    "dials.algorithms.indexing.nave_parameters",
    "dials.algorithms.indexing.basis_vector_search.real_space_grid_search",
    "dials.algorithms.indexing.basis_vector_search.combinations",
    "dials.algorithms.indexing.indexer",
    "dials.algorithms.indexing.lattice_search",
    "dials.algorithms.indexing.lattice_search.low_res_spot_match",
    "dials.algorithms.indexing.lattice_search.ffb_indexer",
]
debug_loggers_to_disable = [
    "dials.algorithms.indexing.symmetry",
    "dials.algorithms.indexing.lattice_search",
    "dials.algorithms.refinement.engine",
]


class manage_loggers:
    """
    A contextmanager for reducing logging levels for the underlying code of
    parallel ssx programs.
    If individual_log_verbosity < 2, only timing info will logged.
    If individual_log_verbosity = 2, then the info logs will be shown, but not debug
    logs
    If individual_log_verbosity > 2, then debug logs will be shown.

    """

    def __init__(
        self, individual_log_verbosity, loggers_to_disable, debug_loggers_to_diable=None
    ):
        self.individual_log_verbosity = individual_log_verbosity
        self.loggers = loggers_to_disable
        if debug_loggers_to_disable:
            self.debug_loggers = debug_loggers_to_disable
        else:
            self.debug_loggers = []

    def __enter__(self):
        if self.individual_log_verbosity < 2:
            for name in self.loggers:
                logging.getLogger(name).disabled = True
        elif self.individual_log_verbosity == 2:
            for name in self.loggers:
                logging.getLogger(name).setLevel(logging.INFO)
            for name in self.debug_loggers:
                logging.getLogger(name).setLevel(logging.INFO)

    def __exit__(self, ex_type, ex_value, ex_traceback):
        # Re-enable the disabled loggers or reset logging levels.
        if self.individual_log_verbosity < 2:
            for logname in self.loggers:
                logging.getLogger(logname).disabled = False
        elif self.individual_log_verbosity == 2:
            for name in self.loggers:
                logging.getLogger(name).setLevel(logging.DEBUG)
            for name in self.debug_loggers:
                logging.getLogger(name).setLevel(logging.DEBUG)


def index_one(
    experiment: Experiment,
    reflection_table: flex.reflection_table,
    params: phil.scope_extract,
    method_list: list[str],
    image_no: int,
    known_crystal_models: list[Crystal] = None,
) -> tuple[ExperimentList, flex.reflection_table] | tuple[bool, bool]:
    elist = ExperimentList([experiment])
    params.indexing.nproc = 1  # make sure none of the processes try to spawn multiprocessing within existing multiprocessing.
    for method in method_list:
        params.indexing.method = method
        idxr = Indexer.from_parameters(
            reflection_table,
            elist,
            params=params,
            known_crystal_models=known_crystal_models,
        )
        try:
            idxr.index()
        except (DialsIndexError, AssertionError) as e:
            logger.info(
                f"Image {image_no+1}: Failed to index with {method} method, error: {e}"
            )
            if method == method_list[-1]:
                return None, None
        else:
            logger.info(
                f"Image {image_no+1}: Indexed {idxr.refined_reflections.size()}/{reflection_table.size()} spots with {method} method."
            )
            return idxr.refined_experiments, idxr.refined_reflections


def wrap_index_one(input_to_index: InputToIndex) -> IndexingResult:
    # First unpack the input and run the function
    expts, table = index_one(
        input_to_index.experiment,
        input_to_index.reflection_table,
        input_to_index.parameters,
        input_to_index.method_list,
        input_to_index.image_no,
        input_to_index.known_crystal_models,
    )

    # Now calculate some useful quantities for the result
    n_strong = input_to_index.reflection_table.get_flags(
        input_to_index.reflection_table.flags.strong
    ).count(True)
    if expts and table:
        result = IndexingResult(
            input_to_index.image_identifier,
            input_to_index.image_no,
            table,
            expts,
            n_strong=n_strong,
            imageset_no=input_to_index.imageset_no,
            unindexed_experiment=input_to_index.experiment,
        )
        for id_, identifier in table.experiment_identifiers():
            selr = table.select(table["id"] == id_)
            calx, caly, _ = selr["xyzcal.px"].parts()
            obsx, obsy, _ = selr["xyzobs.px.value"].parts()
            if "delpsical.rad" in selr:
                delpsi = selr["delpsical.rad"]
                rmsd_z = flex.mean(((delpsi) * RAD2DEG) ** 2) ** 0.5
            else:
                rmsd_z = 0.0
            rmsd_x = flex.mean((calx - obsx) ** 2) ** 0.5
            rmsd_y = flex.mean((caly - obsy) ** 2) ** 0.5
            n_id_ = calx.size()
            result.n_indexed.append(n_id_)
            result.rmsd_x.append(rmsd_x)
            result.rmsd_y.append(rmsd_y)
            result.rmsd_dpsi.append(rmsd_z)
            result.identifiers.append(identifier)
    else:
        result = IndexingResult(
            input_to_index.image_identifier,
            input_to_index.image_no,
            n_strong=n_strong,
            imageset_no=input_to_index.imageset_no,
            unindexed_experiment=input_to_index.experiment,
        )

    # If chosen, output a message to json to show live progress
    if input_to_index.parameters.output.nuggets:
        msg = {
            "image_no": input_to_index.image_no + 1,
            "n_indexed": (len(table) if table else 0),
            "n_strong": n_strong,
            "n_cryst": (len(expts) if expts else 0),
            "image": input_to_index.image_identifier,
        }
        with open(
            input_to_index.parameters.output.nuggets
            / f"nugget_index_{msg['image']}.json",
            "w",
        ) as f:
            f.write(json.dumps(msg))

    return result


def index_all_concurrent(
    experiments: ExperimentList,
    reflections: list[flex.reflection_table],
    params: phil.scope_extract,
    method_list: list[str],
) -> tuple[ExperimentList, flex.reflection_table, dict]:
    input_iterable = []
    results_summary = {
        i: [] for i in range(len(experiments))
    }  # create to give results in order

    # Create a suitable iterable for passing to pool.map
    n = 0
    original_isets = list(experiments.imagesets())
    identifiers_to_scans = {expt.identifier: expt.scan for expt in experiments}
    for n_iset, iset in enumerate(experiments.imagesets()):
        for i in range(len(iset)):
            refl_index = i + n
            if reflections[refl_index]:
                expt = experiments[refl_index]
                input_iterable.append(
                    InputToIndex(
                        reflection_table=reflections[refl_index],
                        experiment=expt,
                        parameters=params,
                        image_identifier=pathlib.Path(
                            iset.get_image_identifier(i)
                        ).name,
                        image_no=refl_index,
                        method_list=method_list,
                        imageset_no=n_iset,
                    )
                )
            else:  # experiments that have already been filtered
                results_summary[refl_index].append(
                    {
                        "Image": pathlib.Path(iset.get_image_identifier(i)).name,
                        "n_indexed": 0,
                        "n_strong": 0,
                    }
                )
        n += len(iset)

    with open(os.devnull, "w") as devnull:
        sys.stdout = devnull  # block printing from rstbx
        with manage_loggers(
            params.individual_log_verbosity,
            loggers_to_disable,
            debug_loggers_to_disable,
        ):
            if params.indexing.nproc > 1:
                with Pool(params.indexing.nproc) as pool:
                    results: list[IndexingResult] = pool.map(
                        wrap_index_one, input_iterable
                    )
            else:
                results: list[IndexingResult] = [
                    wrap_index_one(i) for i in input_iterable
                ]

    sys.stdout = sys.__stdout__
    # prepare tables for output
    indexed_experiments, indexed_reflections = _join_indexing_results(
        results, experiments, original_isets, identifiers_to_scans
    )

    results_summary = _add_results_to_summary_dict(results_summary, results)

    return indexed_experiments, indexed_reflections, results_summary


def _join_indexing_results(
    results: list[IndexingResult],
    experiments,
    original_isets,
    identifiers_to_scans,
) -> tuple[ExperimentList, flex.reflection_table]:
    indexed_experiments = ExperimentList()
    indexed_reflections = flex.reflection_table()

    use_beam = None
    use_gonio = None
    if len(experiments.beams()) == 1:
        use_beam = experiments.beams()[0]
    if len(experiments.goniometers()):  # need a placeholder gonio
        use_gonio = experiments.goniometers()[0]

    n_tot = 0
    for res in results:
        if res.n_indexed:
            identifier = res.unindexed_experiment.identifier
            scan = identifiers_to_scans[identifier]
            for expt in res.experiments:
                expt.scan = scan
                expt.imageset = original_isets[res.imageset_no]
                res.reflection_table["imageset_id"] = flex.int(
                    res.reflection_table.size(), res.imageset_no
                )
                if use_beam:
                    expt.beam = use_beam
                if use_gonio:
                    expt.goniometer = use_gonio

            indexed_experiments.extend(res.experiments)
            table = res.reflection_table
            ids_map = dict(table.experiment_identifiers())
            for k in table.experiment_identifiers().keys():
                del table.experiment_identifiers()[k]
            table["id"] += n_tot
            for k, v in ids_map.items():
                table.experiment_identifiers()[k + n_tot] = v
            n_tot += len(ids_map.keys())
            indexed_reflections.extend(table)

    indexed_reflections.assert_experiment_identifiers_are_consistent(
        indexed_experiments
    )
    return indexed_experiments, indexed_reflections


def _add_results_to_summary_dict(
    results_summary: dict, results: list[IndexingResult]
) -> dict:
    # convert to results_summary dict current format
    for res in results:
        if res.n_indexed:
            for j in range(len(res.n_indexed)):
                results_summary[res.image_no].append(
                    {
                        "Image": res.image,
                        "identifier": res.identifiers[j],
                        "n_indexed": res.n_indexed[j],
                        "n_strong": res.n_strong,
                        "RMSD_X": res.rmsd_x[j],
                        "RMSD_Y": res.rmsd_y[j],
                        "RMSD_dPsi": res.rmsd_dpsi[j],
                    }
                )
        else:
            results_summary[res.image_no].append(
                {
                    "Image": res.image,
                    "n_indexed": 0,
                    "n_strong": res.n_strong,
                }
            )
    return results_summary


def preprocess(
    experiments: ExperimentList,
    observed: flex.reflection_table,
    params: phil.scope_extract,
) -> tuple[list[flex.reflection_table], phil.scope_extract, list[str]]:
    reflections = observed.split_by_experiment_id()

    if len(reflections) != len(experiments):
        # spots may not have been found on every image. In this case, the length
        # of the list of reflection tables will be less than the length of experiments.
        # Add in empty items to the list, so that this can be reported on
        obs = set(observed["id"])
        no_refls = set(range(len(experiments))).difference(obs)
        # need to handle both cases where lots have no refls, or only a few do
        full_list = [None] * len(experiments)
        for i, o in enumerate(obs):
            full_list[o] = reflections[i]
        reflections = full_list
        logger.info(f"Filtered {len(no_refls)} images with no spots found.")
        if len(experiments) != len(reflections):
            raise ValueError(
                f"Unequal number of reflection tables {len(reflections)} and experiments {len(experiments)}"
            )

    # Calculate necessary quantities
    n_filtered_out = 0
    filtered_out_images = []
    for i, (refl, experiment) in enumerate(zip(reflections, experiments)):
        if refl:
            if refl.size() >= params.min_spots:
                elist = ExperimentList([experiment])
                refl["imageset_id"] = flex.int(
                    refl.size(), 0
                )  # needed for centroid_px_to_mm
                refl.centroid_px_to_mm(elist)
                refl.map_centroids_to_reciprocal_space(elist)
            else:
                n_filtered_out += 1
                reflections[i] = None
                filtered_out_images.append(i + 1)
        else:
            filtered_out_images.append(i + 1)
    if n_filtered_out:
        logger.info(
            f"Filtered {n_filtered_out} images with fewer than {params.min_spots} spots"
        )
        if params.output.nuggets:
            jstr = json.dumps({"filtered_images": filtered_out_images})
            with open(
                params.output.nuggets / "nugget_index_filtered_images.json", "w"
            ) as f:
                f.write(jstr)

    if n_filtered_out == len(reflections):
        logger.info("All images filtered out, none left to attempt indexing")
        return reflections, params, []

    # Determine the max cell if applicable
    if (params.indexing.max_cell is Auto) and (
        not params.indexing.known_symmetry.unit_cell
    ):
        if params.individual_log_verbosity <= 2:  # suppress the max cell debug log
            logging.getLogger("dials.algorithms.indexing.max_cell").setLevel(
                logging.INFO
            )
        max_cells = []
        for refl in reflections:
            if refl:
                try:
                    max_cells.append(find_max_cell(refl).max_cell)
                except (DialsIndexError, AssertionError):
                    pass
        if not max_cells:
            raise DialsIndexError("Unable to find a max cell for any images")

        sorted_cells = np.sort(np.array(max_cells))
        n_cells = len(sorted_cells)
        if n_cells > 20:
            centile_95_pos = int(math.floor(0.95 * n_cells))
            limit = sorted_cells[centile_95_pos]
            logger.info(f"Setting max cell to {limit:.1f} " + "\u212b")
            params.indexing.max_cell = limit
        else:
            params.indexing.max_cell = sorted_cells[-1]
            logger.info(f"Setting max cell to {sorted_cells[-1]:.1f} " + "\u212b")

    # Determine which methods to try
    method_list = params.method
    if "real_space_grid_search" in method_list:
        if not params.indexing.known_symmetry.unit_cell:
            logger.info("No unit cell given, real_space_grid_search will not be used")
            method_list.remove("real_space_grid_search")
    methods = ", ".join(method_list)
    pl = "s" if (len(method_list) > 1) else ""
    logger.info(f"Attempting indexing with {methods} method{pl}")

    return reflections, params, method_list


def index(
    experiments: ExperimentList,
    observed: flex.reflection_table,
    params: phil.scope_extract,
) -> tuple[ExperimentList, flex.reflection_table, dict]:
    if params.output.nuggets:
        params.output.nuggets = pathlib.Path(
            params.output.nuggets
        )  # needed if not going via cli
        if not params.output.nuggets.is_dir():
            logger.warning(
                "output.nuggets not recognised as a valid directory path, no nuggets will be output"
            )
            params.output.nuggets = None

    if any(s and not s.is_still() for s in experiments.scans()):
        raise DialsIndexError(
            "Not all experiments are stills. For rotation data, use the dials.index program."
        )

    reflection_tables, params, method_list = preprocess(experiments, observed, params)

    # Do the indexing and generate a summary of results
    indexed_experiments, indexed_reflections, results_summary = index_all_concurrent(
        experiments,
        reflection_tables,
        params,
        method_list,
    )

    # combine detector models if all the same (i.e. haven't been refined in indexing),
    # to enable use case of joint refinement.
    detector_0 = indexed_experiments[0].detector if indexed_experiments else None
    if (len(indexed_experiments.detectors()) > 1) and all(
        d == detector_0 for d in indexed_experiments.detectors()[1:]
    ):
        combine = CombineWithReference(detector=detector_0)
        elist = ExperimentList()
        for expt in indexed_experiments:
            elist.append(combine(expt))
        indexed_experiments = elist

    return indexed_experiments, indexed_reflections, results_summary
