from __future__ import annotations

import logging
import math

import pkg_resources

import iotbx.phil
import libtbx
from cctbx import sgtbx
from dxtbx.model import ExperimentList, ImageSequence, tof_helpers

import dials.util
from dials.algorithms.indexing import (
    DialsIndexError,
    DialsIndexRefineError,
    assign_indices,
)
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from dials.algorithms.indexing.max_cell import find_max_cell
from dials.algorithms.indexing.symmetry import SymmetryHandler
from dials.algorithms.refinement import DialsRefineConfigError, DialsRefineRuntimeError
from dials.array_family import flex
from dials.util.multi_dataset_handling import generate_experiment_identifiers

logger = logging.getLogger(__name__)


max_cell_phil_str = """\
max_cell_estimation
  .expert_level = 1
{
  filter_ice = True
    .type = bool
    .help = "Filter out reflections at typical ice ring resolutions"
            "before max_cell estimation."
  filter_overlaps = True
    .type = bool
    .help = "Filter out reflections with overlapping bounding boxes before"
            "max_cell estimation."
  overlaps_border = 0
    .type = int(value_min=0)
    .help = "Optionally add a border around the bounding boxes before finding"
            "overlaps."
  multiplier = 1.3
    .type = float(value_min=0)
    .help = "Multiply the estimated maximum basis vector length by this value."
    .expert_level = 2
  step_size = 45
    .type = float(value_min=0)
    .help = "Step size, in degrees, of the blocks used to perform the max_cell "
            "estimation."
    .expert_level = 2
  nearest_neighbor_percentile = None
    .type = float(value_min=0, value_max=1)
    .help = "Percentile of NN histogram to use for max cell determination."
    .expert_level = 2
  histogram_binning = linear *log
    .type = choice
    .help = "Choose between linear or logarithmic bins for nearest neighbour"
            "histogram analysis."
  nn_per_bin = 5
    .type = int(value_min=1)
    .help = "Target number of nearest neighbours per histogram bin."
  max_height_fraction = 0.25
    .type = float(value_min=0, value_max=1)
    .expert_level=2
}
"""

phil_str = (
    """\
indexing {
  nproc = 1
    .type = int(value_min=1)
    .help = "The number of processes to use."
  mm_search_scope = 4.0
    .help = "Global radius of origin offset search."
    .type = float(value_min=0)
    .expert_level = 1
  wide_search_binning = 2
    .help = "Modify the coarseness of the wide grid search for the beam centre."
    .type = float(value_min=0)
    .expert_level = 1
  min_cell_volume = 25
    .type = float(value_min=0)
    .help = "Minimum unit cell volume (in Angstrom^3)."
    .expert_level = 1
  min_cell = 3
    .type = float(value_min=0)
    .help = "Minimum length of candidate unit cell basis vectors (in Angstrom)."
    .expert_level = 1
  max_cell = Auto
    .type = float(value_min=0)
    .help = "Maximum length of candidate unit cell basis vectors (in Angstrom)."
    .expert_level = 1
  %s
  sigma_phi_deg = None
    .type = float(value_min=0)
    .help = "Override the phi sigmas for refinement. Mainly intended for single-shot"
            "rotation images where the phi sigma is almost certainly incorrect."
    .expert_level = 2

  known_symmetry {
    space_group = None
      .type = space_group
      .help = "Target space group for indexing."
    unit_cell = None
      .type = unit_cell
      .help = "Target unit cell for indexing."
    relative_length_tolerance = 0.1
      .type = float
      .help = "Relative tolerance for unit cell lengths in unit cell comparison."
      .expert_level = 1
    absolute_angle_tolerance = 5
      .type = float
      .help = "Angular tolerance (in degrees) in unit cell comparison."
      .expert_level = 1
    max_delta = 5
      .type = float(value_min=0)
      .help = "Maximum allowed Le Page delta used in searching for basis vector"
              "combinations that are consistent with the given symmetry."
      .expert_level = 1
  }

  index_assignment {
    method = *simple local
      .type = choice
      .help = "Choose between simple 'global' index assignment and xds-style "
              "'local' index assignment."
      .expert_level = 1
    simple {
      hkl_tolerance = 0.3
        .type = float(value_min=0, value_max=0.5)
        .help = "Maximum allowable deviation from integer-ness for assigning "
                "a miller index to a reciprocal lattice vector."
    }
    local
      .expert_level = 1
    {
      epsilon = 0.05
        .type = float
        .help = "This corresponds to the xds parameter INDEX_ERROR="
      delta = 8
        .type = int
        .help = "This corresponds to the xds parameter INDEX_MAGNITUDE="
      l_min = 0.8
        .type = float
        .help = "This corresponds to the xds parameter INDEX_QUALITY="
      nearest_neighbours = 20
        .type = int(value_min=1)
    }
  }
  check_misindexing {
    grid_search_scope = 0
      .type = int
      .help = "Search scope for testing misindexing on h, k, l."
  }
  debug = False
    .type = bool
    .expert_level = 1
  combine_scans = False
    .type = bool
    .expert_level = 1
  refinement_protocol {
    mode = *refine_shells repredict_only None
      .type = choice
      .expert_level = 1
      .help = "refine_shells: if using sequences indexer, refine in increasing"
              "resolution cutoffs after indexing, if using stills indexer,"
              "refine all data up to d_min_start resolution once only."
              "repredict_only: do not refine after indexing, just update spot"
              "predictions."
              "None: do not refine and do not update spot predictions."
    n_macro_cycles = 5
      .type = int(value_min=1)
      .help = "Maximum number of macro cycles of refinement, reindexing all"
              "reflections using updated geometry at the beginning of each"
              "cycle. Does not apply to stills.indexer=stills."
    d_min_step = Auto
      .type = float(value_min=0.0)
      .help = "Reduction per step in d_min for reflections to include"
              "in refinement. Does not apply to stills.indexer=stills."
    d_min_start = None
      .type = float(value_min=0.0)
      .help = "For sequences/stills indexer, the lower limit of d-spacing"
              "of reflections used in the first/the only round of refinement."
    d_min_final = None
      .type = float(value_min=0.0)
      .help = "Do not ever include reflections below this value in refinement."
              "Does not apply to stills.indexer=stills."
    disable_unit_cell_volume_sanity_check = False
      .type = bool
      .help = "Disable sanity check on unrealistic increases in unit cell volume"
              "during refinement. Does not apply to stills.indexer=stills."
      .expert_level = 1
  }
  multiple_lattice_search
    .expert_level = 1
  {
    recycle_unindexed_reflections_cutoff = 0.1
      .type = float(value_min=0, value_max=1)
      .help = "Attempt another cycle of indexing on the unindexed reflections "
              "if more than the fraction of input reflections are unindexed."
    minimum_angular_separation = 5
      .type = float(value_min=0)
      .help = "The minimum angular separation (in degrees) between two lattices."
    max_lattices = 1
      .type = int
    cluster_analysis {
      method = *dbscan hcluster
        .type = choice
      hcluster {
        linkage {
          method = *ward
            .type = choice
          metric = *euclidean
            .type = choice
        }
        cutoff = 15
          .type = float(value_min=0)
        cutoff_criterion = *distance inconsistent
          .type = choice
      }
      dbscan {
        eps = 0.05
          .type = float(value_min=0.0)
        min_samples = 30
          .type = int(value_min=1)
      }
      min_cluster_size = 20
        .type = int(value_min=0)
      intersection_union_ratio_cutoff = 0.4
        .type = float(value_min=0.0, value_max=1.0)
    }
  }
  stills {
    indexer = *Auto stills sequences
      .type = choice
      .help = Use the stills or sequences indexer.  Auto: choose based on the input \
              imagesets (stills or sequences).
      .expert_level = 1
    ewald_proximity_resolution_cutoff = 2.0
      .type = float
      .help = For still images, this high-resolution cutoff is used to calculate
      .help = the acceptable volume of reciprocal space for spot prediction
    refine_all_candidates = True
      .type = bool
      .help = If False, no attempt is made to refine the model from initial basis \
              vector selection. The indexing solution with the best RMSD is chosen.
    candidate_outlier_rejection = True
      .type = bool
      .expert_level = 1
      .help = If True, while refining candidate basis solutions, also apply Sauter/ \
              Poon (2010) outlier rejection
    refine_candidates_with_known_symmetry = False
      .type = bool
      .expert_level = 2
      .help = If False, when choosing the best set of candidate basis solutions, \
              refine the candidates in the P1 setting. If True, after indexing \
              in P1, convert the candidates to the known symmetry and apply the \
              corresponding change of basis to the indexed reflections.
    rmsd_min_px = 2
      .type = float
      .help = Minimum acceptable RMSD for choosing candidate basis solutions \
              (in pixels)
    ewald_proximal_volume_max = 0.0025
      .type = float
      .help = Maximum acceptable ewald proximal volume when choosing candidate \
              basis solutions
    isoforms
      .help = Constrain the unit cell to specific values during refinement after initial indexing.
      .multiple=True
    {
      name=None
        .type=str
      cell=None
        .type=unit_cell
      lookup_symbol=None
        .type=str
        .help=The sgtbx lookup symbol of the reflections pointgroup
      rmsd_target_mm=None
        .type=float
        .help=Maximum acceptable DIALS positional rmsd, in mm
      beam_restraint=None
        .type=floats(size=2)
        .help=Known beam position in mm X,Y, rmsd_target_mm is reused here as a circle of confusion
        .help=to assure that no images are accepted where the lattice is misindexed by a unit shift.
    }
    set_domain_size_ang_value = None
      .type=float
      .help=If specified, will set the domain size ang value and override the value determined from nave refinement
    set_mosaic_half_deg_value = None
      .type=float
      .help=If specified, will set the mosaic half degree value and override the value determined from nave refinement
  }
}
"""
    % max_cell_phil_str
)

phil_scope = iotbx.phil.parse(phil_str, process_includes=True)


class Indexer:
    def __init__(self, reflections, experiments, params):
        self.reflections = reflections
        self.experiments = experiments

        self.params = params.indexing
        self.all_params = params
        self.refined_experiments = None
        self.hkl_offset = None

        if self.params.index_assignment.method == "local":
            self._assign_indices = assign_indices.AssignIndicesLocal(
                epsilon=self.params.index_assignment.local.epsilon,
                delta=self.params.index_assignment.local.delta,
                l_min=self.params.index_assignment.local.l_min,
                nearest_neighbours=self.params.index_assignment.local.nearest_neighbours,
            )
        else:
            self._assign_indices = assign_indices.AssignIndicesGlobal(
                tolerance=self.params.index_assignment.simple.hkl_tolerance
            )

        if self.all_params.refinement.reflections.outlier.algorithm in (
            "auto",
            libtbx.Auto,
        ):
            if (
                self.experiments[0].goniometer is None
                or self.experiments[0].scan is None
                or self.experiments[0].scan.is_still()
            ):
                self.all_params.refinement.reflections.outlier.algorithm = "sauter_poon"
            else:
                # different default to dials.refine
                # tukey is faster and more appropriate at the indexing step
                self.all_params.refinement.reflections.outlier.algorithm = "tukey"

        for expt in self.experiments[1:]:
            if expt.detector.is_similar_to(self.experiments[0].detector):
                expt.detector = self.experiments[0].detector
            if expt.goniometer is not None and expt.goniometer.is_similar_to(
                self.experiments[0].goniometer
            ):
                expt.goniometer = self.experiments[0].goniometer
                # can only share a beam if we share a goniometer?
                if expt.beam.is_similar_to(self.experiments[0].beam):
                    expt.beam = self.experiments[0].beam
                if self.params.combine_scans and expt.scan == self.experiments[0].scan:
                    expt.scan = self.experiments[0].scan

        if "flags" in self.reflections:
            strong_sel = self.reflections.get_flags(self.reflections.flags.strong)
            if strong_sel.count(True) > 0:
                self.reflections = self.reflections.select(strong_sel)
        if "flags" not in self.reflections or strong_sel.count(True) == 0:
            # backwards compatibility for testing
            self.reflections.set_flags(
                flex.size_t_range(len(self.reflections)), self.reflections.flags.strong
            )

        self._setup_symmetry()
        self.d_min = None

        self.setup_indexing()

    @staticmethod
    def from_parameters(
        reflections, experiments, known_crystal_models=None, params=None
    ):
        if known_crystal_models is not None:
            from dials.algorithms.indexing.known_orientation import (
                IndexerKnownOrientation,
            )

            if params.indexing.known_symmetry.space_group is None:
                params.indexing.known_symmetry.space_group = (
                    known_crystal_models[0].get_space_group().info()
                )
            idxr = IndexerKnownOrientation(
                reflections, experiments, params, known_crystal_models
            )
        else:
            has_stills = False
            has_sequences = False
            for expt in experiments:
                if isinstance(expt.imageset, ImageSequence):
                    has_sequences = True
                else:
                    has_stills = True

            if has_stills and has_sequences:
                raise ValueError(
                    "Please provide only stills or only sequences, not both"
                )

            use_stills_indexer = has_stills

            if not (
                params.indexing.stills.indexer is libtbx.Auto
                or params.indexing.stills.indexer.lower() == "auto"
            ):
                if params.indexing.stills.indexer == "stills":
                    use_stills_indexer = True
                elif params.indexing.stills.indexer == "sequences":
                    use_stills_indexer = False
                else:
                    assert False

            if params.indexing.basis_vector_combinations.max_refine is libtbx.Auto:
                if use_stills_indexer:
                    params.indexing.basis_vector_combinations.max_refine = 5
                else:
                    params.indexing.basis_vector_combinations.max_refine = 50

            if use_stills_indexer:
                # Ensure the indexer and downstream applications treat this as set of stills
                from dxtbx.imageset import ImageSet  # , MemImageSet

                for experiment in experiments:
                    # Elsewhere, checks are made for ImageSequence when picking between algorithms
                    # specific to rotations vs. stills, so here reset any ImageSequences to stills.
                    # Note, dials.stills_process resets ImageSequences to ImageSets already,
                    # and it's not free (the ImageSet cache is dropped), only do it if needed
                    if isinstance(experiment.imageset, ImageSequence):
                        experiment.imageset = ImageSet(
                            experiment.imageset.data(), experiment.imageset.indices()
                        )
                    # if isinstance(imageset, MemImageSet):
                    #   imageset = MemImageSet(imagesequence._images, imagesequence.indices())
                    # else:
                    #   imageset = ImageSet(imagesequence.reader(), imagesequence.indices())
                    #   imageset._models = imagesequence._models
                    experiment.imageset.set_scan(None)
                    experiment.imageset.set_goniometer(None)
                    experiment.scan = None
                    experiment.goniometer = None

            IndexerType = None
            for entry_point in pkg_resources.iter_entry_points(
                "dials.index.basis_vector_search"
            ):
                if params.indexing.method == entry_point.name:
                    if use_stills_indexer:
                        # do something
                        from dials.algorithms.indexing.stills_indexer import (
                            StillsIndexerBasisVectorSearch as IndexerType,
                        )
                    else:
                        from dials.algorithms.indexing.lattice_search import (
                            BasisVectorSearch as IndexerType,
                        )

            if IndexerType is None:
                for entry_point in pkg_resources.iter_entry_points(
                    "dials.index.lattice_search"
                ):
                    if params.indexing.method == entry_point.name:
                        if use_stills_indexer:
                            from dials.algorithms.indexing.stills_indexer import (
                                StillsIndexerLatticeSearch as IndexerType,
                            )
                        else:
                            from dials.algorithms.indexing.lattice_search import (
                                LatticeSearch as IndexerType,
                            )

            assert IndexerType is not None

            idxr = IndexerType(reflections, experiments, params=params)

        return idxr

    def _setup_symmetry(self):
        target_unit_cell = self.params.known_symmetry.unit_cell
        target_space_group = self.params.known_symmetry.space_group
        if target_space_group is not None:
            target_space_group = target_space_group.group()
        else:
            target_space_group = sgtbx.space_group()
            self.params.known_symmetry.space_group = target_space_group.info()
        self._symmetry_handler = SymmetryHandler(
            unit_cell=target_unit_cell,
            space_group=target_space_group,
            max_delta=self.params.known_symmetry.max_delta,
        )
        return

    def setup_indexing(self):
        if len(self.reflections) == 0:
            raise DialsIndexError("No reflections left to index!")

        if "imageset_id" not in self.reflections:
            self.reflections["imageset_id"] = self.reflections["id"]
        self.reflections.centroid_px_to_mm(self.experiments)
        self.reflections.map_centroids_to_reciprocal_space(self.experiments)
        self.reflections.calculate_entering_flags(self.experiments)

        self.find_max_cell()

        if self.params.sigma_phi_deg is not None:
            var_x, var_y, _ = self.reflections["xyzobs.mm.variance"].parts()
            var_phi_rad = flex.double(
                var_x.size(), (math.pi / 180 * self.params.sigma_phi_deg) ** 2
            )
            self.reflections["xyzobs.mm.variance"] = flex.vec3_double(
                var_x, var_y, var_phi_rad
            )

        if self.params.debug:
            self._debug_write_reciprocal_lattice_points_as_pdb()

        self.reflections["id"] = flex.int(len(self.reflections), -1)

    def index(self):
        experiments = ExperimentList()

        had_refinement_error = False
        have_similar_crystal_models = False

        while True:
            if had_refinement_error or have_similar_crystal_models:
                break
            max_lattices = self.params.multiple_lattice_search.max_lattices
            if max_lattices is not None and len(experiments.crystals()) >= max_lattices:
                break
            if len(experiments) > 0:
                cutoff_fraction = self.params.multiple_lattice_search.recycle_unindexed_reflections_cutoff
                d_spacings = 1 / self.reflections["rlp"].norms()
                d_min_indexed = flex.min(d_spacings.select(self.indexed_reflections))
                min_reflections_for_indexing = cutoff_fraction * len(
                    self.reflections.select(d_spacings > d_min_indexed)
                )
                crystal_ids = self.reflections.select(d_spacings > d_min_indexed)["id"]
                if (crystal_ids == -1).count(True) < min_reflections_for_indexing:
                    logger.info(
                        "Finish searching for more lattices: %i unindexed reflections remaining.",
                        (crystal_ids == -1).count(True),
                    )
                    break

            n_lattices_previous_cycle = len(experiments.crystals())

            if self.d_min is None:
                self.d_min = self.params.refinement_protocol.d_min_start

            if len(experiments) == 0:
                new_expts = self.find_lattices()
                generate_experiment_identifiers(new_expts)
                experiments.extend(new_expts)
            else:
                try:
                    new = self.find_lattices()
                    generate_experiment_identifiers(new)
                    experiments.extend(new)
                except DialsIndexError:
                    logger.info("Indexing remaining reflections failed")

            if self.params.refinement_protocol.d_min_step is libtbx.Auto:
                n_cycles = self.params.refinement_protocol.n_macro_cycles
                if self.d_min is None or n_cycles == 1:
                    self.params.refinement_protocol.d_min_step = 0
                else:
                    d_spacings = 1 / self.reflections["rlp"].norms()
                    d_min_all = flex.min(d_spacings)
                    self.params.refinement_protocol.d_min_step = (
                        self.d_min - d_min_all
                    ) / (n_cycles - 1)
                    logger.info(
                        "Using d_min_step %.1f",
                        self.params.refinement_protocol.d_min_step,
                    )

            if len(experiments) == 0:
                raise DialsIndexError("No suitable lattice could be found.")
            elif len(experiments.crystals()) == n_lattices_previous_cycle:
                logger.warning("No more suitable lattices could be found")
                # no more lattices found
                break

            for i_cycle in range(self.params.refinement_protocol.n_macro_cycles):
                if (
                    i_cycle > 0
                    and self.d_min is not None
                    and self.params.refinement_protocol.d_min_step > 0
                ):
                    d_min = self.d_min - self.params.refinement_protocol.d_min_step
                    d_min = max(d_min, 0)
                    if self.params.refinement_protocol.d_min_final is not None:
                        d_min = max(d_min, self.params.refinement_protocol.d_min_final)
                    if d_min >= 0:
                        self.d_min = d_min
                        logger.info("Increasing resolution to %.2f Angstrom", d_min)

                # reset reflection lattice flags
                # the lattice a given reflection belongs to: a value of -1 indicates
                # that a reflection doesn't belong to any lattice so far
                self.reflections["id"] = flex.int(len(self.reflections), -1)

                self.index_reflections(experiments, self.reflections)

                if i_cycle == 0 and self.params.known_symmetry.space_group is not None:
                    self._apply_symmetry_post_indexing(
                        experiments, self.reflections, n_lattices_previous_cycle
                    )

                logger.info("\nIndexed crystal models:")
                self.show_experiments(experiments, self.reflections, d_min=self.d_min)

                if self._remove_similar_crystal_models(experiments):
                    have_similar_crystal_models = True
                    break

                logger.info("")
                logger.info("#" * 80)
                logger.info("Starting refinement (macro-cycle %i)", i_cycle + 1)
                logger.info("#" * 80)
                logger.info("")
                self.indexed_reflections = self.reflections["id"] > -1

                sel = flex.bool(len(self.reflections), False)
                # Note, anything below d_min doesn't get indexed so has an id of -1
                # so code below is redundant?
                lengths = 1 / self.reflections["rlp"].norms()
                if self.d_min is not None:
                    isel = (lengths <= self.d_min).iselection()
                    sel.set_selected(isel, True)

                sel.set_selected(self.reflections["id"] == -1, True)
                self.reflections.unset_flags(sel, self.reflections.flags.indexed)
                # N.B. we don't set self.unindexed_reflections here, as we don't want to overwrite yet
                # if the refinement cycle below fails
                unindexed_reflections = self.reflections.select(
                    sel
                )  # reflections that have an id of -1

                reflections_for_refinement = self.reflections.select(
                    self.indexed_reflections
                )
                if self.params.refinement_protocol.mode == "repredict_only":
                    refined_experiments, refined_reflections = (
                        experiments,
                        reflections_for_refinement,
                    )
                    from dials.algorithms.refinement.prediction.managed_predictors import (
                        ExperimentsPredictorFactory,
                    )

                    ref_predictor = ExperimentsPredictorFactory.from_experiments(
                        experiments,
                        spherical_relp=self.all_params.refinement.parameterisation.spherical_relp_model,
                    )
                    ref_predictor(refined_reflections)
                elif self.params.refinement_protocol.mode is None:
                    refined_experiments, refined_reflections = (
                        experiments,
                        reflections_for_refinement,
                    )
                else:
                    try:
                        refined_experiments, refined_reflections = self.refine(
                            experiments, reflections_for_refinement
                        )
                    except (DialsRefineConfigError, DialsRefineRuntimeError) as e:
                        if len(experiments) == 1:
                            raise DialsIndexRefineError(str(e))
                        had_refinement_error = True
                        logger.info("Refinement failed:")
                        logger.info(e)
                        # need to remove crystals - may be shared!
                        models_to_remove = experiments.where(
                            crystal=experiments[-1].crystal
                        )
                        for model_id in sorted(models_to_remove, reverse=True):
                            del experiments[model_id]
                            # remove experiment id from the reflections associated
                            # with this deleted experiment - indexed flag removed
                            # below
                            # note here we are acting on the table from the last macrocycle
                            # This is guaranteed to exist due to the check if len(experiments) == 1: above
                            # N.B. Need to unset the flags here as the break below means we
                            # don't enter the code after
                            self._remove_id_from_reflections(model_id)
                        break

                self._unit_cell_volume_sanity_check(experiments, refined_experiments)

                self.refined_reflections = refined_reflections
                # id can be -1 if some were determined as outliers in self.refine
                sel = self.refined_reflections["id"] < 0
                if sel.count(True):
                    self.refined_reflections.unset_flags(
                        sel,
                        self.refined_reflections.flags.indexed,
                    )
                    self.refined_reflections["miller_index"].set_selected(
                        sel, (0, 0, 0)
                    )
                    unindexed_reflections.extend(self.refined_reflections.select(sel))
                    self.refined_reflections = self.refined_reflections.select(~sel)
                self.refined_reflections.clean_experiment_identifiers_map()
                self.unindexed_reflections = unindexed_reflections

                for i, expt in enumerate(self.experiments):
                    ref_sel = self.refined_reflections.select(
                        self.refined_reflections["imageset_id"] == i
                    )
                    ref_sel = ref_sel.select(ref_sel["id"] >= 0)
                    for i_expt in set(ref_sel["id"]):
                        refined_expt = refined_experiments[i_expt]
                        expt.detector = refined_expt.detector
                        expt.beam = refined_expt.beam
                        expt.goniometer = refined_expt.goniometer
                        expt.scan = refined_expt.scan
                        refined_expt.imageset = expt.imageset

                if not (
                    self.all_params.refinement.parameterisation.beam.fix == "all"
                    and self.all_params.refinement.parameterisation.detector.fix
                    == "all"
                ):
                    # Experimental geometry may have changed - re-map centroids to
                    # reciprocal space
                    self.reflections.map_centroids_to_reciprocal_space(self.experiments)

                # update for next cycle
                experiments = refined_experiments
                self.refined_experiments = refined_experiments

                logger.info("\nRefined crystal models:")
                self.show_experiments(
                    self.refined_experiments,
                    self.refined_reflections,
                    d_min=self.d_min,
                    unindexed_reflections=self.unindexed_reflections,
                )

                if (
                    i_cycle >= 2
                    and self.d_min == self.params.refinement_protocol.d_min_final
                ):
                    logger.info("Target d_min_final reached: finished with refinement")
                    break

        if self.refined_experiments is None:
            raise DialsIndexRefineError("None of the experiments could refine.")

        if len(self.refined_experiments) > 1:
            from dials.algorithms.indexing.compare_orientation_matrices import (
                rotation_matrix_differences,
            )

            logger.info(
                rotation_matrix_differences(self.refined_experiments.crystals())
            )

        if "xyzcal.mm" in self.refined_reflections:
            self._xyzcal_mm_to_px(self.refined_experiments, self.refined_reflections)

    def _remove_id_from_reflections(self, model_id):
        sel = self.refined_reflections["id"] == model_id
        if sel.count(
            True
        ):  # not the case if failure on first cycle of refinement of new xtal
            logger.info(
                "Removing %d reflections with id %d",
                sel.count(True),
                model_id,
            )
            self.refined_reflections["id"].set_selected(sel, -1)
            del self.refined_reflections.experiment_identifiers()[model_id]
            self.refined_reflections.unset_flags(
                sel, self.refined_reflections.flags.indexed
            )
            self.refined_reflections["miller_index"].set_selected(sel, (0, 0, 0))
            self.unindexed_reflections.extend(self.refined_reflections.select(sel))
            self.refined_reflections = self.refined_reflections.select(~sel)
            self.refined_reflections.clean_experiment_identifiers_map()

    def _unit_cell_volume_sanity_check(self, original_experiments, refined_experiments):
        # sanity check for unrealistic unit cell volume increase during refinement
        # usually this indicates too many parameters are being refined given the
        # number of observations provided.
        if not self.params.refinement_protocol.disable_unit_cell_volume_sanity_check:
            for orig_expt, refined_expt in zip(
                original_experiments, refined_experiments
            ):
                uc1 = orig_expt.crystal.get_unit_cell()
                uc2 = refined_expt.crystal.get_unit_cell()
                volume_change = abs(uc1.volume() - uc2.volume()) / uc1.volume()
                cutoff = 0.5
                if volume_change > cutoff:
                    msg = "\n".join(
                        (
                            "Unrealistic unit cell volume increase during refinement of %.1f%%.",
                            "Please try refining fewer parameters, either by enforcing symmetry",
                            "constraints (space_group=) and/or disabling experimental geometry",
                            "refinement (detector.fix=all and beam.fix=all). To disable this",
                            "sanity check set disable_unit_cell_volume_sanity_check=True.",
                        )
                    ) % (100 * volume_change)
                    raise DialsIndexError(msg)

    def _apply_symmetry_post_indexing(
        self, experiments, reflections, n_lattices_previous_cycle
    ):
        # now apply the space group symmetry only after the first indexing
        # need to make sure that the symmetrized orientation is similar to the P1 model
        for cryst in experiments.crystals()[n_lattices_previous_cycle:]:
            new_cryst, cb_op = self._symmetry_handler.apply_symmetry(cryst)
            new_cryst = new_cryst.change_basis(cb_op)
            cryst.update(new_cryst)
            cryst.set_space_group(self.params.known_symmetry.space_group.group())
            for i_expt, expt in enumerate(experiments):
                if expt.crystal is not cryst:
                    continue
                if not cb_op.is_identity_op():
                    miller_indices = reflections["miller_index"].select(
                        reflections["id"] == i_expt
                    )
                    miller_indices = cb_op.apply(miller_indices)
                    reflections["miller_index"].set_selected(
                        reflections["id"] == i_expt, miller_indices
                    )

    def _remove_similar_crystal_models(self, experiments):
        """
        Checks for too-similar crystal models and removes them.

        Checks whether the most recently added crystal model is similar to previously
        found crystal models, and if so, deletes the last crystal model from the
        experiment list.
        """
        have_similar_crystal_models = False
        cryst_b = experiments.crystals()[-1]
        for i_a, cryst_a in enumerate(experiments.crystals()[:-1]):
            R_ab, axis, angle, cb_op_ab = difference_rotation_matrix_axis_angle(
                cryst_a, cryst_b
            )
            min_angle = self.params.multiple_lattice_search.minimum_angular_separation
            if abs(angle) < min_angle:  # degrees
                # account for potential shared crystal model
                models_to_reject = experiments.where(crystal=cryst_b)
                reject = len(experiments.crystals())
                logger.info(f"Crystal models too similar, rejecting crystal {reject}")
                logger.info(
                    f"Rotation matrix to transform crystal {i_a + 1} to crystal {reject}"
                )
                logger.info(R_ab)
                logger.info(
                    f"Rotation of {angle:.3f} degrees"
                    + " about axis ({:.3f}, {:.3f}, {:.3f})".format(*axis)
                )
                have_similar_crystal_models = True
                for id_ in sorted(models_to_reject, reverse=True):
                    del experiments[id_]
                    # Unset ids in the reflection table
                    self._remove_id_from_reflections(id_)
                break
        return have_similar_crystal_models

    def _xyzcal_mm_to_px(self, experiments, reflections):
        # set xyzcal.px field in reflections
        reflections["xyzcal.px"] = flex.vec3_double(len(reflections))
        for i, expt in enumerate(experiments):
            imgset_sel = reflections["imageset_id"] == i
            refined_reflections = reflections.select(imgset_sel)
            panel_numbers = flex.size_t(refined_reflections["panel"])
            xyzcal_mm = refined_reflections["xyzcal.mm"]
            x_mm, y_mm, z = xyzcal_mm.parts()
            xy_cal_mm = flex.vec2_double(x_mm, y_mm)
            xy_cal_px = flex.vec2_double(len(xy_cal_mm))
            for i_panel in range(len(expt.detector)):
                panel = expt.detector[i_panel]
                sel = panel_numbers == i_panel
                xy_cal_px.set_selected(
                    sel, panel.millimeter_to_pixel(xy_cal_mm.select(sel))
                )
            x_px, y_px = xy_cal_px.parts()
            if expt.scan is not None:
                if expt.scan.has_property("time_of_flight"):
                    tof = expt.scan.get_property("time_of_flight")
                    frames = list(range(len(tof)))
                    tof_to_frame = tof_helpers.tof_to_frame_interpolator(tof, frames)
                    z.set_selected(z < min(tof), min(tof))
                    z.set_selected(z > max(tof), max(tof))
                    z_px = flex.double(tof_to_frame(z))
                else:
                    z_px = expt.scan.get_array_index_from_angle(z, deg=False)
            else:
                # must be a still image, z centroid not meaningful
                z_px = z
            xyzcal_px = flex.vec3_double(x_px, y_px, z_px)
            reflections["xyzcal.px"].set_selected(imgset_sel, xyzcal_px)

    def show_experiments(
        self, experiments, reflections, d_min=None, unindexed_reflections=None
    ):
        if d_min is not None:
            reciprocal_lattice_points = reflections["rlp"]
            d_spacings = 1 / reciprocal_lattice_points.norms()
            reflections = reflections.select(d_spacings > d_min)
            if unindexed_reflections:
                reciprocal_lattice_points = unindexed_reflections["rlp"]
                d_spacings = 1 / reciprocal_lattice_points.norms()
                unindexed_reflections = unindexed_reflections.select(d_spacings > d_min)

        for i_expt, expt in enumerate(experiments):
            logger.info(
                "model %i (%i reflections):",
                i_expt + 1,
                (reflections["id"] == i_expt).count(True),
            )
            logger.info(expt.crystal)

        indexed_flags = reflections.get_flags(reflections.flags.indexed)
        imageset_id = reflections["imageset_id"]
        rows = [["Imageset", "# indexed", "# unindexed", "% indexed"]]
        for i in range(flex.max(imageset_id) + 1):
            imageset_indexed_flags = indexed_flags.select(imageset_id == i)
            indexed_count = imageset_indexed_flags.count(True)
            unindexed_count = imageset_indexed_flags.count(False)
            if unindexed_reflections:
                sel = unindexed_reflections["imageset_id"] == i
                unindexed_count += sel.count(True)
            rows.append(
                [
                    str(i),
                    str(indexed_count),
                    str(unindexed_count),
                    f"{indexed_count / (indexed_count + unindexed_count)*100:.1f}",
                ]
            )
        logger.info(dials.util.tabulate(rows, headers="firstrow"))

    def find_max_cell(self):
        params = self.params.max_cell_estimation
        if self.params.max_cell is libtbx.Auto:
            if self.params.known_symmetry.unit_cell is not None:
                uc_params = self._symmetry_handler.target_symmetry_primitive.unit_cell().parameters()
                self.params.max_cell = params.multiplier * max(uc_params[:3])
                logger.info("Using max_cell: %.1f Angstrom", self.params.max_cell)
            else:
                convert_reflections_z_to_deg = True
                all_tof_experiments = False
                for expt in self.experiments:
                    if expt.scan is not None and expt.scan.has_property(
                        "time_of_flight"
                    ):
                        all_tof_experiments = True
                    elif all_tof_experiments:
                        raise ValueError(
                            "Cannot find max cell for ToF and non-ToF experiments at the same time"
                        )

                if all_tof_experiments:
                    if params.step_size < 100:
                        logger.info("Setting default ToF step size to 500 usec")
                        params.step_size = 500
                        convert_reflections_z_to_deg = False

                self.params.max_cell = find_max_cell(
                    self.reflections,
                    max_cell_multiplier=params.multiplier,
                    step_size=params.step_size,
                    nearest_neighbor_percentile=params.nearest_neighbor_percentile,
                    histogram_binning=params.histogram_binning,
                    nn_per_bin=params.nn_per_bin,
                    max_height_fraction=params.max_height_fraction,
                    filter_ice=params.filter_ice,
                    filter_overlaps=params.filter_overlaps,
                    overlaps_border=params.overlaps_border,
                    convert_reflections_z_to_deg=convert_reflections_z_to_deg,
                ).max_cell
                logger.info("Found max_cell: %.1f Angstrom", self.params.max_cell)

    def index_reflections(self, experiments, reflections):
        self._assign_indices(reflections, experiments, d_min=self.d_min)
        if self.hkl_offset is not None and self.hkl_offset != (0, 0, 0):
            reflections["miller_index"] = apply_hkl_offset(
                reflections["miller_index"], self.hkl_offset
            )
            self.hkl_offset = None

    def refine(self, experiments, reflections):
        from dials.algorithms.indexing.refinement import refine

        properties_to_save = [
            "xyzcal.mm",
            "entering",
            "wavelength_cal",
            "s0_cal",
            "tof_cal",
        ]

        refiner, refined, outliers = refine(self.all_params, reflections, experiments)
        if outliers is not None:
            reflections["id"].set_selected(outliers, -1)

        predicted = refiner.predict_for_indexed()

        for i in properties_to_save:
            if i in predicted:
                reflections[i] = predicted[i]

        reflections.unset_flags(
            flex.bool(len(reflections), True), reflections.flags.centroid_outlier
        )
        assert (
            reflections.get_flags(reflections.flags.centroid_outlier).count(True) == 0
        )
        reflections.set_flags(
            predicted.get_flags(predicted.flags.centroid_outlier),
            reflections.flags.centroid_outlier,
        )
        reflections.set_flags(
            refiner.selection_used_for_refinement(),
            reflections.flags.used_in_refinement,
        )
        return refiner.get_experiments(), reflections

    def _debug_write_reciprocal_lattice_points_as_pdb(
        self, file_name="reciprocal_lattice.pdb"
    ):
        from cctbx import crystal, xray

        cs = crystal.symmetry(
            unit_cell=(1000, 1000, 1000, 90, 90, 90), space_group="P1"
        )
        for i_panel in range(len(self.experiments[0].detector)):
            if len(self.experiments[0].detector) > 1:
                file_name = "reciprocal_lattice_%i.pdb" % i_panel
            with open(file_name, "wb") as f:
                xs = xray.structure(crystal_symmetry=cs)
                reflections = self.reflections.select(
                    self.reflections["panel"] == i_panel
                )
                for site in reflections["rlp"]:
                    xs.add_scatterer(xray.scatterer("C", site=site))
                xs.sites_mod_short()
                f.write(xs.as_pdb_file())

    def export_as_json(
        self, experiments, file_name="indexed_experiments.json", compact=False
    ):
        assert experiments.is_consistent()
        experiments.as_file(file_name)

    def export_reflections(self, reflections, file_name="reflections.pickle"):
        reflections.as_file(file_name)

    def find_lattices(self):
        raise NotImplementedError()


def apply_hkl_offset(indices, offset):
    h, k, l = indices.as_vec3_double().parts()
    h += offset[0]
    k += offset[1]
    l += offset[2]
    return flex.miller_index(h.iround(), k.iround(), l.iround())
