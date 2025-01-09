from __future__ import annotations

import collections
import copy
import glob
import logging
import os
import pickle
import sys
import tarfile
import time
from io import BytesIO

from dxtbx.model.experiment_list import (
    Experiment,
    ExperimentList,
    ExperimentListFactory,
)
from libtbx.phil import parse
from libtbx.utils import Abort, Sorry

import dials.util
from dials.array_family import flex
from dials.util import log

logger = logging.getLogger("dials.command_line.stills_process")

help_message = """
DIALS script for processing still images. Import, index, refine, and integrate are all done for each image
separately.
"""


def _control_phil_str():
    return """
  input {
    file_list = None
      .type = path
      .help = Path to a text file with a list of images
    glob = None
      .type = str
      .help = For large, multi-file datasets, specify the paths using wildcards (e.g. *.cbf)
      .multiple = True
    image_tag = None
      .type = str
      .multiple = True
      .help = Only process images with these tag(s). For single-image files (like CBFs or SMVs), the image \
              tag for each file is the file name. For multi-image files like HDF5, the image tag is        \
              filename_imagenumber (including leading zeros). Use show_image_tags=True to see the list     \
              of image tags that will be used for a dataset.
    show_image_tags = False
      .type = bool
      .help = Show the set of image tags that would be used during processing. To process subsets of image \
              files, use these tags with the image_tag parameter.
    max_images = None
      .type = int
      .help = Limit total number of processed images to max_images
    ignore_gain_mismatch = False
      .type = bool
      .expert_level = 3
      .help = Detector gain should be set on the detector models loaded from the images or in the \
              processing parameters, not both. Override the check that this is true with this flag. \
  }

  dispatch {
    pre_import = False
      .type = bool
      .expert_level = 2
      .help = If True, before processing import all the data. Needed only if processing \
              multiple multi-image files at once (not a recommended use case)
    process_percent = None
      .type = int(value_min=1, value_max=100)
      .help = Percent of events to process
    refine = False
      .expert_level = 2
      .type = bool
      .help = If True, after indexing, refine the experimental models
    squash_errors = True
      .expert_level = 2
      .type = bool
      .help = If True, if an image fails to process, continue to the next image. \
              otherwise, halt processing and show the error.
    find_spots = True
      .expert_level = 2
      .type = bool
      .help = Whether to do spotfinding. Needed for indexing/integration
    index = True
      .expert_level = 2
      .type = bool
      .help = Attempt to index images. find_spots also needs to be True for this to work
    integrate = True
      .expert_level = 2
      .type = bool
      .help = Integrate indexed images. Ignored if index=False or find_spots=False
    coset = False
      .expert_level = 2
      .type = bool
      .help = Within the integrate dispatcher, integrate a sublattice coset intended to represent \
              negative control spots with no Bragg diffraction.
    hit_finder{
      enable = True
        .type = bool
        .help = Whether to do hitfinding. hit_finder=False: process all images
      minimum_number_of_reflections = 16
        .type = int
        .help = If the number of strong reflections on an image is less than this, and \
                 the hitfinder is enabled, discard this image.
      maximum_number_of_reflections = None
       .type = int
       .help = If specified, ignores images with more than this many number of reflections
    }
  }

  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    composite_output = True
      .type = bool
      .help = If True, save one set of experiment/reflection files per process, where each is a \
              concatenated list of all the successful events examined by that process. \
              If False, output a separate experiment/reflection file per image (generates a \
              lot of files).
    logging_dir = None
      .type = str
      .help = Directory output log files will be placed
    logging_option = normal *suppressed disabled
      .type = choice
      .help = normal includes all logging, suppress turns off DIALS refine output
      .help = and disabled removes basically all logging
    experiments_filename = None
      .type = str
      .help = The filename for output experiments. For example, %s_imported.expt
    strong_filename = None
      .type = str
      .help = The filename for strong reflections from spot finder output. For example: \
              %s_strong.refl
    indexed_filename = %s_indexed.refl
      .type = str
      .help = The filename for indexed reflections.
    refined_experiments_filename = %s_refined.expt
      .type = str
      .help = The filename for saving refined experimental models
    integrated_filename = %s_integrated.refl
      .type = str
      .help = The filename for final integrated reflections.
    integrated_experiments_filename = %s_integrated.expt
      .type = str
      .help = The filename for saving final experimental models.
    coset_filename = %s_coset%d.refl
      .type = str
      .help = The filename for final coset reflections.
    coset_experiments_filename = %s_coset%d.expt
      .type = str
      .help = The filename for saving final coset experimental models.
    profile_filename = None
      .type = str
      .help = The filename for output reflection profile parameters
    integration_pickle = None
      .type = str
      .expert_level = 3
      .help = Filename for legacy cxi.merge integration pickle files. Example: int-%d-%s.pickle
  }

  mp {
    method = *multiprocessing sge lsf pbs mpi
      .type = choice
      .help = "The multiprocessing method to use"
    nproc = 1
      .type = int(value_min=1)
      .help = "The number of processes to use."
    composite_stride = None
      .type = int
      .help = For MPI, if using composite mode, specify how many ranks to    \
              aggregate data from.  For example, if you have 100 processes,  \
              composite mode will output N*100 files, where N is the number  \
              of file types (expt, refl, etc). If you specify stride = 25, \
              then each group of 25 process will send their results to 4     \
              processes and only N*4 files will be created. Ideally, match   \
              stride to the number of processors per node.
    debug
      .expert_level = 2
    {
      cProfile = False
        .type = bool
        .help = Enable code profiling. Profiling file will be available in  \
                the debug folder. Use (for example) runsnake to visualize   \
                processing performance
      output_debug_logs = True
        .type = bool
        .help = Whether to write debugging information for every image      \
                processed
    }
  }
"""


def _dials_phil_str():
    return """
  input {
    reference_geometry = None
      .type = str
      .help = Provide an models.expt file with exactly one detector model. Data processing will use \
              that geometry instead of the geometry found in the image headers.
    sync_reference_geom = True
      .type = bool
      .help = ensures the reference hierarchy agrees with the image format
  }

  output {
    shoeboxes = True
      .type = bool
      .help = Save the raw pixel values inside the reflection shoeboxes during spotfinding.
  }

  include scope dials.util.options.geometry_phil_scope
  include scope dials.algorithms.spot_finding.factory.phil_scope
  include scope dials.algorithms.indexing.indexer.phil_scope
  indexing {
      include scope dials.algorithms.indexing.lattice_search.basis_vector_search_phil_scope
  }
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope

  indexing {
    stills {
      method_list = None
        .type = strings
        .help = List of indexing methods. If indexing fails with first method, indexing will be \
                attempted with the next, and so forth
      reflection_subsampling
      {
        enable  = False
          .type = bool
          .help = Enable random subsampling of reflections during indexing. Attempts to index    \
                  will be repeated with random subsampling of spotfinder spots, starting at      \
                  step_start % of total spots, repeating n_attempts_per_step times per step,     \
                  and decreasing by step_size % until step_stop % is reached. With the defaults, \
                  26 total indexing attempts will be made, randomly subsampling from 100% to     \
                  50% by 2 % steps, one attempt per step.
        step_start = 100
          .type = int
          .help = What percent of reflections to start with
        step_stop = 50
          .type = int
          .help = What percent of reflections to end with
        step_size = 2
          .type = int
          .help = What percent of reflections to decrease by per step
        n_attempts_per_step = 1
          .type = int
          .help = How many attempts to make at each step
      }
      known_orientations = None
        .type = path
        .multiple = True
        .expert_level = 2
        .help = Paths to previous processing results including crystal orientations. \
                If specified, images will not be re-indexed, but instead the known \
                orientations will be used. Provide paths to experiment list files, using \
                wildcards as needed.
      require_known_orientation = False
        .type = bool
        .expert_level = 2
        .help = If known_orientations are provided, and an orientation for an image is not \
                found, this is whether or not to attempt to index the image from scratch \
                using indexing.method
    }
  }

  integration {
    include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope
    coset {
      transformation = 6
        .type = int(value_min=0, value_max=6)
        .multiple = False
        .help = The index number(s) of the modulus=2 sublattice transformation(s) used to produce distance coset results. \
                0=Double a, 1=Double b, 2=Double c, 3=C-face centering, 4=B-face centering, 5=A-face centering, 6=Body centering \
                See Sauter and Zwart, Acta D (2009) 65:553
    }

    integration_only_overrides {
      trusted_range = None
        .type = floats(size=2)
        .help = "Override the panel trusted range [min_trusted_value, max_trusted_value] during integration."
        .short_caption = "Panel trusted range"
    }
  }

  profile {
    gaussian_rs {
      parameters {
        sigma_b_cutoff = 0.1
          .type = float
          .help = Maximum sigma_b before the image is rejected
      }
    }
  }
"""


def _program_defaults_phil_str():
    return """
indexing {
  method = fft1d
}
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 1
      action = fix
    }
    beam.fix = all
    detector.fix = all
  }
  reflections {
    weighting_strategy.override = stills
    outlier.algorithm = null
  }
}
integration {
  integrator = stills
  profile.fitting = False
  background {
    algorithm = simple
    simple {
      outlier.algorithm = plane
      model.algorithm = linear2d
    }
  }
}
profile.gaussian_rs.min_spots.overall = 0
"""


control_phil_str = _control_phil_str()
dials_phil_str = _dials_phil_str()
program_defaults_phil_str = _program_defaults_phil_str()

phil_scope = parse(control_phil_str + dials_phil_str, process_includes=True).fetch(
    parse(program_defaults_phil_str)
)


def do_import(filename, load_models=True):
    logger.info("Loading %s", os.path.basename(filename))
    experiments = ExperimentListFactory.from_filenames([filename], load_models=False)
    if len(experiments) == 0:
        try:
            experiments = ExperimentListFactory.from_json_file(filename)
        except ValueError:
            raise Abort(f"Could not load {filename}")

    if len(experiments) == 0:
        raise Abort(f"Could not load {filename}")

    from dxtbx.imageset import ImageSetFactory

    all_experiments = ExperimentList()
    for experiment in experiments:
        # Convert from ImageSequence to ImageSet, if needed
        imageset = ImageSetFactory.imageset_from_anyset(experiment.imageset)
        for i in range(len(imageset)):
            # Preserve original models if they were available (in the case of an image file
            # they will not be, but in the case of a previously processed experiment list,
            # then they may be available
            expt = Experiment(
                imageset=imageset[i : i + 1],
                detector=experiment.detector,
                beam=experiment.beam,
                scan=experiment.scan,
                goniometer=experiment.goniometer,
                crystal=experiment.crystal,
            )
            if load_models:
                expt.load_models()
            all_experiments.append(expt)

    return all_experiments


def sync_geometry(src, dest):
    dest.set_local_frame(
        src.get_local_fast_axis(), src.get_local_slow_axis(), src.get_local_origin()
    )
    if not src.is_panel():
        for src_child, dest_child in zip(src, dest):
            sync_geometry(src_child, dest_child)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import ArgumentParser

        # The script usage
        usage = "usage: dials.stills_process [options] [param.phil] filenames"

        self.tag = None
        self.reference_detector = None

        # Create the parser
        self.parser = ArgumentParser(usage=usage, phil=phil_scope, epilog=help_message)

    def load_reference_geometry(self):
        if self.params.input.reference_geometry is None:
            return

        try:
            ref_experiments = ExperimentListFactory.from_json_file(
                self.params.input.reference_geometry, check_format=False
            )
        except Exception:
            try:
                import dxtbx

                img = dxtbx.load(self.params.input.reference_geometry)
            except Exception:
                raise Sorry(
                    "Couldn't load geometry file %s"
                    % self.params.input.reference_geometry
                )
            else:
                self.reference_detector = img.get_detector()
        else:
            assert len(ref_experiments.detectors()) == 1
            self.reference_detector = ref_experiments.detectors()[0]

    def run(self, args=None):
        """Execute the script."""
        from libtbx import easy_mp

        try:
            from mpi4py import MPI
        except ImportError:
            rank = 0
            size = 1
        else:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
            size = comm.Get_size()  # size: number of processes running in this job

        if rank == 0:
            # Parse the command line
            params, options, all_paths = self.parser.parse_args(
                args, show_diff_phil=False, return_unhandled=True, quick_parse=True
            )

            if params.input.glob:
                all_paths.extend(params.input.glob)
            globbed = []
            for p in all_paths:
                g = glob.glob(p)
                if not g:
                    sys.exit(f"Error: Unhandled path or option: {p}")
                globbed.extend(g)
            all_paths = globbed

            if not all_paths and params.input.file_list is not None:
                all_paths.extend(
                    [path.strip() for path in open(params.input.file_list).readlines()]
                )

            if params.indexing.stills.known_orientations:
                known_orientations = {}
                for path in params.indexing.stills.known_orientations:
                    for g in glob.glob(path):
                        ko_expts = ExperimentList.from_file(g, check_format=False)
                        for expt in ko_expts:
                            assert (
                                len(expt.imageset.indices()) == 1
                                and len(expt.imageset.paths()) == 1
                            )
                            key = (
                                os.path.basename(expt.imageset.paths()[0]),
                                expt.imageset.indices()[0],
                            )
                            if key not in known_orientations:
                                known_orientations[key] = []
                            known_orientations[key].append(expt.crystal)
                if not known_orientations:
                    raise Sorry(
                        "No known_orientations found at the locations specified: %s"
                        % ", ".join(params.indexing.stills.known_orientations)
                    )
                params.indexing.stills.known_orientations = known_orientations
        if size > 1:
            if rank == 0:
                transmitted_info = params, options, all_paths
            else:
                transmitted_info = None
            params, options, all_paths = comm.bcast(transmitted_info, root=0)

        if params.output.logging_option == "suppressed":
            logging.getLogger("dials.algorithms.indexing.nave_parameters").setLevel(
                logging.ERROR
            )
            logging.getLogger("dials.algorithms.indexing.stills_indexer").setLevel(
                logging.ERROR
            )
            logging.getLogger("dials.algorithms.refinement.refiner").setLevel(
                logging.ERROR
            )
            logging.getLogger(
                "dials.algorithms.refinement.reflection_manager"
            ).setLevel(logging.ERROR)

        elif params.output.logging_option == "disabled":
            logging.disable(logging.ERROR)

        # Check we have some filenames
        if not all_paths:
            self.parser.print_help()
            return

        if params.mp.debug.cProfile:
            import cProfile

            self.pr = cProfile.Profile()
            self.pr.enable()

        logger.info(f"Have {len(all_paths)} files")

        # Mask validation
        for mask_path in params.spotfinder.lookup.mask, params.integration.lookup.mask:
            if mask_path is not None and not os.path.isfile(mask_path):
                raise Sorry(f"Mask {mask_path} not found")

        # Save the options
        self.options = options
        self.params = params

        st = time.time()

        if params.mp.method == "mpi":
            # Configure the logging
            if params.output.logging_dir is None:
                logfile = None
            else:
                log_path = os.path.join(
                    params.output.logging_dir, "log_rank%04d.out" % rank
                )
                error_path = os.path.join(
                    params.output.logging_dir, "error_rank%04d.out" % rank
                )
                print(f"Redirecting stdout to {log_path}")
                print(f"Redirecting stderr to {error_path}")
                sys.stdout = open(log_path, "a")
                sys.stderr = open(error_path, "a")
                print("Should be redirected now")

                logfile = os.path.join(
                    params.output.logging_dir, "info_rank%04d.out" % rank
                )

            log.config(verbosity=options.verbose, logfile=logfile)

        else:
            # Configure logging
            log.config(verbosity=options.verbose, logfile="dials.process.log")

        bad_phils = [f for f in all_paths if os.path.splitext(f)[1] == ".phil"]
        if len(bad_phils) > 0:
            self.parser.print_help()
            logger.error(
                "Error: the following phil files were not understood: %s",
                ", ".join(bad_phils),
            )
            return

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if not (
                    self.params.integration.debug.output
                    and not self.params.integration.debug.separate_files
                ):
                    raise Sorry(
                        "Shoeboxes must be saved to integration intermediates to apply an absorption correction. "
                        + "Set integration.debug.output=True, integration.debug.separate_files=False and "
                        + "integration.debug.delete_shoeboxes=True to temporarily store shoeboxes."
                    )

        self.load_reference_geometry()
        from dials.command_line.dials_import import ManualGeometryUpdater

        update_geometry = ManualGeometryUpdater(params)

        # Import stuff
        logger.info("Loading files...")
        pre_import = params.dispatch.pre_import or len(all_paths) == 1
        if pre_import:
            # Handle still imagesets by breaking them apart into multiple experiments
            # Further handle single file still imagesets (like HDF5) by tagging each
            # frame using its index

            experiments = ExperimentList()
            for path in sorted(all_paths):
                experiments.extend(do_import(path, load_models=False))

            indices = []
            basenames = []
            basename_counts = {}
            split_experiments = []
            for i, imageset in enumerate(experiments.imagesets()):
                assert len(imageset) == 1
                paths = imageset.paths()
                indices.append(i)
                basename = os.path.splitext(os.path.basename(paths[0]))[0]
                basenames.append(basename)
                if basename in basename_counts:
                    basename_counts[basename] += 1
                else:
                    basename_counts[basename] = 1
                split_experiments.append(experiments[i : i + 1])
            tags = []
            split_experiments2 = []
            for i, basename in zip(indices, basenames):
                if basename_counts[basename] > 1:
                    tag = "%s_%05d" % (basename, i)
                else:
                    tag = basename
                if (
                    not self.params.input.image_tag
                    or tag in self.params.input.image_tag
                ):
                    tags.append(tag)
                    split_experiments2.append(split_experiments[i])
            split_experiments = split_experiments2

            # Wrapper function
            def do_work(i, item_list, processor=None, finalize=True):
                if not processor:
                    processor = Processor(
                        copy.deepcopy(params), composite_tag="%04d" % i, rank=i
                    )

                for item in item_list:
                    tag = item[0]
                    experiments = split_experiments[item[1]]
                    try:
                        assert len(experiments) == 1
                        experiment = experiments[0]
                        experiment.load_models()
                        imageset = experiment.imageset
                        update_geometry(imageset)
                        experiment.beam = imageset.get_beam()
                        experiment.detector = imageset.get_detector()
                    except RuntimeError as e:
                        logger.warning("Error updating geometry on item %s, %s", tag, e)
                        continue

                    if self.reference_detector is not None:
                        experiment = experiments[0]
                        if self.params.input.sync_reference_geom:
                            imageset = experiment.imageset
                            sync_geometry(
                                self.reference_detector.hierarchy(),
                                imageset.get_detector().hierarchy(),
                            )
                            experiment.detector = imageset.get_detector()
                        else:
                            experiment.detector = copy.deepcopy(self.reference_detector)

                    processor.process_experiments(tag, experiments)
                    imageset.clear_cache()
                if finalize:
                    processor.finalize()
                return processor

            iterable = list(zip(tags, range(len(split_experiments))))

        else:
            basenames = collections.defaultdict(int)
            sorted_paths = sorted(all_paths)
            for filename in sorted_paths:
                basename = os.path.splitext(os.path.basename(filename))[0]
                basenames[basename] += 1
            tags = []
            all_paths2 = []
            for i, (basename, count) in enumerate(basenames.items()):
                if count > 1:
                    tag = "%s_%05d" % (basename, i)
                else:
                    tag = basename
                if (
                    not self.params.input.image_tag
                    or tag in self.params.input.image_tag
                ):
                    tags.append(tag)
                    all_paths2.append(sorted_paths[i])
            all_paths = all_paths2

            # Wrapper function
            def do_work(i, item_list, processor=None, finalize=True):
                if not processor:
                    processor = Processor(
                        copy.deepcopy(params), composite_tag="%04d" % i, rank=i
                    )
                for item in item_list:
                    tag, filename = item

                    experiments = do_import(filename, load_models=True)
                    imagesets = experiments.imagesets()
                    if len(imagesets) == 0 or len(imagesets[0]) == 0:
                        logger.info("Zero length imageset in file: %s", filename)
                        return
                    if len(imagesets) > 1:
                        raise Abort(f"Found more than one imageset in file: {filename}")
                    if len(imagesets[0]) > 1:
                        raise Abort(
                            "Found a multi-image file. Run again with pre_import=True"
                        )

                    try:
                        update_geometry(imagesets[0])
                        experiment = experiments[0]
                        experiment.beam = imagesets[0].get_beam()
                        experiment.detector = imagesets[0].get_detector()
                    except RuntimeError as e:
                        logger.warning("Error updating geometry on item %s, %s", tag, e)
                        continue

                    if self.reference_detector is not None:
                        if self.params.input.sync_reference_geom:
                            imageset = experiments[0].imageset
                            sync_geometry(
                                self.reference_detector.hierarchy(),
                                imageset.get_detector().hierarchy(),
                            )
                            experiments[0].detector = imageset.get_detector()
                        else:
                            experiments[0].detector = copy.deepcopy(
                                self.reference_detector
                            )

                    processor.process_experiments(tag, experiments)
                if finalize:
                    processor.finalize()
                return processor

            iterable = list(zip(tags, all_paths))

        if params.input.max_images:
            iterable = iterable[: params.input.max_images]

        if params.input.show_image_tags:
            print("Showing image tags for this dataset and exiting")
            for tag, item in iterable:
                print(tag)
            return

        # prepare fractions of process_percent, if given
        process_fractions = None
        if params.dispatch.process_percent:
            import fractions

            percent = params.dispatch.process_percent / 100
            process_fractions = fractions.Fraction(percent).limit_denominator(100)

            def process_this_event(nevent):
                # nevent modulo the denominator gives us which fraction we're in
                n_mod_denom = nevent % process_fractions.denominator
                # compare the 0-indexed modulo against the 1-indexed numerator (intentionally not <=)
                n_accept = n_mod_denom < process_fractions.numerator
                return n_accept

        # Process the data
        if params.mp.method == "mpi":
            if size <= 2:  # client/server only makes sense for n>2
                subset = [
                    item for i, item in enumerate(iterable) if (i + rank) % size == 0
                ]
                do_work(rank, subset)
            else:
                processor = Processor(
                    copy.deepcopy(params), composite_tag="%04d" % rank, rank=rank
                )

                if rank == 0:
                    # server process
                    num_iter = len(iterable)
                    for item_num, item in enumerate(iterable):
                        print(
                            "Processing %d / %d shots" % (item_num, num_iter),
                            flush=True,
                        )
                        if process_fractions and not process_this_event(item_num):
                            continue

                        print("Getting next available process")
                        rankreq = comm.recv(source=MPI.ANY_SOURCE)
                        print(f"Process {rankreq} is ready, sending {item[0]}\n")
                        comm.send(item, dest=rankreq)
                    # send a stop command to each process
                    print("MPI DONE, sending stops\n")
                    for rankreq in range(size - 1):
                        rankreq = comm.recv(source=MPI.ANY_SOURCE)
                        print("Sending stop to %d\n" % rankreq)
                        comm.send("endrun", dest=rankreq)
                    print("All stops sent.")

                else:
                    # client process
                    while True:
                        # inform the server this process is ready for an event
                        print("Rank %d getting next task" % rank)
                        comm.send(rank, dest=0)
                        print("Rank %d waiting for response" % rank)
                        item = comm.recv(source=0)
                        if item == "endrun":
                            print("Rank %d received endrun" % rank)
                            break
                        print("Rank %d beginning processing" % rank)
                        try:
                            processor = do_work(rank, [item], processor, finalize=False)
                        except Exception as e:
                            print(
                                "Rank %d unhandled exception processing event" % rank,
                                str(e),
                            )
                        print("Rank %d event processed" % rank)
                processor.finalize()
        else:
            from dxtbx.command_line.image_average import splitit

            if params.mp.nproc == 1:
                do_work(0, iterable)
            else:
                result = list(
                    easy_mp.multi_core_run(
                        myfunction=do_work,
                        argstuples=list(enumerate(splitit(iterable, params.mp.nproc))),
                        nproc=params.mp.nproc,
                    )
                )
                error_list = [r[2] for r in result]
                if error_list.count(None) != len(error_list):
                    print(
                        "Some processes failed execution. Not all images may have processed. Error messages:"
                    )
                    for error in error_list:
                        if error is None:
                            continue
                        print(error)

        # Total Time
        logger.info("")
        logger.info("Total Time Taken = %f seconds", time.time() - st)

        if params.mp.debug.cProfile:
            self.pr.disable()
            self.pr.dump_stats(
                os.path.join(
                    self.params.output.output_dir, "debug", "cpu_%d.prof" % comm.rank
                )
            )


class Processor:
    def __init__(self, params, composite_tag=None, rank=0):
        self.params = params
        self.composite_tag = composite_tag

        # The convention is to put %s in the phil parameter to add a tag to
        # each output datafile. Save the initial templates here.
        self.experiments_filename_template = params.output.experiments_filename
        self.strong_filename_template = params.output.strong_filename
        self.indexed_filename_template = params.output.indexed_filename
        self.refined_experiments_filename_template = (
            params.output.refined_experiments_filename
        )
        self.integrated_filename_template = params.output.integrated_filename
        self.integrated_experiments_filename_template = (
            params.output.integrated_experiments_filename
        )
        if params.dispatch.coset:
            self.coset_filename_template = params.output.coset_filename
            self.coset_experiments_filename_template = (
                params.output.coset_experiments_filename
            )

        debug_dir = os.path.join(params.output.output_dir, "debug")
        if not os.path.exists(debug_dir):
            try:
                os.makedirs(debug_dir)
            except OSError:
                pass  # due to multiprocessing, makedirs can sometimes fail
        assert os.path.exists(debug_dir)
        self.debug_file_path = os.path.join(debug_dir, "debug_%d.txt" % rank)
        write_newline = os.path.exists(self.debug_file_path)
        if write_newline:  # needed if the there was a crash
            self.debug_write("")

        if params.output.composite_output:
            assert composite_tag is not None

            self.all_imported_experiments = ExperimentList()
            self.all_strong_reflections = flex.reflection_table()
            self.all_indexed_experiments = ExperimentList()
            self.all_indexed_reflections = flex.reflection_table()
            self.all_integrated_experiments = ExperimentList()
            self.all_integrated_reflections = flex.reflection_table()
            self.all_int_pickle_filenames = []
            self.all_int_pickles = []
            self.all_coset_experiments = ExperimentList()
            self.all_coset_reflections = flex.reflection_table()

            self.setup_filenames(composite_tag)

    def setup_filenames(self, tag):
        # before processing, set output paths according to the templates
        if (
            self.experiments_filename_template is not None
            and "%s" in self.experiments_filename_template
        ):
            self.params.output.experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.strong_filename_template is not None
            and "%s" in self.strong_filename_template
        ):
            self.params.output.strong_filename = os.path.join(
                self.params.output.output_dir,
                self.strong_filename_template % ("idx-" + tag),
            )
        if (
            self.indexed_filename_template is not None
            and "%s" in self.indexed_filename_template
        ):
            self.params.output.indexed_filename = os.path.join(
                self.params.output.output_dir,
                self.indexed_filename_template % ("idx-" + tag),
            )
        if (
            self.refined_experiments_filename_template is not None
            and "%s" in self.refined_experiments_filename_template
        ):
            self.params.output.refined_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.refined_experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.integrated_filename_template is not None
            and "%s" in self.integrated_filename_template
        ):
            self.params.output.integrated_filename = os.path.join(
                self.params.output.output_dir,
                self.integrated_filename_template % ("idx-" + tag),
            )
        if (
            self.integrated_experiments_filename_template is not None
            and "%s" in self.integrated_experiments_filename_template
        ):
            self.params.output.integrated_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.integrated_experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.params.dispatch.coset
            and self.coset_filename_template is not None
            and "%s" in self.coset_filename_template
        ):
            self.params.output.coset_filename = os.path.join(
                self.params.output.output_dir,
                self.coset_filename_template
                % ("idx-" + tag, self.params.integration.coset.transformation),
            )
        if (
            self.params.dispatch.coset
            and self.coset_experiments_filename_template is not None
            and "%s" in self.coset_experiments_filename_template
        ):
            self.params.output.coset_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.coset_experiments_filename_template
                % ("idx-" + tag, self.params.integration.coset.transformation),
            )

    def debug_start(self, tag):
        if not self.params.mp.debug.output_debug_logs:
            return

        import socket

        self.debug_str = f"{socket.gethostname()},{tag}"
        self.debug_str += ",%s,%s,%s\n"
        self.debug_write("start")

    def debug_write(self, string, state=None):
        if not self.params.mp.debug.output_debug_logs:
            return

        from serialtbx.util.time import timestamp  # XXX move to common timestamp format

        ts = timestamp()  # Now
        debug_file_handle = open(self.debug_file_path, "a")
        if string == "":
            debug_file_handle.write("\n")
        else:
            if state is None:
                state = "    "
            debug_file_handle.write(self.debug_str % (ts, state, string))
        debug_file_handle.close()

    def process_experiments(self, tag, experiments):
        if not self.params.output.composite_output:
            self.setup_filenames(tag)
        self.tag = tag
        self.debug_start(tag)

        if self.params.output.experiments_filename:
            if self.params.output.composite_output:
                self.all_imported_experiments.extend(experiments)
            else:
                experiments.as_json(self.params.output.experiments_filename)

        # Do the processing
        try:
            self.pre_process(experiments)
        except Exception as e:
            print("Error in pre-process", tag, str(e))
            self.debug_write("preprocess_exception", "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.find_spots:
                self.debug_write("spotfind_start")
                observed = self.find_spots(experiments)
            else:
                print("Spot Finding turned off. Exiting")
                self.debug_write("data_loaded", "done")
                return
        except Exception as e:
            print("Error spotfinding", tag, str(e))
            self.debug_write("spotfinding_exception", "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.index:
                if (
                    self.params.dispatch.hit_finder.enable
                    and len(observed)
                    < self.params.dispatch.hit_finder.minimum_number_of_reflections
                ):
                    print("Not enough spots to index", tag)
                    self.debug_write(f"not_enough_spots_{len(observed)}", "stop")
                    return
                if (
                    self.params.dispatch.hit_finder.maximum_number_of_reflections
                    is not None
                ):
                    if (
                        self.params.dispatch.hit_finder.enable
                        and len(observed)
                        > self.params.dispatch.hit_finder.maximum_number_of_reflections
                    ):
                        print("Too many spots to index - Possibly junk", tag)
                        self.debug_write(f"too_many_spots_{len(observed)}", "stop")
                        return
                self.debug_write("index_start")
                experiments, indexed = self.index(experiments, observed)
            else:
                print("Indexing turned off. Exiting")
                self.debug_write(f"spotfinding_ok_{len(observed)}", "done")
                return
        except Exception as e:
            print("Couldn't index", tag, str(e))
            if not self.params.dispatch.squash_errors:
                raise
            self.debug_write(f"indexing_failed_{len(observed)}", "stop")
            return
        self.debug_write("refine_start")
        try:
            experiments, indexed = self.refine(experiments, indexed)
        except Exception as e:
            print("Error refining", tag, str(e))
            self.debug_write(f"refine_failed_{len(indexed)}", "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.integrate:
                self.debug_write("integrate_start")
                integrated = self.integrate(experiments, indexed)
            else:
                print("Integration turned off. Exiting")
                self.debug_write(f"index_ok_{len(indexed)}", "done")
                return
        except Exception as e:
            print("Error integrating", tag, str(e))
            self.debug_write(f"integrate_failed_{len(indexed)}", "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        self.debug_write(f"integrate_ok_{len(integrated)}", "done")

    def pre_process(self, experiments):
        """Add any pre-processing steps here"""
        if (
            self.params.indexing.stills.known_orientations
            and self.params.indexing.stills.require_known_orientation
        ):
            for expt in experiments:
                assert (
                    len(expt.imageset.indices()) == 1
                    and len(expt.imageset.paths()) == 1
                )
                key = (
                    os.path.basename(expt.imageset.paths()[0]),
                    expt.imageset.indices()[0],
                )
                if key not in self.params.indexing.stills.known_orientations:
                    raise RuntimeError("Image not found in set of known orientations")

        if not self.params.input.ignore_gain_mismatch:
            g1 = self.params.spotfinder.threshold.dispersion.gain
            g2 = self.params.integration.summation.detector_gain
            gain = g1 if g1 is not None else g2
            if gain is not None and gain != 1.0:
                for detector in experiments.detectors():
                    for panel in detector:
                        if panel.get_gain() != 1.0 and panel.get_gain() != gain:
                            raise RuntimeError(
                                f"""
The detector is reporting a gain of {panel.get_gain():f} but you have also supplied a gain of {gain:f}. Since the detector gain is not 1.0, your supplied gain will be multiplicatively applied in addition to the detector's gain, which is unlikely to be correct. Please re-run, removing spotfinder.dispersion.gain and integration.summation.detector_gain from your parameters. You can override this exception by setting input.ignore_gain_mismatch=True."""
                            )

    def find_spots(self, experiments):
        st = time.time()

        logger.info("*" * 80)
        logger.info("Finding Strong Spots")
        logger.info("*" * 80)

        # Find the strong spots
        observed = flex.reflection_table.from_observations(
            experiments, self.params, is_stills=True
        )

        # Reset z coordinates for dials.image_viewer; see Issues #226 for details
        xyzobs = observed["xyzobs.px.value"]
        for i in range(len(xyzobs)):
            xyzobs[i] = (xyzobs[i][0], xyzobs[i][1], 0)
        bbox = observed["bbox"]
        for i in range(len(bbox)):
            bbox[i] = (bbox[i][0], bbox[i][1], bbox[i][2], bbox[i][3], 0, 1)

        if self.params.output.composite_output:
            n = len(self.all_strong_reflections.experiment_identifiers())
            for i, experiment in enumerate(experiments):
                refls = observed.select(observed["id"] == i)
                refls["id"] = flex.int(len(refls), n)
                del refls.experiment_identifiers()[i]
                refls.experiment_identifiers()[n] = experiment.identifier
                self.all_strong_reflections.extend(refls)
                n += 1
        else:
            # Save the reflections to file
            logger.info("\n" + "-" * 80)
            if self.params.output.strong_filename:
                self.save_reflections(observed, self.params.output.strong_filename)

        logger.info("")
        logger.info("Time Taken = %f seconds", time.time() - st)
        return observed

    def index(self, experiments, reflections):
        from dials.algorithms.indexing.indexer import Indexer

        st = time.time()

        logger.info("*" * 80)
        logger.info("Indexing Strong Spots")
        logger.info("*" * 80)

        params = copy.deepcopy(self.params)
        # don't do scan-varying refinement during indexing
        params.refinement.parameterisation.scan_varying = False

        if hasattr(self, "known_crystal_models"):
            known_crystal_models = self.known_crystal_models
        elif self.params.indexing.stills.known_orientations:
            known_crystal_models = []
            extended_experiments = ExperimentList()
            for expt in experiments:
                assert (
                    len(expt.imageset.indices()) == 1
                    and len(expt.imageset.paths()) == 1
                )
                key = (
                    os.path.basename(expt.imageset.paths()[0]),
                    expt.imageset.indices()[0],
                )
                if key not in self.params.indexing.stills.known_orientations:
                    if self.params.indexing.stills.require_known_orientation:
                        raise RuntimeError(
                            "Image not found in set of known orientations"
                        )
                    else:
                        oris = [None]
                else:
                    oris = self.params.indexing.stills.known_orientations[key]
                known_crystal_models.extend(oris)
                extended_experiments.extend(ExperimentList([expt] * len(oris)))
            experiments = extended_experiments
        else:
            known_crystal_models = None

        indexing_succeeded = False
        if known_crystal_models:
            try:
                idxr = Indexer.from_parameters(
                    reflections,
                    experiments,
                    known_crystal_models=known_crystal_models,
                    params=params,
                )
                idxr.index()
                logger.info("indexed from known orientation")
                indexing_succeeded = True
            except Exception:
                if self.params.indexing.stills.require_known_orientation:
                    raise

        if not indexing_succeeded:
            if self.params.indexing.stills.reflection_subsampling.enable:
                subsets = range(
                    self.params.indexing.stills.reflection_subsampling.step_start,
                    self.params.indexing.stills.reflection_subsampling.step_stop
                    - self.params.indexing.stills.reflection_subsampling.step_size,
                    -self.params.indexing.stills.reflection_subsampling.step_size,
                )
            else:
                subsets = [100]
            all_reflections = reflections
            subset_indexing_error = None
            for i in subsets:
                if i != 100:
                    reflections = all_reflections.select(
                        flex.random_permutation(len(all_reflections))
                    )[: int(len(all_reflections) * i / 100)]
                try:
                    if params.indexing.stills.method_list is None:
                        idxr = Indexer.from_parameters(
                            reflections,
                            experiments,
                            params=params,
                        )
                        idxr.index()
                    else:
                        ml_indexing_error = None
                        for method in params.indexing.stills.method_list:
                            params.indexing.method = method
                            try:
                                idxr = Indexer.from_parameters(
                                    reflections,
                                    experiments,
                                    params=params,
                                )
                                idxr.index()
                            except Exception as e_ml:
                                logger.info("Couldn't index using method %s", method)
                                ml_indexing_error = e_ml
                            else:
                                ml_indexing_error = None
                                break
                        if ml_indexing_error:
                            raise ml_indexing_error
                except Exception as e_subset:
                    subset_indexing_error = e_subset
                else:
                    logger.info("Indexed using %d%% of the reflections" % i)
                    subset_indexing_error = None
                    break
            if subset_indexing_error:
                raise subset_indexing_error

        indexed = idxr.refined_reflections
        experiments = idxr.refined_experiments

        if known_crystal_models is not None:
            filtered_sel = flex.bool(len(indexed), True)
            for expt_id in range(len(experiments)):
                for idx in set(
                    indexed["miller_index"].select(indexed["id"] == expt_id)
                ):
                    sel = (indexed["miller_index"] == idx) & (indexed["id"] == expt_id)
                    if sel.count(True) > 1:
                        filtered_sel = filtered_sel & ~sel
            filtered = indexed.select(filtered_sel)
            logger.info(
                "Filtered duplicate reflections, %d out of %d remaining",
                len(filtered),
                len(indexed),
            )
            print(
                "Filtered duplicate reflections, %d out of %d remaining"
                % (len(filtered), len(indexed))
            )
            indexed = filtered

        logger.info("")
        logger.info("Time Taken = %f seconds", time.time() - st)
        return experiments, indexed

    def refine(self, experiments, centroids):
        if self.params.dispatch.refine:
            from dials.algorithms.refinement import RefinerFactory

            st = time.time()

            logger.info("*" * 80)
            logger.info("Refining Model")
            logger.info("*" * 80)

            refiner = RefinerFactory.from_parameters_data_experiments(
                self.params, centroids, experiments
            )

            refiner.run()
            experiments = refiner.get_experiments()
            predicted = refiner.predict_for_indexed()
            centroids["xyzcal.mm"] = predicted["xyzcal.mm"]
            centroids["entering"] = predicted["entering"]
            centroids = centroids.select(refiner.selection_used_for_refinement())

            # Re-estimate mosaic estimates
            from dials.algorithms.indexing.nave_parameters import NaveParameters

            nv = NaveParameters(
                params=self.params,
                experiments=experiments,
                reflections=centroids,
                refinery=refiner,
                graph_verbose=False,
            )
            nv()
            acceptance_flags_nv = nv.nv_acceptance_flags
            centroids = centroids.select(acceptance_flags_nv)

        if self.params.output.composite_output:
            if (
                self.params.output.refined_experiments_filename
                or self.params.output.indexed_filename
            ):
                assert (
                    self.params.output.refined_experiments_filename is not None
                    and self.params.output.indexed_filename is not None
                )

                n = len(self.all_indexed_experiments)
                self.all_indexed_experiments.extend(experiments)
                for i, experiment in enumerate(experiments):
                    refls = centroids.select(centroids["id"] == i)
                    refls["id"] = flex.int(len(refls), n)
                    del refls.experiment_identifiers()[i]
                    refls.experiment_identifiers()[n] = experiment.identifier
                    self.all_indexed_reflections.extend(refls)
                    n += 1
        else:
            # Dump experiments to disk
            if self.params.output.refined_experiments_filename:
                experiments.as_json(self.params.output.refined_experiments_filename)

            if self.params.output.indexed_filename:
                self.save_reflections(centroids, self.params.output.indexed_filename)

        if self.params.dispatch.refine:
            logger.info("")
            logger.info("Time Taken = %f seconds", time.time() - st)

        return experiments, centroids

    def integrate(self, experiments, indexed):
        st = time.time()

        logger.info("*" * 80)
        logger.info("Integrating Reflections")
        logger.info("*" * 80)

        if self.params.integration.integration_only_overrides.trusted_range:
            for detector in experiments.detectors():
                for panel in detector:
                    panel.set_trusted_range(
                        self.params.integration.integration_only_overrides.trusted_range
                    )

        if self.params.profile.algorithm == "ellipsoid":
            from dials.command_line.ssx_integrate import (
                EllipsoidIntegrator,
                process_one_image,
            )

            experiment, integrated, _ = process_one_image(
                experiments[0],
                indexed,
                self.params,
                EllipsoidIntegrator,
                collect_data=False,
            )
            experiments = ExperimentList([experiment])

        else:
            indexed, _ = self.process_reference(indexed)

            if self.params.dispatch.coset:
                from dials.algorithms.integration.sublattice_helper import (
                    integrate_coset,
                )

                integrate_coset(self, experiments, indexed)

            # Get the integrator from the input parameters
            logger.info("Configuring integrator from input parameters")
            from dials.algorithms.integration.integrator import create_integrator
            from dials.algorithms.profile_model.factory import ProfileModelFactory

            # Compute the profile model
            # Predict the reflections
            # Match the predictions with the reference
            # Create the integrator
            experiments = ProfileModelFactory.create(self.params, experiments, indexed)
            new_experiments = ExperimentList()
            new_reflections = flex.reflection_table()
            for expt_id, expt in enumerate(experiments):
                if (
                    self.params.profile.gaussian_rs.parameters.sigma_b_cutoff is None
                    or expt.profile.sigma_b()
                    < self.params.profile.gaussian_rs.parameters.sigma_b_cutoff
                ):
                    refls = indexed.select(indexed["id"] == expt_id)
                    refls["id"] = flex.int(len(refls), len(new_experiments))
                    # refls.reset_ids()
                    del refls.experiment_identifiers()[expt_id]
                    refls.experiment_identifiers()[len(new_experiments)] = (
                        expt.identifier
                    )
                    new_reflections.extend(refls)
                    new_experiments.append(expt)
                else:
                    logger.info(
                        "Rejected expt %d with sigma_b %f"
                        % (expt_id, expt.profile.sigma_b())
                    )
            experiments = new_experiments
            indexed = new_reflections
            if len(experiments) == 0:
                raise RuntimeError("No experiments after filtering by sigma_b")
            logger.info("")
            logger.info("=" * 80)
            logger.info("")
            logger.info("Predicting reflections")
            logger.info("")
            predicted = flex.reflection_table.from_predictions_multi(
                experiments,
                dmin=self.params.prediction.d_min,
                dmax=self.params.prediction.d_max,
                margin=self.params.prediction.margin,
                force_static=self.params.prediction.force_static,
            )
            predicted.match_with_reference(indexed)
            logger.info("")
            integrator = create_integrator(self.params, experiments, predicted)

            # Integrate the reflections
            integrated = integrator.integrate()

        # correct integrated intensities for absorption correction, if necessary
        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if abs_params.algorithm == "fuller_kapton":
                    from dials.algorithms.integration.kapton_correction import (
                        multi_kapton_correction,
                    )
                elif abs_params.algorithm == "kapton_2019":
                    from dials.algorithms.integration.kapton_2019_correction import (
                        multi_kapton_correction,
                    )
                elif abs_params.algorithm == "other":
                    continue  # custom abs. corr. implementation should go here
                else:
                    raise ValueError(
                        "absorption_correction.apply=True, "
                        "but no .algorithm has been selected!"
                    )
                experiments, integrated = multi_kapton_correction(
                    experiments, integrated, abs_params.fuller_kapton, logger=logger
                )()

        if self.params.significance_filter.enable:
            from dials.algorithms.integration.stills_significance_filter import (
                SignificanceFilter,
            )

            sig_filter = SignificanceFilter(self.params)
            filtered_refls = sig_filter(experiments, integrated)
            accepted_expts = ExperimentList()
            accepted_refls = flex.reflection_table()
            logger.info(
                "Removed %d reflections out of %d when applying significance filter",
                len(integrated) - len(filtered_refls),
                len(integrated),
            )
            for expt_id, expt in enumerate(experiments):
                refls = filtered_refls.select(filtered_refls["id"] == expt_id)
                if len(refls) > 0:
                    accepted_expts.append(expt)
                    refls["id"] = flex.int(len(refls), len(accepted_expts) - 1)
                    accepted_refls.extend(refls)
                else:
                    logger.info(
                        "Removed experiment %d which has no reflections left after applying significance filter",
                        expt_id,
                    )

            if len(accepted_refls) == 0:
                raise Sorry("No reflections left after applying significance filter")
            experiments = accepted_expts
            integrated = accepted_refls

        # Delete the shoeboxes used for intermediate calculations, if requested
        if self.params.integration.debug.delete_shoeboxes and "shoebox" in integrated:
            del integrated["shoebox"]

        if self.params.output.composite_output:
            if (
                self.params.output.integrated_experiments_filename
                or self.params.output.integrated_filename
            ):
                assert (
                    self.params.output.integrated_experiments_filename is not None
                    and self.params.output.integrated_filename is not None
                )

                n = len(self.all_integrated_experiments)
                self.all_integrated_experiments.extend(experiments)
                for i, experiment in enumerate(experiments):
                    refls = integrated.select(integrated["id"] == i)
                    refls["id"] = flex.int(len(refls), n)
                    del refls.experiment_identifiers()[i]
                    refls.experiment_identifiers()[n] = experiment.identifier
                    self.all_integrated_reflections.extend(refls)
                    n += 1
        else:
            # Dump experiments to disk
            if self.params.output.integrated_experiments_filename:
                experiments.as_json(self.params.output.integrated_experiments_filename)

            if self.params.output.integrated_filename:
                # Save the reflections
                self.save_reflections(
                    integrated, self.params.output.integrated_filename
                )

        self.write_integration_pickles(integrated, experiments)
        from dials.algorithms.indexing.stills_indexer import (
            calc_2D_rmsd_and_displacements,
        )

        rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
        log_str = f"RMSD indexed (px): {rmsd_indexed:f}\n"
        for i in range(6):
            bright_integrated = integrated.select(
                (
                    integrated["intensity.sum.value"]
                    / flex.sqrt(integrated["intensity.sum.variance"])
                )
                >= i
            )
            if len(bright_integrated) > 0:
                rmsd_integrated, _ = calc_2D_rmsd_and_displacements(bright_integrated)
            else:
                rmsd_integrated = 0
            log_str += (
                "N reflections integrated at I/sigI >= %d: % 4d, RMSD (px): %f\n"
                % (i, len(bright_integrated), rmsd_integrated)
            )

        for crystal_model in experiments.crystals():
            if hasattr(crystal_model, "get_domain_size_ang"):
                log_str += f". Final ML model: domain size angstroms: {crystal_model.get_domain_size_ang():f}, half mosaicity degrees: {crystal_model.get_half_mosaicity_deg():f}"

        logger.info(log_str)

        logger.info("")
        logger.info("Time Taken = %f seconds", time.time() - st)
        return integrated

    def write_integration_pickles(self, integrated, experiments, callback=None):
        """
        Write a serialized python dictionary with integrated intensities and other information
        suitible for use by cxi.merge or prime.postrefine.
        @param integrated Reflection table with integrated intensities
        @param experiments Experiment list. One integration pickle for each experiment will be created.
        @param callback Deriving classes can use callback to make further modifications to the dictionary
        before it is serialized. Callback should be a function with this signature:
        def functionname(params, outfile, frame), where params is the phil scope, outfile is the path
        to the pickle that will be saved, and frame is the python dictionary to be serialized.
        """
        if not hasattr(self.params.output, "integration_pickle"):
            return

        if self.params.output.integration_pickle is not None:
            from serialtbx.util.construct_frame import ConstructFrame

            # Split everything into separate experiments for pickling
            for e_number, experiment in enumerate(experiments):
                e_selection = integrated["id"] == e_number
                reflections = integrated.select(e_selection)

                frame = ConstructFrame(reflections, experiment).make_frame()
                frame["pixel_size"] = experiment.detector[0].get_pixel_size()[0]

                if not hasattr(self, "tag") or self.tag is None:
                    try:
                        # if the data was a file on disc, get the path
                        event_timestamp = os.path.splitext(
                            experiments[0].imageset.paths()[0]
                        )[0]
                    except NotImplementedError:
                        # if the data is in memory only, check if the reader set a timestamp on the format object
                        event_timestamp = (
                            experiment.imageset.reader().get_format(0).timestamp
                        )
                    event_timestamp = os.path.basename(event_timestamp)
                    if event_timestamp.find("shot-") == 0:
                        event_timestamp = os.path.splitext(event_timestamp)[
                            0
                        ]  # micromanage the file name
                else:
                    event_timestamp = self.tag
                if hasattr(self.params.output, "output_dir"):
                    outfile = os.path.join(
                        self.params.output.output_dir,
                        self.params.output.integration_pickle
                        % (e_number, event_timestamp),
                    )
                else:
                    outfile = os.path.join(
                        os.path.dirname(self.params.output.integration_pickle),
                        self.params.output.integration_pickle
                        % (e_number, event_timestamp),
                    )

                if callback is not None:
                    callback(self.params, outfile, frame)

                if self.params.output.composite_output:
                    self.all_int_pickle_filenames.append(os.path.basename(outfile))
                    self.all_int_pickles.append(frame)
                else:
                    with open(outfile, "wb") as fh:
                        pickle.dump(frame, fh, protocol=pickle.HIGHEST_PROTOCOL)

    def process_reference(self, reference):
        """Load the reference spots."""
        if reference is None:
            return None, None
        st = time.time()
        assert "miller_index" in reference
        assert "id" in reference
        logger.info("Processing reference reflections")
        logger.info(" read %d strong spots", len(reference))
        mask = reference.get_flags(reference.flags.indexed)
        rubbish = reference.select(~mask)
        if mask.count(False) > 0:
            reference.del_selected(~mask)
            logger.info(" removing %d unindexed reflections", mask.count(True))
        if len(reference) == 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        Expected > %d indexed spots, got %d
      """
                % (0, len(reference))
            )
        mask = reference["miller_index"] == (0, 0, 0)
        if mask.count(True) > 0:
            rubbish.extend(reference.select(mask))
            reference.del_selected(mask)
            logger.info(" removing %d reflections with hkl (0,0,0)", mask.count(True))
        mask = reference["id"] < 0
        if mask.count(True) > 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        %d reference spots have an invalid experiment id
      """
                % mask.count(True)
            )
        logger.info(" using %d indexed reflections", len(reference))
        logger.info(" found %d junk reflections", len(rubbish))
        logger.info(" time taken: %g", time.time() - st)
        return reference, rubbish

    def save_reflections(self, reflections, filename):
        """Save the reflections to file."""
        st = time.time()
        logger.info("Saving %d reflections to %s", len(reflections), filename)
        reflections.as_file(filename)
        logger.info(" time taken: %g", time.time() - st)

    def finalize(self):
        """Perform any final operations"""
        if self.params.output.composite_output:
            if self.params.mp.composite_stride is not None:
                assert self.params.mp.method == "mpi"
                stride = self.params.mp.composite_stride

                from mpi4py import MPI

                comm = MPI.COMM_WORLD
                rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
                size = comm.Get_size()  # size: number of processes running in this job
                comm.barrier()

                if rank % stride == 0:
                    subranks = [rank + i for i in range(1, stride) if rank + i < size]
                    for i in range(len(subranks)):
                        logger.info("Rank %d waiting for sender", rank)
                        (
                            sender,
                            imported_experiments,
                            strong_reflections,
                            indexed_experiments,
                            indexed_reflections,
                            integrated_experiments,
                            integrated_reflections,
                            coset_experiments,
                            coset_reflections,
                            int_pickles,
                            int_pickle_filenames,
                        ) = comm.recv(source=MPI.ANY_SOURCE)
                        logger.info("Rank %d received data from rank %d", rank, sender)

                        def extend_with_bookkeeping(
                            src_expts, src_refls, dest_expts, dest_refls
                        ):
                            n = len(dest_refls.experiment_identifiers())
                            src_refls["id"] += n
                            idents = src_refls.experiment_identifiers()
                            keys = idents.keys()
                            values = idents.values()
                            for key in keys:
                                del idents[key]
                            for i, key in enumerate(keys):
                                idents[key + n] = values[i]
                            dest_expts.extend(src_expts)
                            dest_refls.extend(src_refls)

                        if len(imported_experiments) > 0:
                            extend_with_bookkeeping(
                                imported_experiments,
                                strong_reflections,
                                self.all_imported_experiments,
                                self.all_strong_reflections,
                            )

                        if len(indexed_experiments) > 0:
                            extend_with_bookkeeping(
                                indexed_experiments,
                                indexed_reflections,
                                self.all_indexed_experiments,
                                self.all_indexed_reflections,
                            )

                        if len(integrated_experiments) > 0:
                            extend_with_bookkeeping(
                                integrated_experiments,
                                integrated_reflections,
                                self.all_integrated_experiments,
                                self.all_integrated_reflections,
                            )

                        if len(coset_experiments) > 0:
                            extend_with_bookkeeping(
                                coset_experiments,
                                coset_reflections,
                                self.all_coset_experiments,
                                self.all_coset_reflections,
                            )

                        self.all_int_pickles.extend(int_pickles)
                        self.all_int_pickle_filenames.extend(int_pickle_filenames)

                else:
                    destrank = (rank // stride) * stride
                    logger.info(
                        "Rank %d sending results to rank %d",
                        rank,
                        (rank // stride) * stride,
                    )
                    comm.send(
                        (
                            rank,
                            self.all_imported_experiments,
                            self.all_strong_reflections,
                            self.all_indexed_experiments,
                            self.all_indexed_reflections,
                            self.all_integrated_experiments,
                            self.all_integrated_reflections,
                            self.all_coset_experiments,
                            self.all_coset_reflections,
                            self.all_int_pickles,
                            self.all_int_pickle_filenames,
                        ),
                        dest=destrank,
                    )

                    self.all_imported_experiments = self.all_strong_reflections = (
                        self.all_indexed_experiments
                    ) = self.all_indexed_reflections = (
                        self.all_integrated_experiments
                    ) = self.all_integrated_reflections = self.all_coset_experiments = (
                        self.all_coset_reflections
                    ) = self.all_int_pickles = self.all_integrated_reflections = []

            # Dump composite files to disk
            if (
                len(self.all_imported_experiments) > 0
                and self.params.output.experiments_filename
            ):
                self.all_imported_experiments.as_json(
                    self.params.output.experiments_filename
                )

            if (
                len(self.all_strong_reflections) > 0
                and self.params.output.strong_filename
            ):
                self.save_reflections(
                    self.all_strong_reflections, self.params.output.strong_filename
                )

            if (
                len(self.all_indexed_experiments) > 0
                and self.params.output.refined_experiments_filename
            ):
                self.all_indexed_experiments.as_json(
                    self.params.output.refined_experiments_filename
                )

            if (
                len(self.all_indexed_reflections) > 0
                and self.params.output.indexed_filename
            ):
                self.save_reflections(
                    self.all_indexed_reflections, self.params.output.indexed_filename
                )

            if (
                len(self.all_integrated_experiments) > 0
                and self.params.output.integrated_experiments_filename
            ):
                self.all_integrated_experiments.as_json(
                    self.params.output.integrated_experiments_filename
                )

            if (
                len(self.all_integrated_reflections) > 0
                and self.params.output.integrated_filename
            ):
                self.save_reflections(
                    self.all_integrated_reflections,
                    self.params.output.integrated_filename,
                )

            if self.params.dispatch.coset:
                if (
                    len(self.all_coset_experiments) > 0
                    and self.params.output.coset_experiments_filename
                ):
                    self.all_coset_experiments.as_json(
                        self.params.output.coset_experiments_filename
                    )

                if (
                    len(self.all_coset_reflections) > 0
                    and self.params.output.coset_filename
                ):
                    self.save_reflections(
                        self.all_coset_reflections, self.params.output.coset_filename
                    )

            # Create a tar archive of the integration dictionary pickles
            if len(self.all_int_pickles) > 0 and self.params.output.integration_pickle:
                tar_template_integration_pickle = (
                    self.params.output.integration_pickle.replace("%d", "%s")
                )
                outfile = (
                    os.path.join(
                        self.params.output.output_dir,
                        tar_template_integration_pickle % ("x", self.composite_tag),
                    )
                    + ".tar"
                )
                tar = tarfile.TarFile(outfile, "w")
                for i, (fname, d) in enumerate(
                    zip(self.all_int_pickle_filenames, self.all_int_pickles)
                ):
                    string = BytesIO(pickle.dumps(d, protocol=2))
                    info = tarfile.TarInfo(name=fname)
                    info.size = string.getbuffer().nbytes
                    info.mtime = time.time()
                    tar.addfile(tarinfo=info, fileobj=string)
                tar.close()


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()
