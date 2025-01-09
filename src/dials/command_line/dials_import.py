# LIBTBX_SET_DISPATCHER_NAME dials.import

from __future__ import annotations

import logging
import pickle
import sys
from collections import defaultdict, namedtuple

from ordered_set import OrderedSet

import dxtbx.model.compare as compare
import libtbx.phil
from dxtbx.imageset import ImageSequence
from dxtbx.model.experiment_list import (
    Experiment,
    ExperimentList,
    ExperimentListFactory,
)
from dxtbx.sequence_filenames import template_regex_from_list

from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.multi_dataset_handling import generate_experiment_identifiers
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.import")


def _pickle_load(fh):
    return pickle.load(fh, encoding="bytes")


help_message = """

This program is used to import image data files into a format that can be used
within dials. The program looks at the metadata for each image along with the
filenames to determine the relationship between sets of images. Once all the
images have been analysed, a experiments object is written to file which specifies
the relationship between files. For example if two sets of images which belong
to two rotation scans have been given, two image sequences will be saved. Images to
be processed are specified as command line arguments. Sometimes, there is a
maximum number of arguments that can be given on the command line and the number
of files may exceed this. In this case image filenames can be input on stdin
as shown in the examples below. Alternatively a template can be specified using
the template= parameter where the consecutive digits representing the image
numbers in the filenames are replaced with '#' characters.

The geometry can be set manually, either by using the reference_geometry=
parameter to specify an experiment list .expt file containing
the reference geometry, by using the mosflm_beam_centre= parameter to set
the Mosflm beam centre, or by specifying each variable to be overridden
using various geometry parameters.

Examples::

  dials.import /data/directory-containing-images/

  dials.import image_*.cbf

  dials.import image_1_*.cbf image_2_*.cbf

  dials.import directory/with/images

  dials.import template=image_1_####.cbf

  dials.import directory=directory/with/images

  find . -name "image_*.cbf" | dials.import

  dials.import << EOF
  image_1.cbf
  image_2.cbf
  EOF
"""


# Create the phil parameters
phil_scope = libtbx.phil.parse(
    """

  output {

    experiments = imported.expt
      .type = str
      .help = "The output experiment file"

    log = 'dials.import.log'
      .type = str
      .help = "The log filename"

    compact = False
      .type = bool
      .help = "For experiment output use compact JSON representation"

  }

  identifier_type = *uuid timestamp None
    .type = choice
    .help = "Type of unique identifier to generate."

  input {

    ignore_unhandled = True
      .type = bool
      .help = "Ignore unhandled input (e.g. log files)"

    template = None
      .type = str
      .help = "The image sequence template"
      .multiple = True

    directory = None
      .type = str
      .help = "A directory with images"
      .multiple = True

    split = None
      .type = ints
      .help = "Scan split: either frames_per_block or 1-indexed start,end,frames_per_block"

    reference_geometry = None
      .type = path
      .help = "Experimental geometry from this models.expt "
              "will override the geometry from the "
              "image headers."

    check_reference_geometry = True
      .type = bool
      .expert_level = 2
      .help = "If True, assert the reference geometry is similar to"
              "the image geometry"

    use_beam_reference = True
      .type = bool
      .expert_level = 2
      .help = "If True, the beam from reference_geometry will override "
              "the beam from the image headers."

    use_gonio_reference = True
      .type = bool
      .expert_level = 2
      .help = "If True, the goniometer from reference_geometry will override "
              "the goniometer from the image headers."

    use_detector_reference = True
      .type = bool
      .expert_level = 2
      .help = "If True, the detector from reference_geometry will override "
              "the detector from the image headers."

    allow_multiple_sequences = True
      .type = bool
      .help = "If False, raise an error if multiple sequences are found"

  }

  include scope dials.util.options.format_phil_scope

  include scope dials.util.options.geometry_phil_scope

  lookup {
    mask = None
      .type = str
      .help = "Apply a mask to the imported data"

    gain = None
      .type = str
      .help = "Apply a gain to the imported data"

    pedestal = None
      .type = str
      .help = "Apply a pedestal to the imported data"

    dx = None
      .type = str
      .help = "Apply an x geometry offset"
              "If both dx and dy are set then"
              "OffsetParallaxCorrectedPxMmStrategy will be used"

    dy = None
      .type = str
      .help = "Apply an y geometry offset"
              "If both dx and dy are set then"
              "OffsetParallaxCorrectedPxMmStrategy will be used"
  }
""",
    process_includes=True,
)


def _extract_or_read_imagesets(params):
    """
    Return a list of ImageSets, importing them via alternative means if necessary.

    The "Alternative Means" means via params.input.template or .directory,
    if the images to import haven't been specified directly.

    Args:
        params: The phil.scope_extract from dials.import

    Returns: A list of ImageSet objects
    """

    # Get the experiments
    experiments = flatten_experiments(params.input.experiments)

    # Check we have some filenames
    if len(experiments) == 0:
        # FIXME Should probably make this smarter since it requires editing here
        # and in dials.import phil scope
        try:
            format_kwargs = {
                "dynamic_shadowing": params.format.dynamic_shadowing,
                "multi_panel": params.format.multi_panel,
            }
        except AttributeError:
            format_kwargs = None

        # Check if a template has been set and print help if not, otherwise try to
        # import the images based on the template input
        if len(params.input.template) > 0:
            experiments = ExperimentListFactory.from_templates(
                params.input.template,
                image_range=params.geometry.scan.image_range,
                format_kwargs=format_kwargs,
            )
            if len(experiments) == 0:
                raise Sorry(
                    "No experiments found matching template %s"
                    % params.input.experiments
                )
        elif len(params.input.directory) > 0:
            experiments = ExperimentListFactory.from_filenames(
                params.input.directory, format_kwargs=format_kwargs
            )
            if len(experiments) == 0:
                raise Sorry(
                    "No experiments found in directories %s" % params.input.directory
                )
        else:
            raise Sorry("No experiments found")

    # TODO (Nick):  This looks redundant as the experiments are immediately discarded.
    #               verify this, and remove if it is.
    if params.identifier_type:
        generate_experiment_identifiers(experiments, params.identifier_type)

    # Get a list of all imagesets
    imageset_list = experiments.imagesets()

    # Return the experiments
    return imageset_list


class ReferenceGeometryUpdater:
    """
    A class to replace beam + detector with a reference
    """

    def __init__(self, params):
        """
        Load the reference geometry
        """
        self.params = params
        self.reference = self.load_reference_geometry(params)

    def __call__(self, imageset):
        """
        Replace with the reference geometry
        """
        if self.params.input.check_reference_geometry:
            # Check static detector items are the same
            assert self.reference.detector.is_similar_to(
                imageset.get_detector(), static_only=True
            ), "Reference detector model does not match input detector model"

        # Set beam and detector
        if self.params.input.use_beam_reference:
            imageset.set_beam(self.reference.beam)
        if self.params.input.use_detector_reference:
            imageset.set_detector(self.reference.detector)
        if self.params.input.use_gonio_reference:
            imageset.set_goniometer(self.reference.goniometer)
        return imageset

    def load_reference_geometry(self, params):
        """
        Load a reference geometry file
        """
        # Load reference geometry
        reference_detector = None
        reference_beam = None
        reference_goniometer = None
        if params.input.reference_geometry is not None:
            from dxtbx.serialize import load

            experiments = None
            experiments = load.experiment_list(
                params.input.reference_geometry, check_format=False
            )
            assert experiments, "Could not import reference geometry"
            assert len(experiments.detectors()) >= 1
            if len(experiments.detectors()) > 1:
                raise Sorry(
                    "The reference geometry file contains %d detector definitions, but only a single definition is allowed."
                    % len(experiments.detectors())
                )
            reference_detector = experiments.detectors()[0]
            if self.params.input.use_beam_reference:
                assert len(experiments.beams()) >= 1
                if len(experiments.beams()) > 1:
                    raise Sorry(
                        "The reference geometry file contains %d beam definitions, but only a single definition is allowed."
                        % len(experiments.beams())
                    )
                reference_beam = experiments.beams()[0]
            if self.params.input.use_gonio_reference:
                assert len(experiments.goniometers()) >= 1
                reference_goniometer = experiments.goniometers()[0]
        Reference = namedtuple("Reference", ["detector", "beam", "goniometer"])
        return Reference(
            detector=reference_detector,
            beam=reference_beam,
            goniometer=reference_goniometer,
        )


class ManualGeometryUpdater:
    """
    A class to update the geometry manually
    """

    def __init__(self, params):
        """
        Save the params
        """
        self.params = params
        self.touched = set()

    def __call__(self, imageset):
        """
        Override the parameters
        """
        from copy import deepcopy

        from dxtbx.imageset import ImageSequence, ImageSetFactory
        from dxtbx.model import (
            BeamFactory,
            DetectorFactory,
            GoniometerFactory,
            ScanFactory,
        )

        if self.params.geometry.convert_sequences_to_stills:
            imageset = ImageSetFactory.imageset_from_anyset(imageset)
            for j in range(len(imageset)):
                imageset.set_scan(None, j)
                imageset.set_goniometer(None, j)

        beam = imageset.get_beam()
        detector = imageset.get_detector()
        goniometer = imageset.get_goniometer()
        scan = imageset.get_scan()

        # Create a new model with updated geometry for each model that is not
        # already in the touched set
        if isinstance(imageset, ImageSequence):
            if beam and beam not in self.touched:
                beam = BeamFactory.from_phil(self.params.geometry, imageset.get_beam())
            if detector and detector not in self.touched:
                detector = DetectorFactory.from_phil(
                    self.params.geometry, imageset.get_detector(), beam
                )
            if goniometer and goniometer not in self.touched:
                goniometer = GoniometerFactory.from_phil(
                    self.params.geometry, imageset.get_goniometer()
                )
            if scan and scan not in self.touched:
                scan = ScanFactory.from_phil(
                    self.params.geometry, deepcopy(imageset.get_scan())
                )
            i0, i1 = scan.get_array_range()
            j0, j1 = imageset.get_scan().get_array_range()
            if i0 < j0 or i1 > j1:
                imageset = self.extrapolate_imageset(
                    imageset=imageset,
                    beam=beam,
                    detector=detector,
                    goniometer=goniometer,
                    scan=scan,
                )
            else:
                imageset.set_beam(beam)
                imageset.set_detector(detector)
                imageset.set_goniometer(goniometer)
                imageset.set_scan(scan)
        else:
            for i in range(len(imageset)):
                if beam and beam not in self.touched:
                    beam = BeamFactory.from_phil(
                        self.params.geometry, imageset.get_beam(i)
                    )
                if detector and detector not in self.touched:
                    detector = DetectorFactory.from_phil(
                        self.params.geometry, imageset.get_detector(i), beam
                    )
                if goniometer and goniometer not in self.touched:
                    goniometer = GoniometerFactory.from_phil(
                        self.params.geometry, imageset.get_goniometer(i)
                    )
                if scan and scan not in self.touched:
                    scan = ScanFactory.from_phil(
                        self.params.geometry, imageset.get_scan(i)
                    )
                imageset.set_beam(beam, i)
                imageset.set_detector(detector, i)
                imageset.set_goniometer(goniometer, i)
                imageset.set_scan(scan, i)

        # Add the models from this imageset to the touched set, so they will not
        # have their geometry updated again
        if beam:
            self.touched.add(beam)
        if detector:
            self.touched.add(detector)
        if goniometer:
            self.touched.add(goniometer)
        if scan:
            self.touched.add(scan)
        return imageset

    def extrapolate_imageset(
        self, imageset=None, beam=None, detector=None, goniometer=None, scan=None
    ):
        from dxtbx.imageset import ImageSetFactory

        first, last = scan.get_image_range()
        sequence = ImageSetFactory.make_sequence(
            template=imageset.get_template(),
            indices=list(range(first, last + 1)),
            format_class=imageset.get_format_class(),
            beam=beam,
            detector=detector,
            goniometer=goniometer,
            scan=scan,
            format_kwargs=imageset.params(),
        )
        return sequence


class MetaDataUpdater:
    """
    A class to manage updating the experiments metadata
    """

    def __init__(self, params):
        """
        Init the class
        """
        from dials.util.options import geometry_phil_scope

        self.params = params

        # Create the geometry updater
        self.update_geometry = []
        update_order = []

        # First add reference geometry is present
        if self.params.input.reference_geometry is not None:
            self.update_geometry.append(ReferenceGeometryUpdater(self.params))
            update_order.append("Reference geometry")

        # Then add manual geometry
        working_phil = geometry_phil_scope.format(self.params)
        diff_phil = geometry_phil_scope.fetch_diff(source=working_phil)
        if diff_phil.as_str() != "":
            self.update_geometry.append(ManualGeometryUpdater(self.params))
            update_order.append("Manual geometry")

        if len(update_order) > 0:
            logger.info("")
            logger.info("Applying input geometry in the following order:")
            for i, item in enumerate(update_order, start=1):
                logger.info("  %d. %s", i, item)
            logger.info("")

    def __call__(self, imageset_list):
        """
        Transform the metadata
        """
        # Import the lookup data
        lookup = self.import_lookup_data(self.params)

        # Create the experiments
        experiments = ExperimentList()

        # Loop through imagesets
        for imageset in imageset_list:
            # Set the external lookups
            imageset = self.update_lookup(imageset, lookup)

            # Update the geometry
            for updater in self.update_geometry:
                imageset = updater(imageset)

            # Check beam and detector are present
            if imageset.get_beam() is None or imageset.get_detector() is None:
                raise Sorry(
                    """
          Imageset contains no beam or detector model. This means you will be
          unable to process your data.

          Possible causes of this error are:
             - A problem reading the images with one of the dxtbx format classes
             - A lack of header information in the file itself.

          You can override this by specifying the metadata as geometry parameters
        """
                )

            # Check if dx and dy are set
            if [
                imageset.external_lookup.dx.filename,
                imageset.external_lookup.dy.filename,
            ].count(None) == 0:
                imageset.update_detector_px_mm_data()
            elif [
                imageset.external_lookup.dx.filename,
                imageset.external_lookup.dy.filename,
            ].count(None) == 1:
                raise Sorry(
                    """
          Only 1 offset map is set. Need to set both dx and dy
        """
                )

            # Append to new imageset list
            if isinstance(imageset, ImageSequence):
                if imageset.get_scan().is_still():
                    # make lots of experiments all pointing at one
                    # image set

                    # check if user has overridden the input - if yes, recall
                    # that these are in people numbers (1...) and are inclusive
                    if self.params.geometry.scan.image_range:
                        user_start, user_end = self.params.geometry.scan.image_range
                        start, end = user_start - 1, user_end
                    else:
                        start, end = imageset.get_scan().get_array_range()

                    # offset to get 0-based indexing into the imageset
                    offset = imageset.get_scan().get_array_range()[0]
                    for j in range(start, end):
                        subset = imageset[j - offset : j - offset + 1]
                        experiments.append(
                            Experiment(
                                imageset=imageset,
                                beam=imageset.get_beam(),
                                detector=imageset.get_detector(),
                                goniometer=imageset.get_goniometer(),
                                scan=subset.get_scan(),
                                crystal=None,
                            )
                        )
                else:
                    # have just one experiment
                    experiments.append(
                        Experiment(
                            imageset=imageset,
                            beam=imageset.get_beam(),
                            detector=imageset.get_detector(),
                            goniometer=imageset.get_goniometer(),
                            scan=imageset.get_scan(),
                            crystal=None,
                        )
                    )
            else:
                for i in range(len(imageset)):
                    experiments.append(
                        Experiment(
                            imageset=imageset[i : i + 1],
                            beam=imageset.get_beam(i),
                            detector=imageset.get_detector(i),
                            goniometer=imageset.get_goniometer(i),
                            scan=imageset.get_scan(i),
                            crystal=None,
                        )
                    )

        if self.params.geometry.convert_stills_to_sequences:
            if any(not isinstance(i, ImageSequence) for i in experiments.imagesets()):
                files_to_indiv = defaultdict(int)
                beams = []
                formats = []
                detectors = []
                iset_params = []
                existing_isets_sequences = []
                for i, iset in enumerate(experiments.imagesets()):
                    if not isinstance(iset, ImageSequence):
                        path = iset.get_path(0)
                        if path not in files_to_indiv:
                            beams.append(iset.get_beam())
                            detectors.append(iset.get_detector())
                            formats.append(iset.get_format_class())
                            iset_params.append(iset.params())
                        files_to_indiv[iset.get_path(0)] += 1
                    else:
                        existing_isets_sequences.append(iset)

                from dxtbx.imageset import ImageSetFactory
                from dxtbx.model import GoniometerFactory, Scan

                new_experiments = ExperimentList()
                for i, (file, n) in enumerate(files_to_indiv.items()):
                    if self.params.geometry.scan.image_range:
                        user_start, user_end = self.params.geometry.scan.image_range
                        first, last = user_start, user_end
                    else:
                        first, last = 1, n
                    iset_params[i].update({"lazy": False})
                    sequence = ImageSetFactory.make_sequence(
                        template=file,
                        indices=list(range(first, last + 1)),
                        format_class=formats[i],
                        beam=beams[i],
                        detector=detectors[i],
                        goniometer=GoniometerFactory.make_goniometer(
                            (0.0, 1.0, 0.0), (1, 0, 0, 0, 1, 0, 0, 0, 1)
                        ),
                        scan=Scan(image_range=(first, last), oscillation=(0.0, 0.0)),
                        format_kwargs=iset_params[i],
                    )
                    sequence = self.update_lookup(sequence, lookup)  # for mask etc
                    for j in range(len(sequence)):
                        subset = sequence[j : j + 1]
                        new_experiments.append(
                            Experiment(
                                imageset=sequence,
                                beam=sequence.get_beam(),
                                detector=sequence.get_detector(),
                                goniometer=sequence.get_goniometer(),
                                scan=subset.get_scan(),
                                crystal=None,
                            )
                        )
                if existing_isets_sequences:
                    for expt in experiments:
                        if expt.imageset in existing_isets_sequences:
                            new_experiments.append(expt)
                experiments = new_experiments

        if self.params.identifier_type:
            generate_experiment_identifiers(experiments, self.params.identifier_type)
        # Return the experiments
        return experiments

    def update_lookup(self, imageset, lookup):
        from dxtbx.format.image import ImageBool, ImageDouble

        if lookup.size is not None:
            d = imageset.get_detector()
            assert len(lookup.size) == len(d), "Incompatible size"
            for s, p in zip(lookup.size, d):
                assert s == p.get_image_size()[::-1], "Incompatible size"
            if lookup.mask.filename is not None:
                imageset.external_lookup.mask.filename = lookup.mask.filename
                imageset.external_lookup.mask.data = ImageBool(lookup.mask.data)
            if lookup.gain.filename is not None:
                imageset.external_lookup.gain.filename = lookup.gain.filename
                imageset.external_lookup.gain.data = ImageDouble(lookup.gain.data)
            if lookup.dark.filename is not None:
                imageset.external_lookup.pedestal.filename = lookup.dark.filename
                imageset.external_lookup.pedestal.data = ImageDouble(lookup.dark.data)
            if lookup.dx.filename is not None:
                imageset.external_lookup.dx.filename = lookup.dx.filename
                imageset.external_lookup.dx.data = ImageDouble(lookup.dx.data)
            if lookup.dy.filename is not None:
                imageset.external_lookup.dy.filename = lookup.dy.filename
                imageset.external_lookup.dy.data = ImageDouble(lookup.dy.data)
        return imageset

    def import_lookup_data(self, params):
        """
        Get the lookup data
        """
        # Check the lookup inputs
        mask_filename = None
        gain_filename = None
        dark_filename = None
        dx_filename = None
        dy_filename = None
        mask = None
        gain = None
        dark = None
        dx = None
        dy = None
        lookup_size = None
        if params.lookup.mask is not None:
            mask_filename = params.lookup.mask
            with open(mask_filename, "rb") as fh:
                mask = _pickle_load(fh)
            if not isinstance(mask, tuple):
                mask = (mask,)
            lookup_size = [m.all() for m in mask]
        if params.lookup.gain is not None:
            gain_filename = params.lookup.gain
            with open(gain_filename, "rb") as fh:
                gain = _pickle_load(fh)
            if not isinstance(gain, tuple):
                gain = (gain,)
            if lookup_size is None:
                lookup_size = [g.all() for g in gain]
            else:
                assert len(gain) == len(lookup_size), "Incompatible size"
                for s, g in zip(lookup_size, gain):
                    assert s == g.all(), "Incompatible size"
        if params.lookup.pedestal is not None:
            dark_filename = params.lookup.pedestal
            with open(dark_filename, "rb") as fh:
                dark = _pickle_load(fh)
            if not isinstance(dark, tuple):
                dark = (dark,)
            if lookup_size is None:
                lookup_size = [d.all() for d in dark]
            else:
                assert len(dark) == len(lookup_size), "Incompatible size"
                for s, d in zip(lookup_size, dark):
                    assert s == d.all(), "Incompatible size"
        if params.lookup.dx is not None:
            dx_filename = params.lookup.dx
            with open(dx_filename, "rb") as fh:
                dx = _pickle_load(fh)
            if not isinstance(dx, tuple):
                dx = (dx,)
            if lookup_size is None:
                lookup_size = [d.all() for d in dx]
            else:
                assert len(dx) == len(lookup_size), "Incompatible size"
                for s, d in zip(lookup_size, dx):
                    assert s == d.all(), "Incompatible size"
        if params.lookup.dy is not None:
            dy_filename = params.lookup.dy
            with open(dx_filename, "rb") as fh:
                dy = _pickle_load(fh)
            if not isinstance(dy, tuple):
                dy = (dy,)
            if lookup_size is None:
                lookup_size = [d.all() for d in dy]
            else:
                assert len(dy) == len(lookup_size), "Incompatible size"
                for s, d in zip(lookup_size, dy):
                    assert s == d.all(), "Incompatible size"
        Lookup = namedtuple("Lookup", ["size", "mask", "gain", "dark", "dx", "dy"])
        Item = namedtuple("Item", ["data", "filename"])
        return Lookup(
            size=lookup_size,
            mask=Item(data=mask, filename=mask_filename),
            gain=Item(data=gain, filename=gain_filename),
            dark=Item(data=dark, filename=dark_filename),
            dx=Item(data=dx, filename=dx_filename),
            dy=Item(data=dy, filename=dy_filename),
        )


def print_sequence_diff(sequence1, sequence2, params):
    """
    Print a diff between sequences.
    """
    logger.info(
        compare.sequence_diff(sequence1, sequence2, tolerance=params.input.tolerance)
    )


def diagnose_multiple_sequences(sequences, params):
    """
    Print a diff between sequences.
    """
    logger.info("")
    for i in range(1, len(sequences)):
        logger.info("=" * 80)
        logger.info("Diff between sequence %d and %d", i - 1, i)
        logger.info("")
        print_sequence_diff(sequences[i - 1], sequences[i], params)
    logger.info("=" * 80)
    logger.info("")


def write_experiments(experiments, params):
    """
    Output the experiments to file.
    """
    if params.output.experiments:
        logger.info("-" * 80)
        logger.info("Writing experiments to %s", params.output.experiments)
        experiments.as_file(params.output.experiments, compact=params.output.compact)


def assert_single_sequence(experiments, params):
    """
    Print an error message if more than 1 sequence
    """
    sequences = [
        e.imageset for e in experiments if isinstance(e.imageset, ImageSequence)
    ]

    if len(sequences) > 1:
        # Print some info about multiple sequences
        diagnose_multiple_sequences(sequences, params)

        # Raise exception
        raise Sorry(
            """
    More than 1 sequence was found. Two things may be happening here:

    1. There really is more than 1 sequence. If you expected this to be the
        case, set the parameter allow_multiple_sequences=True. If you don't
        expect this, then check the input to dials.import.

    2. There may be something wrong with your image headers (for example,
        the rotation ranges of each image may not match up). You should
        investigate what went wrong, but you can force dials.import to treat
        your images as a single sequence by using the template=image_####.cbf
        parameter (see help).
    """
        )


def do_import(
    args: list[str] | None = None,
    *,
    phil: libtbx.phil.scope,
    configure_logging: bool = False,
):
    # Create the option parser
    usage = "dials.import [options] /path/to/image/files"
    parser = ArgumentParser(
        usage=usage,
        sort_options=True,
        phil=phil,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    # Parse the command line arguments in two passes
    params, options = parser.parse_args(
        args=args, show_diff_phil=False, quick_parse=True
    )

    if configure_logging:
        log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Parse the command line arguments completely
    if params.input.ignore_unhandled:
        params, options, unhandled = parser.parse_args(
            args=args, show_diff_phil=False, return_unhandled=True
        )
        # Remove any False values from unhandled (eliminate empty strings)
        unhandled = [x for x in unhandled if x]
    else:
        params, options = parser.parse_args(args=args, show_diff_phil=False)
        unhandled = None

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Print a warning if something unhandled
    if unhandled:
        msg = "Unable to handle the following arguments:\n"
        msg += "\n".join(["  %s" % a for a in unhandled])
        msg += "\n"
        logger.warning(msg)

    # Print help if no input
    if len(params.input.experiments) == 0 and not (
        params.input.template or params.input.directory
    ):
        parser.print_help()
        sys.exit(0)

    # Re-extract the imagesets to rebuild experiments from
    imagesets = _extract_or_read_imagesets(params)

    metadata_updater = MetaDataUpdater(params)
    experiments = metadata_updater(imagesets)

    # If the user has requested splitting, split - these are in human numbers
    # so need to subtract 1 from start
    if params.input.split:
        if len(params.input.split) == 3:
            split_start, split_end, step = params.input.split
            split_start -= 1
        elif len(params.input.split) == 1:
            split_start = 0
            split_end = 0
            step = params.input.split[0]
        else:
            sys.exit("split=frames_per_block or split=start,end,frames_per_block")

        new_experiments = ExperimentList()
        for experiment in experiments:
            if split_start == 0 and split_end == 0:
                _split_start, _split_end = experiment.scan.get_image_range()
                _split_start -= 1
            else:
                _split_start = split_start
                _split_end = split_end

            for chunk_start in range(_split_start, _split_end, step):
                end = chunk_start + step
                if end > _split_end:
                    end = _split_end

                tmp = ExperimentListFactory.from_imageset_and_crystal(
                    experiment.imageset, crystal=None
                )[0]
                tmp.scan = experiment.scan[chunk_start:end]
                tmp.imageset = experiment.imageset[chunk_start:end]
                new_experiments.append(tmp)
        experiments = new_experiments

    # Compute some numbers
    num_sweeps = 0
    num_still_sequences = 0
    num_stills = 0
    num_images = 0

    # importing a lot of experiments all pointing at one imageset should
    # work gracefully
    counted_imagesets = []

    for e in experiments:
        if e.imageset in counted_imagesets:
            continue
        if isinstance(e.imageset, ImageSequence):
            if e.imageset.get_scan().is_still():
                num_still_sequences += 1
            else:
                num_sweeps += 1
        else:
            num_stills += 1
        num_images += len(e.imageset)
        counted_imagesets.append(e.imageset)

    unique_formats = OrderedSet()
    unique_templates = OrderedSet()
    for imgset in counted_imagesets:
        unique_formats.add(imgset.get_format_class())
        if scan := imgset.get_scan():
            start, end = scan.get_image_range()
            unique_templates.add(f"{imgset.get_template()}:{start}:{end}")
        else:
            paths = imgset.reader().paths()
            if len(paths) == 1:
                unique_templates.add(paths[0])
            else:
                template, _ = template_regex_from_list(paths)
                unique_templates.add(template)

    # Print out some bulk info
    logger.info("-" * 80)
    for fmt in unique_formats:
        logger.info("  format: %s", fmt)
    for template in unique_templates:
        logger.info("  template: %s", template)
    logger.info("  num images: %d", num_images)
    logger.info("  sequences:")
    logger.info("    still:    %d", num_still_sequences)
    logger.info("    sweep:    %d", num_sweeps)
    logger.info("  num stills: %d", num_stills)

    # Print out info for all experiments
    for experiment in experiments:
        # Print some experiment info - override the output of image range
        # if appropriate
        image_range = params.geometry.scan.image_range
        if isinstance(experiment.imageset, ImageSequence):
            imageset_type = "sequence"
        else:
            imageset_type = "stills"

        logger.debug("-" * 80)
        logger.debug("  format: %s", str(experiment.imageset.get_format_class()))
        logger.debug("  imageset type: %s", imageset_type)
        if image_range is None:
            logger.debug("  num images:    %d", len(experiment.imageset))
        else:
            logger.debug("  num images:    %d", image_range[1] - image_range[0] + 1)

        logger.debug("")
        logger.debug(experiment.imageset.get_beam())
        logger.debug(experiment.imageset.get_goniometer())
        logger.debug(experiment.imageset.get_detector())
        logger.debug(experiment.imageset.get_scan())

    # Only allow a single sequence
    if params.input.allow_multiple_sequences is False:
        assert_single_sequence(experiments, params)

    # Write the experiments to file
    write_experiments(experiments, params)
    return experiments


@show_mail_handle_errors()
def run(args=None, *, phil=phil_scope):
    do_import(args, phil=phil, configure_logging=True)


if __name__ == "__main__":
    run()
