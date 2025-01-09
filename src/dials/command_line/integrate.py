# DIALS_ENABLE_COMMAND_LINE_COMPLETION
"""
This program is used to integrate the reflections on the diffraction images. It
is called with an experiment list outputted from dials.index or dials.refine and
a corresponding set of strong spots from which a profile model is calculated.
The program will output a set of integrated reflections and an experiment list
with additional profile model data. The data can be reintegrated using the same
profile model by inputting this integrated.expt file back into
dials.integrate.

Examples::

  dials.integrate models.expt refined.refl

  dials.integrate models.expt refined.refl output.reflections=integrated.refl

  dials.integrate models.expt refined.refl profile.fitting=False

  dials.integrate models.expt refined.refl background.algorithm=glm
"""

from __future__ import annotations

import logging
import math
import sys
import warnings

from ordered_set import OrderedSet

from dxtbx.model.experiment_list import Experiment, ExperimentList
from libtbx.phil import parse

import dials.util.log
from dials.algorithms.integration.integrator import create_integrator
from dials.algorithms.profile_model.factory import ProfileModelFactory
from dials.array_family import flex
from dials.util import show_mail_handle_errors
from dials.util.command_line import heading
from dials.util.exclude_images import expand_exclude_multiples, set_invalid_images
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.slice import slice_crystal
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.integrate")

# Create the phil scope

phil_scope = parse(
    """

  output {
    experiments = 'integrated.expt'
      .type = str
      .help = "The experiments output filename"

    output_unintegrated_reflections = False
      .type = bool
      .expert_level = 2
      .help = "Include unintegrated reflections in output file"

    reflections = 'integrated.refl'
      .type = str
      .help = "The integrated output filename"

    phil = 'dials.integrate.phil'
      .type = str
      .help = "The output phil file"

    log = 'dials.integrate.log'
      .type = str
      .help = "The log filename"

    report = None
      .type = str
      .help = "The integration report filename (*.xml or *.json)"

    include_bad_reference = False
      .type = bool
      .help = "Include bad reference data including unindexed spots,"
              "and reflections whose predictions are messed up in the"
              "reflection table output. Reflections will have the"
              "'bad_reference' flag set."
  }

  scan_range = None
    .type = ints(size=2)
    .help = "Explicitly specify the images to be processed. Only applicable"
            "when experiment list contains a single imageset."
    .multiple = True

  create_profile_model = True
    .type = bool
    .help = "Create the profile model"

  sampling
    .expert_level = 1
  {

    reflections_per_degree = 50
      .help = "The number of predicted reflections per degree of the sequence "
              "to integrate."
      .type = float(value_min=0.)

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input "
              "prediction list is done"
      .type = int

    maximum_sample_size = None
      .help = "The maximum number of predictions to integrate."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)

    integrate_all_reflections = True
      .help = "Override reflections_per_degree and integrate all predicted"
              "reflections."
      .type = bool

    random_seed = 0
      .help = "Random seed for sampling"
      .type = int

  }

  include scope dials.util.exclude_images.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope
  include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope
""",
    process_includes=True,
)

# Local overrides for dials.integrate
phil_overrides = parse(
    """
integration {
  mp {
    nproc = Auto
  }
}
"""
)
working_phil = phil_scope.fetch(sources=[phil_overrides])


def process_reference(reference):
    """
    Remove bad reflections from the reference.

    Remove unindexed, bad_for_refinement, bad miller index.

    Args:
        reference: A reflection table.

    Returns:
        (tuple): tuple containing:

            reference: A reduction of the input reference reflection table.
            rubbish: A reflection table containing the reflections filtered out of
                the input table.

    Raises:
        ValueError: If no indexed spots, bad id, unmatched panel.
    """

    if reference is None:
        return None, None
    assert "miller_index" in reference
    assert "id" in reference
    logger.info(
        "Processing reference reflections\n read %d strong spots", reference.size()
    )
    mask = reference.get_flags(reference.flags.indexed)
    rubbish = reference.select(~mask)
    n_unindexed = mask.count(False)
    if n_unindexed > 0:
        reference.del_selected(~mask)
        logger.info(" removing %d unindexed reflections", n_unindexed)
    if reference.size() == 0:
        raise ValueError(
            "Invalid input for reference reflections. No indexed spots found."
        )
    mask = reference.get_flags(reference.flags.bad_for_refinement, all=False)
    n_masked = mask.count(True)
    if n_masked:
        rubbish.extend(reference.select(mask))
        reference.del_selected(mask)
        logger.info(" removing %d reflections marked as bad for refinement", n_masked)
    mask = reference["miller_index"] == (0, 0, 0)
    n_masked = mask.count(True)
    if n_masked > 0:
        rubbish.extend(reference.select(mask))
        reference.del_selected(mask)
        logger.info(" removing %d reflections with hkl (0,0,0)", n_masked)
    mask = reference["id"] < 0
    n_masked = mask.count(True)
    if n_masked > 0:
        raise ValueError(
            """
    Invalid input for reference reflections.
    %d reference spots have an invalid experiment id
    """
            % n_masked
        )
    if (reference["panel"] == reference["shoebox"].panels()).count(False) > 0:
        raise ValueError(
            'reflection table "panel" column does not match "shoebox" panel'
        )
    logger.info(" using %d indexed reflections", reference.size())
    logger.info(" found %d junk reflections", rubbish.size())
    return reference, rubbish


def filter_reference_pixels(reference, experiments):
    """
    Set any pixel closer to other reflections to background.

    Args:
        reference: A reflection table
        experiments: The experiment list

    Returns:
        The input reflection table with modified shoeboxes.
    """
    modified_count = 0
    for experiment, indices in reference.iterate_experiments_and_indices(experiments):
        subset = reference.select(indices)
        modified = subset["shoebox"].mask_neighbouring(
            subset["miller_index"],
            experiment.beam,
            experiment.detector,
            experiment.goniometer,
            experiment.scan,
            experiment.crystal,
        )
        modified_count += modified.count(True)
        reference.set_selected(indices, subset)
    logger.info(" masked neighbouring pixels in %d shoeboxes", modified_count)
    return reference


def sample_predictions(experiments, predicted, params):
    """
    Select a random sample of the predicted reflections to integrate.

    Args:
        experiments: The experiment list
        predicted: A reflection table of predicted reflections
        params: The integration phil parameters

    Returns:
        A subset of the original predicted table.
    """

    if params.sampling.random_seed:
        flex.set_random_seed(params.sampling.random_seed)

    nref_per_degree = params.sampling.reflections_per_degree
    min_sample_size = params.sampling.minimum_sample_size
    max_sample_size = params.sampling.maximum_sample_size

    # this code is very similar to David's code in algorithms/refinement/reflection_manager.py!

    working_isel = flex.size_t()
    for iexp, exp in enumerate(experiments):
        sel = predicted["id"] == iexp
        isel = sel.iselection()
        nrefs = sample_size = len(isel)

        # set sample size according to nref_per_degree (per experiment)
        if exp.scan and nref_per_degree:
            sequence_range_rad = exp.scan.get_oscillation_range(deg=False)
            width = math.degrees(abs(sequence_range_rad[1] - sequence_range_rad[0]))
            sample_size = int(nref_per_degree * width)
        else:
            sequence_range_rad = None

        # adjust sample size if below the chosen limit
        sample_size = max(sample_size, min_sample_size)

        # set maximum sample size if requested
        if max_sample_size:
            sample_size = min(sample_size, max_sample_size)

        # determine subset and collect indices
        if sample_size < nrefs:
            isel = isel.select(flex.random_selection(nrefs, sample_size))
        working_isel.extend(isel)

    # create subset
    return predicted.select(working_isel)


def split_for_scan_range(experiments, reference, scan_range):
    """Update experiments when scan range is set.

    Args:
        experiments: An experiment list
        reference: A reflection table of reference reflections
        scan_range (tuple): Range of scan images to be processed

    Returns:
        experiments: A new experiment list with the requested scan ranges
        reference: A reflection table with data from the scan ranges

    Raises:
        ValueError: If bad input for scan range.
    """

    # Only do anything is the scan range is set
    if scan_range is not None and len(scan_range) > 0:
        # Ensure that all experiments have the same imageset and scan
        iset = [e.imageset for e in experiments]
        scan = [e.scan for e in experiments]
        assert all(x == iset[0] for x in iset)
        assert all(x == scan[0] for x in scan)

        # Get the imageset and scan
        iset = experiments[0].imageset
        scan = experiments[0].scan

        # Get the array range
        if scan is not None:
            frames_start, frames_end = scan.get_array_range()
            assert scan.get_num_images() == len(iset)
        else:
            frames_start, frames_end = (0, len(iset))

        # Create the new lists
        new_experiments = ExperimentList()
        new_reference_all = reference.split_by_experiment_id()
        new_reference = flex.reflection_table()
        for i in range(len(new_reference_all) - len(experiments)):
            new_reference_all.append(flex.reflection_table())
        assert len(new_reference_all) == len(experiments)

        # Loop through all the scan ranges and create a new experiment list with
        # the requested scan ranges.
        for scan_start, scan_end in scan_range:
            # Validate the requested scan range
            if scan_end == scan_start:
                raise ValueError(
                    f"Scan range end must be higher than start; pass {scan_start},{scan_start + 1} for single image"
                )
            if scan_end < scan_start:
                raise ValueError("Scan range must be in ascending order")
            elif scan_start < frames_start or scan_end > frames_end:
                raise ValueError(
                    f"Scan range must be within image range {frames_start}..{frames_end}"
                )

            assert scan_end > scan_start
            assert scan_start >= frames_start
            assert scan_end <= frames_end

            index_start = scan_start - frames_start
            index_end = index_start + (scan_end - scan_start)
            assert index_start < index_end
            assert index_start >= 0
            assert index_end <= len(iset)
            new_iset = iset[index_start:index_end]
            if scan is None:
                new_scan = None
            else:
                new_scan = scan[index_start:index_end]

            for i, e1 in enumerate(experiments):
                e2 = Experiment()
                e2.beam = e1.beam
                e2.detector = e1.detector
                e2.goniometer = e1.goniometer
                e2.crystal = slice_crystal(e1.crystal, (index_start, index_end))
                e2.profile = e1.profile
                e2.imageset = new_iset
                e2.scan = new_scan
                new_reference_all[i]["id"] = flex.int(
                    len(new_reference_all[i]), len(new_experiments)
                )
                new_reference.extend(new_reference_all[i])
                new_experiments.append(e2)
        experiments = new_experiments
        reference = new_reference

        # Print some information
        logger.info("Modified experiment list to integrate over requested scan range")
        for scan_start, scan_end in scan_range:
            logger.info(" scan_range = %d -> %d", scan_start, scan_end)

    # Return the experiments
    return experiments, reference


def run_integration(params, experiments, reference=None):
    """Perform the integration.

    Returns:
        experiments: The integrated experiments
        reflections: The integrated reflections
        report(optional): An integration report.

    Raises:
        ValueError: For a number of bad inputs
        RuntimeError: If the profile model creation fails
    """
    predicted = None
    rubbish = None

    for abs_params in params.absorption_correction:
        if abs_params.apply:
            if not (
                params.integration.debug.output
                and not params.integration.debug.separate_files
            ):
                raise ValueError(
                    "Shoeboxes must be saved to integration intermediates to apply an absorption correction. "
                    + "Set integration.debug.output=True, integration.debug.separate_files=False and "
                    + "integration.debug.delete_shoeboxes=True to temporarily store shoeboxes."
                )

    # Print if we're using a mask
    for i, exp in enumerate(experiments):
        mask = exp.imageset.external_lookup.mask
        if mask.filename is not None:
            if mask.data:
                logger.info("Using external mask: %s", mask.filename)
                for tile in mask.data:
                    logger.info(" Mask has %d pixels masked", tile.data().count(False))

    # Print the experimental models
    for i, exp in enumerate(experiments):
        summary = "\n".join(
            (
                "",
                "=" * 80,
                "",
                "Experiments",
                "",
                "Models for experiment %d" % i,
                "",
                str(exp.beam),
                str(exp.detector),
            )
        )
        if exp.goniometer:
            summary += str(exp.goniometer) + "\n"
        if exp.scan:
            summary += str(exp.scan) + "\n"
        summary += str(exp.crystal)
        logger.info(summary)

    logger.info("\n".join(("", "=" * 80, "")))
    logger.info(heading("Initialising"))

    # Load the data
    if reference:
        reference, rubbish = process_reference(reference)

        # Check pixels don't belong to neighbours
        if exp.goniometer is not None and exp.scan is not None:
            reference = filter_reference_pixels(reference, experiments)

        # Modify experiment list if scan range is set.
        experiments, reference = split_for_scan_range(
            experiments, reference, params.scan_range
        )

    # Modify experiment list if exclude_images is set
    if params.exclude_images_multiple:
        params.exclude_images = expand_exclude_multiples(
            experiments,
            params.exclude_images_multiple,
            params.exclude_images,
        )
    if params.exclude_images:
        try:
            experiments = set_invalid_images(experiments, params.exclude_images)
        except ValueError as err:
            # Handle deprecated way of providing exclude_images
            try:
                exclude_images = [
                    int(e)
                    for e in str(params.exclude_images)
                    .replace("[", "")
                    .replace("]", "")
                    .replace("'", "")
                    .split(",")
                ]
            except ValueError:
                raise (err)
            warnings.warn(
                "Providing exclude_images as a single list (e.g. 1,2,3,4,5 etc.) is deprecated.\n"
                + str(err),
                DeprecationWarning,
                stacklevel=2,
            )
            for experiment in experiments:
                for index in exclude_images:
                    experiment.imageset.mark_for_rejection(index, True)

    # Predict the reflections
    logger.info("\n".join(("", "=" * 80, "")))
    logger.info(heading("Predicting reflections"))
    predicted = flex.reflection_table.from_predictions_multi(
        experiments,
        dmin=params.prediction.d_min,
        dmax=params.prediction.d_max,
        margin=params.prediction.margin,
        force_static=params.prediction.force_static,
        padding=params.prediction.padding,
    )
    isets = OrderedSet(e.imageset for e in experiments)
    predicted["imageset_id"] = flex.int(predicted.size(), 0)
    if len(isets) > 1:
        for e in experiments:
            iset_id = isets.index(e.imageset)
            for id_ in predicted.experiment_identifiers().keys():
                identifier = predicted.experiment_identifiers()[id_]
                if identifier == e.identifier:
                    sel = predicted["id"] == id_
                    predicted["imageset_id"].set_selected(sel, iset_id)
                    break

    # Match reference with predicted
    if reference:
        matched, reference, unmatched = predicted.match_with_reference(reference)
        assert len(matched) == len(predicted)
        assert matched.count(True) <= len(reference)
        if matched.count(True) == 0:
            raise ValueError(
                """
        Invalid input for reference reflections.
        Zero reference spots were matched to predictions
    """
            )
        elif unmatched:
            msg = (
                "Warning: %d reference spots were not matched to predictions"
                % unmatched.size()
            )
            border = "\n".join(("", "*" * 80, ""))
            logger.info("".join((border, msg, border)))
            rubbish.extend(unmatched)

        if len(experiments) > 1:
            # filter out any experiments without matched reference reflections
            # f_: filtered

            f_reference = flex.reflection_table()
            f_predicted = flex.reflection_table()
            f_rubbish = flex.reflection_table()
            f_experiments = ExperimentList()
            good_expt_count = 0

            def refl_extend(src, dest, eid):
                old_id = eid
                new_id = good_expt_count
                tmp = src.select(src["id"] == old_id)
                tmp["id"] = flex.int(len(tmp), good_expt_count)
                if old_id in tmp.experiment_identifiers():
                    identifier = tmp.experiment_identifiers()[old_id]
                    del tmp.experiment_identifiers()[old_id]
                    tmp.experiment_identifiers()[new_id] = identifier
                dest.extend(tmp)

            for expt_id, experiment in enumerate(experiments):
                if len(reference.select(reference["id"] == expt_id)) != 0:
                    refl_extend(reference, f_reference, expt_id)
                    refl_extend(predicted, f_predicted, expt_id)
                    refl_extend(rubbish, f_rubbish, expt_id)
                    f_experiments.append(experiment)
                    good_expt_count += 1
                else:
                    logger.info(
                        "Removing experiment %d: no reference reflections matched to predictions",
                        expt_id,
                    )

            reference = f_reference
            predicted = f_predicted
            experiments = f_experiments
            rubbish = f_rubbish

    # Select a random sample of the predicted reflections
    if not params.sampling.integrate_all_reflections:
        predicted = sample_predictions(experiments, predicted, params)

    # Compute the profile model - either load existing or compute
    # can raise RuntimeError
    experiments = ProfileModelFactory.create(params, experiments, reference)
    for expr in experiments:
        if expr.profile is None:
            raise ValueError("No profile information in experiment list")
    del reference

    # Compute the bounding box
    predicted.compute_bbox(experiments)

    # Create the integrator
    integrator = create_integrator(params, experiments, predicted)

    # Integrate the reflections
    reflections = integrator.integrate()

    # Remove unintegrated reflections
    if not params.output.output_unintegrated_reflections:
        keep = reflections.get_flags(reflections.flags.integrated, all=False)
        logger.info(
            "Removing %d unintegrated reflections of %d total",
            keep.count(False),
            keep.size(),
        )

        reflections = reflections.select(keep)

    # Append rubbish data onto the end
    if rubbish is not None and params.output.include_bad_reference:
        mask = flex.bool(len(rubbish), True)
        rubbish.unset_flags(mask, rubbish.flags.integrated_sum)
        rubbish.unset_flags(mask, rubbish.flags.integrated_prf)
        rubbish.set_flags(mask, rubbish.flags.bad_reference)
        reflections.extend(rubbish)

    # Correct integrated intensities for absorption correction, if necessary
    for abs_params in params.absorption_correction:
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
            experiments, reflections = multi_kapton_correction(
                experiments, reflections, abs_params.fuller_kapton, logger=logger
            )()

    if params.significance_filter.enable:
        from dials.algorithms.integration.stills_significance_filter import (
            SignificanceFilter,
        )

        sig_filter = SignificanceFilter(params)
        filtered_refls = sig_filter(experiments, reflections)
        accepted_expts = ExperimentList()
        accepted_refls = flex.reflection_table()
        logger.info(
            "Removed %d reflections out of %d when applying significance filter",
            (reflections.size() - filtered_refls.size()),
            reflections.size(),
        )
        for expt_id, expt in enumerate(experiments):
            refls = filtered_refls.select(filtered_refls["id"] == expt_id)
            if refls:
                accepted_expts.append(expt)
                current_id = expt_id
                new_id = len(accepted_expts) - 1
                refls["id"] = flex.int(len(refls), new_id)
                if expt.identifier:
                    del refls.experiment_identifiers()[current_id]
                    refls.experiment_identifiers()[new_id] = expt.identifier
                accepted_refls.extend(refls)
            else:
                logger.info(
                    "Removed experiment %d which has no reflections left after applying significance filter",
                    expt_id,
                )

        if not accepted_refls:
            raise ValueError("No reflections left after applying significance filter")
        experiments = accepted_expts
        reflections = accepted_refls

    # Write a report if requested
    report = None
    if params.output.report is not None:
        report = integrator.report()

    return experiments, reflections, report


@show_mail_handle_errors()
def run(args=None, phil=working_phil):
    """Run the integration command line script."""
    usage = "usage: dials.integrate [options] models.expt"

    # Create the parser
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging
    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the PHIL diff
    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)
    # Save phil parameters
    if params.output.phil is not None:
        with open(params.output.phil, "w") as outfile:
            outfile.write(parser.diff_phil.as_str())

    reference, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if not reference and not experiments:
        parser.print_help()
        return
    if not experiments:
        sys.exit("No experiment list was specified")
    if not reference:
        reference = None
    elif len(reference) != 1:
        sys.exit("More than 1 reflection file was given")
    else:
        reference = reference[0]

    if reference and "shoebox" not in reference:
        sys.exit("Error: shoebox data missing from reflection table")

    try:
        experiments, reflections, report = run_integration(
            params, experiments, reference
        )
    except (ValueError, RuntimeError) as e:
        sys.exit(e)
    else:
        # Delete the shoeboxes used for intermediate calculations, if requested
        if params.integration.debug.delete_shoeboxes and "shoebox" in reflections:
            del reflections["shoebox"]

        logger.info(
            "Saving %d reflections to %s", reflections.size(), params.output.reflections
        )
        reflections.as_file(params.output.reflections)
        logger.info("Saving the experiments to %s", params.output.experiments)
        experiments.as_file(params.output.experiments)

        if report:
            report.as_file(params.output.report)


if __name__ == "__main__":
    run()
