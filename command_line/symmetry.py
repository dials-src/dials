from __future__ import absolute_import, division, print_function

import copy
import json
import logging
import random
import sys
from dials.util import tabulate

from cctbx import sgtbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from libtbx import Auto
import iotbx.phil
from rstbx.symmetry.constraints import parameter_reduction

from dials.array_family import flex
from dials.util import log, show_mail_on_error
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from dials.util.version import dials_version
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.algorithms.symmetry import resolution_filter_from_reflections_experiments
from dials.algorithms.symmetry.laue_group import LaueGroupAnalysis
from dials.algorithms.merging.merge import prepare_merged_reflection_table
from dials.algorithms.symmetry.absences.screw_axes import ScrewAxisObserver
from dials.algorithms.symmetry.absences.run_absences_checks import (
    run_systematic_absences_checks,
)
from dials.algorithms.symmetry.absences.laue_groups_info import (
    laue_groups as laue_groups_for_absence_analysis,
)

logger = logging.getLogger("dials.command_line.symmetry")

phil_scope = iotbx.phil.parse(
    """\
d_min = Auto
  .type = float(value_min=0)

min_i_mean_over_sigma_mean = 4
  .type = float(value_min=0)

min_cc_half = 0.6
  .type = float(value_min=0, value_max=1)

batch = None
  .type = ints(value_min=0, size=2)
  .help = "Limit batch range for analysis: manually apply results afterwards"

normalisation = kernel quasi ml_iso *ml_aniso
  .type = choice

lattice_group = None
  .type = space_group

seed = 230
  .type = int(value_min=0)

lattice_symmetry_max_delta = 2.0
  .type = float(value_min=0)

relative_length_tolerance = 0.05
  .type = float(value_min=0)

absolute_angle_tolerance = 2
  .type = float(value_min=0)

partiality_threshold = 0.99
  .type = float
  .help = "Use only reflections with a partiality above this threshold."

laue_group = auto
  .type = space_group
  .help = "Optionally specify the Laue group. If set to auto, then test all possible "
          "Laue groups. If set to None, then take the Laue group from the input file."

change_of_basis_op = None
  .type = str

systematic_absences {

  check = True
    .type = bool
    .help = "Check systematic absences for the current laue group."

  significance_level = *0.95 0.975 0.99
    .type = choice
    .help = "Signficance to use when testing whether axial reflections are "
            "different to zero (absences and reflections in reflecting condition)."

}

output {
  log = dials.symmetry.log
    .type = str
  experiments = "symmetrized.expt"
    .type = path
  reflections = "symmetrized.refl"
    .type = path
  json = dials.symmetry.json
    .type = path
  html = "dials-symmetry.html"
    .type = path
    .help = "Filename for html report."
}
""",
    process_includes=True,
)


def map_to_minimum_cell(experiments, reflections, max_delta):
    """
    Map experiments and reflections to the minimum cell

    Map to the minimum cell via the best cell, which appears to guarantee that the
    resulting minimum cells are consistent.

    Args:
        experiments (ExperimentList): a list of experiments.
        reflections (list): a list of reflection tables

    Returns: The experiments and reflections mapped to the minimum cell
    """
    cb_ops = []
    for expt, refl in zip(experiments, reflections):
        groups = metric_subgroups(
            expt.crystal.get_crystal_symmetry(),
            max_delta,
            enforce_max_delta_for_generated_two_folds=True,
        )
        group = groups.result_groups[0]
        cb_op_best_to_min = group["best_subsym"].change_of_basis_op_to_minimum_cell()
        cb_op_inp_min = cb_op_best_to_min * group["cb_op_inp_best"]
        refl["miller_index"] = cb_op_inp_min.apply(refl["miller_index"])
        expt.crystal = expt.crystal.change_basis(cb_op_inp_min)
        expt.crystal.set_space_group(sgtbx.space_group())
        cb_ops.append(cb_op_inp_min)
    return experiments, reflections, cb_ops


def symmetry(experiments, reflection_tables, params=None):
    """
    Run symmetry analysis

    Args:
        experiments: An experiment list.
        reflection_tables: A list of reflection tables.
        params: The dials.symmetry phil scope.
    """
    result = None
    if params is None:
        params = phil_scope.extract()

    if params.laue_group is Auto:
        logger.info("=" * 80)
        logger.info("")
        logger.info("Performing Laue group analysis")
        logger.info("")

        # transform models into miller arrays
        n_datasets = len(experiments)

        experiments, reflection_tables, cb_ops = map_to_minimum_cell(
            experiments, reflection_tables, params.lattice_symmetry_max_delta
        )

        datasets = filtered_arrays_from_experiments_reflections(
            experiments,
            reflection_tables,
            outlier_rejection_after_filter=True,
            partiality_threshold=params.partiality_threshold,
        )
        if len(datasets) != n_datasets:
            raise ValueError(
                """Some datasets have no reflection after prefiltering, please check
    input data and filtering settings e.g partiality_threshold"""
            )

        result = LaueGroupAnalysis(
            datasets,
            normalisation=params.normalisation,
            d_min=params.d_min,
            min_i_mean_over_sigma_mean=params.min_i_mean_over_sigma_mean,
            lattice_symmetry_max_delta=params.lattice_symmetry_max_delta,
            relative_length_tolerance=params.relative_length_tolerance,
            absolute_angle_tolerance=params.absolute_angle_tolerance,
        )
        logger.info("")
        logger.info(result)

        if params.output.json is not None:
            d = result.as_dict()
            d["cb_op_inp_min"] = [str(cb_op) for cb_op in cb_ops]
            # This is not the input symmetry as we have already mapped it to minimum
            # cell, so delete from the output dictionary to avoid confusion
            del d["input_symmetry"]
            json_str = json.dumps(d, indent=2)
            with open(params.output.json, "w") as f:
                f.write(json_str)

        # Change of basis operator from input unit cell to best unit cell
        cb_op_inp_best = result.best_solution.subgroup["cb_op_inp_best"]
        # Get the best space group.
        best_subsym = result.best_solution.subgroup["best_subsym"]
        best_space_group = best_subsym.space_group().build_derived_acentric_group()
        logger.info(
            tabulate(
                [[str(best_subsym.space_group_info()), str(best_space_group.info())]],
                ["Patterson group", "Corresponding MX group"],
            )
        )
        # Reindex the input data
        experiments, reflection_tables = _reindex_experiments_reflections(
            experiments, reflection_tables, best_space_group, cb_op_inp_best
        )

    elif params.laue_group is not None:
        if params.change_of_basis_op is not None:
            cb_op = sgtbx.change_of_basis_op(params.change_of_basis_op)
        else:
            cb_op = sgtbx.change_of_basis_op()
        # Reindex the input data
        experiments, reflection_tables = _reindex_experiments_reflections(
            experiments, reflection_tables, params.laue_group.group(), cb_op
        )

    if params.systematic_absences.check:
        logger.info("=" * 80)
        logger.info("")
        logger.info("Analysing systematic absences")
        logger.info("")

        # Get the laue class from the current space group.
        space_group = experiments[0].crystal.get_space_group()
        laue_group = str(space_group.build_derived_patterson_group().info())
        logger.info("Laue group: %s", laue_group)
        if laue_group not in laue_groups_for_absence_analysis:
            logger.info("No absences to check for this laue group\n")
        else:
            if (params.d_min is Auto) and (result is not None):
                d_min = result.intensities.resolution_range()[1]
            elif params.d_min is Auto:
                d_min = resolution_filter_from_reflections_experiments(
                    reflection_tables,
                    experiments,
                    params.min_i_mean_over_sigma_mean,
                    params.min_cc_half,
                )
            else:
                d_min = params.d_min

            # combine before sys abs test - only triggers if laue_group=None and
            # multiple input files.
            if len(reflection_tables) > 1:
                joint_reflections = flex.reflection_table()
                for table in reflection_tables:
                    joint_reflections.extend(table)
            else:
                joint_reflections = reflection_tables[0]

            merged_reflections = prepare_merged_reflection_table(
                experiments, joint_reflections, d_min
            )
            run_systematic_absences_checks(
                experiments,
                merged_reflections,
                float(params.systematic_absences.significance_level),
            )

    logger.info(
        "Saving reindexed experiments to %s in space group %s",
        params.output.experiments,
        str(experiments[0].crystal.get_space_group().info()),
    )
    experiments.as_file(params.output.experiments)
    if params.output.reflections is not None:
        if len(reflection_tables) > 1:
            joint_reflections = flex.reflection_table()
            for table in reflection_tables:
                joint_reflections.extend(table)
        else:
            joint_reflections = reflection_tables[0]
        logger.info(
            "Saving %s reindexed reflections to %s",
            len(joint_reflections),
            params.output.reflections,
        )
        joint_reflections.as_file(params.output.reflections)

    if params.output.html and params.systematic_absences.check:
        ScrewAxisObserver().generate_html_report(params.output.html)


def _reindex_experiments_reflections(experiments, reflections, space_group, cb_op):
    """Reindex the input data."""
    reindexed_experiments = copy.deepcopy(experiments)
    reindexed_reflections = flex.reflection_table()
    for i, expt in enumerate(reindexed_experiments):
        # Set the space group to the best symmetry and change basis accordingly.
        # Setting the basis to one incompatible with the initial space group is
        # forbidden, so we must first change the space group to P1 to be safe`.
        expt.crystal.set_space_group(sgtbx.space_group("P 1"))
        expt.crystal = expt.crystal.change_basis(cb_op)
        expt.crystal.set_space_group(space_group)

        S = parameter_reduction.symmetrize_reduce_enlarge(
            expt.crystal.get_space_group()
        )
        S.set_orientation(expt.crystal.get_B())
        S.symmetrize()
        expt.crystal.set_B(S.orientation.reciprocal_matrix())
        reindexed_refl = copy.deepcopy(reflections[i])
        reindexed_refl["miller_index"] = cb_op.apply(reindexed_refl["miller_index"])
        reindexed_reflections.extend(reindexed_refl)

    return reindexed_experiments, [reindexed_reflections]


help_message = """
This program implements the methods of
`POINTLESS <http://www.ccp4.ac.uk/html/pointless.html>`_ (
`Evans, P. (2006). Acta Cryst. D62, 72-82. <https://doi.org/10.1107/S0907444905036693>`_ and
`Evans, P. R. (2011). Acta Cryst. D67, 282-292. <https://doi.org/10.1107/S090744491003982X>`_)
for scoring and determination of Laue group symmetry.

The program takes as input a set of one or more integrated experiments and
reflections.

Examples::

  dials.symmetry models.expt observations.refl
"""


def run(args=None):
    """Run symmetry analysis from the command-line."""
    usage = "dials.symmetry [options] models.expt observations.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.seed is not None:
        flex.set_random_seed(params.seed)
        random.seed(params.seed)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    reflections = parse_multiple_datasets(reflections)

    # Cut down reflection lists according to input batch range if set
    if params.batch is not None:
        z0, z1 = map(float, params.batch)
        logger.info("Cutting reflection lists to batch range %d to %d" % (z0, z1))
        trimmed_reflections = []
        for refl in reflections:
            z = refl["xyzcal.px"].parts()[2]
            keep = (z >= z0) & (z <= z1)
            trimmed_reflections.append(refl.select(keep))
        reflections = trimmed_reflections

    if len(experiments) != len(reflections):
        sys.exit(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        symmetry(experiments, reflections, params=params)
    except ValueError as e:
        sys.exit(e)


if __name__ == "__main__":
    with show_mail_on_error():
        run()
