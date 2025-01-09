#!/usr/bin/env python
#
# dials.ssx_integrate.py
#
#  Copyright (C) 2022 Diamond Light Source
#
#  Author: James Beilsten-Edmands
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
"""
This program rums profile modelling and integration on indexed results from a
still sequence i.e. SSX data. This script uses parts of the regular DIALS
integration code, using either the ellipsoid or stills integrator algorithms.

The ellipsoid algorithm refines the unit cell, orientation and a 3D ellipsoidal
mosaicity parameterisation for each crystal, by assessing the pixel-intensity
distribution of the strong spots. The integrated data are saved in batches to
hep with memory management. A html output report is generated, showing integration
and clutering statistics.

Further program documentation can be found at dials.github.io/ssx_processing_guide.html

Usage:
    dials.ssx_integrate indexed.expt indexed.refl
    dials.ssx_integrate refined.expt refined.refl
    dials.ssx_integrate indexed.expt indexed.refl algorithm=stills
"""

from __future__ import annotations

import copy
import json
import logging
import pathlib
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Any

import iotbx.phil
from cctbx import crystal
from dxtbx.model import Experiment, ExperimentList
from libtbx import Auto
from libtbx.utils import Sorry

from dials.algorithms.indexing.ssx.analysis import report_on_crystal_clusters
from dials.algorithms.indexing.ssx.processing import manage_loggers
from dials.algorithms.integration.ssx.ellipsoid_integrate import (
    EllipsoidIntegrator,
    EllipsoidOutputAggregator,
)
from dials.algorithms.integration.ssx.ssx_integrate import (
    OutputAggregator,
    SimpleIntegrator,
    generate_html_report,
)
from dials.algorithms.integration.ssx.stills_integrate import StillsIntegrator
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.combine_experiments import CombineWithReference
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections
from dials.util.system import CPU_COUNT
from dials.util.version import dials_version

logger = logging.getLogger("dials.ssx_integrate")

# Create the phil scope
phil_scope = iotbx.phil.parse(
    """
  algorithm = *ellipsoid stills
    .type = choice
  nproc=Auto
    .type = int
  image_range = None
    .type = str
    .help = "Only process this image range. Input in the format start:end."
            "This range is an inclusive range."
  output {
    batch_size = 50
      .type = int
      .help = "Number of images to save in each output file"
    log = "dials.ssx_integrate.log"
      .type = str
    html = "dials.ssx_integrate.html"
      .type = str
    json = None
      .type = str
    nuggets = None
      .type = path
      .help = "Specify a directory to which a per-image summary json will be saved"
              "during processing, as each image is integrated, to enable live monitoring."
    history = None
      .type = str
      .help = "Output refinement history to json"
  }

  ellipsoid {
    include scope dials.algorithms.profile_model.ellipsoid.algorithm.ellipsoid_algorithm_phil_scope
  }

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope
  include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope


  debug {
    output {
      shoeboxes = False
        .type = bool
    }
  }

""",
    process_includes=True,
)

phil_overrides = phil_scope.fetch(
    source=iotbx.phil.parse(
        """\
profile {
    gaussian_rs {
        min_spots {
            overall=0
        }
    }
    fitting = False
    ellipsoid {
        refinement {
            n_cycles = 1
        }
    }
}
integration {
    background {
        simple {
            outlier {
                algorithm = Null
            }
        }
    }
}
"""
    )
)
working_phil = phil_scope.fetch(sources=[phil_overrides])

working_phil.adopt_scope(
    iotbx.phil.parse(
        """
    individual_log_verbosity = 1
    .type =int
"""
    )
)

loggers_to_disable = ["dials", "dials.array_family.flex_ext"]

loggers_to_disable_for_stills = loggers_to_disable + [
    "dials.algorithms.integration.integrator",
    "dials.algorithms.profile_model.gaussian_rs.calculator",
    "dials.command_line.integrate",
    "dials.algorithms.spot_prediction.reflection_predictor",
]


def process_one_image(experiment, table, params, integrator_class, collect_data=True):
    integrator = integrator_class(params, collect_data)
    try:
        experiment, table, collector = integrator.run(experiment, table)
    except RuntimeError as e:
        logger.info(f"Processing failed due to error: {e}")
        return (None, None, None)
    else:
        return experiment, table, collector


def setup(reflections, params):
    # calculate the batches for processing
    if params.image_range:
        if not len(params.image_range.split(":")) == 2:
            raise ValueError("Image range must be given in the form first:last")
        first, last = params.image_range.split(":")
        try:
            first = int(first)
        except ValueError:
            raise ValueError(f"Issue interpreting {first} as an integer")
        try:
            last = int(last)
        except ValueError:
            raise ValueError(f"Issue interpreting {last} as an integer")
        lowest_index = max(0, first - 1)
        highest_image = min(last, len(reflections))
        batches = list(range(lowest_index, highest_image, params.output.batch_size))
        batches.append(highest_image)
    else:
        batches = list(range(0, len(reflections), params.output.batch_size))
        batches.append(len(reflections))

    # Note, memory processing logic can go here
    if params.nproc is Auto:
        params.nproc = CPU_COUNT
    logger.info(f"Using {params.nproc} processes for integration")

    # aggregate some output for json, html etc
    if params.algorithm == "ellipsoid":
        process = EllipsoidIntegrator
        aggregator = EllipsoidOutputAggregator()
    elif params.algorithm == "stills":
        process = StillsIntegrator
        aggregator = OutputAggregator()
    else:
        raise ValueError("Invalid algorithm choice")

    configuration = {
        "process": process,
        "aggregator": aggregator,
        "params": params,
        "loggers_to_disable": (
            loggers_to_disable
            if params.algorithm == "ellipsoid"
            else loggers_to_disable_for_stills
        ),
    }

    return batches, configuration


@dataclass
class InputToIntegrate:
    integrator_class: type[SimpleIntegrator]
    experiment: Experiment
    table: flex.reflection_table
    params: Any
    crystalno: int
    imageset_index: int = 0


@dataclass
class IntegrationResult:
    experiment: Experiment
    table: flex.reflection_table
    collector: Any
    crystalno: int
    imageset_index: int = 0


def wrap_integrate_one(input_to_integrate: InputToIntegrate):
    collect_data = (
        input_to_integrate.params.output.html or input_to_integrate.params.output.json
    )
    expt, refls, collector = process_one_image(
        input_to_integrate.experiment,
        input_to_integrate.table,
        input_to_integrate.params,
        input_to_integrate.integrator_class,
        collect_data,
    )

    result = IntegrationResult(
        expt,
        refls,
        collector,
        input_to_integrate.crystalno,
        input_to_integrate.imageset_index,
    )
    if expt and refls:
        if not input_to_integrate.params.debug.output.shoeboxes:
            del result.table["shoebox"]
        logger.info(f"Processed crystal {input_to_integrate.crystalno}")
    if input_to_integrate.params.output.nuggets:
        img = input_to_integrate.experiment.imageset.get_image_identifier(0).split("/")[
            -1
        ]
        msg = {
            "crystal_no": input_to_integrate.crystalno,
            "n_integrated": 0,
            "i_over_sigma_overall": 0,
            "image": img,
        }
        if expt and refls:
            msg["n_integrated"] = collector.data["n_integrated"]
            msg["i_over_sigma_overall"] = round(
                collector.data["i_over_sigma_overall"], 2
            )

        with open(
            input_to_integrate.params.output.nuggets
            / f"nugget_integrated_{input_to_integrate.crystalno}.json",
            "w",
        ) as f:
            f.write(json.dumps(msg))
    return result


def process_batch(sub_tables, sub_expts, configuration, batch_offset=0):
    # create iterable
    input_iterable: list[InputToIntegrate] = []
    from dxtbx.imageset import ImageSequence, ImageSet

    original_isets = list(sub_expts.imagesets())
    identifiers_to_scans = {}
    n_iset = 0
    if any(sub_expts.scans()):
        identifiers_to_scans = {expt.identifier: expt.scan for expt in sub_expts}
    for i, (table, expt) in enumerate(zip(sub_tables, sub_expts)):
        if expt.scan:  # backcompatibility for indexed.expt without scans
            iset = expt.imageset
            if isinstance(iset, ImageSequence):
                # note, need to subtract imageset offset, as imageset may no longer start from same
                # index as the imported scan
                idx = (
                    expt.scan.get_array_range()[0] - iset.get_scan().get_batch_offset()
                )
                n_iset = original_isets.index(iset)
                subset = iset[idx : idx + 1]
                expt.imageset = ImageSet(
                    subset.data(), subset.indices()
                )  # Needed for stills integration code
                expt.scan = None  # Needed for some aspect of integration code, unclear what exactly.
        input_iterable.append(
            InputToIntegrate(
                configuration["process"],
                expt,
                table,
                configuration["params"],
                i + 1 + batch_offset,
                imageset_index=n_iset,
            )
        )
    input_iterable = sorted(input_iterable, key=lambda i: i.table.size(), reverse=True)
    with manage_loggers(
        configuration["params"].individual_log_verbosity,
        configuration["loggers_to_disable"],
    ):
        if configuration["params"].nproc > 1:
            with Pool(configuration["params"].nproc) as pool:
                results: list[IntegrationResult] = pool.map(
                    wrap_integrate_one, input_iterable
                )
        else:
            results: list[IntegrationResult] = [
                wrap_integrate_one(i) for i in input_iterable
            ]

    # then join
    integrated_reflections = flex.reflection_table()
    integrated_experiments = []

    use_beam = None
    use_gonio = None
    use_detector = None
    if len(sub_expts.beams()) == 1:
        use_beam = sub_expts.beams()[0]
    if len(sub_expts.goniometers()) == 1:
        use_gonio = sub_expts.goniometers()[0]
    if len(sub_expts.detectors()) == 1:
        use_detector = sub_expts.detectors()[0]

    n_integrated = 0
    for result in sorted(results, key=lambda result: result.crystalno):
        if result.table:
            if identifiers_to_scans:
                result.experiment.scan = identifiers_to_scans[
                    result.experiment.identifier
                ]
                result.experiment.imageset = original_isets[result.imageset_index]
                result.table["imageset_id"] = flex.int(
                    result.table.size(), result.imageset_index
                )
                if use_beam:
                    result.experiment.beam = use_beam
                if use_gonio:
                    result.experiment.goniometer = use_gonio
                if use_detector:
                    result.experiment.detector = use_detector
            ids_map = dict(result.table.experiment_identifiers())
            del result.table.experiment_identifiers()[list(ids_map.keys())[0]]
            result.table["id"] = flex.int(result.table.size(), n_integrated)
            result.table.experiment_identifiers()[n_integrated] = list(
                ids_map.values()
            )[0]
            n_integrated += 1
            integrated_reflections.extend(result.table)
            integrated_experiments.append(result.experiment)
            configuration["aggregator"].add_dataset(result.collector, result.crystalno)

    integrated_experiments = ExperimentList(integrated_experiments)
    integrated_reflections.assert_experiment_identifiers_are_consistent(
        integrated_experiments
    )
    return integrated_experiments, integrated_reflections


def run_integration(reflections, experiments, params):
    assert len(reflections) == len(experiments)
    if params.output.nuggets:
        params.output.nuggets = pathlib.Path(params.output.nuggets)
        if not params.output.nuggets.is_dir():
            logger.warning(
                "output.nuggets not recognised as a valid directory path, no nuggets will be output"
            )
            params.output.nuggets = None
    batches, configuration = setup(reflections, params)

    # now process each batch, and do parallel processing within a batch
    for i, b in enumerate(batches[:-1]):
        end_ = batches[i + 1]
        logger.info(f"Processing images {b+1} to {end_}")
        sub_tables = reflections[b:end_]
        sub_expts = experiments[b:end_]

        integrated_experiments, integrated_reflections = process_batch(
            sub_tables, sub_expts, configuration, batch_offset=b
        )
        yield (
            integrated_experiments,
            integrated_reflections,
            configuration["aggregator"],
        )


@show_mail_handle_errors()
def run(args: list[str] = None, phil=working_phil) -> None:
    """
    Run dials.ssx_integrate from the command-line.

    This program takes an indexed experiment list and reflection table and
    performs parallelised integration for synchrotron serial crystallography
    experiments. The programs acts as a wrapper to run one of two algorithms,
    the stills integrator or the 'ellipsoid' integrator (which uses a generalised
    ellipsoidal profile model). Analysis statistics are captured and output as
    a html report, while the output data are saved in batches for memory
    management.
    """

    parser = ArgumentParser(
        usage="dials.ssx_integrate indexed.expt indexed.refl [options]",
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
    )
    # Check the number of arguments is correct

    # Parse the command line
    params, options = parser.parse_args(args=args, show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    if len(reflections) == 0 or len(experiments) == 0:
        parser.print_help()
        return

    if len(reflections) != 1:
        raise Sorry(
            "Only a single reflection table file can be input (this can be a multi-still table)"
        )

    # Configure logging
    log.config(verbosity=options.verbose, logfile=params.output.log)
    params.individual_log_verbosity = options.verbose
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    ## FIXME - experiment identifiers approach won't work if input strong.refl and refined.expt
    # for now - check image path and update identifiers to that of refined.expt?
    if len(set(reflections[0]["id"]).difference({-1})) > 1:
        logger.info("Attempting to split multi-still reflection table")
        reflections = reflections[0].split_by_experiment_id()
        if not (len(reflections) == len(experiments)):
            raise Sorry(
                "Unequal number of reflection tables and experiments after splitting"
            )

    integrated_crystal_symmetries = []

    for i, (int_expt, int_refl, aggregator) in enumerate(
        run_integration(reflections, experiments, params)
    ):
        # combine beam and detector models if not already
        if len(int_expt.detectors()) > 1 or len(int_expt.beams()) > 1:
            combine = CombineWithReference(
                detector=int_expt[0].detector, beam=int_expt[0].beam
            )
            elist = ExperimentList()
            for expt in int_expt:
                elist.append(combine(expt))
            int_expt = elist
        reflections_filename = f"integrated_{i+1}.refl"
        experiments_filename = f"integrated_{i+1}.expt"
        logger.info(f"Saving {int_refl.size()} reflections to {reflections_filename}")
        int_refl.as_file(reflections_filename)
        logger.info(f"Saving the experiments to {experiments_filename}")
        int_expt.as_file(experiments_filename)

        integrated_crystal_symmetries.extend(
            [
                crystal.symmetry(
                    unit_cell=copy.deepcopy(cryst.get_unit_cell()),
                    space_group=copy.deepcopy(cryst.get_space_group()),
                )
                for cryst in int_expt.crystals()
            ]
        )

    plots, cluster_plots = ({}, {})
    if integrated_crystal_symmetries:
        cluster_plots, _ = report_on_crystal_clusters(
            integrated_crystal_symmetries,
            make_plots=(params.output.html or params.output.json),
        )

    if params.output.html or params.output.json:
        # now generate plots using the aggregated data.
        plots = aggregator.make_plots()
        plots.update(cluster_plots)

    if params.output.history:
        history = aggregator.make_history_json()
        with open(params.output.history, "w") as outfile:
            json.dump(history, outfile, indent=2)

    if params.output.html and plots:
        logger.info(f"Writing html report to {params.output.html}")
        generate_html_report(plots, params.output.html)
    if params.output.json and plots:
        logger.info(f"Saving plot data in json format to {params.output.json}")
        with open(params.output.json, "w") as outfile:
            json.dump(plots, outfile, indent=2)

    logger.info(
        "Further program documentation can be found at dials.github.io/ssx_processing_guide.html"
    )


if __name__ == "__main__":
    run()
