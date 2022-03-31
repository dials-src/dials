# LIBTBX_SET_DISPATCHER_NAME dev.dials.ssx_index
"""
This program runs indexing on the spotfinding results from a
still sequence i.e. SSX data. This wraps a call to the regular
indexing code, and so all parameters from dials.index can be set.

If a unit cell is given, indexing of each image will be attempted with the
fft1d algorithm, followed by the real_space_grid_search algorithm if indexing
with the fft1d algorithm was unsuccessful. If no unit cell is given, only fft1d
indexing can be attempted.

Indexing statistics are reported in a table and a unit cell clustering analysis
is performed, which can be useful for assessing the crystal symmetry if the unit
cell is unknown. An extensive html report is generated showing the indexing
and clustering statistics. The indexed data are saved into a single reflection
file and a single experiment list file, with a joint detector and beam model.

Further program documentation can be found at dials.github.io/ssx_processing_guide.html

Usage:
    dev.dials.ssx_index imported.expt strong.refl
    dev.dials.ssx_index imported.expt strong.refl unit_cell=x space_group=y
"""

from __future__ import annotations

import json
import logging
import sys
import time

from cctbx import crystal
from libtbx import Auto, phil
from libtbx.introspection import number_of_processors

from dials.algorithms.indexing.ssx.analysis import (
    generate_html_report,
    generate_plots,
    make_summary_table,
    report_on_crystal_clusters,
)
from dials.algorithms.indexing.ssx.processing import index
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

try:
    from typing import List
except ImportError:
    pass

logger = logging.getLogger("dials")

program_defaults_phil_str = """
indexing {
  method = fft1d
  stills {
    indexer = stills
  }
}
output.log = dials.ssx_index.log
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 1
      action = fix
    }
    beam.fix = all
    detector.fix = all
    scan_varying = False
  }
  reflections {
    weighting_strategy.override = stills
    outlier.algorithm = null
  }
}
"""

phil_scope = phil.parse(
    """
method = *fft1d *real_space_grid_search
    .type = choice(multi=True)
nproc = Auto
    .type = int
    .expert_level = 1
    .help = "Set the number of processors to use in indexing"
min_spots = 10
    .type = int
    .expert_level = 2
    .help = "Images with fewer than this number of strong spots will not be indexed"
output.html = dials.ssx_index.html
    .type = str
output.json = None
    .type = str
include scope dials.command_line.index.phil_scope
""",
    process_includes=True,
).fetch(phil.parse(program_defaults_phil_str))

phil_scope.adopt_scope(
    phil.parse(
        """
    individual_log_verbosity = 1
    .type =int
"""
    )
)


@show_mail_handle_errors()
def run(args: List[str] = None, phil: phil.scope = phil_scope) -> None:
    """
    Run dev.dials.ssx_index as from the command line.

    This program takes an imported experiment list and a reflection table
    of strong spots and performs parallelised indexing for synchrotron
    serial crystallography experiments. This is done by calling the regular
    dials indexing code and capturing output to provide a html report, and
    outputs a multi-image indexed.expt and indexed.refl file containing the
    indexed data.
    """

    parser = ArgumentParser(
        usage="dev.dials.ssx_index imported.expt strong.refl [options]",
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=__doc__,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    log.config(verbosity=options.verbose, logfile=params.output.log)
    params.individual_log_verbosity = options.verbose
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    if params.nproc is Auto:
        params.nproc = number_of_processors(return_value_if_unknown=1)

    if params.nproc > 1:
        params.indexing.nproc = params.nproc

    logger.info(f"Using {params.indexing.nproc} processes for indexing")

    st = time.time()
    indexed_experiments, indexed_reflections, summary_data = index(
        experiments, reflections[0], params
    )

    summary_table = make_summary_table(summary_data)
    logger.info("\nSummary of images sucessfully indexed\n" + summary_table)

    n_images = len([1 for v in summary_data.values() if v[0]["n_indexed"]])
    logger.info(f"{indexed_reflections.size()} spots indexed on {n_images} images\n")

    crystal_symmetries = [
        crystal.symmetry(
            unit_cell=expt.crystal.get_unit_cell(),
            space_group=expt.crystal.get_space_group(),
        )
        for expt in indexed_experiments
    ]
    if crystal_symmetries:
        cluster_plots, _ = report_on_crystal_clusters(
            crystal_symmetries,
            make_plots=(params.output.html or params.output.json),
        )

    logger.info(f"Saving indexed experiments to {params.output.experiments}")
    indexed_experiments.as_file(params.output.experiments)
    logger.info(f"Saving indexed reflections to {params.output.reflections}")
    indexed_reflections.as_file(params.output.reflections)

    if (params.output.html or params.output.json) and indexed_experiments:
        summary_plots = generate_plots(summary_data)
        if cluster_plots:
            summary_plots.update(cluster_plots)
        if params.output.html:
            generate_html_report(summary_plots, params.output.html)
        if params.output.json:
            with open(params.output.json, "w") as outfile:
                json.dump(summary_plots, outfile)

    logger.info(f"Total time: {time.time() - st:.2f}s")
    logger.info(
        "Further program documentation can be found at dials.github.io/ssx_processing_guide.html"
    )


if __name__ == "__main__":
    run()
