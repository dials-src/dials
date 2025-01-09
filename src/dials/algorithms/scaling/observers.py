"""
Helper functions/classes for making HTML report and scaling summary output.
"""

from __future__ import annotations

import json
import logging
import math

from jinja2 import ChoiceLoader, Environment, PackageLoader
from ordered_set import OrderedSet

from cctbx import uctbx
from dxtbx import flumpy
from scitbx.array_family import flex

from dials.algorithms.scaling.error_model.error_model import (
    calc_deltahl,
    calc_sigmaprime,
)
from dials.algorithms.scaling.error_model.error_model_target import (
    calculate_regression_x_y,
)
from dials.algorithms.scaling.model.model import (
    make_combined_plots,
    plot_scaling_models,
)
from dials.algorithms.scaling.plots import (
    error_model_variance_plot,
    error_regression_plot,
    normal_probability_plot,
    plot_outliers,
)
from dials.algorithms.scaling.scale_and_filter import make_scaling_filtering_plots
from dials.algorithms.scaling.scaling_library import (
    DialsMergingStatisticsError,
    merging_stats_from_scaled_array,
)
from dials.report.analysis import (
    make_merging_statistics_summary,
    reflection_tables_to_batch_dependent_properties,
    table_1_summary,
)
from dials.report.plots import (
    AnomalousPlotter,
    IntensityStatisticsPlots,
    ResolutionPlotsAndStats,
    i_over_sig_i_vs_batch_plot,
    i_over_sig_i_vs_i_plot,
    make_image_range_table,
    scale_rmerge_vs_batch_plot,
)
from dials.util import tabulate
from dials.util.batch_handling import batch_manager, get_image_ranges
from dials.util.exclude_images import get_valid_image_ranges
from dials.util.resolution_analysis import resolution_cc_half

logger = logging.getLogger("dials")


def assert_is_json_serialisable(thing, name, path=None):
    path = path or []
    if isinstance(thing, list):
        for n, element in enumerate(thing):
            assert_is_json_serialisable(element, name, path + [n])
    elif isinstance(thing, dict):
        for key, value in thing.items():
            assert_is_json_serialisable(value, name, path + [repr(key)])
    else:
        try:
            json.dumps(thing)
        except TypeError as e:
            raise TypeError(
                "JSON serialisation error '{}' for value '{}' type {} in {}{}".format(
                    e,
                    str(thing),
                    type(thing),
                    name,
                    "".join("[%s]" % step for step in path),
                )
            )


class ScalingSummaryContextManager:
    def __init__(self, script):
        self.script = script

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        print_scaling_summary(self.script)


def print_scaling_summary(script):
    """Log summary information after scaling."""
    logger.info(print_scaling_model_error_summary(script.experiments))
    valid_ranges = get_valid_image_ranges(script.experiments)
    image_ranges = get_image_ranges(script.experiments)
    msg = []
    for img, valid, refl in zip(image_ranges, valid_ranges, script.reflections):
        if valid:
            if len(valid) > 1 or valid[0][0] != img[0] or valid[-1][1] != img[1]:
                msg.append(
                    f"Excluded images for experiment id: {refl.experiment_identifiers().keys()[0]}, image range: {list(img)}, limited range: {list(valid)}"
                )
    if msg:
        msg = ["Summary of image ranges removed:"] + msg
        logger.info("\n".join(msg))

    # report on partiality of dataset
    partials = flex.double()
    for r in script.reflections:
        if "partiality" in r:
            partials.extend(r["partiality"])
    not_full_sel = partials < 0.99
    not_zero_sel = partials > 0.01
    gt_half = partials > 0.5
    lt_half = partials < 0.5
    partial_gt_half_sel = not_full_sel & gt_half
    partial_lt_half_sel = not_zero_sel & lt_half
    logger.info("Summary of dataset partialities")
    header = ["Partiality (p)", "n_refl"]
    rows = [
        ["all reflections", str(partials.size())],
        ["p > 0.99", str(not_full_sel.count(False))],
        ["0.5 < p < 0.99", str(partial_gt_half_sel.count(True))],
        ["0.01 < p < 0.5", str(partial_lt_half_sel.count(True))],
        ["p < 0.01", str(not_zero_sel.count(False))],
    ]
    logger.info(tabulate(rows, header))
    logger.info(
        """
Reflections below a partiality_cutoff of %s are not considered for any
part of the scaling analysis or for the reporting of merging statistics.
Additionally, if applicable, only reflections with a min_partiality > %s
were considered for use when refining the scaling model.
""",
        script.params.cut_data.partiality_cutoff,
        script.params.reflection_selection.min_partiality,
    )
    stats = script.merging_statistics_result
    if stats:
        anom_stats, cut_stats, cut_anom_stats = (None, None, None)
        if not script.scaled_miller_array.space_group().is_centric():
            anom_stats = script.anom_merging_statistics_result
        logger.info(make_merging_statistics_summary(stats))
        try:
            if script.params.cut_data.d_min is None:
                d_min = resolution_cc_half(stats, limit=0.3).d_min
            else:
                d_min = script.params.cut_data.d_min
        except RuntimeError as e:
            logger.debug(f"Resolution fit failed: {e}")
        else:
            max_current_res = stats.bins[-1].d_min
            if d_min and d_min - max_current_res > 0.005:
                logger.info(
                    "Resolution limit suggested from CC"
                    + "\u00bd"
                    + " fit (limit CC"
                    + "\u00bd"
                    + "=0.3): %.2f",
                    d_min,
                )
                try:
                    cut_stats, cut_anom_stats = merging_stats_from_scaled_array(
                        script.scaled_miller_array.resolution_filter(d_min=d_min),
                        script.params.output.merging.nbins,
                        script.params.output.use_internal_variance,
                        additional_stats=script.params.output.additional_stats,
                    )
                except DialsMergingStatisticsError:
                    pass
                else:
                    if script.scaled_miller_array.space_group().is_centric():
                        cut_anom_stats = None
        logger.info(table_1_summary(stats, anom_stats, cut_stats, cut_anom_stats))


class ScalingHTMLContextManager:
    def __init__(self, script):
        self.script = script

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        _make_scaling_html(self.script)


def _make_scaling_html(scaling_script):
    """Collect data from the individual observers and write the html."""
    html_file = scaling_script.params.output.html
    json_file = scaling_script.params.output.json
    if not (html_file or json_file):
        return
    data = {}
    data.update(
        make_scaling_model_plots(scaling_script.experiments, scaling_script.reflections)
    )
    data.update(
        make_outlier_plots(scaling_script.reflections, scaling_script.experiments)
    )
    data.update(
        make_error_model_plots(scaling_script.params, scaling_script.experiments)
    )
    data.update(make_merging_stats_plots(scaling_script))
    data.update(make_filtering_plots(scaling_script))
    if html_file:
        logger.info("Writing html report to %s", html_file)
        loader = ChoiceLoader(
            [
                PackageLoader("dials", "templates"),
                PackageLoader("dials", "static", encoding="utf-8"),
            ]
        )
        env = Environment(loader=loader)
        template = env.get_template("scaling_report.html")
        assert_is_json_serialisable(data, "self.data")
        html = template.render(
            page_title="DIALS scaling report",
            scaling_model_graphs=data["scaling_model"],
            scaling_tables=data["scaling_tables"],
            error_model_summary=data["error_model_summary"],
            resolution_plots=data["resolution_plots"],
            scaling_outlier_graphs=data["outlier_plots"],
            error_model_plots=data["error_model_plots"],
            anom_plots=data["anom_plots"],
            batch_plots=data["batch_plots"],
            image_range_tables=data["image_range_tables"],
            misc_plots=data["misc_plots"],
            filter_plots=data["filter_plots"],
        )
        with open(html_file, "wb") as f:
            f.write(html.encode("utf-8", "xmlcharrefreplace"))
    if json_file:
        logger.info("Writing html report data to %s", json_file)
        with open(json_file, "w", encoding="utf-8") as outfile:
            json.dump(data, outfile)


def make_scaling_model_plots(experiments, reflection_tables):
    """Collect scaling model plots for html report."""
    data = {i: e.scaling_model for i, e in enumerate(experiments)}
    d = {}
    combined_plots = make_combined_plots(data)
    if combined_plots:
        d.update(combined_plots)
    for key in sorted(data.keys()):
        scaling_model_plots = plot_scaling_models(data[key], reflection_tables[key])
        for plot in scaling_model_plots.values():
            plot["layout"]["title"] += f" (dataset {key})"
        for name, plot in scaling_model_plots.items():
            d[name + "_" + str(key)] = plot
    graphs = {"scaling_model": d}
    return graphs


def print_scaling_model_error_summary(experiments):
    """Get a summary of the error distribution of the models."""
    models = [e.scaling_model.to_dict() for e in experiments]
    first_model = models[0]
    component = first_model["configuration_parameters"]["corrections"][0]
    msg = ""
    if "est_standard_devs" in first_model[component]:
        p_sigmas = flex.double()
        for model in models:
            for component in model["configuration_parameters"]["corrections"]:
                if "est_standard_devs" in model[component]:
                    params = flex.double(model[component]["parameters"])
                    sigmas = flex.double(model[component]["est_standard_devs"])
                    null_value = flex.double(
                        len(params), model[component]["null_parameter_value"]
                    )
                    p_sigmas.extend(flex.abs(params - null_value) / sigmas)
        log_p_sigmas = flex.log(p_sigmas)
        frac_high_uncertainty = (log_p_sigmas < math.log(2)).count(True) / len(
            log_p_sigmas
        )
        if frac_high_uncertainty > 0.5:
            msg = (
                f"Warning: Over half ({frac_high_uncertainty * 100:.2f}%) of model parameters have significant\n"
                "uncertainty (sigma/abs(parameter) > 0.5), which could indicate a\n"
                "poorly-determined scaling problem or overparameterisation.\n"
            )
        else:
            msg = (
                f"{frac_high_uncertainty * 100:.2f}% of model parameters have significant uncertainty\n"
                "(sigma/abs(parameter) > 0.5)\n"
            )
    return msg


def make_outlier_plots(reflection_tables, experiments):
    """Make outlier plots for the HTML report."""
    data = {}
    for j, (table, expt) in enumerate(zip(reflection_tables, experiments)):
        outliers = table.get_flags(table.flags.outlier_in_scaling)
        x, y, z = table["xyzobs.px.value"].select(outliers).parts()
        if expt.scan and (expt.scan.get_oscillation()[1] != 0.0):
            zrange = [
                i / expt.scan.get_oscillation()[1]
                for i in expt.scan.get_oscillation_range()
            ]
        else:
            zrange = [0, 0]
        data[j] = {
            "x": list(x),
            "y": list(y),
            "z": list(z),
            "image_size": expt.detector[0].get_image_size(),
            "z_range": zrange,
        }
    d = {}
    for key in sorted(data):
        outlier_plots = plot_outliers(data[key])
        for plot in outlier_plots.values():
            if plot:  # may be null if no outliers
                plot["layout"]["title"] += f" (dataset {key})"
        d["outlier_plot_" + str(key)] = outlier_plots["outlier_xy_positions"]
        d["outlier_plot_z" + str(key)] = outlier_plots["outliers_vs_z"]
    graphs = {"outlier_plots": d}
    return graphs


def make_error_model_plots(params, experiments):
    """Generate normal probability plot data."""
    d = {"error_model_plots": {}, "error_model_summary": "No error model applied"}
    error_model_data = []  # a list of dicts of error model data
    if experiments[0].scaling_model.error_model:
        error_models = [e.scaling_model.error_model for e in experiments]
        unique_error_models = OrderedSet(error_models)
        if len(unique_error_models) == 1:
            d["error_model_summary"] = str(error_models[0])
        else:
            d["error_model_summary"] = ""
            for i, e in enumerate(unique_error_models):
                indices = [str(j + 1) for j, x in enumerate(error_models) if e is x]
                d["error_model_summary"] += (
                    f"\nError model {i+1}, applied to sweeps {', '.join(indices)}:"
                    + str(e)
                )
        for em in unique_error_models:
            if em.filtered_Ih_table:
                data_i = {}
                table = em.filtered_Ih_table
                data_i["intensity"] = table.intensities
                sigmaprime = calc_sigmaprime(em.parameters, table)
                data_i["delta_hl"] = calc_deltahl(table, table.calc_nh(), sigmaprime)
                data_i["inv_scale"] = table.inverse_scale_factors
                data_i["sigma"] = sigmaprime * data_i["inv_scale"]
                data_i["binning_info"] = em.binner.binning_info
                em.clear_Ih_table()
                if params.weighting.error_model.basic.minimisation == "regression":
                    x, y = calculate_regression_x_y(em.filtered_Ih_table)
                    data_i["regression_x"] = x
                    data_i["regression_y"] = y
                    data_i["model_a"] = em.parameters[0]
                    data_i["model_b"] = em.parameters[1]
                error_model_data.append(data_i)

    if error_model_data:
        for i, emd in enumerate(error_model_data):
            d["error_model_plots"].update(normal_probability_plot(emd, label=i + 1))
            d["error_model_plots"].update(
                i_over_sig_i_vs_i_plot(
                    flumpy.from_numpy(emd["intensity"]),
                    flumpy.from_numpy(emd["sigma"]),
                    label=i + 1,
                )
            )
            d["error_model_plots"].update(error_model_variance_plot(emd, label=i + 1))
            if "regression_x" in emd:
                d["error_model_plots"].update(error_regression_plot(emd, label=i + 1))
    return d


def make_filtering_plots(script):
    """Make filtering plots for HTML report"""
    if script.filtering_results:
        data = {
            "merging_stats": script.filtering_results.get_merging_stats(),
            "initial_expids_and_image_ranges": script.filtering_results.initial_expids_and_image_ranges,
            "cycle_results": script.filtering_results.get_cycle_results(),
            "expids_and_image_ranges": script.filtering_results.expids_and_image_ranges,
            "mode": script.params.filtering.deltacchalf.mode,
        }
        d = make_scaling_filtering_plots(data)
        return {"filter_plots": d}
    return {"filter_plots": {}}


def make_merging_stats_plots(script, run_xtriage_analysis=False, make_batch_plots=True):
    """Make merging stats plots for HTML report"""
    d = {
        "scaling_tables": {},
        "resolution_plots": {},
        "batch_plots": {},
        "misc_plots": {},
        "anom_plots": {},
        "image_range_tables": [],
    }
    if make_batch_plots:
        (
            batches,
            rvb,
            isigivb,
            svb,
            batch_data,
        ) = reflection_tables_to_batch_dependent_properties(  # pylint: disable=unbalanced-tuple-unpacking
            script.reflections,
            script.experiments,
            script.scaled_miller_array,
        )
        bm = batch_manager(batches, batch_data)
        image_range_tables = make_image_range_table(script.experiments, bm)
        d["image_range_tables"] = [image_range_tables]
        d["batch_plots"].update(scale_rmerge_vs_batch_plot(bm, rvb, svb))
        d["batch_plots"].update(i_over_sig_i_vs_batch_plot(bm, isigivb))

    if script.merging_statistics_result:
        stats = script.merging_statistics_result
        anom_stats = script.anom_merging_statistics_result
        is_centric = script.scaled_miller_array.space_group().is_centric()
        # Now calculate batch data
        plotter = ResolutionPlotsAndStats(stats, anom_stats, is_centric)
        d["resolution_plots"].update(plotter.make_all_plots())
        d["scaling_tables"] = plotter.statistics_tables()
        plotter = IntensityStatisticsPlots(
            script.scaled_miller_array, run_xtriage_analysis=run_xtriage_analysis
        )
        d["xtriage_output"] = plotter.generate_xtriage_output()
        d["resolution_plots"].update(plotter.generate_resolution_dependent_plots())
        if d["resolution_plots"]["cc_one_half"]["data"][2]:
            cc_anom = d["resolution_plots"]["cc_one_half"]["data"][2]["y"]
            significance = d["resolution_plots"]["cc_one_half"]["data"][3]["y"]
            sig = flex.double(cc_anom) > flex.double(significance)
            max_anom = 0
            for i, v in enumerate(sig):
                if v:
                    max_anom = i
                else:
                    break
            d_min = uctbx.d_star_sq_as_d(plotter.binner.limits())[max_anom + 1]
        else:
            d_min = 0.0
        d["misc_plots"].update(plotter.generate_miscellanous_plots())
        intensities_anom = script.scaled_miller_array.as_anomalous_array()
        intensities_anom = intensities_anom.map_to_asu().customized_copy(
            info=script.scaled_miller_array.info()
        )
        anom_plotter = AnomalousPlotter(intensities_anom, strong_cutoff=d_min)
        d["anom_plots"].update(anom_plotter.make_plots())
    return d
