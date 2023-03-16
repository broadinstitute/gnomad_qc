import base64
import io
import logging
import math
import random
import uuid
from collections import defaultdict
from os import path

import hail as hl
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bokeh.layouts import gridplot
from bokeh.models import DataRange1d, Legend, LegendItem, Range1d, Span, Title
from bokeh.models.mappers import CategoricalColorMapper
from bokeh.models.widgets import Panel, Tabs
from bokeh.palettes import Category10, viridis
from bokeh.plotting import figure
from hail.utils.misc import new_temp_file
from IPython.core.display import display
from IPython.display import Image, display_html
from upsetplot import from_indicators, plot

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)

outlier_metrics = {
    "n_snp",
    "n_singleton",
    "r_ti_tv",
    "r_insertion_deletion",
    "n_insertion",
    "n_deletion",
    "r_snp_indel",
    "r_het_hom_var",
    "n_transition",
    "n_transversion",
    "r_ti_tv_singleton",
}

tmp_prefix = "gs://gnomad-tmp/julia/outlier_filter/"


def get_hist_with_cutoff(
    hist_expr,
    bins=100,
    cutoff_locations=[],
    title=None,
    range=None,
    plot_width=500,
    plot_height=400,
    color=None,
):
    hist = hl.plot.histogram(hist_expr, bins=bins, range=range)
    if color is not None:
        for r in hist.rendeoutrers:
            r.glyph.fill_color = r.glyph.line_color = color

    hist.renderers.extend(
        [
            Span(location=x, dimension="height", line_dash="dashed", line_color="red")
            for x in cutoff_locations
        ]
    )
    hist.legend.visible = False
    hist.plot_width = plot_width
    hist.plot_height = plot_height

    title_lines = title.split("\n")
    title_lines.reverse()
    for line_text in title_lines:
        hist.add_layout(Title(text=line_text), "above")

    return hist


def get_group_count_frac(ht, groupings, metric, cutoffs):
    return ht.group_by(*groupings).aggregate(
        **{
            f"n_{cutoff}": hl.agg.count_where(ht[metric] < cutoff) for cutoff in cutoffs
        },
        **{
            f"frac_{cutoff}": hl.agg.fraction(ht[metric] < cutoff) for cutoff in cutoffs
        },
    )


def get_pc_plots(
    pcs_ht,
    pc_name,
    title,
    label,
    hover_fields=None,
    n_pcs=10,
    plot_height=800,
    plot_width=1300,
    range_ht=None,
    colors=None,
):
    plots = []
    for pc in range(1, n_pcs, 2):
        p = hl.plot.scatter(
            pcs_ht[f"{pc_name}{pc}"],
            pcs_ht[f"{pc_name}{pc + 1}"],
            label=label,
            hover_fields=hover_fields,
            width=plot_width,
            height=plot_height,
            colors=colors,
            title=title,
            xlabel=f"PC{pc}",
            ylabel=f"PC{pc + 1}",
        )
        if range_ht is not None:
            p.x_range = Range1d(
                range_ht.aggregate(hl.agg.min(range_ht[f"{pc_name}{pc}"])),
                range_ht.aggregate(hl.agg.max(range_ht[f"{pc_name}{pc}"])),
            )
            p.y_range = Range1d(
                range_ht.aggregate(hl.agg.min(range_ht[f"{pc_name}{pc + 1}"])),
                range_ht.aggregate(hl.agg.max(range_ht[f"{pc_name}{pc + 1}"])),
            )
        plots.append(Panel(child=p, title=f"PC{pc} vs PC{pc + 1}"))
    return Tabs(tabs=plots)


def format_bokeh_plot(
    p,
    plot_width=600,
    plot_height=400,
    title_font_size="14pt",
    label_font_size="12pt",
    axis_font_size="12pt",
    legend_position=None,
    alpha=None,
):
    p.plot_width = plot_width
    p.plot_height = plot_height
    if p.title is not None:
        p.title.text_font_size = title_font_size
    p.xaxis.axis_label_text_font_size = label_font_size
    p.yaxis.axis_label_text_font_size = label_font_size
    p.xaxis.major_label_text_font_size = axis_font_size
    p.yaxis.major_label_text_font_size = axis_font_size
    p.legend.label_text_font_size = label_font_size
    p.xgrid.visible = False
    p.ygrid.visible = False
    p.background_fill_color = None
    p.outline_line_color = None

    if p.legend and legend_position is not None:
        p.add_layout(p.legend[0], legend_position)

    if alpha is not None:
        for r in p.renderers:
            r.glyph.fill_alpha = alpha
            r.glyph.line_alpha = alpha

    return p


def pair_plot(
    ht,
    metrics=None,
    label_col=None,
    per_metric_label_col=None,
    colors=None,
    hover_fields=None,
    cutoffs=None,
    **kwargs,
):
    """
    Plot each column of `ht` against each other and returns a grid of plots.

    The diagonal contains a histogram of each column, or a density plot if labels are provided.
    The lower diagonal contains scatter plots of each column against each other.
    The upper diagonal is empty.

    All columns should be numerical with the exception of the `label_col` if provided.
    If a color dict containing provided mapping labels to specific colors can be specified using `color_dict`
    """
    if label_col is None and per_metric_label_col is None and colors is not None:
        logger.warning(
            "`colors_dict` ignored since no `label_col` or `per_metric_label_col` is"
            " specified"
        )

    colors_col = "__pair_plot_color"

    colors_dict = {}
    if label_col is None and per_metric_label_col is None:
        ht = ht.annotate(**{colors_col: viridis(1)})
        labels = None
    elif per_metric_label_col is not None:
        ht = ht.annotate(
            **{
                f"{label_col}{m}": hl.str(ht[f"{label_col}{m}"])
                for m, m_label in metrics
            }
        )
        labels = [
            "x outlier",
            "y outlier",
            "x and y outlier",
        ]  # ht.aggregate(hl.agg.collect_as_set(ht[
        # label_col])) - {None}
        if not isinstance(colors, dict):
            color_palette = viridis(len(labels)) if colors is None else colors
            colors_dict = {l: color_palette[i] for i, l in enumerate(labels)}
        else:
            colors_dict = colors
        # ht = ht.annotate(**{colors_col: hl.literal(colors_dict).get(ht[label_col], "gray")})
    else:
        labels = ht.aggregate(hl.agg.collect_as_set(ht[label_col])) - {None}
        if not isinstance(colors, dict):
            color_palette = viridis(len(labels)) if colors is None else colors
            colors_dict = {l: color_palette[i] for i, l in enumerate(labels)}
        else:
            colors_dict = colors
        # ht = ht.annotate(**{colors_col: hl.literal(colors_dict).get(ht[label_col], "gray")})

    factors = []
    palette = []
    for f, c in colors_dict.items():
        factors.append(f)
        palette.append(c)
    colors_mapper = CategoricalColorMapper(
        factors=factors,
        palette=palette,
    )
    if metrics is None:
        metrics = [(m, m) for m in set(ht.row.keys()) - set(ht.key.keys())]
    metrics = [
        (m, m_label) for m, m_label in metrics if m not in [colors_col, label_col]
    ]
    range_agg = {}
    for m in metrics:
        m = m[0]
        range_agg[f"{m}_min"] = hl.agg.min(ht[m])
        range_agg[f"{m}_max"] = hl.agg.max(ht[m])

    min_max = ht.aggregate(hl.struct(**range_agg))
    data_ranges = []
    for m in metrics:
        rmin = min_max[f"{m[0]}_min"]
        rmax = min_max[f"{m[0]}_max"]
        data_ranges.append(
            DataRange1d(
                start=rmin - (abs(rmin - rmax) * 0.05),
                end=rmax + (abs(rmin - rmax) * 0.05),
            )
        )

    n_cols = len(metrics)
    plot_grid = []
    renderer_list = []
    for i in range(n_cols):
        row = [None] * n_cols
        for j in range(i + 1):
            if cutoffs is not None:
                i_slopes = []
                for intercept in [
                    cutoffs[metrics[i][0]].lower,
                    cutoffs[metrics[i][0]].upper,
                ]:
                    if (intercept != math.inf) and (intercept != -math.inf):
                        i_slopes.append(
                            Span(
                                location=intercept,
                                dimension="width",
                                line_color="red",
                                line_dash="dashed",
                                line_width=3.5,
                            )
                        )
                j_slopes = []
                for intercept in [
                    cutoffs[metrics[j][0]].lower,
                    cutoffs[metrics[j][0]].upper,
                ]:
                    if (intercept != math.inf) and (intercept != -math.inf):
                        j_slopes.append(
                            Span(
                                location=intercept,
                                dimension="height",
                                line_color="red",
                                line_dash="dashed",
                                line_width=3.5,
                            )
                        )
            else:
                i_slopes = None
                j_slopes = None
            x_axis_label = metrics[j][1] if i == n_cols - 1 else ""
            y_axis_label = metrics[i][1] if j == 0 else ""
            x_range = data_ranges[j]

            if i == j:
                kwargs2 = kwargs.copy()
                p = figure()
                try:
                    p = hl.plot.histogram(ht[metrics[j][0]], bins=50)
                    if j_slopes is not None:
                        for slope in j_slopes:
                            p.add_layout(slope)
                    p.legend.items = p.legend.items[1:]
                    kwargs2["alpha"] = None
                    p = format_bokeh_plot(p, **kwargs2)
                except BaseException:
                    errors = True
            else:
                y_range = data_ranges[i]
                if label_col is not None:
                    p = hl.plot.scatter(
                        ht[metrics[j][0]],
                        ht[metrics[i][0]],
                        label={label_col: ht[label_col]},
                        colors=colors_mapper,
                        hover_fields={l: ht[m] for l, m in hover_fields.items()},
                    )
                elif per_metric_label_col is not None:
                    ht = ht.annotate(
                        color_col=hl.case()
                        .when(
                            (ht[f"{label_col}{metrics[j][0]}"] == "true")
                            & (ht[f"{label_col}{metrics[i][0]}"] == "true"),
                            "x and y outlier",
                        )
                        .when(ht[f"{label_col}{metrics[j][0]}"] == "true", "x outlier")
                        .when(ht[f"{label_col}{metrics[i][0]}"] == "true", "y outlier")
                        .or_missing()
                    )
                    p = hl.plot.scatter(
                        ht[metrics[j][0]],
                        ht[metrics[i][0]],
                        label={"color_col": ht["color_col"]},
                        colors=colors_mapper,
                        hover_fields={l: ht[m] for l, m in hover_fields.items()},
                    )
                else:
                    p = hl.plot.scatter(
                        ht[metrics[j][0]],
                        ht[metrics[i][0]],
                        hover_fields={l: ht[m] for l, m in hover_fields.items()},
                        # label={colors_col: ht[colors_col]},
                    )
                    p.legend.items = p.legend.items[1:]
                p.y_range = y_range
                p = format_bokeh_plot(p, legend_position=None, **kwargs)
                for l in p.legend.items:
                    c = colors_dict.get(l.label["value"], "gray")
                    for r in l.renderers:
                        r.glyph.fill_color = r.glyph.line_color = c

                if i_slopes is not None:
                    for slope in i_slopes:
                        p.add_layout(slope)
                if j_slopes is not None:
                    for slope in j_slopes:
                        p.add_layout(slope)
            p.x_range = x_range
            p.xaxis.axis_label = x_axis_label
            p.yaxis.axis_label = y_axis_label
            if p.legend is not None:
                p.legend.visible = False
            renderer_list.extend(p.renderers)
            row[j] = p
        plot_grid.append(row)
    legend_items = [
        LegendItem(
            label=label,
            renderers=[
                renderer
                for renderer in renderer_list
                if renderer.glyph.fill_color == color
            ],
        )
        for label, color in colors_dict.items()
    ]
    legend_fig = figure(
        plot_width=2000, plot_height=75, outline_line_alpha=0, toolbar_location=None
    )
    for c in [
        legend_fig.grid[0],
        legend_fig.ygrid[0],
        legend_fig.xaxis[0],
        legend_fig.yaxis[0],
    ]:
        c.visible = False
    legend_fig.renderers += renderer_list
    legend_fig.x_range.end = 4005
    legend_fig.x_range.start = 4000
    legend_fig.add_layout(
        Legend(
            click_policy="hide",
            location="top_left",
            border_line_alpha=0,
            items=legend_items,
            orientation="horizontal",
        )
    )
    legend_fig.legend.visible = True
    legend_fig.legend.margin = 10
    legend_fig.legend.spacing = 10
    legend_fig.legend.label_text_font_size = "16pt"
    legend_fig.margin = (40, 0, 10, 100)
    legend_fig.legend.glyph_width = 30
    legend_fig.legend.glyph_height = 30
    legend_fig.background_fill_color = None
    legend_fig.border_fill_color = None

    return gridplot(
        [[legend_fig], [gridplot(plot_grid, toolbar_location="left")]],
        toolbar_location=None,
    )


def my_pair_plot(
    ht,
    cutoffs,
    hue=None,
    output_tab=None,
    title=None,
    file_name=None,
    use_fig_if_exists=False,
    add_kde_plot=True,
):
    if (
        use_fig_if_exists
        and file_name is not None
        and path.exists(f"plot_pngs/pair_plot_{file_name}.png")
    ):
        if output_tab is not None:
            output_tab.append_display_data(
                Image(f"plot_pngs/pair_plot_{file_name}.png")
            )
        else:
            display(Image(f"plot_pngs/pair_plot_{file_name}.png"))

        return

    df = ht.to_pandas()

    def my_hist(x, label, color, cutoffs):
        ax0 = plt.gca()
        ax = ax0.twinx()
        ax.yaxis.tick_left()
        ax.yaxis.label_position = "left"
        ax.hist(x.dropna(), label=label, color=color)
        if x.name in cutoffs:
            for intercept, c in [
                (cutoffs[x.name].lower, "red"),
                (cutoffs[x.name].upper, "purple"),
            ]:
                if (
                    (intercept != math.inf)
                    and (intercept != -math.inf)
                    and (intercept is not None)
                    and not math.isnan(intercept)
                ):
                    ax.axvline(intercept, color=c, linestyle="--")

    def col_nan_scatter(x, y, **kwargs):
        df_no_na = pd.DataFrame({"x": x[:], "y": y[:]})
        df_no_na = df_no_na.dropna()
        x = df_no_na["x"]
        y = df_no_na["y"]
        plt.gca()
        sns.scatterplot(x=x, y=y, **kwargs)

    def kdeplot_catch_error(x, y, **kwargs):
        plt.gca()
        try:
            sns.kdeplot(x=x, y=y, **kwargs)
        except ValueError:
            pass
        except np.linalg.LinAlgError:
            pass

    if hue is not None:
        df[hue] = df[hue].map({True: "True", False: "False"})
    g = sns.PairGrid(df, diag_sharey=False, dropna=True)
    if hue is None:
        g = g.map_lower(col_nan_scatter)
    else:
        g = g.map_lower(col_nan_scatter, hue=df[hue], hue_order=["False", "True"])

    if add_kde_plot:
        g = g.map_upper(kdeplot_catch_error)
    g = g.map_diag(my_hist, cutoffs=cutoffs)

    for i, ax_r in enumerate(g.axes):
        ylab = ax_r[0].properties()["ylabel"]
        if "|" in ylab:
            ylab = ylab.split("|")[1]
        for j, ax in enumerate(ax_r):
            xlab = g.axes[-1][j].properties()["xlabel"]
            if "|" in xlab:
                xlab = xlab.split("|")[1]
            if i == j:
                ax.tick_params(
                    axis="both",
                    left=False,
                    labelleft=False,
                    bottom=True,
                    labelbottom=True,
                )
                ax.xaxis.set_label_text(ylab, visible=True)
                ax.yaxis.set_label_text("Counts", visible=True)
                ax.yaxis.labelpad = 30
            else:
                if xlab in cutoffs:
                    for intercept, color in [
                        (cutoffs[xlab].lower, "red"),
                        (cutoffs[xlab].upper, "purple"),
                    ]:
                        if (
                            (intercept != math.inf)
                            and (intercept != -math.inf)
                            and (intercept is not None)
                            and not math.isnan(intercept)
                        ):
                            ax.axvline(intercept, color=color, linestyle="--")
                if ylab in cutoffs:
                    for intercept, color in [
                        (cutoffs[ylab].lower, "red"),
                        (cutoffs[ylab].upper, "purple"),
                    ]:
                        if (
                            (intercept != math.inf)
                            and (intercept != -math.inf)
                            and (intercept is not None)
                            and not math.isnan(intercept)
                        ):
                            ax.axhline(intercept, color=color, linestyle="--")
                ax.tick_params(axis="both", labelleft=True, labelbottom=True)
                ax.xaxis.set_label_text(xlab, visible=True)
                ax.yaxis.set_label_text(ylab, visible=True)

    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    if title is not None:
        plt.suptitle(title, y=1.02)
    # plt.show()
    if file_name is None:
        file_name = str(uuid.uuid4())
    plt.savefig(f"plot_pngs/pair_plot_{file_name}.png", bbox_inches="tight")
    fig = plt.gcf()
    plt.close(fig)
    if output_tab is not None:
        output_tab.append_display_data(Image(f"plot_pngs/pair_plot_{file_name}.png"))
    else:
        display(Image(f"plot_pngs/pair_plot_{file_name}.png"))


def box_plot_strat_pop_platform(
    outlier_ht,
    strat_pop_platform_ht,
    metric_labels,
    pair_plot_metrics,
    figsize=(30, 10),
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize

    tabs = []
    strat_cutoff_dict = hl.eval(strat_pop_platform_ht.qc_metrics_stats)
    strat = list(strat_cutoff_dict.keys())
    strat_pops = sorted(set([strata[0] for strata in strat]) - {None})
    strat_platform = sorted(set([strata[1] for strata in strat]) - {None})

    order = []
    for platform in strat_platform:
        for pop in strat_pops:
            if (pop, platform) in strat:
                order.append(pop + " / " + platform)

    outputs = []
    outputs.append(widgets.Output())
    for m in [m[0] for m in pair_plot_metrics]:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)
    for i, m in enumerate([m[0] for m in pair_plot_metrics]):
        plot_grid = []
        ht_strata_pd = outlier_ht.select(
            m, "pop", pop_platform=outlier_ht.pop + " / " + outlier_ht.platform
        ).to_pandas()
        o = outputs[i + 1]
        fig, ax = plt.subplots(figsize=(8, 60))
        g = sns.boxplot(
            data=ht_strata_pd,
            x=m,
            y="pop_platform",
            showfliers=False,
            ax=ax,
            hue="pop",
            order=order,
            hue_order=strat_pops,
            dodge=False,
        )
        g.xaxis.set_label_text(metric_labels[m], visible=True)
        g.set(title=metric_labels[m])
        # plt.show()
        file_name = str(uuid.uuid4())
        plt.savefig(
            f"plot_pngs/{m}_box_plot_strat_pop_platform_{file_name}.png",
            bbox_inches="tight",
        )
        plt.close(fig)
        o.append_display_data(
            Image(f"plot_pngs/{m}_box_plot_strat_pop_platform_{file_name}.png")
        )
        tab.set_title(i + 1, m)

    display(tab)


def box_plot_regress_pop_strat_platform(
    outlier_ht,
    strat_pop_platform_ht,
    metric_labels,
    pair_plot_metrics,
    figsize=(30, 10),
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize
    tabs = []
    strat_cutoff_dict = hl.eval(strat_pop_platform_ht.qc_metrics_stats)
    strat = list(strat_cutoff_dict.keys())
    strat_pops = sorted(set([strata[0] for strata in strat]) - {None})
    strat_platform = sorted(set([strata[1] for strata in strat]) - {None})

    order = []
    for platform in strat_platform:
        for pop in strat_pops:
            if (pop, platform) in strat:
                order.append(pop + " / " + platform)

    outputs = []
    outputs.append(widgets.Output())
    for m in [m[0] for m in pair_plot_metrics]:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)
    for i, m in enumerate([m[0] for m in pair_plot_metrics]):
        plot_grid = []
        ht_strata_pd = outlier_ht.select(
            f"regress_pop_strat_platform|{m}_residual", "pop", "platform"
        ).to_pandas()  # ,
        # pop_platform=outlier_ht.pop + " / " + outlier_ht.platform).to_pandas()
        o = outputs[i + 1]
        # fig, ax = plt.subplots(figsize=(8, 60))
        g = sns.catplot(
            data=ht_strata_pd,
            x=f"regress_pop_strat_platform|{m}_residual",
            y="pop",
            col="platform",
            kind="box",
            col_wrap=2,
            showfliers=False,
            hue="pop",
            order=strat_pops,
            hue_order=strat_pops,
            dodge=False,
            sharex=False,
            sharey=False,
        )
        # g = sns.boxplot(data=ht_strata_pd, x=f"regress_pop_strat_platform|{
        # m}_residual", y="pop_platform", showfliers = False, ax=ax, hue="pop",
        # order=order, hue_order=strat_pops, dodge = False)
        for i, ax_r in enumerate(g.axes.flatten()):
            ax_r.xaxis.set_label_text(f"{m}_residual", visible=True)
        # g.set(title=f"{m}_residual")
        # plt.show()
        file_name = str(uuid.uuid4())
        plt.savefig(
            f"plot_pngs/{m}_box_plot_regress_pop_strat_platform_{file_name}.png",
            bbox_inches="tight",
        )
        fig = plt.gcf()
        plt.close(fig)
        o.append_display_data(
            Image(f"plot_pngs/{m}_box_plot_regress_pop_strat_platform_{file_name}.png")
        )
        # tab.set_title(i+1, f"{m}_residual")
        tab.set_title(i + 1, m)

    display(tab)


def box_plot_regress_pop_platform(
    outlier_ht,
    strat_pop_platform_ht,
    metric_labels,
    pair_plot_metrics,
    figsize=(30, 10),
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize
    tabs = []
    strat_cutoff_dict = hl.eval(strat_pop_platform_ht.qc_metrics_stats)
    strat = list(strat_cutoff_dict.keys())
    strat_pops = sorted(set([strata[0] for strata in strat]) - {None})
    strat_platform = sorted(set([strata[1] for strata in strat]) - {None})

    order = []
    for platform in strat_platform:
        for pop in strat_pops:
            if (pop, platform) in strat:
                order.append(pop + " / " + platform)

    outputs = []
    outputs.append(widgets.Output())
    for m in [m[0] for m in pair_plot_metrics]:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)
    for i, m in enumerate([m[0] for m in pair_plot_metrics]):
        plot_grid = []
        ht_strata_pd = outlier_ht.select(
            f"regress_pop_platform|{m}_residual",
            "pop",
            pop_platform=outlier_ht.pop + " / " + outlier_ht.platform,
        ).to_pandas()
        o = outputs[i + 1]
        fig, ax = plt.subplots(figsize=(8, 60))
        g = sns.boxplot(
            data=ht_strata_pd,
            x=f"regress_pop_platform|{m}_residual",
            y="pop_platform",
            showfliers=False,
            ax=ax,
            hue="pop",
            order=order,
            hue_order=strat_pops,
            dodge=False,
        )
        g.xaxis.set_label_text(f"{m}_residual", visible=True)
        g.set(title=f"{m}_residual")
        # plt.show()
        file_name = str(uuid.uuid4())
        plt.savefig(
            f"plot_pngs/{m}_box_plot_regress_pop_platform_{file_name}.png",
            bbox_inches="tight",
        )
        plt.close(fig)
        o.append_display_data(
            Image(f"plot_pngs/{m}_box_plot_regress_pop_platform_{file_name}.png")
        )
        tab.set_title(i + 1, f"{m}_residual")

    display(tab)


def get_hist_plots_strat_pop_platform(
    outlier_sample_qc,
    cutoffs,
    fail_prefix="",
    qc_metrics=list(outlier_metrics),
    figsize=(30, 10),
    bins=100,
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize
    cols = ["pop", "platform"]
    key = "s"
    # colnames = cols + [key] + [f'{metric}' for metric in qc_metrics]
    # fail_colnames = cols + [key] + [f'{fail_prefix}{metric}' for metric in qc_metrics]
    colnames = cols + [f"{metric}" for metric in qc_metrics]
    fail_colnames = cols + [f"{fail_prefix}{metric}" for metric in qc_metrics]
    sample_qc_pd = outlier_sample_qc.select(*colnames).to_pandas()
    sample_qc_fail_pd = (
        outlier_sample_qc.select(*fail_colnames)
        .rename(dict(zip(fail_colnames, colnames)))
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)
    plots = None
    tables_pop_platform = []
    tables_pop = []
    tables_platform = []

    def my_hist(data, x, color, cutoffs):
        ax = plt.gca()
        # ax = ax0.twinx()
        # ax.yaxis.tick_left()
        # ax.yaxis.label_position = "left"
        pop = data["pop"][0]
        platform = data["platform"][0]
        ax.hist(data[x].dropna(), color=color, bins=bins)
        if (pop, platform) in cutoffs and x in cutoffs[(pop, platform)]:
            cutoffs = cutoffs[(pop, platform)][x]
            for intercept, c in [(cutoffs.lower, "red"), (cutoffs.upper, "purple")]:
                if (
                    (intercept != math.inf)
                    and (intercept != -math.inf)
                    and (intercept is not None)
                    and not math.isnan(intercept)
                ):
                    ax.axvline(intercept, color=c, linestyle="--")

    outputs = [widgets.Output()]
    for metric in qc_metrics:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)

    for i, metric in enumerate(qc_metrics):
        curve_dict = {}
        fail_table = (
            sample_qc_fail_pd.groupby(cols)[metric].value_counts().unstack().fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_pop_platform.append((metric, fail_table))

        fail_table = (
            sample_qc_fail_pd.groupby(["pop"])[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_pop.append((metric, fail_table))

        fail_table = (
            sample_qc_fail_pd.groupby(["platform"])[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_platform.append((metric, fail_table))

        o = outputs[i + 1]
        g = sns.FacetGrid(
            sample_qc_pd, col="pop", row="platform", sharex=False, sharey=False
        )
        g.map_dataframe(my_hist, x=metric, cutoffs=cutoffs)
        # plt.show()
        file_name = str(uuid.uuid4())
        plt.savefig(
            f"plot_pngs/{metric}_hist_plots_strat_pop_platform_{file_name}.png",
            bbox_inches="tight",
        )
        fig = plt.gcf()
        plt.close(fig)
        o.append_display_data(
            Image(f"plot_pngs/{metric}_hist_plots_strat_pop_platform_{file_name}.png")
        )
        tab.set_title(i + 1, metric)

    display(tab)

    return tables_pop_platform, tables_pop, tables_platform


def display_tables(table_list):
    head = """<html>
<body>
<table style="width:100%; display: block; overflow-x: auto; white-space: nowrap;">
<thead align="center">
<tr>
{}
</tr>
</thead>
<tr>
"""
    titles = ""
    row = ""
    for title, serie in table_list:
        s = serie.copy()
        titles += f'<th style="text-align: center;">{title}</th>\n'
        s.name = ""
        row += "<td>{}</td>".format(s.to_html())

    head = head.format(titles)
    head += row
    head += """
</tr>
</table>
</body>
</html>"""
    display_html(head, raw=True)


def pair_plot_strat_pop_platform(
    outlier_ht,
    strat_pop_platform_ht,
    pair_plot_metrics,
    read_if_exists=True,
    figsize=(30, 30),
    tmp_dir_prefix=tmp_prefix,
    strat_pops=None,
    plot_dir_prefix="",
    use_fig_if_exists=False,
    add_kde_plot=True,
    fail_filter_name="fail_strat_pop_platform",
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["figure.figsize"] = figsize
    strat_cutoff_dict = hl.eval(strat_pop_platform_ht.qc_metrics_stats)
    # print(strat_cutoff_dict)
    strat = list(strat_cutoff_dict.keys())
    # print(strat)
    if strat_pops is None:
        strat_pops = list(set([strata[0] for strata in strat]) - {None})
    strat_pops.sort()
    strat_platform = sorted(set([strata[1] for strata in strat]) - {None})

    # print(strat_pops)
    # print(strat_platform)

    for pop in strat_pops:
        outputs = [widgets.Output()]
        hts = []
        for platform in strat_platform:
            strata = (pop, platform)
            # print(strata)
            if strata not in strat_cutoff_dict:
                # print("Not in cutoff dict")
                continue
            cutoffs = strat_cutoff_dict[strata]
            ht_strata = outlier_ht.filter(
                (outlier_ht.pop == strata[0])
                & (outlier_ht.platform == strata[1])
                & (outlier_ht.r_insertion_deletion < 1)
            )
            ht_strata = ht_strata.repartition(10).checkpoint(
                f"{tmp_dir_prefix}strat_pop_platform_{strata[0]}_{strata[1]}.ht",
                _read_if_exists=read_if_exists,
                overwrite=not read_if_exists,
            )
            if ht_strata.count() > 1:
                outputs.append(widgets.Output())
                hts.append((ht_strata, cutoffs, strata))

        tab = widgets.Tab(children=outputs)
        tab.set_title(0, strata[0])

        for i, (strat_ht_cutoff, output_tab) in enumerate(zip(hts, outputs[1:])):
            ht_strata, cutoffs, strata = strat_ht_cutoff
            ht_strata = ht_strata.select(
                fail_filter_name, *[m[0] for m in pair_plot_metrics]
            )
            my_pair_plot(
                ht_strata,
                cutoffs,
                hue=fail_filter_name,
                output_tab=output_tab,
                title=f"Pop: {strata[0]}, Platform: {strata[1]}",
                file_name=(
                    f"{plot_dir_prefix}.strat_pop_platform.{strata[0]}_{strata[1]}"
                ),
                use_fig_if_exists=use_fig_if_exists,
                add_kde_plot=add_kde_plot
            )

            tab.set_title(i + 1, strata[1])

        display(tab)


def pair_plot_regress_pop_strat_platform(
    outlier_ht,
    regress_pop_strat_platform_ht,
    pair_plot_metrics,
    read_if_exists=True,
    figsize=(30, 30),
    tmp_dir_prefix=tmp_prefix,
    plot_dir_prefix="",
    use_fig_if_exists=False,
    add_kde_plot=True,
    fail_filter_name="fail_regress_pop_strat_platform",
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["figure.figsize"] = figsize
    strat_cutoff_dict = hl.eval(regress_pop_strat_platform_ht.qc_metrics_stats)
    strat = list(strat_cutoff_dict.keys())
    strat_platform = sorted(set([strata[0] for strata in strat]) - {None})

    renamed_metrics = [
        (
            "regress_pop_strat_platform|" + m[0] + "_residual",
            "regress_pop_strat_platform|" + m[1] + "_residual",
        )
        for m in pair_plot_metrics
    ]

    outputs = [widgets.Output()]
    hts = []
    for platform in strat_platform:
        strata = (platform,)
        if strata not in strat_cutoff_dict:
            continue
        cutoffs = strat_cutoff_dict[strata]
        cutoffs = {"regress_pop_strat_platform|" + m: cutoffs[m] for m in cutoffs}
        ht_strata = outlier_ht.filter(
            (outlier_ht.platform == strata[0]) & (outlier_ht.r_insertion_deletion < 1)
        )
        ht_strata = ht_strata.repartition(50).checkpoint(
            f"{tmp_dir_prefix}regress_pop_strat_platform_ht_{strata[0]}.ht",
            _read_if_exists=read_if_exists,
            overwrite=not read_if_exists,
        )
        if ht_strata.count() > 1:
            outputs.append(widgets.Output())
            hts.append((ht_strata, cutoffs, strata))

    tab = widgets.Tab(children=outputs)

    for i, (strat_ht_cutoff, output_tab) in enumerate(zip(hts, outputs[1:])):
        ht_strata, cutoffs, strata = strat_ht_cutoff
        ht_strata = ht_strata.select(
            fail_filter_name, *[m[0] for m in renamed_metrics]
        )
        my_pair_plot(
            ht_strata,
            cutoffs,
            hue=fail_filter_name,
            output_tab=output_tab,
            title=f"Platform: {strata[0]}",
            file_name=f"{plot_dir_prefix}.regress_pop_strat_platform.{strata[0]}",
            use_fig_if_exists=use_fig_if_exists,
            add_kde_plot=add_kde_plot,
        )

        tab.set_title(i + 1, f"{strata[0]}")

    display(tab)

def pair_plot_regress_pop_strat_platform_strat_pop(
    outlier_ht,
    regress_pop_strat_platform_ht,
    pair_plot_metrics,
    strat_pops=None,
    read_if_exists=True,
    figsize=(30, 30),
    tmp_dir_prefix=tmp_prefix,
    plot_dir_prefix="",
    use_fig_if_exists=False,
    add_kde_plot=True,
    fail_filter_name="regress_pop_strat_platform",
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["figure.figsize"] = figsize
    strat_cutoff_dict = hl.eval(regress_pop_strat_platform_ht.qc_metrics_stats)
    strat = list(strat_cutoff_dict.keys())

    if strat_pops is None:
        strat_pops = list(set([strata[0] for strata in strat]) - {None})
    strat_pops.sort()

    strat_platform = sorted(set([strata[0] for strata in strat]) - {None})

    renamed_metrics = [
        (
            "regress_pop_strat_platform_redo_regression|" + m[0] + "_residual",
            "regress_pop_strat_platform_redo_regression|" + m[1] + "_residual",
        )
        for m in pair_plot_metrics
    ]

    for pop in strat_pops:
        outputs = [widgets.Output()]
        hts = []
        for platform in strat_platform:
            strata = (pop, platform)
            cutoff_strata = (platform,)
            if cutoff_strata not in strat_cutoff_dict:
                continue
            cutoffs = strat_cutoff_dict[cutoff_strata]
            cutoffs = {"regress_pop_strat_platform_redo_regression|" + m: cutoffs[m] for m in cutoffs}
            ht_strata = outlier_ht.filter(
                (outlier_ht.pop == strata[0]) & (outlier_ht.platform == strata[1]) & (outlier_ht.r_insertion_deletion < 1)
            )
            ht_strata = ht_strata.repartition(50).checkpoint(
                f"{tmp_dir_prefix}regress_pop_strat_platform_ht_{strata[0]}.{strata[1]}.ht",
                _read_if_exists=read_if_exists,
                overwrite=not read_if_exists,
            )
            if ht_strata.count() > 1:
                outputs.append(widgets.Output())
                hts.append((ht_strata, cutoffs, strata))

        tab = widgets.Tab(children=outputs)

        for i, (strat_ht_cutoff, output_tab) in enumerate(zip(hts, outputs[1:])):
            ht_strata, cutoffs, strata = strat_ht_cutoff
            ht_strata = ht_strata.select(
                fail_filter_name, *[m[0] for m in renamed_metrics]
            )
            my_pair_plot(
                ht_strata,
                cutoffs,
                hue=fail_filter_name,
                output_tab=output_tab,
                title=f"Pop: {strata[0]}, Platform: {strata[1]}",
                file_name=f"{plot_dir_prefix}.regress_pop_strat_platform.{strata[0]}.{strata[1]}",
                use_fig_if_exists=use_fig_if_exists,
                add_kde_plot=add_kde_plot,
            )

            tab.set_title(i + 1, f"{strata[1]}")

        display(tab)




def pair_plot_regress_pop_platform(
    outlier_ht,
    regress_pop_platform_ht,
    pair_plot_metrics,
    read_if_exists=True,
    figsize=(30, 30),
    tmp_dir_prefix=tmp_prefix,
    plot_dir_prefix="",
    use_fig_if_exists=False,
    add_kde_plot=True,
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["figure.figsize"] = figsize
    strat_cutoff_dict = hl.eval(regress_pop_platform_ht.qc_metrics_stats)
    strat_cutoff_dict = {
        "regress_pop_platform|" + m: strat_cutoff_dict[m] for m in strat_cutoff_dict
    }

    renamed_metrics = [
        (
            "regress_pop_platform|" + m[0] + "_residual",
            "regress_pop_platform|" + m[1] + "_residual",
        )
        for m in pair_plot_metrics
    ]

    outputs = []
    hts = []

    cutoffs = {
        "regress_pop_platform|" + m: strat_cutoff_dict[m] for m in strat_cutoff_dict
    }
    ht_strata = outlier_ht.filter(outlier_ht.r_insertion_deletion < 1)
    ht_strata = ht_strata.checkpoint(
        f"{tmp_dir_prefix}regress_pop_platform_ht.ht",
        _read_if_exists=read_if_exists,
        overwrite=not read_if_exists,
    )

    ht_strata = ht_strata.select(
        "fail_regress_pop_platform", *[m[0] for m in renamed_metrics]
    )
    my_pair_plot(
        ht_strata,
        cutoffs,
        hue="fail_regress_pop_platform",
        file_name=f"{plot_dir_prefix}.regress_pop_platform",
        use_fig_if_exists=use_fig_if_exists,
        add_kde_plot=add_kde_plot,
    )


def get_hist_plots_regress_pop_strat_platform(
    outlier_sample_qc,
    cutoffs,
    prefix="",
    qc_metrics=list(outlier_metrics),
    figsize=(30, 10),
    bins=100,
    hist_only=False,
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize
    cols = ["pop", "platform"]
    key = "s"
    # colnames = cols + [key] + [f'{metric}' for metric in qc_metrics]
    # fail_colnames = cols + [key] + [f'{fail_prefix}{metric}' for metric in qc_metrics]
    colnames = cols + [f"{prefix}{metric}_residual" for metric in qc_metrics]
    colnames_rename = cols + [f"{metric}_residual" for metric in qc_metrics]
    fail_colnames = cols + [f"{prefix}fail_{metric}_residual" for metric in qc_metrics]
    sample_qc_pd = (
        outlier_sample_qc.select(*colnames)
        .rename(dict(zip(colnames, colnames_rename)))
        .to_pandas()
    )
    sample_qc_fail_pd = (
        outlier_sample_qc.select(*fail_colnames)
        .rename(dict(zip(fail_colnames, colnames_rename)))
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)
    plots = None
    tables_pop_platform = []
    tables_pop = []
    tables_platform = []

    def my_hist(data, x, color, cutoffs):
        ax = plt.gca()
        # ax = ax0.twinx()
        # ax.yaxis.tick_left()
        # ax.yaxis.label_position = "left"
        platform = data["platform"][0]
        ax.hist(data[x].dropna(), color=color, bins=bins)
        if (platform,) in cutoffs and x in cutoffs[(platform,)]:
            cutoffs = cutoffs[(platform,)][x]
            for intercept, c in [(cutoffs.lower, "red"), (cutoffs.upper, "purple")]:
                if (
                    (intercept != math.inf)
                    and (intercept != -math.inf)
                    and (intercept is not None)
                    and (intercept is not None)
                    and not math.isnan(intercept)
                ):
                    ax.axvline(intercept, color=c, linestyle="--")

    outputs = [widgets.Output()]
    for metric in qc_metrics:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)

    for i, metric in enumerate([f"{metric}_residual" for metric in qc_metrics]):
        curve_dict = {}
        if not hist_only:
            fail_table = (
                sample_qc_fail_pd.groupby(cols)[metric].value_counts().unstack().fillna(0)
            )
            # fail_table = fail_table.rename_axis(mapper="None")
            fail_table.columns = ["Pass", "Fail"]
            fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
            decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
            fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
            fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
            fail_table = fail_table.round(2)  # (decimals)
            tables_pop_platform.append((metric, fail_table))

            fail_table = (
                sample_qc_fail_pd.groupby(["pop"])[metric]
                .value_counts()
                .unstack()
                .fillna(0)
            )
            # fail_table = fail_table.rename_axis(mapper="None")
            fail_table.columns = ["Pass", "Fail"]
            fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
            decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
            fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
            fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
            fail_table = fail_table.round(2)  # (decimals)
            tables_pop.append((metric, fail_table))

            fail_table = (
                sample_qc_fail_pd.groupby(["platform"])[metric]
                .value_counts()
                .unstack()
                .fillna(0)
            )
            # fail_table = fail_table.rename_axis(mapper="None")
            fail_table.columns = ["Pass", "Fail"]
            fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
            decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
            fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
            fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
            fail_table = fail_table.round(2)  # (decimals)
            tables_platform.append((metric, fail_table))

        o = outputs[i + 1]
        g = sns.FacetGrid(
            sample_qc_pd, col="platform", sharex=False, sharey=False, col_wrap=4
        )
        g.map_dataframe(my_hist, x=metric, cutoffs=cutoffs)
        # plt.show()
        file_name = str(uuid.uuid4())
        plt.savefig(
            f"plot_pngs/{metric}_hist_plots_regress_pop_strat_platform_{file_name}.png",
            bbox_inches="tight",
        )
        fig = plt.gcf()
        plt.close(fig)
        o.append_display_data(
            Image(
                f"plot_pngs/{metric}_hist_plots_regress_pop_strat_platform_{file_name}.png"
            )
        )

        tab.set_title(i + 1, metric)

    display(tab)

    if not hist_only:
        return tables_pop_platform, tables_pop, tables_platform


def get_hist_plots_regress_pop_platform(
    outlier_sample_qc,
    cutoffs,
    prefix="",
    qc_metrics=list(outlier_metrics),
    figsize=(30, 10),
    bins=100,
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize
    cols = ["pop", "platform"]
    key = "s"
    # colnames = cols + [key] + [f'{metric}' for metric in qc_metrics]
    # fail_colnames = cols + [key] + [f'{fail_prefix}{metric}' for metric in qc_metrics]
    colnames = cols + [f"{prefix}{metric}_residual" for metric in qc_metrics]
    colnames_rename = cols + [f"{metric}_residual" for metric in qc_metrics]
    fail_colnames = cols + [f"{prefix}fail_{metric}_residual" for metric in qc_metrics]
    sample_qc_pd = (
        outlier_sample_qc.select(*colnames)
        .rename(dict(zip(colnames, colnames_rename)))
        .to_pandas()
    )
    sample_qc_fail_pd = (
        outlier_sample_qc.select(*fail_colnames)
        .rename(dict(zip(fail_colnames, colnames_rename)))
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)
    plots = None
    tables_pop_platform = []
    tables_pop = []
    tables_platform = []

    def my_hist(data, x, color, cutoffs):
        ax = plt.gca()
        # ax = ax0.twinx()
        # ax.yaxis.tick_left()
        # ax.yaxis.label_position = "left"
        m = data.iloc[0]["metric"]
        ax.hist(data[x].dropna(), bins=bins)  # , color=color)
        if m in cutoffs:
            cutoffs = cutoffs[m]
            for intercept, c in [(cutoffs.lower, "red"), (cutoffs.upper, "purple")]:
                if (
                    (intercept != math.inf)
                    and (intercept != -math.inf)
                    and (intercept is not None)
                    and not math.isnan(intercept)
                ):
                    ax.axvline(intercept, color=c, linestyle="--")

    for i, metric in enumerate([f"{metric}_residual" for metric in qc_metrics]):
        curve_dict = {}
        fail_table = (
            sample_qc_fail_pd.groupby(cols)[metric].value_counts().unstack().fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_pop_platform.append((metric, fail_table))

        fail_table = (
            sample_qc_fail_pd.groupby(["pop"])[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_pop.append((metric, fail_table))

        fail_table = (
            sample_qc_fail_pd.groupby(["platform"])[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_platform.append((metric, fail_table))

    colnames = [f"{prefix}{metric}_residual" for metric in qc_metrics]
    colnames_rename = [f"{metric}_residual" for metric in qc_metrics]
    sample_qc_pd = (
        outlier_sample_qc.select(*colnames)
        .rename(dict(zip(colnames, colnames_rename)))
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)
    sample_qc_pd_melted = sample_qc_pd.melt(var_name="metric")
    g = sns.FacetGrid(
        sample_qc_pd_melted, col="metric", sharex=False, sharey=False, col_wrap=4
    )
    g.map_dataframe(my_hist, x="value", cutoffs=cutoffs)
    # plt.show()
    file_name = str(uuid.uuid4())
    plt.savefig(
        f"plot_pngs/{metric}_hist_plots_regress_pop_platform_{file_name}.png",
        bbox_inches="tight",
    )
    fig = plt.gcf()
    plt.close(fig)
    Image(f"plot_pngs/{metric}_hist_plots_regress_pop_platform_{file_name}.png")

    return tables_pop_platform, tables_pop, tables_platform


# TODO: Fix residuals to work like above
def get_tables_only(
    outlier_sample_qc, fail_prefix="", residual=False, qc_metrics=list(outlier_metrics)
):
    cols = ["pop", "platform"]
    key = "s"
    # colnames = cols + [key] + [f'{metric}' for metric in qc_metrics]
    # fail_colnames = cols + [key] + [f'{fail_prefix}{metric}' for metric in qc_metrics]
    colnames = cols + [f"{metric}" for metric in qc_metrics]
    colnames_rename = cols + [f"{metric}_residual" for metric in qc_metrics]
    fail_colnames = cols + [
        f"{fail_prefix}{metric}{'_residual' if residual else ''}"
        for metric in qc_metrics
    ]
    sample_qc_pd = outlier_sample_qc.select(*colnames).to_pandas()
    sample_qc_fail_pd = (
        outlier_sample_qc.select(*fail_colnames)
        .rename(dict(zip(fail_colnames, colnames)))
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)
    plots = None
    tables_pop_platform = []
    tables_pop = []
    tables_platform = []

    for i, metric in enumerate(qc_metrics):
        curve_dict = {}
        fail_table = (
            sample_qc_fail_pd.groupby(cols)[metric].value_counts().unstack().fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_pop_platform.append((metric, fail_table))

        fail_table = (
            sample_qc_fail_pd.groupby(["pop"])[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_pop.append((metric, fail_table))

        fail_table = (
            sample_qc_fail_pd.groupby(["platform"])[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables_platform.append((metric, fail_table))

    return tables_pop_platform, tables_pop, tables_platform


def plot_pcs(ht: hl.Table, pca_scores_expr, pcs, label: str, pop_colormap,
             collect_all: bool):
    tabs = []
    for pc_i, pc_j in pcs:
        plot = hl.plot.scatter(
            pca_scores_expr[pc_i - 1],
            pca_scores_expr[pc_j - 1],
            xlabel=f'PC{pc_i}',
            ylabel=f'PC{pc_j}',
            label={f'{label}': ht[f'{label}'], },
            hover_fields={
                's': ht.s,
                'pop': ht.pop,
            },
            collect_all=collect_all,
            colors={f'{label}': pop_colormap, },
        )

        plot.xaxis.axis_label_text_font_size = "20pt"
        plot.yaxis.axis_label_text_font_size = "20pt"
        plot.legend.label_text_font_size = '20pt'

        tabs.append(
            Panel(
                title=f'PC{pc_i} vs PC{pc_j}',
                child=plot
            )
        )
    return Tabs(tabs=tabs)


def upset_plot(
    ht,
    cols_to_keep=None,
    filter_expr=None,
    min_subset_size=50,
    strip_prefix=True,
    figsize=(45, 10),
    output_tab=None,
    file_name="",
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    fig = plt.figure(figsize=figsize)
    # fig.canvas.width = '12in'
    # fig.canvas.height = '5in'
    plt.rcParams["figure.dpi"] = 300
    # plt.rcParams['figure.figsize'] = (80, 20)
    if cols_to_keep is None:
        cols_to_keep = list(ht.row.keys())
    if filter_expr is not None:
        cols_to_keep = list(filter(filter_expr, cols_to_keep))
    cols_to_keep_rename = cols_to_keep
    if strip_prefix:
        cols_to_keep_rename = []
        for k in cols_to_keep:
            cols_to_keep_rename.append(k.split("|")[-1])
    ht2 = ht.select(*cols_to_keep)
    ht2 = ht2.filter(hl.any(hl.array([ht2[x] for x in cols_to_keep])))
    meta_ht2_pd = ht2.rename(dict(zip(cols_to_keep, cols_to_keep_rename))).to_pandas()
    meta_ht2_pd_counts = from_indicators(
        cols_to_keep_rename, meta_ht2_pd.drop("s", axis=1)
    )

    plot(
        meta_ht2_pd_counts,
        show_counts=True,
        min_subset_size=min_subset_size,
        fig=fig,
        element_size=None,
    )
    # plt.show()
    if output_tab is None:
        my_stringIObytes = io.BytesIO()
        plt.savefig(my_stringIObytes, format="png", bbox_inches="tight")
        my_stringIObytes.seek(0)
        img_data = base64.b64encode(my_stringIObytes.read()).decode()
        plt.close(fig)
        head = """<html>
<style>
#upsetplot img{
    max-width: unset;
    align: left;
    display: block;
    padding-top: 50px;
    padding-right: 10px;
    padding-bottom: 50px;
    padding-left: 0px;
}
</style>
<div id="upsetplot">""" + """<img src="data:image/png;base64,{}"></div>
</html>
""".format(
        img_data
    )
        display_html(head, raw=True)
    else:
        plt.savefig(f"plot_pngs/upset_{file_name}.png", bbox_inches="tight")
        fig = plt.gcf()
        plt.close(fig)
        output_tab.append_display_data(
            Image(f"plot_pngs/upset_{file_name}.png")
        )

def get_nn_hists(
    outlier_ht,
    prefix1="nn_pop_strat_platform|",
    prefix2="regress_pop_strat_platform|",
    qc_metrics=list(outlier_metrics),
    figsize=(30, 10),
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    plt.rcParams["figure.figsize"] = figsize
    samples_keep_pass_both = {}
    samples_keep_fail_both = {}
    samples_keep_fail_residual = {}
    samples_keep_fail_nn = {}
    samples = outlier_ht.aggregate(
        hl.struct(
            pass_both={
                metric: hl.agg.filter(
                    ~outlier_ht[f"{prefix1}fail_{metric}"]
                    & ~outlier_ht[f"{prefix2}fail_{metric}_residual"],
                    hl.agg.collect(outlier_ht.s),
                )
                for metric in qc_metrics
            },
            fail_both={
                metric: hl.agg.filter(
                    outlier_ht[f"{prefix1}fail_{metric}"]
                    & outlier_ht[f"{prefix2}fail_{metric}_residual"],
                    hl.agg.collect(outlier_ht.s),
                )
                for metric in qc_metrics
            },
            fail_residual={
                metric: hl.agg.filter(
                    ~outlier_ht[f"{prefix1}fail_{metric}"]
                    & outlier_ht[f"{prefix2}fail_{metric}_residual"],
                    hl.agg.collect(outlier_ht.s),
                )
                for metric in qc_metrics
            },
            fail_nn={
                metric: hl.agg.filter(
                    outlier_ht[f"{prefix1}fail_{metric}"]
                    & ~outlier_ht[f"{prefix2}fail_{metric}_residual"],
                    hl.agg.collect(outlier_ht.s),
                )
                for metric in qc_metrics
            },
        )
    )

    for metric in qc_metrics:
        samples_keep_pass_both[metric] = random.sample(samples["pass_both"][metric], 5)
        samples_keep_fail_both[metric] = random.sample(samples["fail_both"][metric], 5)
        samples_keep_fail_residual[metric] = random.sample(
            samples["fail_residual"][metric], 5
        )
        samples_keep_fail_nn[metric] = random.sample(samples["fail_nn"][metric], 5)

    all_samples_keep = set([])
    by_sample = {}
    for k, v in samples_keep_pass_both.items():
        for s in v:
            by_sample[s] = (k, "pass_both")
            all_samples_keep.add(s)

    for k, v in samples_keep_fail_both.items():
        for s in v:
            by_sample[s] = (k, "fail_both")
            all_samples_keep.add(s)

    for k, v in samples_keep_fail_residual.items():
        for s in v:
            by_sample[s] = (k, "fail_residual")
            all_samples_keep.add(s)

    for k, v in samples_keep_fail_nn.items():
        for s in v:
            by_sample[s] = (k, "fail_nn")
            all_samples_keep.add(s)

    outlier_ht = outlier_ht.filter(hl.literal(all_samples_keep).contains(outlier_ht.s))
    outlier_ht = outlier_ht.repartition(1).checkpoint(
        f"{tmp_prefix}all_random_samples_temp.ht", overwrite=True
    )

    def get_hist(
        ht_pd,
        lower_cutoff,
        upper_cutoff,
        sample_metric_val,
        fail_val,
        old_fail_val,
        ax,
        metric,
    ):
        sns.histplot(data=ht_pd, x="comparison_qc_metrics", ax=ax)
        for intercept, c in [
            (sample_metric_val, "red"),
            (lower_cutoff, "blue"),
            (upper_cutoff, "purple"),
        ]:
            if (
                not pd.isna(intercept)
                and (intercept != math.inf)
                and (intercept != -math.inf)
                and (intercept is not None)
                and (intercept is not None)
                and not math.isnan(intercept)
            ):
                ax.axvline(intercept, color=c, linestyle="--")
        ax.xaxis.set_label_text(metric, visible=True)
        ax.title.set_text(
            f"{metric}{' FAIL NN' if fail_val else ''}{' FAIL Regression' if old_fail_val else ''}"
        )

    outputs = []
    outputs.append(widgets.Output())
    for metric in qc_metrics:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)

    hist_data_all = defaultdict(list)
    by_sample = hl.literal(by_sample)
    outlier_ht = outlier_ht.annotate(sample_info=by_sample[outlier_ht.s])
    outlier_ht = outlier_ht.annotate(
        metric=outlier_ht.sample_info[0],
        pass_fail_type=outlier_ht.sample_info[1],
    )
    # outlier_ht.describe()
    if "r_ti_tv_singleton" in qc_metrics:
        outlier_ht = outlier_ht.annotate(
            r_ti_tv_singleton=outlier_ht.r_ti_tv_singleton_nn
        )
    outlier_ht = outlier_ht.select(
        *qc_metrics,
        qc_metrics_stats=hl.coalesce(
            *[
                hl.or_missing(
                    m == outlier_ht.metric, outlier_ht[f"{prefix1}qc_metrics_stats"][m]
                )
                for m in qc_metrics
            ]
        ),
        comparison_qc_metrics=outlier_ht[f"{prefix1}comparison_qc_metrics"].map(
            lambda x: x[outlier_ht.metric]
        ),
        fail1=hl.coalesce(
            *[
                hl.or_missing(m == outlier_ht.metric, outlier_ht[f"{prefix1}fail_{m}"])
                for m in qc_metrics
            ]
        ),
        fail2=hl.coalesce(
            *[
                hl.or_missing(
                    m == outlier_ht.metric, outlier_ht[f"{prefix2}fail_{m}_residual"]
                )
                for m in qc_metrics
            ]
        ),
    )
    # outlier_ht.describe()
    outlier_ht = outlier_ht.explode(outlier_ht.comparison_qc_metrics)
    outlier_ht = outlier_ht.checkpoint(
        "gs://gnomad-tmp/julia/nn_example_plots_explode.ht", overwrite=True
    )
    ht_pd = outlier_ht.to_pandas()

    for i, metric in enumerate(qc_metrics):
        for k in range(5):
            hist_data_all[metric].append(
                [
                    ht_pd[ht_pd["s"] == samples_keep_pass_both[metric][k]],
                    ht_pd[ht_pd["s"] == samples_keep_fail_both[metric][k]],
                    ht_pd[ht_pd["s"] == samples_keep_fail_residual[metric][k]],
                    ht_pd[ht_pd["s"] == samples_keep_fail_nn[metric][k]],
                ]
            )

    for i, metric in enumerate(qc_metrics):
        hist_data = hist_data_all[metric]
        o = outputs[i + 1]
        fig, axs = plt.subplots(ncols=4, nrows=5, figsize=(25, 25))
        for k, df_row in enumerate(hist_data):
            for j, df in enumerate(df_row):
                print(df)
                upper_cutoff = df["qc_metrics_stats.upper"].iat[0]
                lower_cutoff = df["qc_metrics_stats.lower"].iat[0]
                sample_metric_val = df[metric].iat[0]
                fail_val = df["fail1"].iat[0]
                old_fail_val = df["fail2"].iat[0]
                get_hist(
                    df,
                    lower_cutoff,
                    upper_cutoff,
                    sample_metric_val,
                    fail_val,
                    old_fail_val,
                    ax=axs[k][j],
                    metric=metric,
                )
        plt.subplots_adjust(wspace=0.5, hspace=0.5)
        file_name = str(uuid.uuid4())
        plt.savefig(
            f"plot_pngs/{metric}_nn_example_{file_name}.png", bbox_inches="tight"
        )
        plt.close(fig)
        o.append_display_data(Image(f"plot_pngs/{metric}_nn_example_{file_name}.png"))
        tab.set_title(i + 1, metric)

    display(tab)


def filter_num_bar_plots(
    tables_pop_platform, figsize=(8, 60), order_platform_first=False
):
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    plt.rcParams["figure.dpi"] = 200
    if order_platform_first:
        label_name = "platform_pop"
    else:
        label_name = "pop_platform"

    def get_bar(df_mod, metric, o):
        fig, ax = plt.subplots(figsize=figsize)
        ax2 = sns.barplot(
            data=df_mod,
            x="Pct_fail",
            y=label_name,
            orient="h",
            ax=ax,
            hue="pop",
            dodge=False,
        )
        ax2.set(title=metric)
        all_labels = []
        for x in ax2.properties()["yticklabels"]:
            fail = df_mod.loc[df_mod[label_name] == x.get_text(), "Fail"].values[0]
            pas = df_mod.loc[df_mod[label_name] == x.get_text(), "Pass"].values[0]
            label = f"{fail} / {pas + fail}"
            all_labels.append(label)
        for c in ax2.containers:
            ax2.bar_label(c, labels=all_labels)

        # plt.show()
        file_name = str(uuid.uuid4())
        plt.savefig(f"plot_pngs/{metric}_bar_plot_{file_name}.png", bbox_inches="tight")
        plt.close(fig)
        o.append_display_data(Image(f"plot_pngs/{metric}_bar_plot_{file_name}.png"))

    outputs = []
    outputs.append(widgets.Output())
    for m in tables_pop_platform:
        outputs.append(widgets.Output())

    tab = widgets.Tab(children=outputs)

    for i, m in enumerate(tables_pop_platform):
        metric = m[0]
        df_mod = m[1].reset_index()
        if order_platform_first:
            df_mod[label_name] = df_mod["platform"] + " / " + df_mod["pop"]
        else:
            df_mod[label_name] = df_mod["pop"] + " / " + df_mod["platform"]
        df_mod = df_mod.sort_values(by=[label_name])
        o = outputs[i + 1]
        get_bar(df_mod, metric, o)
        tab.set_title(i + 1, metric)

    display(tab)
