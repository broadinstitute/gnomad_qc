import hail as hl
import numpy as np
from bokeh.plotting import figure, show
from bokeh.models import (HoverTool, ColumnDataSource, Span, DataTable,
                          TableColumn, HTMLTemplateFormatter, Div)
from bokeh.layouts import column
from bokeh.transform import factor_cmap
from bokeh.palettes import Category10


ANCESTRIES = ["nfe", "afr", "amr", "eas", "asj", "fin", "sas"]


def _hom_depletion_field(ht, ancestry=None):
    """Return the hom-depletion boolean field to filter on.

    :param ancestry: genetic-ancestry group prefix (e.g. "nfe"). When None,
        uses the overall depletion flag ``overall_hom_depletion``; otherwise
        uses ``{ancestry}_hom_depletion``.
    :raises ValueError: if the resolved field is not present on the table.
    """
    field = "overall_hom_depletion" if ancestry is None else f"{ancestry}_hom_depletion"
    if field not in ht.row:
        raise ValueError(
            f"Field {field!r} not found on the table. "
            f"Pass ancestry=None for overall, or one of {ANCESTRIES}."
        )
    return ht[field]


def _ancestry_suffix(ancestry=None):
    """Human-readable title suffix for the selected ancestry (empty for overall)."""
    return "" if ancestry is None else f" ({ancestry})"


def plot_gerp_vs_expected_hom(ht, exp_hom_col="group_total_expected_hom", ancestry=None):
    """
    Scatter plot of GERP score vs expected homozygotes.
    Points are coloured by GERP tier. Hover shows gene name.

    :param exp_hom_col: expected-homozygote field to plot (default group-total).
    """
    df = (
        ht
        .filter(_hom_depletion_field(ht, ancestry))
        .select("most_severe_gene", exp_hom_col, "gerp")
        .to_pandas()
        .rename(columns={"most_severe_gene": "gene"})
    )

    def gerp_color(g):
        if g >= 4:  return "firebrick"
        if g >= 2:  return "darkorange"
        if g >= 0:  return "olivedrab"
        return "gray"

    df["color"] = df["gerp"].apply(gerp_color)
    df["gerp_tier"] = df["gerp"].apply(lambda g:
        "GERP ≥ 4" if g >= 4 else
        "GERP 2–4" if g >= 2 else
        "GERP 0–2" if g >= 0 else "GERP < 0"
    )

    hover = HoverTool(tooltips=[
        ("Gene",     "@gene"),
        ("Exp hom",  f"@{exp_hom_col}{{0.2f}}"),
        ("GERP",     "@gerp{0.2f}"),
        ("Tier",     "@gerp_tier"),
    ])

    p = figure(
        title=f"GERP score vs expected homozygotes (depleted variants){_ancestry_suffix(ancestry)}",
        x_axis_label=f"Expected homozygotes ({exp_hom_col})",
        y_axis_label="GERP score",
        tools=[hover, "pan", "wheel_zoom", "box_zoom", "reset", "save"],
        width=750, height=450,
    )

    for y_val, color in [(2, "olivedrab"), (4, "firebrick")]:
        p.add_layout(Span(location=y_val, dimension="width",
                          line_color=color, line_dash="dashed",
                          line_width=1, line_alpha=0.5))

    for tier, color in [
        ("GERP ≥ 4", "firebrick"),
        ("GERP 2–4", "darkorange"),
        ("GERP 0–2", "olivedrab"),
        ("GERP < 0", "gray"),
    ]:
        sub = df[df["gerp_tier"] == tier]
        src = ColumnDataSource(sub)
        p.scatter(
            x=exp_hom_col, y="gerp",
            source=src,
            size=8, alpha=0.75,
            color=color,
            legend_label=tier,
        )

    p.legend.click_policy = "hide"
    p.legend.location = "top_right"
    return p


def plot_top_candidates(ht, n=20, exp_hom_col="group_total_expected_hom", ancestry=None):
    """
    Horizontal bar chart of top N depleted variants by expected hom.
    Bars are coloured by GERP tier.

    :param exp_hom_col: expected-homozygote field to rank/plot (default group-total).
    """
    _ht = ht.filter(_hom_depletion_field(ht, ancestry)).key_by()
    _ht = _ht.select(
        "most_severe_gene", exp_hom_col, "gerp",
        locus_str=hl.str(_ht.locus),
        alleles_str=hl.delimit(_ht.alleles, "/"),
    )
    df = (
        _ht.to_pandas()
        .rename(columns={"most_severe_gene": "gene", "locus_str": "locus", "alleles_str": "alleles"})
        .dropna(subset=["gene"])  # variants with no gene can't form a y-axis label
        .sort_values(exp_hom_col, ascending=False)
        .head(n)
        .reset_index(drop=True)
    )

    def gerp_color(g):
        if g >= 4:  return "firebrick"
        if g >= 2:  return "darkorange"
        if g >= 0:  return "olivedrab"
        return "gray"

    df["color"] = df["gerp"].apply(gerp_color)
    df["label"] = df["gene"] + "  " + df["locus"].astype(str)
    df = df.iloc[::-1].reset_index(drop=True)

    src = ColumnDataSource(df)

    hover = HoverTool(tooltips=[
        ("Gene",    "@gene"),
        ("Locus",   "@locus"),
        ("Alleles", "@alleles"),
        ("Exp hom", f"@{exp_hom_col}{{0.2f}}"),
        ("GERP",    "@gerp{0.2f}"),
    ])

    p = figure(
        title=f"Top {n} candidates by expected homozygotes{_ancestry_suffix(ancestry)}",
        x_axis_label=f"Expected homozygotes ({exp_hom_col})",
        y_range=df["label"].tolist(),
        tools=[hover, "pan", "wheel_zoom", "reset", "save"],
        width=750, height=max(350, n * 26),
    )

    p.hbar(
        y="label", right=exp_hom_col,
        height=0.7, source=src,
        color="color", alpha=0.8,
        line_color="white", line_width=0.5,
    )

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    return p


def plot_gerp_distribution(ht, bins=50, ancestry=None):
    """
    Histogram of GERP scores for all depleted variants.
    Reference lines mark the conserved (2) and highly conserved (4) thresholds.
    """
    df = (
        ht
        .filter(_hom_depletion_field(ht, ancestry))
        .select("gerp")
        .to_pandas()
    )

    bin_edges = np.linspace(df["gerp"].min(), df["gerp"].max(), bins + 1)
    hist, edges = np.histogram(df["gerp"].values, bins=bin_edges)

    def gerp_color(left):
        if left >= 4:  return "firebrick"
        if left >= 2:  return "darkorange"
        if left >= 0:  return "olivedrab"
        return "gray"

    src = ColumnDataSource(dict(
        top=hist,
        left=edges[:-1],
        right=edges[1:],
        color=[gerp_color(e) for e in edges[:-1]],
        bin_left=[f"{e:.2f}" for e in edges[:-1]],
        bin_right=[f"{e:.2f}" for e in edges[1:]],
    ))

    hover = HoverTool(tooltips=[
        ("Range", "@bin_left – @bin_right"),
        ("Count", "@top"),
    ])

    p = figure(
        title=f"GERP score distribution (depleted variants){_ancestry_suffix(ancestry)}",
        x_axis_label="GERP score",
        y_axis_label="Count",
        tools=[hover, "pan", "wheel_zoom", "reset", "save"],
        width=750, height=400,
    )

    p.quad(
        source=src,
        top="top", bottom=0,
        left="left", right="right",
        color="color", alpha=0.8,
        line_color="white", line_width=0.5,
    )

    for x_val, label, color in [
        (2, "GERP = 2", "olivedrab"),
        (4, "GERP = 4", "firebrick"),
    ]:
        p.add_layout(Span(location=x_val, dimension="height",
                          line_color=color, line_dash="dashed",
                          line_width=1.5, line_alpha=0.7))

    return p


def plot_ancestry_breakdown(ht, n=20, exp_hom_col="group_total_expected_hom", ancestry=None):
    """
    Stacked horizontal bar chart of per-ancestry expected homozygotes
    for the top N depleted variants. Each bar segment is one ancestry.
    Hover tooltips work because every segment uses a ColumnDataSource
    that carries all tooltip columns.

    :param exp_hom_col: expected-homozygote field used to rank candidates
        (default group-total, which is the sum of the per-ancestry segments).
    """
    ancestries = ["nfe", "afr", "amr", "eas", "asj", "fin", "sas"]
    anc_colors = {
        "nfe": "#1D9E75", "afr": "#D85A30", "amr": "#BA7517",
        "eas": "#185FA5", "asj": "#993556", "fin": "#639922", "sas": "#888780",
    }
    anc_labels = {
        "nfe": "NFE", "afr": "AFR", "amr": "AMR",
        "eas": "EAS", "asj": "ASJ", "fin": "FIN", "sas": "SAS",
    }

    select_fields = ["most_severe_gene"] + [a + "_expected_hom" for a in ancestries]
    # exp_hom_col may itself be one of the per-ancestry fields; include it once.
    if exp_hom_col not in select_fields:
        select_fields.append(exp_hom_col)
    _ht = ht.filter(_hom_depletion_field(ht, ancestry)).key_by()
    _ht = _ht.select(
        *select_fields,
        locus_str=hl.str(_ht.locus),
    )
    df = (
        _ht.to_pandas()
        .rename(columns={"most_severe_gene": "gene", "locus_str": "locus"})
        .dropna(subset=["gene"])  # variants with no gene can't form a y-axis label
        .sort_values(exp_hom_col, ascending=False)
        .head(n)
        .reset_index(drop=True)
    )

    df["label"] = df["gene"] + "  " + df["locus"]
    df = df.iloc[::-1].reset_index(drop=True)

    labels = df["label"].tolist()

    hover = HoverTool(tooltips=[
        ("Gene",         "@gene"),
        ("Locus",        "@locus"),
        ("Total exp",    f"@{exp_hom_col}{{0.2f}}"),
        ("NFE",          "@nfe_expected_hom{0.2f}"),
        ("AFR",          "@afr_expected_hom{0.2f}"),
        ("AMR",          "@amr_expected_hom{0.2f}"),
        ("EAS",          "@eas_expected_hom{0.2f}"),
        ("ASJ",          "@asj_expected_hom{0.2f}"),
        ("FIN",          "@fin_expected_hom{0.2f}"),
        ("SAS",          "@sas_expected_hom{0.2f}"),
    ])

    p = figure(
        title=f"Ancestry breakdown — top {n} candidates (expected homozygotes){_ancestry_suffix(ancestry)}",
        x_axis_label="Expected homozygotes",
        y_range=labels,
        tools=[hover, "pan", "wheel_zoom", "reset", "save"],
        width=800, height=max(350, n * 26),
    )

    # FIX: use numpy arrays for arithmetic and wrap each segment in its own
    # ColumnDataSource so that @column_name references resolve in the tooltip.
    left_arr = np.zeros(len(df))
    for anc in ancestries:
        col_name = anc + "_expected_hom"
        right_arr = left_arr + df[col_name].values

        seg_src = ColumnDataSource(dict(
            y=labels,
            left=left_arr.tolist(),
            right=right_arr.tolist(),
            # carry every tooltip column on every segment source
            gene=df["gene"].tolist(),
            locus=df["locus"].tolist(),
            **{exp_hom_col: df[exp_hom_col].tolist()},
            **{a + "_expected_hom": df[a + "_expected_hom"].tolist() for a in ancestries},
        ))

        p.hbar(
            y="y", left="left", right="right",
            height=0.7,
            source=seg_src,
            color=anc_colors[anc],
            alpha=0.85,
            legend_label=anc_labels[anc],
            line_color="white", line_width=0.5,
        )
        left_arr = right_arr

    p.legend.click_policy = "hide"
    p.legend.location = "bottom_right"
    p.xgrid.grid_line_color = None
    return p


def plot_omim_scatter(ht, exp_hom_col="group_total_expected_hom", ancestry=None):
    """
    Scatter of GERP vs expected hom, with OMIM recessive genes
    highlighted as diamonds and non-OMIM as circles.

    :param exp_hom_col: expected-homozygote field to plot (default group-total).
    """
    _ht = ht.filter(_hom_depletion_field(ht, ancestry)).key_by()
    _ht = _ht.select(
        "most_severe_gene", exp_hom_col, "gerp", "is_omim_recessive_gene",
        locus_str=hl.str(_ht.locus),
    )
    df = (
        _ht.to_pandas()
        .rename(columns={"most_severe_gene": "gene", "locus_str": "locus"})
    )

    omim_df  = df[df["is_omim_recessive_gene"]].reset_index(drop=True)
    non_omim = df[~df["is_omim_recessive_gene"]].reset_index(drop=True)

    hover = HoverTool(tooltips=[
        ("Gene",    "@gene"),
        ("Locus",   "@locus"),
        ("Exp hom", f"@{exp_hom_col}{{0.2f}}"),
        ("GERP",    "@gerp{0.2f}"),
        ("OMIM AR", "@is_omim_recessive_gene"),
    ])

    p = figure(
        title=f"GERP vs expected hom — OMIM recessive genes highlighted{_ancestry_suffix(ancestry)}",
        x_axis_label=f"Expected homozygotes ({exp_hom_col})",
        y_axis_label="GERP score",
        tools=[hover, "pan", "wheel_zoom", "box_zoom", "reset", "save"],
        width=750, height=450,
    )

    for y_val, color in [(2, "olivedrab"), (4, "firebrick")]:
        p.add_layout(Span(location=y_val, dimension="width",
                          line_color=color, line_dash="dashed",
                          line_width=1, line_alpha=0.5))

    p.scatter(
        x=exp_hom_col, y="gerp",
        source=ColumnDataSource(non_omim),
        size=7, alpha=0.5, color="steelblue",
        legend_label="Non-OMIM",
    )
    p.scatter(
        x=exp_hom_col, y="gerp",
        marker="diamond",
        source=ColumnDataSource(omim_df),
        size=14, alpha=0.9, color="firebrick",
        line_color="darkred", line_width=0.8,
        legend_label="OMIM AR gene",
    )

    p.legend.click_policy = "hide"
    p.legend.location = "top_right"
    return p


def priority_table_df(ht, n=50, exp_hom_col=None, ancestry=None):
    """
    Build the ranked priority DataFrame (data only, no plotting).

    Same data as :func:`plot_priority_table`, returned as a pandas DataFrame so
    it can be exported. To load into Google Sheets, write a CSV and import it::

        df = priority_table_df(ht, n=50)
        df.to_csv("priority_table.csv", index=False)
        # then in Google Sheets: File → Import → Upload priority_table.csv

    Priority = weighted sum of rank-normalized depletion strength (0.35),
    consequence severity (0.25), gene constraint / LOEUF (0.20), and conservation
    (0.10), plus additive bonuses for ClinVar pathogenic (+0.15), mouse-KO lethal
    (+0.10), and OMIM recessive gene (+0.10), minus a penalty for likely SV/CNV
    artifacts (-0.30, applied only when likely_sv_artifact is present). Ranks span
    the full depleted set.

    :param exp_hom_col: expected-homozygote field used for the depletion-strength
        axis and the "Exp hom" column. Defaults to the selected ancestry's
        ``{ancestry}_expected_hom`` (or ``group_total_expected_hom`` for overall).
    :param ancestry: genetic-ancestry group to filter on (see
        :func:`_hom_depletion_field`); None uses overall depletion.
    """
    base_fields = [
        "most_severe_gene", "most_severe_csq", "gerp",
        "is_omim_recessive_gene", "phylop_score",
        "nfe_expected_hom", "afr_expected_hom", "amr_expected_hom",
        "eas_expected_hom", "asj_expected_hom", "fin_expected_hom", "sas_expected_hom",
    ]
    # Default the depletion metric to the selected ancestry's expected-hom field
    # (or the across-ancestry total for overall) unless the caller set it.
    if exp_hom_col is None:
        exp_hom_col = "group_total_expected_hom" if ancestry is None else f"{ancestry}_expected_hom"

    # exp_hom_col may itself be one of the per-ancestry fields; include it once.
    if exp_hom_col not in base_fields:
        base_fields.append(exp_hom_col)

    # likely_sv_artifact is added downstream by filter_sv_cnv_artifacts(remove=False);
    # include the column only when present so the export still works without it.
    has_sv_flag = "likely_sv_artifact" in ht.row
    if has_sv_flag:
        base_fields.append("likely_sv_artifact")

    _ht = ht.filter(_hom_depletion_field(ht, ancestry)).key_by()
    _ht = _ht.select(
        *base_fields,
        locus_str=hl.str(_ht.locus),
        alleles_str=hl.delimit(_ht.alleles, "/"),
        # Flatten nested annotations into table-friendly scalars/strings.
        loeuf=_ht.constraint.lof.oe_ci_upper,
        clinvar_sig=hl.delimit(_ht.clinvar.clinical_significance, "; "),
        clinvar_disease=hl.delimit(_ht.clinvar.disease, "; "),
        impc_viability=hl.delimit(hl.array(_ht.impc.viability), "; "),
    )
    df = (
        _ht.to_pandas()
        .rename(columns={"most_severe_gene": "gene", "locus_str": "locus", "alleles_str": "alleles"})
    )

    # --- Priority score: weighted sum of evidence axes on a common [0, 1] scale,
    #     plus additive bonuses for corroborating biological flags. Ranks are
    #     computed across the full depleted set (before the top-n slice below). ---
    lof_csqs = {
        "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
        "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    }
    moderate_csqs = {
        "missense_variant", "inframe_insertion", "inframe_deletion",
        "protein_altering_variant",
    }

    def csq_weight(c):
        if c in lof_csqs:       return 1.0
        if c in moderate_csqs:  return 0.6
        return 0.2

    df["s_csq"] = df["most_severe_csq"].map(csq_weight)
    df["s_expected"] = df[exp_hom_col].rank(pct=True)
    df["s_conserved"] = df[["gerp", "phylop_score"]].rank(pct=True).mean(axis=1)
    df["s_constraint"] = (1 - df["loeuf"].rank(pct=True)).fillna(0.0)  # low LOEUF = more constrained

    df["is_clinvar_path"] = df["clinvar_sig"].str.contains(
        "pathogenic", case=False, na=False
    ) & ~df["clinvar_sig"].str.contains("conflicting", case=False, na=False)
    df["is_ko_lethal"] = df["impc_viability"].str.contains("lethal", case=False, na=False)

    df["priority"] = (
        0.35 * df["s_expected"]
        + 0.25 * df["s_csq"]
        + 0.20 * df["s_constraint"]
        + 0.10 * df["s_conserved"]
        + 0.15 * df["is_clinvar_path"]
        + 0.10 * df["is_ko_lethal"]
        + 0.10 * df["is_omim_recessive_gene"].fillna(False)
    )
    # Demote likely SV/CNV artifacts: their depletion is probably technical rather
    # than biological, so subtract a penalty on the order of the depletion weight.
    if has_sv_flag:
        df["priority"] = df["priority"] - 0.30 * df["likely_sv_artifact"].fillna(False)

    # Global priority standing across the full depleted set (before any top-n
    # slice), so it stays meaningful even when rows are shown in locus order.
    df["priority_rank"] = df["priority"].rank(ascending=False, method="first").astype(int)

    df["omim_flag"] = df["is_omim_recessive_gene"].map({True: "OMIM AR", False: ""})
    if has_sv_flag:
        df["sv_flag"] = df["likely_sv_artifact"].map({True: "True", False: ""})
    df["driving_ancestry"] = df[
        ["nfe_expected_hom","afr_expected_hom","amr_expected_hom",
         "eas_expected_hom","asj_expected_hom","fin_expected_hom","sas_expected_hom"]
    ].idxmax(axis=1).str.replace("_expected_hom","").str.upper()

    df = df.sort_values("priority", ascending=False).head(n).reset_index(drop=True)

    # Round numeric columns to match the 2-decimal display in plot_priority_table.
    for col in [exp_hom_col, "gerp", "phylop_score", "loeuf", "priority"]:
        df[col] = df[col].round(2)

    # Select, reorder, and rename to exactly the columns shown in the plotted table.
    display_columns = [
        ("priority_rank", "Priority rank"),
        ("gene", "Gene"),
        ("locus", "Locus"),
        ("alleles", "Alleles"),
        ("most_severe_csq", "Consequence"),
        (exp_hom_col, "Exp hom"),
        ("gerp", "GERP"),
        ("phylop_score", "phyloP"),
        ("omim_flag", "OMIM"),
        ("clinvar_sig", "ClinVar"),
        ("clinvar_disease", "ClinVar disease"),
        ("impc_viability", "Mouse KO"),
        ("loeuf", "LOEUF"),
        ("priority", "Priority"),
        ("driving_ancestry", "Driving anc."),
    ]
    if has_sv_flag:
        display_columns.insert(-1, ("sv_flag", "likely SV artifact"))
    df = df[[c for c, _ in display_columns]].rename(
        columns={c: title for c, title in display_columns}
    )
    return df


def plot_priority_table(ht, n=50, exp_hom_col=None,
                        ancestry=None,
                        gnomad_dataset="gnomad_r4",
                        sort_by="locus",            # "locus" (enables shading) or "priority"
                        neighborhood_gap=100_000,   # bp gap that starts a new neighborhood
                        orange_thresh=100, gray_thresh=50):
    """
    DataTable of depleted variants. With sort_by="locus" (default), rows are
    ordered genomically and shaded by neighborhood: a run of variants whose loci
    are within neighborhood_gap bp is one neighborhood, and if its strongest
    variant has Exp hom > orange_thresh the whole block is sand-colored, > gray_thresh
    gray. A tight cluster of hom-depleted variants is the paralog/segdup/mismapping
    signature, so shaded blocks flag suspect regions and the rarer variants riding
    along in them are guilty by association. sort_by="priority" restores the old
    weighted-priority ranking (no shading).

    Priority = weighted sum of rank-normalized depletion strength (0.35),
    consequence severity (0.25), gene constraint / LOEUF (0.20), and conservation
    (0.10), plus additive bonuses for ClinVar pathogenic (+0.15), mouse-KO lethal
    (+0.10), and OMIM recessive gene (+0.10), minus a penalty for likely SV/CNV
    artifacts (-0.30, applied only when likely_sv_artifact is present). Ranks span
    the full depleted set.

    :param exp_hom_col: expected-homozygote field for the depletion axis, the
        "Exp hom" column, and the neighborhood shading threshold. Defaults to the
        selected ancestry's ``{ancestry}_expected_hom`` (or ``group_total_expected_hom``).
    :param ancestry: genetic-ancestry group to filter on (see
        :func:`_hom_depletion_field`); None uses overall depletion.
    :param n: rows to show. With sort_by="locus", None shows all depleted variants
        (recommended, so interspersed rare variants aren't dropped from blocks).
    """
    base_fields = [
        "most_severe_gene", "most_severe_csq", "gerp",
        "is_omim_recessive_gene", "phylop_score",
        "nfe_expected_hom", "afr_expected_hom", "amr_expected_hom",
        "eas_expected_hom", "asj_expected_hom", "fin_expected_hom", "sas_expected_hom",
    ]
    # Default the depletion metric to the selected ancestry's expected-hom field
    # (or the across-ancestry total for overall) unless the caller set it.
    if exp_hom_col is None:
        exp_hom_col = "group_total_expected_hom" if ancestry is None else f"{ancestry}_expected_hom"

    if exp_hom_col not in base_fields:
        base_fields.append(exp_hom_col)

    # likely_sv_artifact is added downstream by filter_sv_cnv_artifacts(remove=False);
    # include the column only when present so the table still works without it.
    has_sv_flag = "likely_sv_artifact" in ht.row
    if has_sv_flag:
        base_fields.append("likely_sv_artifact")

    _ht = ht.filter(_hom_depletion_field(ht, ancestry)).key_by()
    _ht = _ht.select(
        *base_fields,
        locus_str=hl.str(_ht.locus),
        alleles_str=hl.delimit(_ht.alleles, "/"),
        loeuf=_ht.constraint.lof.oe_ci_upper,
        clinvar_sig=hl.delimit(_ht.clinvar.clinical_significance, "; "),
        clinvar_disease=hl.delimit(_ht.clinvar.disease, "; "),
        impc_viability=hl.delimit(hl.array(_ht.impc.viability), "; "),
    )
    df = (
        _ht.to_pandas()
        .rename(columns={"most_severe_gene": "gene", "locus_str": "locus", "alleles_str": "alleles"})
    )

    df["variant_id"] = (
        df["locus"].str.replace("chr", "", regex=False).str.replace(":", "-")
        + "-" + df["alleles"].str.replace("/", "-")
    )

    lof_csqs = {
        "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
        "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    }
    moderate_csqs = {
        "missense_variant", "inframe_insertion", "inframe_deletion",
        "protein_altering_variant",
    }

    def csq_weight(c):
        if c in lof_csqs:       return 1.0
        if c in moderate_csqs:  return 0.6
        return 0.2

    df["s_csq"] = df["most_severe_csq"].map(csq_weight)
    df["s_expected"] = df[exp_hom_col].rank(pct=True)
    df["s_conserved"] = df[["gerp", "phylop_score"]].rank(pct=True).mean(axis=1)
    df["s_constraint"] = (1 - df["loeuf"].rank(pct=True)).fillna(0.0)

    df["is_clinvar_path"] = df["clinvar_sig"].str.contains(
        "pathogenic", case=False, na=False
    ) & ~df["clinvar_sig"].str.contains("conflicting", case=False, na=False)
    df["is_ko_lethal"] = df["impc_viability"].str.contains("lethal", case=False, na=False)

    df["priority"] = (
        0.35 * df["s_expected"]
        + 0.25 * df["s_csq"]
        + 0.20 * df["s_constraint"]
        + 0.10 * df["s_conserved"]
        + 0.15 * df["is_clinvar_path"]
        + 0.10 * df["is_ko_lethal"]
        + 0.10 * df["is_omim_recessive_gene"].fillna(False)
    )
    # Demote likely SV/CNV artifacts: their depletion is probably technical rather
    # than biological, so subtract a penalty on the order of the depletion weight.
    if has_sv_flag:
        df["priority"] = df["priority"] - 0.30 * df["likely_sv_artifact"].fillna(False)

    # Global priority standing across the full depleted set (before any top-n
    # slice), so it stays meaningful even when rows are shown in locus order.
    df["priority_rank"] = df["priority"].rank(ascending=False, method="first").astype(int)

    df["omim_flag"] = df["is_omim_recessive_gene"].map({True: "OMIM AR", False: ""})
    if has_sv_flag:
        df["sv_flag"] = df["likely_sv_artifact"].map({True: "True", False: ""})
    df["driving_ancestry"] = df[
        ["nfe_expected_hom","afr_expected_hom","amr_expected_hom",
         "eas_expected_hom","asj_expected_hom","fin_expected_hom","sas_expected_hom"]
    ].idxmax(axis=1).str.replace("_expected_hom","").str.upper()

    # --- order + neighborhood shading ---
    if sort_by == "locus":
        import numpy as np
        p = df["locus"].str.split(":", n=1, expand=True)
        df["_chrom"], df["_pos"] = p[0], pd.to_numeric(p[1], errors="coerce")

        def contig_key(c):
            c = str(c).replace("chr", "")
            return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(c, int(c) if c.isdigit() else 99)

        df["_ck"] = df["_chrom"].map(contig_key)
        df = df.sort_values(["_ck", "_pos"], kind="mergesort").reset_index(drop=True)

        gap = df["_pos"] - df["_pos"].shift()
        new_nbhd = (df["_chrom"] != df["_chrom"].shift()) | (gap > neighborhood_gap) | gap.isna()
        nbhd = new_nbhd.cumsum()
        nbhd_max = df.groupby(nbhd)[exp_hom_col].transform("max")
        df["row_color"] = np.select(
            [nbhd_max > orange_thresh, nbhd_max > gray_thresh],
            ["#FBEADB", "#F0F0F0"], default="white")   # soft sand / faint gray
        if n is not None:
            df = df.head(n).reset_index(drop=True)
        df = df.drop(columns=["_chrom", "_pos", "_ck"])
    else:  # priority
        df = df.sort_values("priority", ascending=False)
        if n is not None:
            df = df.head(n)
        df = df.reset_index(drop=True)
        df["row_color"] = "white"

    # Bokeh's HTMLTemplateFormatter skips the template entirely when a cell value
    # is null/NaN (`if (value == null) return ""`), so the row-color wrapper below
    # is never emitted and the cell renders bare white. Text columns are missing
    # for most variants (no ClinVar/OMIM/Mouse-KO entry), which broke the shading
    # across those columns. Fill blanks with a non-breaking space: it's non-null so
    # the template runs, and unlike "" it doesn't collapse, so the cell gets a
    # full-height line box and the fill covers it edge to edge.
    text_cols = ["gene", "locus", "alleles", "most_severe_csq", "omim_flag",
                 "clinvar_sig", "clinvar_disease", "impc_viability", "driving_ancestry"]
    if has_sv_flag:
        text_cols.append("sv_flag")
    for _c in text_cols:
        df[_c] = df[_c].fillna(" ").replace("", " ")

    src = ColumnDataSource(df)

    # cell background = row_color; wrap every column's template with this shell.
    # negative margins eat the cell's default padding so fills touch edge-to-edge,
    # box-shadow bleeds the color over the row border so no hairline gap shows.
    def shaded(inner):
        return ('<div style="background:<%= row_color %>;'
                'margin:-6px -8px;padding:6px 8px;height:100%;'
                'box-shadow:0 0 0 1px <%= row_color %>;">' + inner + '</div>')

    gerp_inner = ('color:<%= gerp >= 4 ? \'#A32D2D\' : gerp >= 2 ? \'#BA7517\' : '
                  'gerp >= 0 ? \'#3B6D11\' : \'#888780\' %>;font-weight:500;')
    gerp_template = shaded('<span style="' + gerp_inner + '">'
                           '<%= (gerp == null || isNaN(gerp)) ? "&#160;" : gerp.toFixed(2) %></span>')
    omim_template = shaded('<span style="<%= is_omim_recessive_gene ? '
                           '\'font-weight:600;color:#8B0000;\' : \'\' %>"><%= omim_flag %></span>')
    sv_template = shaded('<span style="<%= likely_sv_artifact ? '
                         '\'font-weight:600;color:#6A1B9A;\' : \'\' %>"><%= sv_flag %></span>')

    lof_js = "[" + ",".join(f"'{c}'" for c in sorted(lof_csqs)) + "]"
    csq_template = shaded('<span style="<%= ' + lof_js + '.indexOf(value) > -1 '
                          '? \'color:#A32D2D;font-weight:600;\' : \'\' %>"><%= value %></span>')

    gnomad_template = shaded(
        '<a href="https://gnomad.broadinstitute.org/variant/'
        '<%= variant_id %>?dataset=' + gnomad_dataset + '" target="_blank"><%= value %></a>')

    def num_template(field, fmt="toFixed(2)"):   # shaded numeric cell
        # missing numerics arrive as NaN (not null); show a non-breaking space so
        # the cell is blank instead of "NaN" yet still keeps its row-color fill.
        return shaded(f'<%= (value == null || isNaN(value)) ? "&#160;" : value.{fmt} %>')

    def txt_template():                          # shaded plain-text cell
        return shaded('<%= value == null ? "" : value %>')

    def col(field, title, width, tmpl):
        return TableColumn(field=field, title=title, width=width,
                           formatter=HTMLTemplateFormatter(template=tmpl))

    columns = [
        col("priority_rank",   "Priority rank",   90,  num_template("priority_rank", "toFixed(0)")),
        col("gene",            "Gene",            90,  txt_template()),
        col("locus",           "Locus",           140, gnomad_template),
        col("alleles",         "Alleles",         80,  txt_template()),
        col("most_severe_csq", "Consequence",     160, csq_template),
        col(exp_hom_col,       "Exp hom",         80,  num_template(exp_hom_col)),
        col("gerp",            "GERP",             70,  gerp_template),
        col("phylop_score",    "phyloP",           70,  num_template("phylop_score")),
        col("omim_flag",       "OMIM",            75,  omim_template),
        col("clinvar_sig",     "ClinVar",         150, txt_template()),
        col("clinvar_disease", "ClinVar disease", 200, txt_template()),
        col("impc_viability",  "Mouse KO",        110, txt_template()),
        col("loeuf",           "LOEUF",            70,  num_template("loeuf")),
        col("priority",        "Priority",        75,  num_template("priority")),
        col("driving_ancestry","Driving anc.",    90,  txt_template()),
    ]
    if has_sv_flag:
        columns.insert(-1, col("sv_flag", "likely SV artifact", 115, sv_template))

    table = DataTable(
        source=src, columns=columns,
        width=1400, height=min(800, (len(df) + 1) * 25 + 40),
        sortable=True, selectable=True, index_position=None,
    )

    if sort_by == "locus":
        note = (f"sorted by locus &nbsp;|&nbsp; neighborhoods (&le;{neighborhood_gap//1000}kb "
                f"apart) shaded by max Exp hom: <span style='background:#FBEADB;'>&gt;{orange_thresh}</span> "
                f"<span style='background:#F0F0F0;'>&gt;{gray_thresh}</span>")
    else:
        note = "sorted by weighted priority (depletion + consequence + constraint + conservation)"
    sv_legend = (" &nbsp;|&nbsp; <span style='color:#6A1B9A;font-weight:600;'>likely SV artifact</span>"
                 if has_sv_flag else "")
    header = Div(
        text=(f"<b>Priority ranking</b>{_ancestry_suffix(ancestry)} &mdash; {len(df)} depleted variants &nbsp;&nbsp;"
              "<span style='font-size:12px;color:#666;'>" + note +
              " &nbsp;|&nbsp; click a locus to open gnomAD &nbsp;|&nbsp; "
              "LoF consequences &amp; OMIM AR in red" + sv_legend + "</span>"),
        width=1400,
    )
    return column(header, table)


def export_priority_xlsx(ht, gs_path, n=50, exp_hom_col=None, ancestry=None,
                         gnomad_dataset="gnomad_r4", neighborhood_gap=100_000,
                         orange_thresh=100, gray_thresh=50):
    """
    Write the locus-shaded priority table to a styled .xlsx in a GCS bucket.

    Reuses :func:`plot_priority_table` (``sort_by="locus"``) for the fully-computed,
    genomically-ordered frame, then reproduces its coloring in Excel: neighborhood
    row shading (sand/gray) plus per-cell text colors -- LoF consequence & OMIM AR
    in red, GERP on its conservation scale, likely SV artifacts in purple. Locus
    cells are clickable gnomAD links. Requires ``openpyxl``.

    :param gs_path: destination, e.g. "gs://bucket/priority_table.xlsx".
    :param n: rows to include (None = all depleted variants).
    :param exp_hom_col, ancestry, gnomad_dataset, neighborhood_gap, orange_thresh,
        gray_thresh: forwarded to :func:`plot_priority_table`.
    :return: the ``gs_path`` written.
    """
    import io

    # Reuse the plot's exact locus-sorted, shaded, flagged frame.
    layout = plot_priority_table(
        ht, n=n, exp_hom_col=exp_hom_col, ancestry=ancestry,
        gnomad_dataset=gnomad_dataset, sort_by="locus",
        neighborhood_gap=neighborhood_gap,
        orange_thresh=orange_thresh, gray_thresh=gray_thresh,
    )
    table = next(m for m in layout.references() if isinstance(m, DataTable))
    raw = pd.DataFrame(dict(table.source.data))

    # Same columns/titles/order the plot shows.
    ordered = [(c.field, c.title) for c in table.columns]
    display = raw[[f for f, _ in ordered]].copy()

    # Clickable gnomAD locus links (Excel HYPERLINK formula).
    base_url = "https://gnomad.broadinstitute.org/variant/"
    display["locus"] = [
        f'=HYPERLINK("{base_url}{vid}?dataset={gnomad_dataset}", "{loc}")'
        for vid, loc in zip(raw["variant_id"], raw["locus"])
    ]

    # Round numeric columns; blank the non-breaking-space fills.
    num_fields = [f for f in ("gerp", "phylop_score", "loeuf", "priority") if f in display]
    exp_field = next((f for f, t in ordered if t == "Exp hom"), None)
    if exp_field:
        num_fields.append(exp_field)
    for f in num_fields:
        display[f] = pd.to_numeric(display[f], errors="coerce").round(2)
    display = display.replace("\xa0", "")

    display = display.rename(columns={f: t for f, t in ordered})

    lof_csqs = {
        "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
        "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    }

    def _row_bg(row):
        color = raw.loc[row.name, "row_color"]
        css = "" if color == "white" else f"background-color: {color}"
        return [css] * len(row)

    def _gerp_css(v):
        if pd.isna(v):  return ""
        if v >= 4:      return "color: #A32D2D; font-weight: 600"
        if v >= 2:      return "color: #BA7517; font-weight: 600"
        if v >= 0:      return "color: #3B6D11; font-weight: 600"
        return "color: #888780"

    styler = display.style.apply(_row_bg, axis=1)
    styler = styler.apply(lambda s: [_gerp_css(v) for v in s], subset=["GERP"])
    styler = styler.apply(
        lambda s: ["color: #A32D2D; font-weight: 600" if v in lof_csqs else "" for v in s],
        subset=["Consequence"],
    )
    styler = styler.apply(
        lambda s: ["color: #8B0000; font-weight: 600" if v == "OMIM AR" else "" for v in s],
        subset=["OMIM"],
    )
    if "likely SV artifact" in display.columns:
        styler = styler.apply(
            lambda s: ["color: #6A1B9A; font-weight: 600" if v == "True" else "" for v in s],
            subset=["likely SV artifact"],
        )

    buf = io.BytesIO()
    styler.to_excel(buf, engine="openpyxl", index=False)
    with hl.hadoop_open(gs_path, "wb") as f:
        f.write(buf.getvalue())
    return gs_path