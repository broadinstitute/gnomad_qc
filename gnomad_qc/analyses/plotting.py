import hail as hl
import numpy as np
from bokeh.plotting import figure, show
from bokeh.models import HoverTool, ColumnDataSource, Span
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


def priority_table_df(ht, n=50, exp_hom_col="group_total_expected_hom", ancestry=None):
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
    (+0.10), and OMIM recessive gene (+0.10). Ranks span the full depleted set.

    :param exp_hom_col: expected-homozygote field used for the depletion-strength
        axis and the "Exp hom" column (default group-total).
    :param ancestry: genetic-ancestry group to filter on (see
        :func:`_hom_depletion_field`); None uses overall depletion.
    """
    base_fields = [
        "most_severe_gene", "most_severe_csq", "gerp",
        "is_omim_recessive_gene", "phylop_score",
        "nfe_expected_hom", "afr_expected_hom", "amr_expected_hom",
        "eas_expected_hom", "asj_expected_hom", "fin_expected_hom", "sas_expected_hom",
    ]
    # exp_hom_col may itself be one of the per-ancestry fields; include it once.
    if exp_hom_col not in base_fields:
        base_fields.append(exp_hom_col)

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

    df["omim_flag"] = df["is_omim_recessive_gene"].map({True: "OMIM AR", False: ""})
    df["driving_ancestry"] = df[
        ["nfe_expected_hom","afr_expected_hom","amr_expected_hom",
         "eas_expected_hom","asj_expected_hom","fin_expected_hom","sas_expected_hom"]
    ].idxmax(axis=1).str.replace("_expected_hom","").str.upper()

    df = df.sort_values("priority", ascending=False).head(n).reset_index(drop=True)
    df["rank"] = df.index + 1

    # Round numeric columns to match the 2-decimal display in plot_priority_table.
    for col in [exp_hom_col, "gerp", "phylop_score", "loeuf", "priority"]:
        df[col] = df[col].round(2)

    # Select, reorder, and rename to exactly the columns shown in the plotted table.
    display_columns = [
        ("rank", "#"),
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
    df = df[[c for c, _ in display_columns]].rename(
        columns={c: title for c, title in display_columns}
    )
    return df


def plot_priority_table(ht, n=50, exp_hom_col="group_total_expected_hom", ancestry=None):
    """
    DataTable of top N depleted variants ranked by a weighted priority score.

    Priority = weighted sum of rank-normalized depletion strength (0.35),
    consequence severity (0.25), gene constraint / LOEUF (0.20), and conservation
    (0.10), plus additive bonuses for ClinVar pathogenic (+0.15), mouse-KO lethal
    (+0.10), and OMIM recessive gene (+0.10). Ranks span the full depleted set.
    Sortable interactively in the notebook.

    :param exp_hom_col: expected-homozygote field used for the depletion-strength
        axis and the "Exp hom" column (default group-total).
    """
    from bokeh.models import DataTable, TableColumn, NumberFormatter, StringFormatter, HTMLTemplateFormatter
    from bokeh.models import Div

    base_fields = [
        "most_severe_gene", "most_severe_csq", "gerp",
        "is_omim_recessive_gene", "phylop_score",
        "nfe_expected_hom", "afr_expected_hom", "amr_expected_hom",
        "eas_expected_hom", "asj_expected_hom", "fin_expected_hom", "sas_expected_hom",
    ]
    # exp_hom_col may itself be one of the per-ancestry fields; include it once.
    if exp_hom_col not in base_fields:
        base_fields.append(exp_hom_col)

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

    df["omim_flag"] = df["is_omim_recessive_gene"].map({True: "OMIM AR", False: ""})
    df["driving_ancestry"] = df[
        ["nfe_expected_hom","afr_expected_hom","amr_expected_hom",
         "eas_expected_hom","asj_expected_hom","fin_expected_hom","sas_expected_hom"]
    ].idxmax(axis=1).str.replace("_expected_hom","").str.upper()

    df = df.sort_values("priority", ascending=False).head(n).reset_index(drop=True)
    df["rank"] = df.index + 1

    src = ColumnDataSource(df)

    gerp_template = '<div style="color: <%= gerp >= 4 ? \'#A32D2D\' : gerp >= 2 ? \'#BA7517\' : gerp >= 0 ? \'#3B6D11\' : \'#888780\' %>; font-weight:500;"><%= gerp.toFixed(2) %></div>'
    omim_template = '<div style="<%= is_omim_recessive_gene ? \'font-weight:600;color:#8B0000;\' : \'\' %>"><%= omim_flag %></div>'

    columns = [
        TableColumn(field="rank",                  title="#",            width=40),
        TableColumn(field="gene",                  title="Gene",         width=90),
        TableColumn(field="locus",                 title="Locus",        width=140),
        TableColumn(field="alleles",               title="Alleles",      width=80),
        TableColumn(field="most_severe_csq",       title="Consequence",  width=160),
        TableColumn(field=exp_hom_col,             title="Exp hom",      width=80,
                    formatter=NumberFormatter(format="0.00")),
        TableColumn(field="gerp",                  title="GERP",         width=70,
                    formatter=HTMLTemplateFormatter(template=gerp_template)),
        TableColumn(field="phylop_score",          title="phyloP",       width=70,
                    formatter=NumberFormatter(format="0.00")),
        TableColumn(field="omim_flag",             title="OMIM",         width=75,
                    formatter=HTMLTemplateFormatter(template=omim_template)),
        TableColumn(field="clinvar_sig",           title="ClinVar",      width=150),
        TableColumn(field="clinvar_disease",       title="ClinVar disease", width=200),
        TableColumn(field="impc_viability",        title="Mouse KO",     width=110),
        TableColumn(field="loeuf",                 title="LOEUF",        width=70,
                    formatter=NumberFormatter(format="0.00")),
        TableColumn(field="priority",              title="Priority",     width=75,
                    formatter=NumberFormatter(format="0.00")),
        TableColumn(field="driving_ancestry",      title="Driving anc.", width=90),
    ]

    table = DataTable(
        source=src, columns=columns,
        width=1400, height=min(800, (n + 1) * 25 + 40),
        sortable=True, selectable=True,
        index_position=None,
    )

    header = Div(
        text=(
            f"<b>Priority ranking</b>{_ancestry_suffix(ancestry)} &mdash; top {n} depleted variants &nbsp;&nbsp;"
            "<span style='font-size:12px;color:#666;'>"
            "sorted by weighted priority (depletion + consequence + constraint + "
            "conservation, boosted by ClinVar/mouse-KO/OMIM) &nbsp;|&nbsp; "
            "click column headers to re-sort &nbsp;|&nbsp; "
            "OMIM AR genes in red</span>"
        ),
        width=1400,
    )

    return column(header, table)


# ── Cell 1: GERP plots ────────────────────────────────────────────────────────
show(column(
    plot_gerp_vs_expected_hom(ht),
    plot_top_candidates(ht, n=20),
    plot_gerp_distribution(ht),
))

# ── Cell 2 (run separately): Ancestry & OMIM ─────────────────────────────────
# show(column(
#     plot_ancestry_breakdown(ht, n=20),
#     plot_omim_scatter(ht),
# ))

# ── Cell 3 (run separately): Priority table ───────────────────────────────────
# show(plot_priority_table(ht, n=50))
