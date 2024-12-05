"""Script to create release de novo HT for v4.0 exomes."""

import argparse
import logging
from datetime import datetime
from typing import Dict

import hail as hl
from gnomad.resources.grch38.gnomad import all_sites_an, public_release
from gnomad.sample_qc.relatedness import filter_to_trios
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import process_consequences

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.create_release.create_release_sites_ht import (
    custom_region_flags_select,
    get_config,
    get_final_ht_fields,
    get_freq_array_readme,
    join_hts,
    qc_temp_prefix,
)
from gnomad_qc.v4.create_release.create_release_utils import (
    DBSNP_VERSION,
    GENCODE_VERSION,
    MANE_SELECT_VERSION,
    POLYPHEN_VERSION,
    SIFT_VERSION,
)
from gnomad_qc.v4.resources.annotations import get_info, get_trio_stats
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.sample_qc import pedigree
from gnomad_qc.v4.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_de_novo_release_ht")
logger.setLevel(logging.INFO)

FINALIZED_SCHEMA = {
    "globals": [
        "exomes_freq_globals",
        "genomes_freq_globals",
        "joint_freq_globals",
        "exomes_an_globals",
        "filtering_model",
        "inbreeding_coeff_cutoff",
        "tool_versions",
        "vep_globals",
        "frequency_README",
        "date",
        "version",
    ],
    "rows": [
        "de_novo_stats",
        "de_novo_stats_raw",
        "exomes_freq",
        "genomes_freq",
        "joint_freq",
        "exomes_an",
        "a_index",
        "was_split",
        "rsid",
        "filters",
        "vep",
        "region_flags",
        "allele_info",
        "exomes_qual_hists",
        "in_silico_predictors",
    ],
}

TABLES_FOR_RELEASE = [
    "de_novos",
    "dbsnp",
    "filters",
    "info",
    "in_silico",
    "vep",
    "joint_freq",
    "exomes_an",
]


def filter_de_novos(ht: hl.Table) -> hl.Table:
    """Filter de novo stats Table."""
    ht = ht.filter(
        (ht.n_de_novos_raw > 0)
        & ~get_info().ht()[ht.key].AS_lowqual
        & (ht.alleles[1] != "*")
    )
    return ht


def densify_releasable_trios_mt(
    vds: hl.vds.VariantDataset, fam_ht: hl.Table, var_ht: hl.Table
) -> hl.matrixtable:
    """
    Densify the releasable trios MT with a selected set of variants.

    :param vds: Variant Dataset.
    :param fam_ht: Table containing family information.
    :param var_ht: Table containing variants to be densified.
    :return: Densified MT.
    """
    # Filter the family table to only the releasable trios.
    meta = vds.variant_data.cols()
    fam_ht = fam_ht.annotate(
        id_releasable=meta[fam_ht.key].meta.project_meta.releasable,
        pat_releasable=meta[fam_ht.pat_id].meta.project_meta.releasable,
        mat_releasable=meta[fam_ht.mat_id].meta.project_meta.releasable,
    )
    fam_ht = fam_ht.filter(
        fam_ht.id_releasable & fam_ht.pat_releasable & fam_ht.mat_releasable
    )

    # Filter the variant data and reference data to only the trios.
    vds = filter_to_trios(vds, fam_ht)

    mt = vds.variant_data
    mt = mt.filter_rows(hl.is_defined(var_ht.key_by("locus")[mt.locus]))
    mt = mt.checkpoint(hl.utils.new_temp_file("de_novo", "mt"))
    mt = hl.experimental.sparse_split_multi(mt)
    mt = mt.semi_join_rows(var_ht)
    mt = mt.checkpoint(hl.utils.new_temp_file("de_novo.split", "mt"))
    vds = hl.vds.VariantDataset(vds.reference_data, mt)
    mt = hl.vds.to_dense_mt(vds).naive_coalesce(100)

    return mt


def custom_de_novo_select(ht, **_):
    """Select fields from de novo stats Table."""
    # Filter out sites with n_de_novo_raw == 0.
    ht = ht.filter(ht.n_de_novo_raw > 0)
    return {
        f"de_novo_stats{'' if g == 'adj' else '_raw'}": hl.struct(
            n_de_novo=ht[f"n_de_novos_{g}"],
            parents=hl.struct(
                AC=ht[f"ac_parents_{g}"],
                AN=ht[f"an_parents_{g}"],
            ),
            children=hl.struct(
                AC=ht[f"ac_children_{g}"],
                AN=ht[f"an_children_{g}"],
            ),
        )
        for g in ["adj", "raw"]
    }


def custom_joint_freq_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
    """
    Select joint freq fields for release.

    :param ht: Joint freq Hail Table.
    :return: Select expression dict.
    """
    selects = {
        f"{data_type}_freq": ht[data_type].freq
        for data_type in ["exomes", "genomes", "joint"]
    }
    selects["exomes_qual_hists"] = ht.exomes.histograms.qual_hists

    return selects


def custom_joint_freq_select_globals(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select joint freq globals for release.

    :param ht: Joint freq Hail Table.
    :return: Select expression dict
    """
    return {
        f"{data_type}_freq_globals": ht[f"{data_type}_globals"].select(
            "freq_meta", "freq_index_dict", "freq_meta_sample_count"
        )
        for data_type in ["exomes", "genomes", "joint"]
    }


def custom_info_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
    """
    Select fields for info Hail Table annotation in de novo release.

    :param ht: Info Hail Table.
    :return: Select expression dict.
    """
    return {"allele_info": ht.allele_info.drop("nonsplit_alleles")}


def custom_an_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
    """
    Select AN for release.

    :param ht: Hail Table with AN.
    :return: Select expression dict.
    """
    selects = {"exomes_an": ht.AN, **custom_region_flags_select(ht, data_type="")}
    selects["region_flags"] = selects["region_flags"].annotate(
        **{r: ht[r] for r in ht.row_value if r.endswith("_region")}
    )

    return selects


def custom_an_select_globals(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select AN Table globals for release.

    :param ht: Hail Table with AN.
    :return: Select expression dict.
    """
    return {
        "exomes_an_globals": hl.struct(
            meta=ht.strata_meta,
            sample_count=ht.strata_sample_count,
        )
    }


def restructure_for_tsv(ht: hl.Table) -> hl.Table:
    """
    Restructure the de novo release HT for TSV export.

    :param ht: De novo release HT.
    :return: Restructured HT.
    """
    ht = ht.filter((ht.filters.length() == 0) & (ht.de_novo_stats.n_de_novo > 0))
    ht = process_consequences(ht, penalize_flags=False, has_polyphen=False)
    return (
        ht.select(
            "de_novo_stats",
            exomes=ht.exomes_freq[0],
            genomes=ht.genomes_freq[0],
            joint=ht.joint_freq[0],
            worst_csq_for_variant_canonical=ht.vep.worst_csq_for_variant_canonical.select(
                "most_severe_consequence",
                "transcript_id",
                "gene_id",
                "gene_symbol",
                "mane_select",
                "lof",
                "hgnc_id",
            ),
            **ht.in_silico_predictors,
        )
        .flatten()
        .rename(
            {
                "de_novo_stats.n_de_novo": "n_de_novo",
                "de_novo_stats.parents.AC": "trio_parents.AC",
                "de_novo_stats.parents.AN": "trio_parents.AN",
                "de_novo_stats.children.AC": "trio_children.AC",
                "de_novo_stats.children.AN": "trio_children.AN",
            }
        )
    )


def main(args):
    """Create de novo release ht."""
    hl.init(
        log=f"/create_de_novo_release_ht.log",
        tmp_dir="gs://gnomad-tmp-4day",
        default_reference="GRCh38",
    )
    data_type = "exomes"
    config = get_config(data_type="exomes")
    config["info"].update(
        {
            "select": ["was_split", "a_index"],
            "custom_select": custom_info_select,
        }
    )
    config.update(
        {
            "filters": {
                "ht": final_filter(all_variants=True, only_filters=True).ht(),
                "path": final_filter(all_variants=True, only_filters=True).path,
                "select": ["filters"],
            },
            "de_novos": {
                "ht": get_trio_stats(releasable_only=True).ht(),
                "path": get_trio_stats(releasable_only=True).path,
                "custom_select": custom_de_novo_select,
            },
            "joint_freq": {
                "ht": public_release("joint").ht(),
                "path": public_release("joint").path,
                "custom_select": custom_joint_freq_select,
                "custom_globals_select": custom_joint_freq_select_globals,
            },
            "exomes_an": {
                "ht": all_sites_an("exomes").ht(),
                "path": all_sites_an("exomes").path,
                "custom_select": custom_an_select,
                "custom_globals_select": custom_an_select_globals,
            },
        }
    )
    ht = join_hts(
        "de_novos",
        TABLES_FOR_RELEASE,
        "exomes",
        config,
        version=args.version,
        new_n_partitions=5 if args.test else 1000,
        checkpoint_tables=True,
        track_included_datasets=False,
        use_annotate=True,
        test=args.test,
    )
    ht.describe()

    # Filter out chrM, AS_lowqual, alt is '*' sites (these sites are dropped in the
    # final_filters HT so will not have information in `filters`)
    logger.info("Filtering out chrM, AS_lowqual and alt is '*'...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

    logger.info("Finalizing the release HT global and row fields...")
    # Add additional globals that were not present on the joined HTs.
    ht = ht.annotate_globals(
        vep_globals=ht.vep_globals.annotate(
            gencode_version=GENCODE_VERSION,
            mane_select_version=MANE_SELECT_VERSION,
        ),
        tool_versions=ht.tool_versions.annotate(
            dbsnp_version=DBSNP_VERSION,
            sift_version=SIFT_VERSION,
            polyphen_version=POLYPHEN_VERSION,
        ),
        date=datetime.now().isoformat(),
        version=CURRENT_RELEASE,
        frequency_README=get_freq_array_readme(data_type="exomes"),
    )

    # Organize the fields in the release HT to match the order of FINALIZED_SCHEMA when
    # the fields are present in the HT.
    final_fields = get_final_ht_fields(ht, schema=FINALIZED_SCHEMA)
    ht = ht.select(*final_fields["rows"]).select_globals(*final_fields["globals"])

    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/julia/denovo/denovo_releasable_only.release.intermediate.ht",
        # _read_if_exists=True,
        overwrite=True,
    )

    output_path = f"{qc_temp_prefix(data_type=data_type)}release/gnomad.exomes.de_novo.test.{datetime.today().strftime('%Y-%m-%d')}.ht"
    logger.info(f"Writing out de novo release HT to %s", output_path)
    ht = ht.naive_coalesce(args.n_partitions).checkpoint(
        output_path,
        args.overwrite,
    )

    logger.info("Final de novo release HT schema:")
    ht.describe()

    logger.info("Final de novo count: %d", ht.count())

    ht = restructure_for_tsv(ht)
    ht.describe()
    ht.export(
        f"{output_path.replace('.ht', '.tsv.bgz')}",
        header=True,
        delimiter="\t",
    )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "-v",
        "--version",
        help="The version of gnomAD.",
        default=CURRENT_RELEASE,
    )
    parser.add_argument(
        "-t",
        "--test",
        help="Runs a test on chr20",
        action="store_true",
    )
    parser.add_argument(
        "--generate-dense-mt",
        help="Generate a dense MT for releasable trios with the de novo variants.",
        action="store_true",
    )
    parser.add_argument(
        "-j",
        "--tables-for-join",
        help="Tables to join for release",
        default=TABLES_FOR_RELEASE,
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--base-table",
        help="Base table for interval partition calculation.",
        default="de_novos",
        choices=TABLES_FOR_RELEASE,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of partitions to naive coalesce the release Table to.",
        type=int,
        default=100,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
