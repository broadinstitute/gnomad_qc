"""Script to create release de novo HT for v4 exomes."""

import argparse
import logging
from datetime import datetime
from typing import Dict

import hail as hl
from gnomad.resources.grch38.gnomad import all_sites_an, public_release
from gnomad.sample_qc.relatedness import default_get_de_novo_expr
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import process_consequences
from hail.utils import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.create_release.create_release_sites_ht import (
    custom_region_flags_select,
    get_config,
    get_final_ht_fields,
    get_freq_array_readme,
    join_hts,
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
from gnomad_qc.v4.resources.release import release_de_novo
from gnomad_qc.v4.resources.sample_qc import dense_trio_mt, pedigree, trio_denovo_ht
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
    """
    Filter the trio stats Table to de novo variants.

    Removes variants that are not de novo, AS_lowqual, and '*' alleles.

    :param ht: Trio stats Table.
    :return: Filtered Table.
    """
    ht = ht.filter(
        (ht.n_de_novos_raw > 0)
        & ~get_info().ht()[ht.key].AS_lowqual
        & (ht.alleles[1] != "*")
    )
    return ht


def get_releasable_de_novo_calls_ht(
    mt: hl.MatrixTable,
    priors_ht: hl.Table,
    ped: hl.Pedigree,
    test: bool = False,
) -> hl.Table:
    """
    Get de novo calls Hail Table.

    :param mt: Dense MatrixTable of the releasable trios.
    :param priors_ht: Table with AFs used as population frequency priors in de novo calculations.
    :param ped: Pedigree.
    :param test: Whether to filter to chr20 for testing. Default is False.
    :return: Hail Table with de novo calls.
    """
    if test:
        logger.info("Filtering to chr20 for testing...")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])

    # The split MT was made by adding these extra annotations .
    mt = mt.annotate_rows(
        n_unsplit_alleles=hl.len(mt.alleles),
        mixed_site=(
            (hl.len(mt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:])
        ),
    )
    mt = hl.experimental.sparse_split_multi(mt).checkpoint(
        dense_trio_mt(releasable=True, split=True, test=test).path, _read_if_exists=True
    )
    mt = annotate_adj(mt)

    # V3 and v4 VDS have the PL and AD fields for homref genotypes intentionally
    # removed to save storage space and costs. We need to approximate the
    # AD and PL fields when missing.
    mt = mt.annotate_entries(
        AD=hl.or_else(mt.AD, [mt.DP, 0]),
        PL=hl.or_else(mt.PL, [0, mt.GQ, 2 * mt.GQ]),
    )
    mt = mt.annotate_rows(prior=priors_ht[mt.row_key].freq[0].AF)
    mt = mt.select_entries("GT", "AD", "DP", "GQ", "PL", "adj")

    tm = hl.trio_matrix(mt, ped, complete_trios=True)
    tm = tm.transmute_cols(is_xx=tm.is_female)
    tm = tm.checkpoint(new_temp_file("trio_matrix", "mt"))

    ht = tm.entries()

    ht = ht.annotate(
        de_novo_call_info=default_get_de_novo_expr(
            locus_expr=ht.locus,
            alleles_expr=ht.alleles,
            proband_expr=ht.proband_entry,
            father_expr=ht.father_entry,
            mother_expr=ht.mother_entry,
            is_xx_expr=ht.is_xx,
            freq_prior_expr=ht.prior,
        )
    )

    ht = ht.filter(ht.de_novo_call_info.is_de_novo).naive_coalesce(1000)
    return ht


# TODO: Change when we decide on the final de novo data to release.
def custom_de_novo_select(ht, **_):
    """Select fields from de novo stats Table."""
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


def get_de_novo_config() -> Dict[str, dict]:
    """
    Get de novo release config.

    :return: Config dict.
    """
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
            # TODO: Change when we decide on the final de novo data to release.
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

    return config


def get_de_novo_release_ht(
    tables_for_join: list[str], version: str, test: bool = False
) -> hl.Table:
    """
    Get de novo release HT.

    :param tables_for_join: Tables to join for release.
    :param version: Version of gnomAD.
    :param test: Run test on chr20. Default is False.
    """
    config = get_de_novo_config()
    ht = join_hts(
        "de_novos",
        tables_for_join,
        "exomes",
        config,
        version=version,
        new_n_partitions=5 if test else 1000,
        checkpoint_tables=True,
        track_included_datasets=False,
        use_annotate=True,
        test=test,
    )

    # Filter out chrM, AS_lowqual, alt is '*' sites (these sites are dropped in the
    # final_filters HT so will not have information in `filters`).
    logger.info("Filtering out chrM, AS_lowqual and alt is '*'...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

    # TODO: Change when we decide on the final de novo data to release.
    ht = ht.filter(ht.de_novo_stats_raw.n_de_novo > 0)

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
    return ht.select(*final_fields["rows"]).select_globals(*final_fields["globals"])


def restructure_for_tsv(ht: hl.Table) -> hl.Table:
    """
    Restructure the de novo release HT for TSV export.

    :param ht: De novo release HT.
    :return: Restructured HT.
    """
    # TODO: Change when we decide on the final de novo data to release.
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
        # TODO: Change to be set outside init, used this to run with hail 0.2.122.
        default_reference="GRCh38",
    )

    test = args.test
    overwrite = args.overwrite

    if args.generate_de_novo_calls:
        mt = dense_trio_mt(releasable=True).mt()
        priors_ht = public_release("exomes").ht()
        ped = hl.Pedigree.read(pedigree().path, delimiter="\t")

        ht = get_releasable_de_novo_calls_ht(mt, priors_ht, ped, test)
        ht.write(trio_denovo_ht(releasable=True, test=test).path, overwrite=overwrite)

    if args.generate_final_ht:
        ht = get_de_novo_release_ht(args.tables_for_join, args.version, test=test)

        logger.info("Final de novo release HT schema:")
        ht.describe()

        ht = ht.checkpoint(new_temp_file("de_novo_release", "ht"))
        output_path = release_de_novo(test=test).path
        logger.info(f"Writing out de novo release HT to %s", output_path)
        ht = ht.naive_coalesce(args.n_partitions).checkpoint(output_path, overwrite)

        logger.info("Final de novo count: %d", ht.count())

    if args.generate_final_tsv:
        ht = release_de_novo(test=test).ht()
        ht = restructure_for_tsv(ht)
        ht.export(
            release_de_novo(test=test).path.replace(".ht", ".tsv.bgz"),
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
        "--generate-de-novo-calls",
        help="Generate de novo calls HT.",
        action="store_true",
    )
    parser.add_argument(
        "--de-novo-n-partitions",
        help="Number of partitions to naive coalesce the de novo dense MT to.",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--generate-final-ht",
        help="Generate the final HT for de novo release.",
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
        "--n-partitions",
        help="Number of partitions to naive coalesce the release Table to.",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--generate-final-tsv",
        help="Generate the final TSV for de novo release.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
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
