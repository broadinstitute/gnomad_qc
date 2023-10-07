"""Script to create release sites HT for v4 genomes."""

import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.basics import meta as v3_meta
from gnomad_qc.v3.resources.release import (
    hgdp_tgp_subset,
    hgdp_tgp_subset_annotations,
    release_sites,
)
from gnomad_qc.v3.resources.sample_qc import relatedness as v3_relatedness
from gnomad_qc.v4.resources.basics import meta as v4_meta
from gnomad_qc.v4.resources.sample_qc import (
    hgdp_tgp_duplicated_to_exomes,
    hgdp_tgp_related_to_nonsubset,
)
from gnomad_qc.v4.resources.sample_qc import relatedness as v4_relatedness

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(
    "Create v4 genomes release sites HT with updated HGDP/TGP "
    "metadata and new annotations"
)
logger.setLevel(logging.INFO)


def remove_missing_vep_fields(vep_expr: hl.StructExpression) -> hl.StructExpression:
    """
    Remove fields from VEP 105 annotations that have been excluded in past releases or are missing in all rows.

    :param vep_expr: StructExpression containing VEP 105 annotations.
    :return: StructExpression containing VEP 105 annotations with missing fields removed.
    """
    vep_expr = vep_expr.drop("colocated_variants", "context")

    vep_expr = vep_expr.annotate(
        transcript_consequences=vep_expr.transcript_consequences.map(
            lambda x: x.drop("minimised", "swissprot", "trembl", "uniparc")
        )
    )

    for consequence in [
        "intergenic_consequences",
        "motif_feature_consequences",
        "regulatory_feature_consequences",
    ]:
        vep_expr = vep_expr.annotate(
            **{consequence: vep_expr[consequence].map(lambda x: x.drop("minimised"))}
        )
    return vep_expr


def get_hgdp_tgp_related_to_nonsubset(
    v3_meta_ht: hl.Table, rel_ht: hl.Table
) -> hl.Table:
    """
    Get the samples in the HGDP + 1KG subset that were related to samples outside the subset in v3 release and were not included in the v3 release.

    :param v3_meta_ht: Table with the v3.1 release metadata
    :param rel_ht: Table with the v3.1 release relatedness, here we use the
       results from pc_relate of the full dataset.
    :return: Table with the samples in the HGDP + 1KG subset that are related
       to samples outside the subset in v3 release.
    """
    v3_meta_ht = v3_meta_ht.select(
        hgdp_tgp=v3_meta_ht.subsets.tgp | v3_meta_ht.subsets.hgdp,
        release=v3_meta_ht.release,
    ).select_globals()

    # Samples that are either in the HGDP + 1KG subset or in the v3 release.
    v3_meta_release_ht = v3_meta_ht.filter(v3_meta_ht.release | v3_meta_ht.hgdp_tgp)

    rel_ht = rel_ht.annotate(
        **{
            f"{x}_meta": hl.struct(
                v3_release=hl.coalesce(v3_meta_release_ht[rel_ht[x].s].release, False),
                hgdp_tgp=hl.coalesce(v3_meta_release_ht[rel_ht[x].s].hgdp_tgp, False),
            )
            for x in ["i", "j"]
        }
    )

    rel_ht = rel_ht.filter(
        # Filter to pairs where at least one of the samples was in the v3 release.
        (rel_ht.i_meta.v3_release | rel_ht.j_meta.v3_release)
        # Filter to pairs with 2nd degree or closer relatedness.
        & (rel_ht.relationship != "unrelated")
        # Filter to pairs where one and only one of the samples is in the HGDP +
        # 1KG subset.
        & ((hl.int(rel_ht.i_meta.hgdp_tgp) + hl.int(rel_ht.j_meta.hgdp_tgp)) == 1)
        # Exclude pairs where the HGDP/1KG sample in the pair is also a v3 release
        # sample.
        & ~(rel_ht.i_meta.hgdp_tgp & rel_ht.i_meta.v3_release)
        & ~(rel_ht.j_meta.hgdp_tgp & rel_ht.j_meta.v3_release)
    )

    rel_ht = rel_ht.naive_coalesce(1)

    rel_ht = rel_ht.checkpoint(
        new_temp_file("hgdp_tgp_related_temp", extension="ht"), overwrite=True
    )

    ht1 = rel_ht.filter(rel_ht.i_meta.hgdp_tgp).key_by().select("i")
    ht1 = ht1.select(s=ht1.i.s)

    ht2 = rel_ht.filter(rel_ht.j_meta.hgdp_tgp).key_by().select("j")
    ht2 = ht2.select(s=ht2.j.s)

    ht = ht1.union(ht2)
    ht = ht.select(s=ht.s.replace("v3.1::", "")).key_by("s").distinct()
    return ht


def get_hgdp_tgp_duplicated_to_exomes(
    v3_meta: hl.Table,
    v4_meta: hl.Table,
    rel_ht: hl.Table,
) -> hl.Table:
    """
    Get the samples in the HGDP + 1KG subset that were duplicated in the v4 exomes release.

    .. note::

        The duplicated samples are defined as samples that were in the v3.1.2
        subset release and are also in the v4 exomes release. The duplicated samples
        have to be removed because we will have combined frequencies from v4 exomes
        and genomes.

    :param v3_meta: Table with the v3.1.2 release metadata
    :param v4_meta: Table with the v4.0 exomes release metadata
    :param rel_ht: Table with the v3.1.2 and v4 joint relatedness, it's based on
        cuKing relatedness results.
    :return: Table with the samples in the HGDP + 1KG subset that are duplicated in
        the v4 exomes release.
    """
    # get hgdp + 1kg subset samples in v3 meta
    v3_meta = (
        v3_meta.filter(v3_meta.subsets.hgdp | v3_meta.subsets.tgp)
        .select()
        .select_globals()
    )

    # get release samples in v4.0 exomes
    v4_meta = v4_meta.filter(v4_meta.release).select().select_globals()

    # Get samples that are in the v3 genomes and are also in the v4 exomes
    rel_ht = rel_ht.filter(rel_ht.gnomad_v3_duplicate)

    # Check if the duplicates are still included in v4 exomes release and these
    # samples belong to the HGDP + 1KG subset.
    ht = rel_ht.annotate(
        in_v4_exomes_release=hl.is_defined(v4_meta[rel_ht.j.s]),
        in_hgdp_tgp_subset=hl.is_defined(v3_meta[rel_ht.i.s]),
    )

    ht = ht.filter(ht.in_v4_exomes_release & ht.in_hgdp_tgp_subset)
    ht = ht.key_by()
    ht = ht.select(s=ht.i.s)
    # in case some samples had a prefix in v3
    ht = ht.select(s=ht.s.replace("v3.1::", "")).key_by("s").distinct()
    logger.info(
        "%d HGDP/TGP samples are duplicated in the v4 exomes release", ht.count()
    )
    return ht


def get_v4_genomes_release_resources(overwrite: bool) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed to create the gnomAD v4 genomes release.

    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        gnomAD v4 genomes release pipeline.
    """
    # Initialize gnomAD v4 genomes release pipeline resource collection.
    hgdp_tgp_res = {
        "meta_ht": hgdp_tgp_subset_annotations(sample=True).versions["3.1.2"],
        "dense_mt": hgdp_tgp_subset(dense=True, public=True).versions["3.1.2"],
        "sites_ht": release_sites(public=True).versions["3.1.2"],
    }
    v4_genome_release_pipeline = PipelineResourceCollection(
        pipeline_name="gnomad_v4_genomes_release",
        overwrite=overwrite,
        pipeline_resources={"Released HGDP + 1KG resources": hgdp_tgp_res},
    )

    # Create resource collection for each step of the v4 genomes release pipeline.

    # Add all steps to the gnomAD v4 genomes release pipeline resource collection.
    v4_genome_release_pipeline.add_steps({})

    return v4_genome_release_pipeline


def main(args):
    """
    Script to update call stats for HGDP + 1KG subset for v4.

    This code is specifically designed for the update HGDP + 1KG subset in
    v4 release HT. There are a few major changes compared to the v3 release:
        - the following new sample filters are applied: hard filters, pop PC outliers,
        relatedness within the subset and relatedness to the rest of the release.
      - the new pop labels.
      - the new split of the `Han` and `Papuan` samples.

    In order to avoid re-calculating the callstats for the whole subset / whole
    release, we will calculate the callstats for the samples that will be added
    and subtracted, then merge the callstats with the old callstats in the
    release HT.
    """
    overwrite = args.overwrite

    hl.init(
        log="/create_release_v4_genomes.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    v4_genome_release_resources = get_v4_genomes_release_resources(overwrite=overwrite)

    if args.get_related_to_nonsubset:
        ht = get_hgdp_tgp_related_to_nonsubset(v3_meta.ht(), v3_relatedness.ht())
        ht.write(hgdp_tgp_related_to_nonsubset.path, overwrite=overwrite)

    if args.get_duplicated_to_exomes:
        ht = get_hgdp_tgp_duplicated_to_exomes(
            v3_meta.ht(), v4_meta.ht(), v4_relatedness().ht()
        )
        ht.write(hgdp_tgp_duplicated_to_exomes.path, overwrite=overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script updates AFs for HGDP + 1KG subset for v4 release HT."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--get-related-to-nonsubset",
        help="Get the relatedness to nonsubset samples.",
        action="store_true",
    )
    parser.add_argument(
        "--get-duplicated-to-exomes",
        help="Get the duplicated samples to exomes.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
