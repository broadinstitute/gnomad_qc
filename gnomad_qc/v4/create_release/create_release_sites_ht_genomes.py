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
    Get the samples in the HGDP + 1KG subset that were related to samples outside the subset in the v3.1 release and were not included in the v3.1 release.

    :param v3_meta_ht: Table with the v3.1 release metadata.
    :param rel_ht: Table with the v3.1 release relatedness, here we use the results
        from pc_relate of the full dataset.
    :return: Table with the samples in the HGDP + 1KG subset that are related to
        samples outside the subset in v3 release.
    """
    # Filter to samples that are either in the HGDP + 1KG subset or in the v3.1 release.
    v3_meta_ht = v3_meta_ht.annotate(
        hgdp_tgp=v3_meta_ht.subsets.tgp | v3_meta_ht.subsets.hgdp
    )
    v3_meta_ht = v3_meta_ht.filter(v3_meta_ht.release | v3_meta_ht.hgdp_tgp)

    rel_ht = rel_ht.annotate(
        **{
            f"{x}_meta": hl.struct(
                v3_release=hl.coalesce(v3_meta_ht[rel_ht[x].s].release, False),
                hgdp_tgp=hl.coalesce(v3_meta_ht[rel_ht[x].s].hgdp_tgp, False),
            )
            for x in ["i", "j"]
        }
    )

    rel_ht = rel_ht.filter(
        # Filter to pairs where at least one of the samples was in the v3.1 release.
        (rel_ht.i_meta.v3_release | rel_ht.j_meta.v3_release)
        # Filter to pairs with 2nd degree or closer relatedness.
        & (rel_ht.relationship != "unrelated")
        # Filter to pairs where one and only one of the samples is in the HGDP +
        # 1KG subset.
        & ((hl.int(rel_ht.i_meta.hgdp_tgp) + hl.int(rel_ht.j_meta.hgdp_tgp)) == 1)
        # Exclude pairs where the HGDP/1KG sample in the pair is also a v3.1 release
        # sample.
        & ~(rel_ht.i_meta.hgdp_tgp & rel_ht.i_meta.v3_release)
        & ~(rel_ht.j_meta.hgdp_tgp & rel_ht.j_meta.v3_release)
    )

    # Filter the relatedness HT twice to get the HGDP/1KG sample in each pair.
    ht1 = rel_ht.filter(rel_ht.i_meta.hgdp_tgp).key_by().select("i")
    ht1 = ht1.select(s=ht1.i.s)
    ht2 = rel_ht.filter(rel_ht.j_meta.hgdp_tgp).key_by().select("j")
    ht2 = ht2.select(s=ht2.j.s)

    # Merge the two HTs and remove the v3.1 prefix from the sample IDs.
    ht = ht1.union(ht2)
    ht = ht.select(s=ht.s.replace("v3.1::", "")).key_by("s").distinct()

    logger.info(
        "%d HGDP/1KG samples are related to samples outside the subset", ht.count()
    )

    return ht


def get_hgdp_tgp_v4_exome_duplicates(
    v3_meta_ht: hl.Table,
    v4_meta_ht: hl.Table,
    rel_ht: hl.Table,
) -> hl.Table:
    """
    Get the samples in the HGDP + 1KG subset that are duplicates of an exome in the v4 release.

    .. note::

        The duplicated samples are defined as samples that were in the v3.1.2
        subset release and are also in the v4 exomes release. The duplicated samples
        have to be removed because we will have combined frequencies from v4 exomes
        and genomes.

    :param v3_meta_ht: Table with the v3.1.2 release metadata.
    :param v4_meta_ht: Table with the v4.0 exomes release metadata.
    :param rel_ht: Table with the v3.1.2 and v4 joint relatedness, it's based on
        cuKING relatedness results.
    :return: Table with the samples in the HGDP + 1KG subset that are duplicates of an
        exome in the v4 exomes release.
    """
    # Get HGDP + 1KG subset samples in v3.1 meta.
    v3_meta_ht = v3_meta_ht.filter(v3_meta_ht.subsets.hgdp | v3_meta_ht.subsets.tgp)

    # Get release samples in v4.0 exomes.
    v4_meta_ht = v4_meta_ht.filter(v4_meta_ht.release)

    # Get samples that are in the v3.1 genomes and are also in the v4 exomes.
    rel_ht = rel_ht.filter(rel_ht.gnomad_v3_duplicate)

    # Check if the duplicates are still included in v4.0 exomes release and these
    # samples belong to the HGDP + 1KG subset.
    ht = rel_ht.annotate(
        **{
            f"{x}_meta": hl.struct(
                in_v4_exomes_release=hl.is_defined(v4_meta_ht[rel_ht[x].s]),
                in_hgdp_tgp_subset=hl.is_defined(v3_meta_ht[rel_ht[x].s]),
            )
            for x in ["i", "j"]
        }
    )

    ht = ht.filter(
        (ht.i_meta.in_v4_exomes_release & ht.j_meta.in_hgdp_tgp_subset)
        | (ht.j_meta.in_v4_exomes_release & ht.i_meta.in_hgdp_tgp_subset)
    )

    # Filter the relatedness HT twice to get the HGDP/1KG sample in each pair.
    ht1 = ht.filter(ht.i_meta.hgdp_tgp).key_by().select("i")
    ht1 = ht1.select(s=ht1.i.s)
    ht2 = ht.filter(ht.j_meta.hgdp_tgp).key_by().select("j")
    ht2 = ht2.select(s=ht2.j.s)

    # Merge the two HTs and remove the v3.1 prefix from the sample IDs.
    ht = ht1.union(ht2)
    ht = ht.select(s=ht.s.replace("v3.1::", "")).key_by("s").distinct()

    logger.info(
        "%d HGDP/1KG samples are duplicated in the v4 exomes release", ht.count()
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
    get_related_to_nonsubset = PipelineStepResourceCollection(
        "--get-related-to-nonsubset",
        input_resources={
            "v3.1 metadata": {"v3_meta_ht": v3_meta},
            "v3.1 relatedness": {"v3_relatedness_ht": v3_relatedness},
        },
        output_resources={"related_to_nonsubset_ht": hgdp_tgp_related_to_nonsubset},
    )
    hgdp_tgp_v4_exome_duplicates = PipelineStepResourceCollection(
        "--get-duplicated-to-exomes",
        input_resources={
            "v3.1 metadata": {"v3_meta_ht": v3_meta},
            "v4.0 metadata": {"v4_meta_ht": v4_meta},
            "v3.1 and v4.0 relatedness": {"v4_relatedness_ht": v4_relatedness()},
        },
        output_resources={"duplicated_to_exomes": hgdp_tgp_duplicated_to_exomes},
    )

    # Add all steps to the gnomAD v4 genomes release pipeline resource collection.
    v4_genome_release_pipeline.add_steps(
        {
            "get_related_to_nonsubset": get_related_to_nonsubset,
            "get_hgdp_tgp_v4_exome_duplicates": hgdp_tgp_v4_exome_duplicates,
        }
    )

    return v4_genome_release_pipeline


def main(args):
    """
    Create the v4.0 genomes release Table.

    Update call stats to include samples from the HGDP + 1KG subset that were
    unintentionally excluded whole populations within the HGDP + 1KG subset that were
    the most genetically unique and had small sample sizes (more specifically: San,
    Mbuti Pygmy, Biaka Pygmy, Bougainville, and Papuan) compared to other populations
    within the gnomAD v3.1 callset.

    In order to avoid re-calculating the call stats for the whole subset / whole
    release, we calculate the call stats for the samples that will be added
    and subtracted, then merge the call stats with the old call stats in the
    release HT.

    Changes compared to the v3 release:
      - Some small updates to samples that are hard filtered.
      - Use a population PC outlier approach to filter the HGDP + 1KG samples instead
        of the sample QC metric outlier filtering approach used on the full dataset
        that caused some samples to be unintentionally excluded.
      - HGDP + 1KG release samples are determined using relatedness (pc_relate) run on
        only samples within the subset as well as relatedness to the rest of the
        release.
      - Some population labels were updated.
      - The `Han` and `Papuan` populations were split into more specific groupings.
    """
    overwrite = args.overwrite

    hl.init(
        log="/create_release_v4_genomes.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    v4_genome_release_resources = get_v4_genomes_release_resources(overwrite=overwrite)

    if args.get_related_to_nonsubset:
        res = v4_genome_release_resources.get_related_to_nonsubset
        res.check_resource_existence()
        ht = get_hgdp_tgp_related_to_nonsubset(
            res.v3_meta_ht.ht(), res.v3_relatedness_ht.ht()
        )
        ht.write(res.related_to_nonsubset_ht.path, overwrite=overwrite)

    if args.get_hgdp_tgp_v4_exome_duplicates:
        res = v4_genome_release_resources.get_duplicated_to_exomes
        res.check_resource_existence()
        ht = get_hgdp_tgp_v4_exome_duplicates(
            res.v3_meta_ht.ht(), res.v4_meta_ht.ht(), res.v4_relatedness_ht.ht()
        )
        ht.write(res.hgdp_tgp_duplicated_to_exomes.path, overwrite=overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script creates v4 genomes release HT."
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
        "--get-hgdp-tgp-v4-exome-duplicates",
        help="Get the duplicated samples to exomes.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
