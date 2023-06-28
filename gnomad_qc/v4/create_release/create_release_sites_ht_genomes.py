"""Script to create release sites HT for v4 genomes."""
import hail as hl


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


def replace_oth_with_remaining(ht: hl.Table) -> hl.Table:
    """
    Replace 'oth' with 'remaining' in Global fields of a Table.

    .. note::
        This function is used to rename ancestry group fields to match v4 exomes.
        Which means that all "oth" will be changed to "remaining".
        2 Global fields are renamed: freq_meta (in the value of the key "pop") and freq_index_dict (in the keys).

    :param ht: release sites Table to be modified
    :return: release sites Table with 'oth' replaced with 'remaining'
    """
    # in freq_meta, type: array<dict<str, str>>, the value of key "pop" is
    # changed from 'oth' to 'remaining'.
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: x.replace("oth", "remaining"))
        )
    )

    # in freq_index_dict, type: dict<str, int>, the keys containing 'oth'
    # inside are changed to 'remaining'.
    oldkeys = hl.eval(ht.freq_index_dict.keys())
    newkeys = [s.replace("oth", "remaining") for s in oldkeys]
    vals = hl.eval(ht.freq_index_dict.values())
    freq_index_dict_new = {k: v for k, v in zip(newkeys, vals)}
    ht = ht.annotate_globals(freq_index_dict=freq_index_dict_new)
    return ht


# TODO: drop old_locus, old_alleles, seems to be already dropped in v3.1.4 release table
# TODO: rename ancestry group fields to match v4 exomes
# TODO: drop missing fields from VEP 105 annotations
# TODO: merge vep_ht with v3.1.4 release table
# TODO: rerun inbreeding coefficient with callstats instead of GT calls
# TODO: annotate with newest dbsnp release
# TODO: add Pangolin prediction
# TODO: add zoonomia constraint scores


def main(args):
    """Script to generate release sites ht for v4 genomes."""
    logger.info("Loading release table of v3.1.4...")
    # ht = release_sites(public=True).versions["3.1.4"].ht()
    ht = hl.read_table(
        "gs://gnomad/release/3.1.4/ht/genomes/gnomad.genomes.v3.1.4.sites.ht"
    )

    logger.info("Replacing ancestry group 'oth' with 'remaining'...")
    ht = replace_oth_with_remaining(ht)

    logger.info("Loading VEP-annotated table...")
    vep_ht = get_vep(version=CURRENT_VERSION, data_type="genomes").ht()

    logger.info("Removing missing VEP fields...")
    vep_ht = ht.annotate(vep=remove_missing_vep_fields(vep_ht.vep))

    logger.info("Loading dbsnp table...")
    dbsnp_ht = dbsnp.ht().select("rsid")

    logger.info("Loading pangolin table...")
    pangolin_ht = hl.read_table(
        "gs://gnomad/v4.0/annotations/genomes/gnomad.genomes.v4.0.pangolin.ht"
    )  # TODO: make Pangolin a versioned resource
