# noqa: D100

import argparse
import logging
from datetime import datetime
from functools import reduce

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: Update gnomAD_methods to dict: key version, value subsets
)
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.utils.annotations import region_flag_expr
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import AS_FIELDS, SITE_FIELDS

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import (
    analyst_annotations,
    get_freq,
    get_info,
    vep,
)
from gnomad_qc.v4.resources.basics import qc_temp_prefix
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.release import FIELD_DESCRIPTIONS, release_sites
from gnomad_qc.v4.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)

# Remove InbreedingCoeff from allele-specific fields and SOR from SITE
# (processed separately from other fields)
AS_FIELDS.remove("InbreedingCoeff")
SITE_FIELDS.remove("SOR")
<<<<<<< HEAD

"""
Configurations of FIELD DESCRIPTION dictionary.
Format:
'Globals fields': {
    '<Name of global annotation>': {
            'definition': 'Explanation of global annotation.',
            'subfields': '<Dict of subfields and their definitions.>',
        },
'Row fields': {
    '<Name of row annotation>': {
            'definition': 'Explanation of row annotation.',
            'subfields': '<Dict of subfields and their definitions.>',
        },
"""

FIELD_DESCRIPTIONS = {
    "Global fields": {
        "freq_meta": {
            "definition": (
                "Allele frequency metadata. An ordered list containing the frequency"
                " aggregation group for each element of the ‘freq’ array row"
                " annotation."
            )
        },
        "freq_index_dict": {
            "defintion": (
                "Dictionary keyed by specified label grouping combinations (group:"
                " adj/raw, ancestry_group: gnomAD inferred global ancestry group, sex:"
                " sex karyotype), with values describing the corresponding index of"
                " each grouping entry in the ‘freq’ array row annotation. Please visit"
                " our FAQ, 'How do I access the gnomAD Hail Table frequency"
                " annotation?' on our site"
                " https://gnomad.broadinstitute.org/help#technical-details for"
                " instructions on using this annotation."
            )
        },
        "faf_index_dict": {
            "definition": (
                "Dictionary keyed by specified label grouping combinations (group:"
                " adj/raw, pop: gnomAD inferred global population, sex: sex karyotype),"
                " with values describing the corresponding index of each grouping entry"
                " in the filtering allele frequency (‘faf’) row annotation."
            )
        },
        "faf_meta": {
            "definition": (
                "Filtering allele frequency metadata. An ordered list containing the"
                " frequency aggregation group for each element of the ‘faf’ array row"
                " annotation."
            )
        },
        "vep_version": {"definition": "VEP version that was run on the callset."},
        "vep_csq_header": {"definition": " VEP header for VCF export."},
        "dbsnp_version": {"definition": "dbSNP version used in the callset."},
        "cadd_version": {"definition": "CADD version used in the callset."},
        "revel_version": {"definition": "REVEL version used in the callset."},
        "spliceai_version": {"definition": "SpliceAI version used in the callset."},
        "primateai_version": {"definition": "PrimateAI version used in the callset."},
        "filtering_model": {
            "definition": "The variant filtering model used and its specific cutoffs.",
            "subfields": {
                "model_name": {
                    "definition": (
                        "Variant filtering model name used in the 'filters' row"
                        " annotation, indicating the variant was filtered by this model"
                        " during variant QC."
                    )
                },
                "score_name": {
                    "definition": (
                        "Annotation name of the score used for variant filtering."
                    )
                },
                "snv_cutoff": {
                    "definition": "SNV filtering cutoff information.",
                    "subfields": {
                        "bin": {"definition": "Filtering percentile cutoff for SNVs."},
                        "min_score": {
                            "definition": (
                                "Minimum score at SNV filtering percentile cutoff."
                            )
                        },
                    },
                },
                "indel_cutoff": {
                    "definition": "Indel filtering cutoff information.",
                    "subfields": {
                        "bin": {
                            "definition": "Filtering percentile cutoff for indels."
                        },
                        "min_score": {
                            "definition": (
                                "Minimum score at indel filtering percentile cutoff."
                            )
                        },
                    },
                },
                "model_id": {
                    "definition": (
                        "Variant filtering model ID for score data (used for"
                        " internal specification of the model)."
                    )
                },
                "snv_training_variables": (
                    "Variant annotations used as features in the SNV filtering model."
                ),
                "indel_training_variables": (
                    "Variant annotations used as features in the indel filtering model."
                ),
            },
        },
        "age_distribution": {
            "defintion": "Callset-wide age histogram calculated on release samples.",
            "subfields": {
                "bin_edges": {"defintion": "Bin edges for the age histogram."},
                "bin_freq": {
                    "definition": (
                        "Bin frequencies for the age histogram. This is the number of"
                        " records found in each bin."
                    )
                },
                "n_smaller": {
                    "definition": (
                        "Count of age values falling below lowest histogram bin edge."
                    )
                },
                "n_larger": {
                    "definition": (
                        "Count of age values falling above highest histogram bin edge."
                    )
                },
            },
        },
        "age_index_dict": {
            "definition": (
                "Dictionary keyed by specified subset with values describing the"
                " corresponding index of each subset entry in the 'age_hist_het' and"
                " 'age_hist_hom' array row annotation."
            )
        },
        "age_meta": {
            "definition": (
                "Age metadata. An ordered list containing the groups for each"
                " element of the 'age_hist_het' and 'age_hist_hom' array row"
                " annotation."
            )
        },
        "freq_sample_count": {
            "definition": (
                "A sample count per sample grouping defined in the 'freq_meta' global"
                " annotation."
            )
        },
    },
    "Row fields": {
        "locus": {
            "definition": "Variant locus. Contains contig and position information."
        },
        "alleles": {"definition": "Variant Alleles."},
        "freq": {
            "definition": (
                "Array of allele frequency information (AC, AN, AF, homozygote count)"
                " for each frequency aggregation group in the gnomAD release."
            ),
            "subfields": {
                "AC": {"definition": "Alternate allele count in release."},
                "AF": {
                    "definition": "Alternate allele frequency, (AC/AN), in release."
                },
                "AN": {"definition": "Total number of alleles in release."},
                "homozygote_count": {
                    "definition": (
                        " Count of homozygous alternate individuals in release."
                    )
                },
            },
        },
        "raw_qual_hists": {
            "definition": (
                "Genotype quality metric histograms for all genotypes as opposed to"
                " high quality genotypes."
            ),
            "subfields": {
                "gq_hist_all": {
                    "definition": "Histogram for GQ calculated on all genotypes.",
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the GQ histogram calculated on all"
                                " genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                            ),
                            "bin_freq": {
                                "definition": (
                                    "Bin frequencies for the GQ histogram calculated on"
                                    " all genotypes. The number of records found in"
                                    " each bin."
                                )
                            },
                            "n_smaller": {
                                "definition": (
                                    "Count of GQ values falling below lowest histogram"
                                    " bin edge, for GQ calculated on all genotypes."
                                )
                            },
                            "n_larger": {
                                "definition": (
                                    "Count of GQ values falling above highest histogram"
                                    " bin edge, for GQ calculated on all genotypes."
                                )
                            },
                        },
                    },
                },
                "dp_hist_all": {
                    "definition": "Histogram for DP calculated on all genotypes.",
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the DP histogram calculated on all"
                                " genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the DP histogram calculated on all"
                                " genotypes. The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "Count of DP values falling below lowest histogram bin"
                                " edge, for DP calculated on all genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of DP values falling above highest histogram bin"
                                " edge, for DP calculated on all genotypes."
                            )
                        },
                    },
                },
                "gq_hist_alt": {
                    "definition": (
                        " Histogram for GQ in heterozygous individuals calculated on"
                        " all genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the histogram of GQ in heterozygous"
                                " individuals calculated on all genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the histogram of GQ in"
                                " heterozygous individuals calculated on all genotypes."
                                " The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "Count of GQ values in heterozygous individuals falling"
                                " below lowest histogram bin edge, calculated on all"
                                " genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of GQ values in heterozygous individuals falling"
                                " above highest histogram bin edge, calculated on all"
                                " genotypes."
                            )
                        },
                    },
                },
                "dp_hist_alt": {
                    "definition": (
                        " Histogram for DP in heterozygous individuals calculated on"
                        " all genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the histogram of DP in heterozygous"
                                " individuals calculated on all genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the histogram of DP in"
                                " heterozygous individuals calculated on all genotypes."
                                " The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "Count of DP values in heterozygous individuals falling"
                                " below lowest histogram bin edge, calculated on all"
                                " genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of DP values in heterozygous individuals falling"
                                " above highest histogram bin edge, calculated on all"
                                " genotypes."
                            )
                        },
                    },
                },
                "ab_hist_alt": {
                    "definition": (
                        "Histogram for AB in heterozygous individuals calculated on all"
                        " genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the histogram of AB in heterozygous"
                                " individuals calculated on all genotypes are:"
                                " 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the histogram of AB in"
                                " heterozygous individuals calculated on all genotypes."
                                " The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "CCount of AB values in heterozygous individuals"
                                " falling below lowest histogram bin edge, calculated"
                                " on all genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of AB values in heterozygous individuals falling"
                                " above highest histogram bin edge, calculated on all"
                                " genotypes."
                            )
                        },
                    },
                },
            },
        },
        "popmax": {
            "definition": (
                "Allele frequency information (AC, AN, AF, homozygote count) for the"
                " population with maximum allele frequency."
            ),
            "subfields": {
                "AC": {
                    "definition": (
                        "Alternate allele count in the population with the maximum"
                        " allele frequency."
                    )
                },
                "AF": {
                    "definition": (
                        "Maximum alternate allele frequency, (AC/AN), across"
                        " populations in gnomAD."
                    )
                },
                "AN": {
                    "definition": (
                        "Total number of alleles in the population with the maximum"
                        " allele frequency."
                    )
                },
                "homozygote_count": {
                    "definition": (
                        "Count of homozygous individuals in the population with the"
                        " maximum allele frequency."
                    )
                },
                "pop": {"definition": "Population with the maximum allele frequency."},
                "faf95": {
                    "definition": (
                        "Filtering allele frequency (using Poisson 95% CI) for the"
                        " population with the maximum allele frequency."
                    )
                },
            },
        },
        "qual_hists": {
            "definition": (
                "Genotype quality metric histograms for high quality genotypes."
            ),
            "subfields": {
                "gq_hist_all": {
                    "definition": (
                        "Histogram for GQ calculated on high quality genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the GQ histogram calculated on high"
                                " quality"
                                " genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                            ),
                            "bin_freq": {
                                "definition": (
                                    "Bin frequencies for the GQ histogram calculated on"
                                    " high quality genotypes. The number of records"
                                    " found in each bin."
                                )
                            },
                            "n_smaller": {
                                "definition": (
                                    "Count of GQ values falling below lowest histogram"
                                    " bin edge, for GQ calculated on high quality"
                                    " genotypes."
                                )
                            },
                            "n_larger": {
                                "definition": (
                                    "Count of GQ values falling above highest histogram"
                                    " bin edge, for GQ calculated on high quality"
                                    " genotypes."
                                )
                            },
                        },
                    },
                },
                "dp_hist_all": {
                    "definition": (
                        "Histogram for DP calculated on high quality genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the DP histogram calculated on all"
                                " genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the DP histogram calculated on"
                                " high quality genotypes. The number of records found"
                                " in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "Count of DP values falling below lowest histogram bin"
                                " edge, for DP calculated on high quality genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of DP values falling above highest histogram bin"
                                " edge, for DP calculated on high quality genotypes."
                            )
                        },
                    },
                },
                "gq_hist_alt": {
                    "definition": (
                        " Histogram for GQ in heterozygous individuals calculated on"
                        " high quality genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the histogram of GQ in heterozygous"
                                " individuals calculated on high quality genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the histogram of GQ in"
                                " heterozygous individuals calculated on high quality"
                                " genotypes. The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "Count of GQ values in heterozygous individuals falling"
                                " below lowest histogram bin edge, calculated on high"
                                " quality genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of GQ values in heterozygous individuals falling"
                                " above highest histogram bin edge, calculated on high"
                                " quality genotypes."
                            )
                        },
                    },
                },
                "dp_hist_alt": {
                    "definition": (
                        " Histogram for DP in heterozygous individuals calculated on"
                        " high quality genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the histogram of DP in heterozygous"
                                " individuals calculated on high quality genotypes are:"
                                " 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the histogram of DP in"
                                " heterozygous individuals calculated on high quality"
                                " genotypes. The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "Count of DP values in heterozygous individuals falling"
                                " below lowest histogram bin edge, calculated on high"
                                " quality genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of DP values in heterozygous individuals falling"
                                " above highest histogram bin edge, calculated on high"
                                " quality genotypes."
                            )
                        },
                    },
                },
                "ab_hist_alt": {
                    "definition": (
                        "Histogram for AB in heterozygous individuals calculated on"
                        " high quality genotypes."
                    ),
                    "subfields": {
                        "bin_edges": {
                            "definition": (
                                "Bin edges for the histogram of AB in heterozygous"
                                " individuals calculated on high quality genotypes are:"
                                " 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                            )
                        },
                        "bin_freq": {
                            "definition": (
                                "Bin frequencies for the histogram of AB in"
                                " heterozygous individuals calculated on high quality"
                                " genotypes. The number of records found in each bin."
                            )
                        },
                        "n_smaller": {
                            "definition": (
                                "CCount of AB values in heterozygous individuals"
                                " falling below lowest histogram bin edge, calculated"
                                " on high quality genotypes."
                            )
                        },
                        "n_larger": {
                            "definition": (
                                "Count of AB values in heterozygous individuals falling"
                                " above highest histogram bin edge, calculated on high"
                                " quality genotypes."
                            )
                        },
                    },
                },
            },
        },
        "faf": {
            "definition": "Filtering allele frequency.",
            "subfields": {
                "faf95": {
                    "definition": (
                        "Filtering allele frequency (using Poisson 95% confidence"
                        " interval)"
                    )
                },
                "faf99": {
                    "definition": (
                        "Filtering allele frequency (using Poisson 99% confidence"
                        " interval)"
                    )
                },
                "faf_data_type": {
                    "definition": (
                        "The data type, exome or genome, of the filtering allele"
                        " frequency."
                    )
                },
            },
        },
        "a_index": {
            "defintion": (
                "The original index of this alternate allele in the multiallelic"
                " representation (1 is the first alternate allele or the only alternate"
                " allele in a biallelic variant)."
            )
        },
        "was_split": {
            "defintion": (
                "True if this variant was originally multiallelic, otherwise False."
            )
        },
        "rsid": {
            "definition": (
                "dbSNP reference SNP identification (rsID) numbers. See dbsnp_version"
                " global annotation for version."
            )
        },
        "filters": {
            "definition": (  # TODO: This list may need to be updated depending on Variant QC
                "Variant filters; AC0: Allele count is zero after filtering out"
                " low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het"
                " calls), AS_VQSR: Failed allele-specific VQSR filtering thresholds of"
                " -2.7739 for SNPs and -1.0606 for indels, InbreedingCoeff: GATK"
                " InbreedingCoeff < -0.3, PASS: Passed all variant filters."
            )
        },
        "info": {
            "definition": (
                "Struct containing typical GATK allele-specific (AS) info fields and"
                " additional variant QC fields."
            ),
            "subfields": {
                "QUALapprox": {
                    "definition": (
                        "Sum of PL[0] values; used to approximate the QUAL score."
                    )
                },
                "SB": {
                    "definition": (
                        "Per-sample component statistics which comprise the Fisher's"
                        " exact test to detect strand bias. Values are: depth of"
                        " reference allele on forward strand, depth of reference allele"
                        " on reverse strand, depth of alternate allele on forward"
                        " strand, depth of alternate allele on reverse strand."
                    )
                },
                "MQ": {
                    "definition": (
                        "Root mean square of the mapping quality of reads across all"
                        " samples."
                    )
                },
                "MQRankSum": {
                    "definition": (
                        "Z-score from Wilcoxon rank sum test of alternate vs. reference"
                        " read mapping qualities."
                    )
                },
                "VarDP": {
                    "definition": (
                        "Depth over variant genotypes (does not include depth of"
                        " reference samples)."
                    )
                },
                "AS_ReadPosRankSum": {
                    "definition": (
                        "Allele-specific z-score from Wilcoxon rank sum test of"
                        " alternate vs. reference read position bias."
                    )
                },
                "AS_pab_max": {
                    "definition": (
                        "Maximum p-value over callset for binomial test of observed"
                        " allele balance for a heterozygous genotype, given expectation"
                        " of 0.5."
                    )
                },
                "AS_QD": {
                    "definition": (
                        "Allele-specific variant call confidence normalized by depth of"
                        " sample reads supporting a variant."
                    )
                },
                "AS_MQ": {
                    "definition": (
                        "Allele-specific root mean square of the mapping quality of"
                        " reads across all samples."
                    )
                },
                "QD": {
                    "definition": (
                        "Variant call confidence normalized by depth of sample reads"
                        " supporting a variant."
                    )
                },
                "AS_MQRankSum": {
                    "definition": (
                        "Allele-specific z-score from Wilcoxon rank sum test of"
                        " alternate vs. reference read mapping qualities."
                    )
                },
                "FS": {
                    "definition": (
                        "Phred-scaled p-value of Fisher's exact test for strand bias."
                    )
                },
                "AS_FS": {
                    "definition": (
                        "Allele-specific phred-scaled p-value of Fisher's exact test"
                        " for strand bias."
                    )
                },
                "ReadPosRankSum": {
                    "definition": (
                        "Z-score from Wilcoxon rank sum test of alternate vs. reference"
                        " read position bias."
                    )
                },
                "AS_QUALapprox": {
                    "definition": (
                        "Allele-specific sum of PL[0] values; used to approximate the"
                        " QUAL score."
                    )
                },
                "AS_SB_TABLE": {
                    "definition": (
                        "Allele-specific forward/reverse read counts for strand bias"
                        " tests."
                    )
                },
                "AS_VarDP": {
                    "definition": (
                        "Allele-specific depth over variant genotypes (does not include"
                        " depth of reference samples)."
                    )
                },
                "AS_SOR": {
                    "definition": (
                        "Allele-specific strand bias estimated by the symmetric odds"
                        " ratio test."
                    )
                },
                "SOR": {
                    "definition": (
                        "Strand bias estimated by the symmetric odds ratio test."
                    )
                },
                "singleton": {"definition": "Variant is seen once in the callset."},
                "transmitted_singleton": {
                    "definition": (
                        "Variant was a callset-wide doubleton that was transmitted"
                        " within a family from a parent to a child (i.e., a singleton"
                        " amongst unrelated samples in cohort)."
                    )
                },
                "omni": {
                    "definition": (
                        "Variant is present on the Omni 2.5 genotyping array and found"
                        " in 1000 Genomes data."
                    )
                },
                "mills": {
                    "definition": "Indel is present in the Mills and Devine data."
                },
                "monoallelic": {
                    "definition": (
                        "All samples are homozygous alternate for the variant."
                    )
                },
                "AS_VQSLOD": {
                    "definition": (
                        "Allele-specific log-odds ratio of being a true variant versus"
                        " being a false positive under the trained VQSR Gaussian"
                        " mixture model."
                    )
                },
                "InbreedingCoeff": {
                    "definition": (
                        "Inbreeding coefficient, the excess heterozygosity at a variant"
                        " site, computed as 1 - (the number of heterozygous genotypes)"
                        " / (the number of heterozygous genotypes expected under"
                        " Hardy-Weinberg equilibrium)."
                    )
                },
            },
        },
        "vep": {
            "definition": (
                "Consequence annotations from Ensembl VEP. More details about VEP"
                " output is described at"
                " https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#output."
                " VEP was run using the LOFTEE plugin and information about the"
                " additional LOFTEE annotations can be found at"
                " https://github.com/konradjk/loftee."
            )
        },
        "vqsr": {  # TODO: Update to RF struct if Variant QC uses it
            "definition": "VQSR related variant annotations.",
            "subfields": {
                "AS_VQSLOD": {
                    "definition": (
                        "Allele-specific log-odds ratio of being a true variant versus"
                        " being a false positive under the trained VQSR Gaussian"
                        " mixture model."
                    )
                },
                "AS_culprit": {
                    "definition": (
                        "Allele-specific worst-performing annotation in the VQSR"
                        " Gaussian mixture model."
                    )
                },
                "NEGATIVE_TRAIN_SITE": {
                    "definition": (
                        "Variant was used to build the negative training set of"
                        " low-quality variants for VQSR."
                    )
                },
                "POSITIVE_TRAIN_SITE": {
                    "definition": (
                        "Variant was used to build the positive training set of"
                        " high-quality variants for VQSR."
                    )
                },
            },
        },
        "region_flag": {
            "definition": "Struct containing flags for problematic regions.",
            "subfields": {
                "lcr": {"definition": "Variant falls within a low complexity region."},
                "segdup": {
                    "definition": "Variant falls within a segmental duplication region."
                },
                "nonpar": {
                    "definiton": "Variant falls within a non-pseudoautosomal region."
                },
            },
        },
        "allele_info": {
            "definition": "Allele information.",
            "subfields": {
                "variant_type": {
                    "defintion": (
                        "Variant type (snv, indel, multi-snv, multi-indel, or mixed)."
                    )
                },
                "allele_type": {
                    "defintion": "Allele type (snv, insertion, deletion, or mixed)."
                },
                "n_alt_alleles": {
                    "defintion": (
                        "Total number of alternate alleles observed at variant locus."
                    )
                },
                "was_mixed": {"defintion": "Variant type was mixed."},
            },
        },
        "age_hist_het": {
            "definition": (
                "A list of age histogram information for all heterozygous release"
                " samples calculated on high quality genotypes. One entry for the"
                " callset and each subset."
            ),
            "subfields": {
                "bin_edges": {"definition": "Bin edges for the age histogram."},
                "bin_freq": {
                    "definition": (
                        "Bin frequencies for the age histogram. This is the number of"
                        " records found in each bin."
                    )
                },
                "n_smaller": {
                    "definition": (
                        "Count of age values falling below lowest histogram bin edge."
                    )
                },
                "n_larger": {
                    "definition": (
                        "Count of age values falling above highest histogram bin edge."
                    )
                },
            },
        },
        "age_hist_hom": {
            "definition": (
                "A list of age histogram information for all homozygous (and hemizygous"
                " on sex chromosomes) release samples calculated on high quality"
                " genotypes. One entry for the callset and each subset."
            ),
            "subfields": {
                "bin_edges": {"definition": "Bin edges for the age histogram."},
                "bin_freq": {
                    "definition": (
                        "Bin frequencies for the age histogram. This is the number of"
                        " records found in each bin."
                    )
                },
                "n_smaller": {
                    "definition": (
                        "Count of age values falling below lowest histogram bin edge."
                    )
                },
                "n_larger": {
                    "definition": (
                        "Count of age values falling above highest histogram bin edge."
                    )
                },
            },
        },
        "cadd": {
            "definition": (
                "Combined Annotation Depenedent Depletion tool for scoring the"
                " deleteriousness of single nucleotide variants as well as"
                " insertion/deletions variants in the human genome."
            ),
            "subfields": {
                "phred": {
                    "definition": (
                        "Cadd Phred-like scores ('scaled C-scores') ranging from 1 to"
                        " 99, based on the rank of each variant relative to all"
                        " possible 8.6 billion substitutions in the human reference"
                        " genome. Larger values are more deleterious."
                    )
                },
                "raw_score": {
                    "definition": (
                        "Raw CADD scores are interpretable as the extent to which the"
                        " annotation profile for a given variant suggests that the"
                        " variant is likely to be 'observed' (negative values) vs"
                        " 'simulated' (positive values). Larger values are more"
                        " deleterious."
                    )
                },
                "has_duplicate": {  # TODO: Confirm this is still in schema, v3.1 removed because 100% missing
                    "definition": (
                        "Whether the variant has more than one CADD score associated"
                        " with it."
                    )
                },
            },
        },
        "revel": {
            "definition": (
                "An ensemble method for predicting the pathogenicity of rare missense"
                " variants."
            ),
            "subfields": {
                "revel_score": {
                    "defintion": (
                        "dbNSFP's Revel score from 0 to 1. Variants with higher scores"
                        " are predicted to be more likely to be deleterious."
                    )
                },
                "has_duplicate": {
                    "defintion": (
                        "Whether the variant has more than one revel_score associated"
                        " with it."
                    )
                },
            },
        },
        "splice_ai": {  # TODO: Confirm which splice tool we're using and add appropriate fields, i.e. are we using both spliceAI and pangolin?
            "definition": "A deep learning-based tool to identify splice variants.",
            "subfeilds": {
                "splice_ai": {
                    "definition": (
                        "The maximum delta score, interpreted as the probability of the"
                        " variant being splice-altering."
                    )
                },
                "splice_consequence": {
                    "definition": (
                        "The consequence term associated with the max delta score in"
                        " 'splice_ai_max_ds'."
                    )
                },
                "has_duplicate": {
                    "definition": (
                        "Whether the variant has more than one splice_ai score"
                        " associated with it."
                    )
                },
            },
        },
        "primate_ai": {
            "definition": (
                "A deep residual neural network for classifying the pathogenicity of"
                " missense mutations."
            ),
            "subfields": {
                "primate_ai_score": {
                    "definition": (
                        "PrimateAI's deleteriousness score from 0 (less deleterious) to"
                        " 1 (more deleterious)."
                    )
                },
                "has_duplicate": {
                    "definition": (
                        "Whether the variant has more than one primate_ai_score"
                        " associated with it."
                    )
                },
            },
        },
        "subsets": {
            "definition": "Subsets a variant is found in.",
            "subfields": {
                "ukb": {"definition": "If a variant is in the UK Biobank subset."},
                "non-ukb": {
                    "definition": "If a variant is in the non-UK Biobank subset."
                },
                "non-topmed": {
                    "definition": "If a variant is found in the non-TOPMed subset."
                },
            },
        },
    },
}


=======
>>>>>>> eb58740 (Move field description dict to resources)
TABLES_FOR_RELEASE = [
    "dbsnp",
    "filters",
    "freq",
    "info",
    "subsets",
    "region_flags",
    "in_silico",
    "vep",
]


def custom_region_flags_select(ht):
    """
    Select region flags for release.

    :param ht: hail Table.
    :return: select expression dict
    """
    selects = {}
    selects = {
        "region_flags": region_flag_expr(
            ht,
            non_par=True,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
        )
    }
    return selects


def custom_filters_select(ht):
    """
    Select gnomad filter HT fields for release dataset.

    Extract fields like 'filters', 'vqsr', and generates 'alelle_info' struct.
    :param ht: hail table.
    :return: select expression dict.
    """
    selects = {}
    selects["allele_info"] = hl.struct(
        variant_type=ht.variant_type,
        allele_type=ht.allele_type,
        n_alt_alleles=ht.n_alt_alleles,
        was_mixed=ht.was_mixed,
    )
    return selects


def custom_subset_select(ht):
    """
    Select release subset field using freq HT AN value.

    :param ht: hail table
    :return: select expression dict
    """
    selects = {
        subset: hl.if_else(
            ht.freq[ht.freq_index_dict[f"{subset}-adj"]].AN > 0, True, False
        )
        for subset in SUBSETS
    }
    return selects


def custom_info_select(ht):
    """
    Select fields for info hail Table annotation in release.

    The info field requires fields from the freq HT and the filters HT
    so those are pulled in here along with all info HT fields.

    :param ht: hail table
    :return: select expression dict
    """
    freq_ht = CONFIG.get("freq")["ht"]
    freq_info_dict = {"InbreedingCoeff": freq_ht[ht.key]["InbreedingCoeff"]}

    filters_ht = CONFIG.get("filters")["ht"]
    filters = filters_ht[ht.key]
    filters_info_fields = [
        "singleton",
        "transmitted_singleton",
        "omni",
        "mills",
        "monoallelic",
        "SOR",
    ]
    filters_info_dict = {field: filters[field] for field in filters_info_fields}
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    filters_info_dict.update({**{f"{score_name}": filters[f"{score_name}"]}})

    info_dict = {field: ht.info[field] for field in SITE_FIELDS + AS_FIELDS}
    info_dict.update(filters_info_dict)
    info_dict.update(freq_info_dict)

    selects = {"info": hl.struct(**info_dict)}

    return selects


def get_select_global_fields(ht):
    """
    Generate a dictionary of globals to select by checking the configs of all tables joined.

    :param ht: Final joined HT with globals.
    """
    t_globals = [
        get_select_fields(CONFIG.get(t)["select_globals"], ht)
        for t in args.tables_for_join
        if "select_globals" in CONFIG.get(t)
    ]
    t_globals = reduce(lambda a, b: dict(a, **b), t_globals)

    return t_globals


def get_select_fields(selects, base_ht):
    """
    Generate a select dict from traversing the base_ht and extracting annotations.

    :param selects: mapping or list of selections
    :param base_ht: base_ht to traverse
    :return: select mapping from annotation name to base_ht annotation
    """
    select_fields = {}
    if selects is not None:
        if isinstance(selects, list):
            select_fields = {selection: base_ht[selection] for selection in selects}
        elif isinstance(selects, dict):
            for key, val in selects.items():
                # Grab the field and continually select it from the hail table.
                ht = base_ht
                for attr in val.split("."):
                    ht = ht[attr]
                select_fields[key] = ht
    return select_fields


def get_ht(dataset, _intervals, test) -> hl.Table:
    """
    Return the appropriate hail table with selects applied.

    :param dataset: Hail Table to join.
    :param _intervals: Intervals for reading in hail Table.
    :param test: Whether call is for a test run.
    :return: Hail Table with fields to select.
    """
    config = CONFIG[dataset]
    ht_path = config["path"]
    logger.info("Reading in %s", dataset)
    base_ht = hl.read_table(ht_path, _intervals=_intervals)

    if test:
        base_ht = hl.filter_intervals(
            base_ht,
            [hl.parse_locus_interval("chr1:1-1000000", reference_genome="GRCh38")],
        )

    select_fields = get_select_fields(config.get("select"), base_ht)

    if "custom_select" in config:
        custom_select_fn = config["custom_select"]
        select_fields = {**select_fields, **custom_select_fn(base_ht)}

    if "field_name" in config:
        field_name = config.get("field_name")
        select_query = {field_name: hl.struct(**select_fields)}
    else:
        select_query = select_fields

    return base_ht.select(**select_query)


def join_hts(base_table, tables, new_partition_percent, test):
    """
    Outer join a list of hail tables.

    :param base_table: Dataset to use for interval partitioning.
    :param tables: List of tables to join.
    :param new_partition_percent: Percent of base_table partitions used for final release hail Table.
    :param test: Whether this is for a test run.
    """
    logger.info(
        "Reading in %s to determine partition intervals for efficient join",
        base_table,
    )
    base_ht_path = CONFIG[base_table]["path"]
    base_ht = hl.read_table(base_ht_path)
    if test:
        base_ht = hl.filter_intervals(
            base_ht,
            [hl.parse_locus_interval("chr1:1-1000000", reference_genome="GRCh38")],
        )
    partition_intervals = base_ht._calculate_new_partitions(
        base_ht.n_partitions() * new_partition_percent
    )

    hts = [get_ht(table, _intervals=partition_intervals, test=test) for table in tables]
    # TODO: Check with hail if an intermediate checkpoint be helpful here?
    joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "left")), hts)

    # Track the dataset we've added as well as the source path.
    included_dataset = {k: v["path"] for k, v in CONFIG.items() if k in tables}

    joined_ht = joined_ht.annotate_globals(
        date=datetime.now().isoformat(),
        datasets=hl.dict(included_dataset),
    )
    joined_ht.describe()
    return joined_ht


"""
Configurations of dataset to combine.
Format:
'<Name of dataset>': {
        'path': 'gs://path/to/hailtable.ht',
        'select': '<Optional list of fields to select or dict of new field name to location of old fieldin the dataset.>',
        'field_name': '<Optional name of root annotation in combined dataset, defaults to name of dataset.>',
        'custom_select': '<Optional function name of custom select function that is needed for more advanced logic>',
        'select_globals': '<Optional list of globals to select or dict of new global field name to old global field name. If not specified, all globals are selected.>
    },
"""

CONFIG = {
    "dbsnp": {
        "ht": dbsnp.ht(),
        "path": dbsnp.path,
        "select": ["rsid"],
        "select_globals": {
            "dbsnp_version": "version",  # TODO: Need to add global to this table with version
        },
    },
    "filters": {
        "ht": final_filter().ht(),
        "path": final_filter().path,
        "select": ["filters", "vqsr"],
        "custom_select": custom_filters_select,
        "select_globals": ["filtering_model", "inbreeding_coeff_cutoff"],
    },
    "in_silico": {
        "ht": analyst_annotations.ht(),
        "path": analyst_annotations.path,
        # TODO: Update these once we knew which tools we will be usings
        "select": ["cadd", "revel", "splice_ai", "primate_ai"],
        # TODO: Update these once we knew which tools we will be usings
        "select_globals": ["cadd_version", "revel_version"],
    },
    "info": {
        "ht": get_info().ht(),
        "path": get_info().path,
        "select": ["was_split", "a_index"],
        "custom_select": custom_info_select,
    },
    "freq": {
        "ht": get_freq().ht(),
        "path": get_freq().path,
        "select": [
            "freq",
            "faf",
            "popmax",
            "qual_hists",
            "raw_qual_hists",
            "age_hist_het",
            "age_hist_hom",
        ],
        "select_globals": [
            "freq_meta",
            "freq_index_dict",
            "faf_meta",
            "faf_index_dict",
            "age_meta",
            "age_index_dict",
            "age_distribution",
        ],
    },
    "subsets": {
        "ht": get_freq().ht(),
        "path": get_freq().path,
        "custom_select": custom_subset_select,
        "field_name": "subsets",
    },
    "vep": {
        "ht": vep.ht(),
        # TODO: drop 100% missing? Module to do this after all annotations added?
        "select": ["vep"],
        "select_globals": ["vep_version"],  # TODO: Confirm or add this is a global
    },
    "region_flags": {
        "ht": get_freq().ht(),
        "path": get_freq().path,
        "custom_select": custom_region_flags_select,
    },
    "release": {
        "ht": release_sites().ht(),
        "path": release_sites().path,
        "select": [r for r in release_sites().ht()._row],
        "select_globals": [g for g in release_sites().ht()._global],
    },
}


def main(args):
    """Create release ht."""
    hl.init(
        log="/create_release_ht.log",
        tmp_dir="gs://gnomad-tmp-4day",
        default_reference="GRCh38",
    )

    ht = join_hts(
        args.base_table,
        args.tables_for_join,
        args.new_partition_percent,
        args.test,
    )

    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

    t_globals = get_select_global_fields(ht)

    ht = ht.annotate_globals(
        **t_globals,
        gnomad_qc_version=args.gnomad_qc_version,
        # TODO: See if we can pull this from the cluster
        gnomad_methods_version=args.gnomad_methods_version,
        README=FIELD_DESCRIPTIONS,  # TODO: Make version dict for this and have it live in methods?
        version=args.version,
    )

    logger.info("Writing out release HT to %s", release_sites().path)
    ht = ht.checkpoint(
        qc_temp_prefix() + "/release / gnomad.genomes.sites.test.ht"
        if args.test
        else release_sites().path,
        args.overwrite,
    )

    logger.info("Final variant count: %d", ht.count())
    ht.describe()
    ht.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--new-partition-percent",
        help=(
            "Percent of start dataset partitions to use for release HT. Default is 1.1"
            " (110%)"
        ),
        default=1.1,
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "-v",
        "--version",
        help="The version of gnomAD.",
        default=CURRENT_RELEASE,
        # TODO: Consider using this for iterative releases - requires restructing
        # of config select fields where there are changes
        required=True,
    )
    parser.add_argument(
        "-t",
        "--test",
        help="Runs a test on the first two partitions of the HT.",
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
        default="freq",
        choices=TABLES_FOR_RELEASE,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite existing HT.",
        action="store_true",
    )
    parser.add_argument(
        "--gnomad-qc-version",
        help="Version of gnomAD QC repo",
        required=True,
    )
    parser.add_argument(
        "--gnomad-methods-version",
        help="Version of gnomAD methods repo",
        required=True,
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
