"""Script containing release related resources."""
import logging
from typing import Optional

from gnomad.resources.grch38.gnomad import coverage, public_release
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v4.resources.constants import (
    COVERAGE_RELEASES,
    CURRENT_COVERAGE_RELEASE,
    CURRENT_RELEASE,
    RELEASES,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("release_resources")
logger.setLevel(logging.INFO)

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


def _release_root(
    version: str = CURRENT_RELEASE,
    test: bool = False,
    data_type: str = "exomes",
    extension: str = "ht",
) -> str:
    """
    Get root path to the release files.

    :param version: Version of release path to return.
    :param test: Whether to use a tmp path for testing.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
        Default is "exomes".
    :param extension: File extension of release file. Default is "ht".
    :return: Root path of the release files.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/release/{extension}/{data_type}"
        if test
        else f"gs://gnomad/release/{version}/{extension}/{data_type}"
    )


def annotation_hists_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Return path to file containing ANNOTATIONS_HISTS dictionary.

    Dictionary contains histogram values for each metric.
    For example, "InbreedingCoeff": [-0.25, 0.25, 50].

    :return: Path to file with annotation histograms
    """
    return f"gs://gnomad/release/{release_version}/json/annotation_hists.json"


def qual_hists_json_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Fetch filepath for qual histograms JSON.

    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: File path for histogram JSON
    """
    return f"gs://gnomad/release/{release_version}/json/gnomad.exomes.v{release_version}.json"


def release_ht_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables.

    :param data_type: 'exomes' or 'genomes'
    :param release_version: Release version
    :param public: Determines whether release sites Table is read from public or private bucket. Defaults to private
    :return: File path for desired Hail Table
    """
    if public:
        if file_exists(public_release(data_type).versions[release_version].path):
            return public_release(data_type).versions[release_version].path
        else:
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.ht"


def release_sites(public: bool = False) -> VersionedTableResource:
    """
    Retrieve versioned resource for sites-only release Table.

    :param public: Determines whether release sites Table is read from public or private bucket. Defaults to private
    :return: Sites-only release Table
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_ht_path(release_version=release, public=public)
            )
            for release in RELEASES
        },
    )


def release_vcf_path(
    release_version: Optional[str] = None,
    contig: Optional[str] = None,
) -> str:
    """
    Fetch bucket for release (sites-only) VCFs.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used.
    :param contig: String containing the name of the desired reference contig. Defaults to the full (all contigs) sites VCF path
        sites VCF path
    :return: Filepath for the desired VCF
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    if contig:
        return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}.sites.{contig}.vcf.bgz"
    else:
        # If contig is None, return path to sharded vcf bucket.
        # NOTE: need to add .bgz or else hail will not bgzip shards.
        return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}.sites.vcf.bgz"


def release_header_path(release_version: Optional[str] = None) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used
    :return: Filepath for header dictionary pickle
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}_header_dict.pickle"


def append_to_vcf_header_path(
    subset: str, release_version: str = CURRENT_RELEASE
) -> str:
    """
    Fetch path to TSV file containing extra fields to append to VCF header.

    Extra fields are VEP and dbSNP versions.

    :param subset: One of the possible release subsets
    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: Filepath for extra fields TSV file
    """
    return (
        f"gs://gnomad/release/{release_version}/vcf/exomes/extra_fields_for_header{f'_{subset}' if subset else ''}.tsv"
    )


def release_coverage_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
    test: bool = False,
) -> str:
    """
    Fetch filepath for coverage release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version.
    :param public: Determines whether release coverage Table is read from public or
        private bucket. Default is public.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: File path for desired coverage Hail Table.
    """
    if public:
        if test:
            raise ValueError("Cannot use test=True with public=True!")
        try:
            cov = coverage(data_type)
            if release_version in cov.versions:
                path = cov.versions[release_version].path
            else:
                path = None
        except DataException:
            path = None
        if path is None:
            logger.warning(
                "No public coverage Table found for data_type %s and release %s. "
                "Using 'gs://gnomad-public-requester-pays' path.",
                data_type,
                release_version,
            )
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.coverage.ht"
        else:
            return path
    else:
        return (
            f"{_release_root(release_version, test=test, data_type=data_type)}/gnomad.{data_type}.v{release_version}.coverage.ht"
        )


def release_coverage_tsv_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_COVERAGE_RELEASE["exomes"],
    test: bool = False,
) -> str:
    """
    Fetch path to coverage TSV file.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version. Default is CURRENT_COVERAGE_RELEASE["exomes"].
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Coverage TSV path.
    """
    return (
        f"{_release_root(release_version, test=test, data_type=data_type, extension='tsv')}/gnomad.{data_type}.v{release_version}.coverage.tsv.bgz"
    )


def release_coverage(
    data_type: str = "exomes", public: bool = False, test: bool = False
) -> VersionedTableResource:
    """
    Retrieve versioned resource for coverage release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param public: Determines whether release coverage Table is read from public or
        private bucket. Default is private.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Coverage release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_COVERAGE_RELEASE[data_type],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    data_type=data_type,
                    release_version=release,
                    public=public,
                    test=test,
                )
            )
            for release in COVERAGE_RELEASES[data_type]
        },
    )
