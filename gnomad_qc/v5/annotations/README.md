# gnomAD v5 Frequency Calculation

This directory contains scripts for calculating variant frequencies and generating age histograms for gnomAD v5.

## Overview

The v5 frequency calculation processes consent withdrawal samples and All of Us (AoU) dataset samples that need to be added to gnomAD v5. The script:

1. **Processes gnomAD consent withdrawal samples** by subtracting their frequencies and age histograms from the v4 frequency HT
2. **Processes All of Us samples** by calculating their frequencies using imported AN values from pre-calculated consent data
3. **Merges both datasets** to create the final v5 frequency table with updated FAF/grpmax annotations

## Key Features

### Consent Withdrawal Processing

The script handles gnomAD samples that need to be removed from v5 due to consent withdrawal:

- **Loads v4 frequency HT** as the base dataset containing both frequencies and age histograms
- **Identifies consent withdrawal samples** using the `consent_samples_to_drop` resource
- **Calculates frequencies and age histograms** for consent samples using corrected genotypes (after hom alt depletion fix and sex ploidy adjustment)
- **Subtracts consent data** from v4 frequencies and age histograms using `merge_freq_arrays` with `operation="diff"`
- **Applies FAF and grpmax annotations** to the updated frequency data

### All of Us Integration

The script processes AoU samples to be added to gnomAD v5:

- **Resolves sample ID collisions** between gnomAD and AoU datasets
- **Filters related samples** that should be excluded
- **Calculates complete frequency structs** using AC/homozygote counts from variant data and AN values imported from `get_consent_ans`
- **Generates age histograms** for AoU samples
- **Uses efficient `agg_by_strata` approach** for frequency calculations

### Quality Corrections Applied

Both datasets undergo consistent quality corrections:

- **Adj filtering** using standard gnomAD quality criteria
- **Hom alt depletion fix** to correct for systematic genotype calling issues
- **Sex ploidy adjustment** for proper handling of X/Y chromosomes
- **Age histogram calculation** using final corrected genotypes to ensure consistency

## Usage

### Process gnomAD Consent Withdrawals

```bash
python gnomad_qc/v5/annotations/generate_frequency.py \
    --process-gnomad \
    --environment production
```

### Process All of Us Dataset

```bash
python gnomad_qc/v5/annotations/generate_frequency.py \
    --process-aou \
    --environment production
```

### Merge Both Datasets

```bash
python gnomad_qc/v5/annotations/generate_frequency.py \
    --merge-datasets \
    --environment production
```

### Testing Mode

```bash
python gnomad_qc/v5/annotations/generate_frequency.py \
    --process-gnomad \
    --data-test \
    --environment development
```

## Output Files

The script generates frequency tables with embedded age histograms:

- **gnomAD frequency table**: v4 frequencies with consent withdrawals removed, includes embedded age histograms
- **AoU frequency table**: Complete frequencies for AoU samples with separate age histogram table
- **Merged frequency table**: Combined dataset with complete frequency data and age histograms

## Methodology

### gnomAD Consent Withdrawal Processing

1. **Sample Identification**: Uses `consent_samples_to_drop` resource
2. **VDS Preparation**: Loads v4 genomes VDS and filters to consent samples
3. **Quality Corrections**: Applies adj filtering, hom alt depletion fix, and sex ploidy adjustment
4. **Frequency Calculation**: Uses `agg_by_strata` for efficient frequency computation
5. **Age Histograms**: Calculated using corrected genotypes to match frequency calculations
6. **Subtraction**: Removes consent data from v4 frequencies using `merge_freq_arrays`
7. **Post-processing**: Adds FAF, grpmax, and inbreeding coefficient annotations

### All of Us Processing

1. **Sample Collision Resolution**: Handles overlapping sample IDs with gnomAD
2. **Relatedness Filtering**: Removes related samples using v5 relatedness criteria
3. **Variant Data Processing**: Calculates AC and homozygote counts from sparse variant data
4. **AN Import**: Uses pre-calculated AN values from `get_consent_ans` (calculated by separate script)
5. **Complete Frequency Struct**: Builds proper frequency arrays with AC, AF, AN, homozygote_count
6. **Age Histograms**: Calculated separately for AoU samples

### Dataset Merging

1. **Frequency Merging**: Combines gnomAD and AoU frequencies using `merge_freq_arrays` with `operation="sum"`
2. **Age Histogram Merging**: Combines age histograms using `merge_histograms` with `operation="sum"`
3. **Global Metadata**: Merges frequency metadata and sample counts from both datasets

## Key Dependencies

- **v4 genomes VDS**: Source data for consent withdrawal processing
- **v4 frequency HT**: Base frequency data with embedded age histograms
- **AoU VDS**: All of Us variant dataset
- **Consent AN data**: Pre-calculated allele numbers from `get_consent_ans`
- **Sample resources**: `consent_samples_to_drop`, `related_samples_to_drop`, `sample_id_collisions`
- **Group membership tables**: For frequency stratification

## Technical Details

### Genotype Processing Order

For proper consistency between frequencies and age histograms:

1. **Adj annotation**: Quality filtering applied first
2. **Hom alt depletion fix**: Corrects systematic calling issues (requires call expression)
3. **Sex ploidy adjustment**: Handles X/Y chromosome ploidy (converts to integer for frequency calc)
4. **Age histogram calculation**: Uses corrected, sex-adjusted genotypes

### Frequency Stratification

- **Sex karyotype**: XX/XY stratification
- **Genetic ancestry**: Population-based stratification
- **Dataset source**: gnomAD vs AoU identification
- **Quality filtering**: adj vs raw call stratification

## Notes

- **AN Calculation**: AoU frequencies use AN values calculated by a separate script and imported via `get_consent_ans`
- **Age Histogram Consistency**: Calculated after all genotype corrections to match frequency calculations
- **Sex Chromosome Handling**: Proper ploidy adjustment ensures correct het/hom classification
- **Memory Efficiency**: Uses sparse matrix operations and checkpointing to handle large datasets
