# Additional annotations

## CADD

CADD (Rentzsch et al., Nucleic Acids Research, 2018), is a score predicting deleteriousness for both SNVs and indels. The gnomAD browser displays the CADD Phred scores, which range from 1 to 99, based on the rank of each variant relative to all possible 8.6 billion substitutions in the human reference genome. Higher scores are predicted to be more deleterious/damaging. Variants in gnomAD v3.1.1 were annotated with CADD v1.6. Pre-computed CADD scores are available [here](https://cadd.gs.washington.edu/download).

For the gnomAD v4 release, we used the latest CADD v1.6 post-release 1 (released on Mar 22, 2021) to compute a score for indels new to the release (~32,561,253). The following Hail commands were used to generate 1000 small VCF files: 

```commandline
import hail as hl

from gnomad_qc.v4.resources.basics import gnomad_v4_genotypes
from gnomad.resources.grch38.gnomad import public_release

vds = gnomad_v4_genotypes.vds()
ht = vds.variant_data.rows()
split_ht = hl.split_multi_hts(ht)
indels = split_ht.filter(hl.is_indel(split_ht.alleles[0], split_ht.alleles[1]))

v3_ht = public_release("genomes").ht()

v4_new_indels = indels.anti_join(v3_ht)
v4_new_indels = v4_new_indels.checkpoint("gs://gnomad-tmp/qin/gnomad_v4_new_indels.split.ht", overwrite=True) # Checkpoint the matrix table to disk by writing and reading using a fast, but less space-efficient codec.
v4_new_indels = hl.read_table("gs://gnomad-tmp/qin/gnomad_v4_new_indels.split.ht", _n_partitions=1000)
hl.export_vcf(v4_new_indels,'gs://gnomad-tmp/qin/gnomad_v4_new_indels_small.vcf.bgz', parallel="header_per_shard") 
```

The VCF files were unzipped to remove the 'chr' for each line, otherwise it will run into "Encountered uncovered chromosome
Possible precedence issue with control flow operator" (refer to this [issue](https://github.com/kircherlab/CADD-scripts/issues/37)): 

```commandline
# rename *.bgz as *.gz
for file in *.bgz; do 
    mv -- "$file" "${file%.bgz}.pre.gz"
done

for i in *.pre.gz; do gunzip $i; done

for i in *.pre; do awk '{gsub(/^chr/,""); print}' $i > $i.vcf ; done

for i in *.pre.vcf; do gzip $i; done

for file in *.pre.vcf.gz; do 
    mv -- "$file" "${file%.pre.vcf.gz}.vcf.gz"
done
```

CADD was installed with the 216G annotation file downloaded on an attached disk to a 224-core n2d-highcpu VM on google cloud (note: in order to avoid the loss of the annotation file when switching between machine types for scaling and debugging purposes, installation on an attached disk is recommended; also we needed to ask for a reasonable size of bootable disk for the intermediate files created during the run, in /tmp/tmp.* by default). Since CADD requires only one core to run for each VCF file, we parallelized it for 215 VCF files at the same time: 

```commandline
nohup cat parallel.list | parallel -j215 time /mnt/disks/qin/cadd/CADD-scripts-1.6.post1/CADD.sh {} > ../run_parallel215.log 2>&1 &
```

The job itself finished in ~6 hours for 999 VCF files, only 1 file ran into an unknown error, but the problem disappeared by splitting this VCF to smaller ones (technical details described [here](https://github.com/broadinstitute/gnomad_production/issues/782)). 







