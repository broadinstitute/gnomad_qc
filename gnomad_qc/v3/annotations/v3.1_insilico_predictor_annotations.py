#!/usr/bin/env python
# coding: utf-8

# In[1]:


import hail as hl
import gnomad


# In[2]:


def create_grch38_revel_info(path, ow = False, give_object = False):
    revel_information = hl.import_table("gs://gnomad-wphu/revel_grch38_all_chromosomes.csv", delimiter=",", types={'hg19_pos':hl.tint,'grch38_pos':hl.tstr,'REVEL': hl.tfloat64})
    revel_information_grch38 = revel_information.drop("hg19_pos")
    revel_information_grch38 = revel_information_grch38.filter(revel_information_grch38.grch38_pos.contains("."),keep=False)
    revel_information_grch38 = revel_information_grch38.annotate(grch38_pos_int = hl.int(revel_information_grch38.grch38_pos))
    revel_information_grch38 = revel_information_grch38.transmute(grch38_pos = revel_information_grch38.grch38_pos_int)
    revel_information_grch38 = revel_information_grch38.transmute(chr = "chr" + revel_information_grch38.chr)
    revel_information_grch38 = revel_information_grch38.annotate(locus = hl.locus(revel_information_grch38.chr, revel_information_grch38.grch38_pos, reference_genome="GRCh38"))
    revel_information_grch38 = revel_information_grch38.annotate(alleles = hl.array([revel_information_grch38.ref, revel_information_grch38.alt]))
    revel_information_grch38 = revel_information_grch38.select("locus","alleles","REVEL", "aaref", "aaalt")
    revel_information_grch38 = revel_information_grch38.key_by("locus", "alleles")
    if ow:
        revel_information_grch38.write(path, overwrite=ow)
    if give_object:
        return revel_information_grch38
    
def load_grch38_revel_info():
    return(hl.read_table("gs://gnomad-wphu/revel_annotations_grch38.ht"))


# In[3]:


from gnomad.utils.vcf import ht_to_vcf_mt
def get_chr():
    info_ht = gnomad_vars.select()
    info_ht = info_ht.filter(info_ht.locus.contig=="chr1")
    info_ht = info_ht.head(100000000)
    hl.export_vcf(info_ht, "gs://gnomad-wphu/info_split_large_chr1.vcf")


# In[4]:


def create_cadd_info():
    cadd = hl.import_table("gs://gnomad-wphu/gnomad.genomes.r3.0.indel.tsv", comment = "#", types={'Pos':hl.tint32, "RawScore":hl.tfloat, "PHRED":hl.tfloat})
    cadd = cadd.transmute(Chrom = "chr" + cadd.Chrom)
    cadd = cadd.annotate(locus = hl.locus(cadd.Chrom, cadd.Pos, reference_genome="GRCh38"))
    cadd = cadd.annotate(alleles = hl.array([cadd.Ref, cadd.Alt]))
    cadd = cadd.select("locus","alleles","RawScore", "PHRED")
    cadd = cadd.key_by("locus", "alleles")
    return cadd


# In[5]:


def load_CADD(path, n_partitions, force_bgz = False):
    column_names = {'f0': 'chrom', 'f1': 'pos', 'f2': 'ref', 'f3': 'alt', 'f4': 'RawScore', 'f5': 'PHRED'}
    types = {'f0': hl.tstr, 'f1': hl.tint32, 'f4': hl.tfloat32, 'f5': hl.tfloat32}
    cadd_ht = hl.import_table(path, comment="#", no_header=True, types=types, min_partitions=n_partitions, force_bgz = force_bgz)
    cadd_ht = cadd_ht.rename(column_names)
    chrom = hl.format("chr%s", cadd_ht.chrom)
    locus = hl.locus(chrom, cadd_ht.pos, reference_genome="GRCh38")
    alleles = hl.array([cadd_ht.ref, cadd_ht.alt])
    cadd_ht = cadd_ht.transmute(locus=locus, alleles=alleles)
    cadd_union_ht = cadd_ht.head(0)
    for contigs in (range(1, 10), list(range(10, 23)) + ["X", "Y", "MT"]):
        contigs = ["chr%s" % contig for contig in contigs]
        cadd_ht_subset = cadd_ht.filter(hl.array(list(map(str, contigs))).contains(cadd_ht.locus.contig))
        cadd_union_ht = cadd_union_ht.union(cadd_ht_subset)

    cadd_union_ht = cadd_union_ht.select("locus", "alleles", "RawScore", "PHRED")
    cadd_union_ht = cadd_union_ht.key_by("locus", "alleles")
    cadd_union_ht = cadd_union_ht.annotate_globals(source_file_path = path)

    cadd_union_ht.describe()

    return cadd_union_ht

def make_unified_CADD():
    snvs = hl.read_table("gs://gnomad-wphu/CADD-v1.6-SNVs.ht")
    release3_indels = hl.read_table("gs://gnomad-wphu/CADD-v1.6-indels-updated.ht")
    raw31_indels = hl.read_table("gs://gnomad-wphu/CADD-indels-gnomad.3.1.ht")
    raw31_complex = hl.read_table("gs://gnomad-wphu/CADD-1.6-gnomad-complex-variants.ht")
    unified = snvs.head(0)
    unified = unified.union(snvs,release3_indels,raw31_indels, raw31_complex)
    unified = unified.annotate_globals(source_file_path = {"snvs":"gs://gnomad-wphu/CADD-v1.6-SNVs.ht",
                                                         "v3-indels":"gs://gnomad-wphu/CADD-v1.6-indels-updated.ht",
                                                         "v3.1-indels":"gs://gnomad-wphu/CADD-indels-gnomad.3.1.ht",
                                                         "v3.1-complex":"gs://gnomad-wphu/CADD-1.6-gnomad-complex-variants.ht"
                                                        })
    #unified.describe()
    return unified

def convert_CADD_indels_64_32():
    release3_indels = hl.read_table("gs://gnomad-wphu/CADD-v1.6-indels.ht")
    release3_indels = release3_indels.transmute(RawScore = hl.float32(release3_indels.RawScore))
    release3_indels = release3_indels.transmute(PHRED = hl.float32(release3_indels.PHRED))
    release3_indels.describe()
    release3_indels.write("gs://gnomad-wphu/CADD-v1.6-indels-updated.ht", overwrite = True)
    
#unified_CADD = make_unified_CADD()
#unified_CADD = unified_CADD.write("gs://gnomad-wphu/complete-CADD-v1.6-annotations.ht", overwrite=True)


# In[6]:


def export_for_CADD_analysis(hl_tbl,path):
    export = hl_tbl.select()
    export = export.filter(export.locus.contig=="chrM",keep=False)
    hl.methods.export_vcf(export, path)


# In[7]:


def combine_splice_ai():
    recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    splice_snps_skip_invalid = hl.import_vcf("gs://gnomad-wphu/splice_ai_data/splice_ai_data/genome_scores_v1.3_ds.20a701bc58ab45b59de2576db79ac8d0/spliceai_scores.masked.snv.hg38.vcf.gz",
                                force_bgz= True,
                                min_partitions=3000,
                                reference_genome='GRCh38', contig_recoding=recode, skip_invalid_loci= True
                               )
    splice_snps_skip_invalid.annotate_globals(source_file_path="gs://gnomad-wphu/splice_ai_data/splice_ai_data/genome_scores_v1.3_ds.20a701bc58ab45b59de2576db79ac8d0/spliceai_scores.masked.snv.hg38.vcf.gz")
    
    splice_indels_skip_invalid = hl.import_vcf("gs://gnomad-wphu/splice_ai_data/gnomAD_v3.1_SpliceAI_scores-selected/spliceai_scores.masked.gnomad_indel.hg38.vcf.gz",
                                   force_bgz=True,
                                   reference_genome='GRCh38', contig_recoding=recode,skip_invalid_loci=True,
                                   min_partitions=1000)
    splice_indels_skip_invalid.annotate_globals(source_file_path = "gs://gnomad-wphu/splice_ai_data/gnomAD_v3.1_SpliceAI_scores-selected/spliceai_scores.masked.gnomad_indel.hg38.vcf.gz")
    
    spliceAi_info_skip_invalid = splice_snps_skip_invalid.union_rows(splice_indels_skip_invalid)
    return spliceAi_info_skip_invalid

def annotate_spliceAi(mt):
    delta_scores = mt.info.SpliceAI[0].split(delim="\\|")[2:6]
    splice_split = mt.info.annotate(
        SpliceAI=hl.map(lambda x: hl.float32(x), delta_scores)
    ).rename({"SpliceAI":"splice_ai"})
    mt = mt.annotate_rows(info=splice_split)

    # Annotate info.max_DS with the max of DS_AG, DS_AL, DS_DG, DS_DL in info.
    # delta_score array is |DS_AG|DS_AL|DS_DG|DS_DL
    consequences = hl.literal(
        ["acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss"]
    )
    mt = mt.annotate_rows(info=mt.info.annotate(max_ds=hl.max(mt.info.splice_ai)))
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            splice_consequence=hl.if_else(
                mt.info.max_ds > 0,
                consequences[mt.info.splice_ai.index(mt.info.max_ds)],
                "no_consequence",
            )
        )
    )

    return mt
#spliceAi_info_skip_invalid = combine_splice_ai()
#spliceAi_info_skip_invalid = annotate_spliceAi(spliceAi_info_skip_invalid)
#spliceAi_info_skip_invalid.write("gs://gnomad-wphu/spliceai-scores-updated.ht", overwrite=True)


# In[8]:


#chr, pos, ref, alt, refAA, altAA, strand_1pos_0neg, trinucleotide_context, UCSC_gene, ExAC_coverage,primateDL_score
def create_primate_ai_info(write = False, path = "", rewrite = False):
    primate_ai= hl.import_table("gs://gnomad-wphu/PrimateAI_scores_v0.2_hg38.tsv.gz",
                                          force=True, comment="#", skip_blank_lines=True,
                                          types={"pos":hl.tint32, 'primateDL_score':hl.tfloat32, 'ExAC_coverage':hl.tfloat32}
                                         )
    primate_ai = primate_ai.annotate_globals(source_file_path="gs://gnomad-wphu/PrimateAI_scores_v0.2_hg38.tsv.gz")
    primate_ai = primate_ai.annotate(locus = hl.locus(primate_ai.chr,primate_ai.pos, reference_genome="GRCh38"))
    primate_ai = primate_ai.annotate(alleles = hl.array([primate_ai.ref, primate_ai.alt]))
    primate_ai = primate_ai.select("locus", "alleles", "primateDL_score", "ExAC_coverage","strand_1pos_0neg","refAA", "altAA", "trinucleotide_context", "UCSC_gene")
    primate_ai = primate_ai.key_by("locus","alleles")
    if primate_ai.n_partitions() < 10:
        primate_ai = primate_ai.repartition(500)
    if write:
        if overwrite:
            primate_ai.write(path, overwrite = rewrite)
        else:
            try:
                primate_ai.write(path)
            except Exception as e: print(e)
    primate_ai.describe()
    return primate_ai


# In[9]:


gnomad_variants = hl.read_table("gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1/gnomad_genomes_v3.1_info.split.ht")
gnomad_indels = hl.read_table("gs://gnomad-wphu/gnomad_indels.ht")
gnomad_variants_v3 = hl.read_table("gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_info.split.ht")
gnomad3_release = hl.read_table("gs://gnomad-public-requester-pays/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht")

#seqr_annotations_tbl = hl.read_table("gs://seqr-reference-data/GRCh38/all_reference_data/v2/combined_reference_data_grch38-2.0.3.ht")
#gnomad_vars_in_seqr = gnomad_variants.semi_join(seqr_annotations_tbl)
#gnomad_vars_not_in_seqr = gnomad_variants.anti_join(seqr_annotations_tbl)

#cadd_tbl = hl.read_table("gs://seqr-reference-data/GRCh38/CADD/CADD_snvs_and_indels.v1.4.ht")
#CADD_indels = create_cadd_indel_info()
CADD_indels = hl.read_table("gs://gnomad-wphu/CADD-v1.6-indels-updated.ht")
gnomad_indels_not_in_CADD = hl.read_table("gs://gnomad-wphu/gnomad-indels-anti-CADD.ht")
gnomad_indels_in_CADD = hl.read_table("gs://gnomad-wphu/gnomad-indels-semi-CADD.ht")
CADD_snps = hl.read_table("gs://gnomad-wphu/CADD-v1.6-SNVs.ht")
cadd = hl.read_table("gs://gnomad-wphu/complete-CADD-v1.6-annotations.ht")

primate_ai_info_old = hl.read_table("gs://seqr-reference-data/GRCh38/primate_ai/PrimateAI_scores_v0.2.liftover_grch38.ht")
primate_ai_info = hl.read_table("gs://gnomad-wphu/primate-ai-info.ht")
gnomad_vars_not_in_primate_ai = hl.read_table("gs://gnomad-wphu/gnomad-vars-anti-primate.ht")


#revel_information = hl.import_table("gs://gnomad-wphu/revel_grch38_all_chromosomes.csv", delimiter=",", types={'hg19_pos':hl.tint,'grch38_pos':hl.tstr,'REVEL': hl.tfloat64})
revel_information_grch38 = hl.read_table("gs://gnomad-wphu/revel_annotations_grch38.ht")
gnomad_vars_not_in_revel = hl.read_table("gs://gnomad-wphu/gnomad-vars-anti-revel.ht")
gnomad_vars_in_revel = hl.read_table("gs://gnomad-wphu/gnomad-vars-semi-revel.ht")

spliceAi_info_skip_invalid = hl.read_matrix_table("gs://gnomad-wphu/spliceai-scores-updated.ht")
spliceAi_snps_skip_invalid = hl.read_matrix_table("gs://gnomad-wphu/splice_snps_skip_invalid.ht")
spliceAi_indels_skip_invalid = hl.read_matrix_table("gs://gnomad-wphu/splice_indels_skip_invalid.ht")

full_annotations = hl.read_table("gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1/gnomad_genomes_v3.1.analyst_annotations.ht")


# In[50]:


#gnomad_vars_in_revel = gnomad_variants.semi_join(revel_information_grch38)
#gnomad_vars_not_in_revel = gnomad_variants.anti_join(revel_information_grch38)

#gnomad_vars_not_in_revel = gnomad_vars_not_in_revel.checkpoint("gs://gnomad-wphu/gnomad-vars-anti-revel.ht")
#gnomad_indels_in_revel = gnomad_vars_in_revel.checkpoint("gs://gnomad-wphu/gnomad-vars-semi-revel.ht")

#gnomad_indels = gnomad_indels.checkpoint("gs://gnomad-wphu/gnomad-indels-v3.1.ht")
#gnomad_indels_not_in_CADD = gnomad_indels.anti_join(CADD_indels)
#gnomad_indels_in_CADD = gnomad_indels.semi_join(CADD_indels)

#gnomad_indels_not_in_CADD = gnomad_indels_not_in_CADD.checkpoint("gs://gnomad-wphu/gnomad-indels-anti-CADD.ht")
#gnomad_indels_in_CADD = gnomad_indels_in_CADD.checkpoint("gs://gnomad-wphu/gnomad-indels-semi-CADD.ht")


#gnomad_indels = gnomad_indels.checkpoint("gs://gnomad-wphu/gnomad_indels.ht", overwrite = True)
#gnomad_indels.count()
#gnomad_indels.write("gs://gnomad-wphu/gnomad_indels.ht", overwrite = True)

#gnomad3_release = hl.read_table("gs://gnomad-public-requester-pays/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht")
#gnomad3_release.describe()

#gnomad3_release_indels = gnomad3_release.filter(
                            #hl.is_indel(gnomad3_release.alleles[0], gnomad3_release.alleles[1]))

#gnomad3_release_indels.count()


# In[17]:


def annotate_gnomad_v31(list_of_annotations, gnomad_v31_tbl):
    import re
    global_annotations = {}
    for annotation in list_of_annotations:
        annotation = re.sub(r'(_| |-)',"",annotation.lower())
        print(annotation)
        if annotation == "cadd":
            cadd = hl.read_table("gs://gnomad-wphu/complete-CADD-v1.6-annotations.ht")
            cadd = cadd.transmute(cadd = hl.struct(raw_score = cadd.RawScore,
                                                   phred = cadd.PHRED
                                                  ))
            #CADD = hl.read_table()
            #gnomad_v31_tbl = gnomad_v31_tbl.join(CADD.transmute(CADD = hl.struct(RawScore = CADD.RawScore, PHRED = CADD.PHRED)), how="left")
            gnomad_v31_tbl = gnomad_v31_tbl.join(cadd, how = "left")
            global_annotations["cadd"]="gs://gnomad-wphu/complete-CADD-v1.6-annotations.ht"
            
        elif annotation == "revel":
            revel = hl.read_table("gs://gnomad-wphu/revel_annotations_grch38.ht")
            revel = revel.transmute(revel = hl.struct(revel_score = revel.REVEL,
                                                      ref_aa = revel.aaref,
                                                      alt_aa = revel.aaalt
                                                     ))
            gnomad_v31_tbl = gnomad_v31_tbl.join(revel, how = "left")
            global_annotations["revel"] = "gs://gnomad-wphu/revel_annotations_grch38.ht"
            
        elif annotation == "spliceai":
            spliceai = hl.read_matrix_table("gs://gnomad-wphu/spliceai-scores-updated.ht")
            #spliceai = spliceai.annotate_rows(info = spliceai.info.annotate(
                                                #risd = spliceai.rsid,
                                                #qual = spliceai.qual,
                                                #filters = spliceai.filters
                                                #)
                                             #)
            spliceai = spliceai.rename({"info" : "splice_ai"})
            gnomad_v31_tbl = gnomad_v31_tbl.join(spliceai.make_table(), how = "left")
            global_annotations["splice_ai"] = "gs://gnomad-wphu/spliceai-scores-updated.ht"
            
        elif annotation == "primateai":
            primateai = hl.read_table("gs://gnomad-wphu/primate-ai-info.ht")
            primateai = primateai.transmute(primate_ai = hl.struct(primate_ai_score = primateai.primateDL_score))
            gnomad_v31_tbl = gnomad_v31_tbl.join(primateai, how = "left")
            global_annotations["primate_ai"] = "gs://gnomad-wphu/primate-ai-info.ht"
    gnomad_v31_tbl = gnomad_v31_tbl.select_globals()
    gnomad_v31_tbl = gnomad_v31_tbl.annotate_globals(annotation_file_path=global_annotations)
    gnomad_v31_tbl.describe()
    gnomad_v31_tbl = gnomad_v31_tbl.select("cadd", "revel", "splice_ai", "primate_ai")
    return gnomad_v31_tbl
        
result = annotate_gnomad_v31(["CADD", "REVEL", "SPLICE-AI", "PRIMATE-AI"], gnomad_variants)
result.describe()


# In[18]:


#r = gnomad_variants.anti_join(full_annotations)
#r.show()


# In[15]:


#full_annotations.count()


# In[16]:


#gnomad_variants.count()


# In[46]:


#test = gnomad_variants
#print(test.count())
#revel = hl.read_table("gs://gnomad-wphu/revel_annotations_grch38.ht")
#revel = revel.transmute(revel = hl.struct(revel_score = revel.REVEL,
#                                                      ref_aa = revel.aaref,
#                                                      alt_aa = revel.aaalt
#                                                     ))
#test = test.join(revel, how = "left")
#test.count()


# In[47]:


#r = test.anti_join(gnomad_variants)
#r.count()


# In[50]:


#test.distinct().count()


# In[77]:


#a = full_annotations.collect_by_key()
#b = a.filter(hl.len(a.values)>1)
#full_annotations.semi_join(b).show()


# In[89]:


#print(spliceAi_info_skip_invalid.make_table().count())
#print(spliceAi_info_skip_invalid.make_table().distinct().count())
#a = spliceAi_info_skip_invalid.make_table().collect_by_key()
#b = a.filter(hl.len(a.values)>1)
#spliceAi_info_skip_invalid.make_table().semi_join(b).show()


# In[113]:


#a = spliceAi_snps_skip_invalid.make_table()
#a.filter(hl.is_defined(a.info.SpliceAI)).semi_join(b).show()


# In[86]:


#print(revel_information_grch38.count())
#print(revel_information_grch38.distinct().count())
#a = revel_information_grch38.collect_by_key()
#b = a.filter(hl.len(a.values)>1)
#revel_information_grch38.semi_join(b).show()


# In[82]:


#print(cadd.count())
#print(cadd.distinct().count())


# In[84]:


#print(primate_ai_info.count())
#print(primate_ai_info.distinct().count())
#a = primate_ai_info.collect_by_key()
#b = a.filter(hl.len(a.values)>1)
#primate_ai_info.semi_join(b).show()


# # Splice AI

# In[65]:


spliceAi_info_skip_invalid.describe()


# In[63]:


spliceAi_info_skip_invalid.describe()


# In[22]:


#splice_indels = hl.import_vcf("gs://gnomad-wphu/splice_ai_data/gnomAD_v3.1_SpliceAI_scores-selected/spliceai_scores.masked.gnomad_indel.hg38.vcf.gz",
                                   #force_bgz=True,
                                   #reference_genome='GRCh38', contig_recoding=recode,skip_invalid_loci=True,
                                   #min_partitions=1000)

#splice_snps = hl.import_vcf("gs://gnomad-wphu/splice_ai_data/splice_ai_data/genome_scores_v1.3_ds.20a701bc58ab45b59de2576db79ac8d0/spliceai_scores.raw.snv.hg38.vcf.gz",
                                #force_bgz= True,
                                #min_partitions=10000,
                                #reference_genome='GRCh38', contig_recoding=recode
                               #)
#splice_snps_withAnno = splice_snps.filter_rows(hl.len(splice_snps.info.SpliceAI)>0, keep=True)

#splice_indels = hl.import_vcf("gs://gnomad-wphu/splice_ai_data/gnomAD_v3.1_SpliceAI_scores-selected/spliceai_scores.masked.gnomad_indel.hg38.vcf.gz",
                                   #force_bgz=True,
                                   #reference_genome='GRCh38', contig_recoding=recode, skip_invalid_loci=True,
                                   #min_partitions=1000)


# In[23]:


#print(spliceAi_indels_skip_invalid.count())
#print(spliceAi_snps_skip_invalid.count())
#print(spliceAi_info_skip_invalid.count())


# In[24]:


#spliceAi_info_skipInvalid_tbl = spliceAi_info_skipInvalid.make_table()
#spliceAi_info_skipInvalid_tbl.show(5)
#spliceAi_info_skipInvalid_tbl.describe()
#print(gnomad_variants.anti_join(spliceAi_info_skipInvalid_tbl).count())
#gnomad_variants.show(5)
#spliceAi_snps_skipInvalid.show()
#print(gnomad_variants.count())

#numbers do not match because spliceAi only does intergenic variants; should pull up list of genes and compare missing variants against gene list.


# In[25]:


#spliceAi_genes = hl.import_table("gs://gnomad-wphu/grch38.tsv", types = {"TX_START":hl.tint64, "TX_END":hl.tint64})
#spliceAi_genes = hl.import_table("gs://gnomad-wphu/grch38.tsv")
#spliceAi_genes = spliceAi_genes.transmute(CHROM = "chr" + spliceAi_genes.CHROM)
#spliceAi_genes = spliceAi_genes.annotate(INTERVAL = hl.parse_locus_interval(spliceAi_genes.CHROM + ":"+spliceAi_genes.TX_START + "-" +spliceAi_genes.TX_END, reference_genome="GRCh38"))


# In[26]:


#gnomad_variants_in_genomic_regions = hl.filter_intervals(gnomad_variants, spliceAi_genes.INTERVAL.collect())
#gnomad_variants_in_genomic_regions.count()


# In[27]:


#vars_in_genes_not_found = hl.filter_intervals(gnomad_variants.anti_join(spliceAi_info_skipInvalid_tbl),spliceAi_genes.INTERVAL.collect())


# In[28]:


#gnomad_indels_not_in_CADD.aggregate(hl.agg.collect_as_set(gnomad_indels_not_in_CADD.locus.contig))


#gnomad_indels_not_in_CADD_mitoStripped.aggregate(
    #hl.agg.collect_as_set(gnomad_indels_not_in_CADD_mitoStripped.locus.contig))
#gnomad_indels_not_in_CADD_mitoStripped.count()


#export_indels_for_CADD_analysis(gnomad_indels_not_in_CADD_mitoStripped, "gs://gnomad-wphu/CADD_indels_for_upload_mitoStripped.vcf.bgz")


# In[29]:


#test = gnomad_indels_not_in_CADD.head(150000)
#export_indels_for_CADD_analysis(test, "gs://gnomad-wphu/CADD_indels_for_upload_test_90000.vcf")


# In[30]:


#CADD_snps_test = load_CADD(path = "gs://gnomad-wphu/whole_genome_SNVs.tsv", n_partitions=5000)
#CADD_snps_test.show()

#CADD_snps = CADD_snps.checkpoint("gs://gnomad-wphu/CADD-v1.6-SNVs.ht")


# In[31]:


#cadd_ht = hl.import_table("gs://gnomad-wphu/whole_genome_SNVs.tsv", comment="#", no_header=True, min_partitions=5000)


# # Primate AI

# In[32]:


#primate_ai_info.count()


# In[34]:


#primate_ai_info_new = create_primate_ai_info()


# In[35]:


#primate_ai_info_new = primate_ai_info_new.repartition(500)
#primate_ai_info_new = primate_ai_info_new.checkpoint("gs://gnomad-wphu/primate-ai-info.ht", overwrite=False)


# # Combined Annotations

# In[36]:


#gnomad_variants.count()
#gnomad_variants.describe()
#gnomad_variants.show()


# In[37]:


#seqr_annotations_tbl.show()
#seqr_annotations_tbl.describe()


# In[11]:


#result = result.checkpoint("gs://gnomad-wphu/gnomad-3.1-all-variants-annotations.ht")
result.write("gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1/gnomad_genomes_v3.1.analyst_annotations.ht", overwrite=True)


# In[73]:


result.show()


# # Other Things

# In[ ]:


#missing_cadd_scores = result.filter(hl.is_defined(result.cadd.phred),keep=False)
#missing_cadd_scores.count()


# In[ ]:


#missing_cadd_scores.filter(missing_cadd_scores.locus.contig=="chrM", keep=False).count()


# In[ ]:


#primate_ai_info.ExAC_coverage.summarize()


# In[ ]:


#print(gnomad_indels_in_CADD.count())
#print(gnomad_indels_not_in_CADD.count())
#print(gnomad_indels.count())

#misisng_cadd_scores = missing_cadd_scores.write("gs://gnomad-tmp/missing_cadd_scores.ht", overwrite=True)
#missing_cadd_scores = hl.read_table("gs://gnomad-tmp/missing_cadd_scores.ht")


# In[ ]:


#missing_cadd_scores_no_M = missing_cadd_scores.filter(missing_cadd_scores.locus.contig=="chrM",keep=False)
#print(missing_cadd_scores.count())
#print(missing_cadd_scores_no_M.count())

#print(missing_cadd_scores.filter(missing_cadd_scores.locus.contig=="chrM",keep=True).count())
#print(gnomad_variants.filter(gnomad_variants.locus.contig=="chrM", keep=True).count())
#missing_cadd_scores_no_M.alleles.collect()


# In[ ]:


#hl.len(missing_cadd_scores_no_M.old_alleles).summarize()
#missing_cadd_scores_no_M.count()
#hl.len(missing_cadd_scores_no_M.old_alleles).show()
#missing_cadd_scores_no_M.select(missing_cadd_scores_no_M.old_alleles).write("gs://gnomad-tmp/gnomad-31-complex-variants.ht")
#missing_cadd_scores_no_M.select(missing_cadd_scores_no_M.old_alleles).export("gs://gnomad-tmp/gnomad-31-complex-variants.tsv")


# In[ ]:


#missing_cadd_scores_diff_key = missing_cadd_scores_no_M.select(missing_cadd_scores_no_M.old_alleles)
#missing_cadd_scores_diff_key = missing_cadd_scores_diff_key.key_by(missing_cadd_scores_diff_key.old_alleles)
#missing_cadd_scores_diff_key = missing_cadd_scores_diff_key.collect_by_key()


# In[ ]:


#hl.len(missing_cadd_scores_diff_key.values)>1


# In[ ]:


#len(set([frozenset(x) for x in missing_cadd_scores_no_M.old_alleles.collect()]))


# In[ ]:


#missing_cadd_scores_diff_key.export("gs://gnomad-tmp/gnomad-variants-complex-by-old-alleles.tsv")


# In[ ]:


#missing_cadd_scores_diff_key.count()


# In[ ]:


#export_for_CADD_analysis(missing_cadd_scores_no_M, "gs://gnomad-wphu/gnomad-large-deletions-cadd.vcf")


# In[ ]:


#example = ['AGGCTGACCTCTGTCCGCGTGGGAGGGGCCGGTGTGAGGCAAGGGGCTCAGGCTGACCTCTGTCCGCGTGGGAGGGGCCGGTGTGAGGCAAGGGGCTCAGGCTGACCTCTGTCCGCGTGGGAGGGGCCGGGGTGAGGCAAGGGCTCACACTGACCTCTCTCAGCGTGGGAGGGGCCGGTGTGAGGCAAGGGGCTCGGGCTGACCTCTCTCAGCGTGGGAGGGGCCGGTGTGAGGCAAGGGGCTCGGGCTGACCTCTCTCAGCGTGGGAGGGGCCGGTGTGAGGCAAGGGGCTCG', 'G']
#(~(hl.is_snp(example[0], example[1]))).show()
#hl.is_indel(example[0], example[1]).show()


# In[ ]:


#not_snps = gnomad_variants.filter(hl.is_snp(gnomad_variants.alleles[0],gnomad_variants.alleles[1]),keep=False)
#not_snps.count()


# In[ ]:


#gnomad_indels.count()


# In[ ]:


#questionable = not_snps.anti_join(gnomad_indels)
#questionable.count()


# In[ ]:


#print(missing_cadd_scores_no_M.semi_join(questionable).count())
#print(missing_cadd_scores_no_M.filter(hl.is_complex(missing_cadd_scores_no_M.alleles[0], missing_cadd_scores_no_M.alleles[1])).count())


# In[ ]:


#hl.filter_intervals(questionable,spliceAi_genes.INTERVAL.collect()).select().export("gs://gnomad-wphu/missing-splice-ai-variants.tsv")


# In[ ]:


#gnomad_variants.join(CADD_snps.transmute(CADD = hl.struct(RawScore = CADD_snps.RawScore, PHRED = CADD_snps.PHRED)), how="left").show()
#CADD_snps.transmute(CADD = hl.struct(RawScore = CADD_snps.RawScore, PHRED = CADD_snps.PHRED))


# In[ ]:


#spliceAi_info_skip_invalid.describe()
#spliceAi_info_skip_invalid.rsid.show()
#spliceAi_info_skip_invalid.qual.show()
#spliceAi_info_skip_invalid.filters.show()
#spliceAi_info_skip_invalid.info.show()

#reannotate = spliceAi_info_skip_invalid.annotate_rows(info = spliceAi_info_skip_invalid.info.annotate(
#                                                risd = spliceAi_info_skip_invalid.rsid,
#                                                qual = spliceAi_info_skip_invalid.qual,
#                                                filters = spliceAi_info_skip_invalid.filters
#                                                )
#                                        )

#reannotate = reannotate.rename({"info" : "SPLICE_AI"})
#reannotate.describe()


# In[ ]:


#gnomadgnomad3_release.filter(hl.is_indel(gnomad3_release.alleles[0], gnomad3_release.alleles[1]))


# In[ ]:


#gnomad_31_CADD_indels = hl.import_table("gs://gnomad-julia/gnomad_v3.1/cadd_indel_output/CADD_gnomad3.1_scores_*.tsv.gz", comment="#",no_header=True, force_bgz=True)
#gnomad_31_CADD_indels = load_CADD("gs://gnomad-julia/gnomad_v3.1/cadd_indel_output/CADD_gnomad3.1_scores_*.tsv.gz", n_partitions=3000, force_bgz=True)


# In[ ]:


#gnomad_31_CADD_complex = load_CADD("gs://gnomad-julia/gnomad_v3.1/cadd_indel_output_extra/CADD_gnomad3.1_scores_*.tsv.gz",n_partitions=3000, force_bgz=True)
#gnomad_31_CADD_complex.describe()
#gnomad_31_CADD_complex.checkpoint("gs://gnomad-wphu/CADD-1.6-gnomad-complex-variants.ht")


# In[ ]:


#gnomad_31_CADD_complex.count()


# In[ ]:


#gnomad_31_CADD_indels.count()


# In[ ]:


#gnomad_indels_not_in_CADD.filter(gnomad_indels_not_in_CADD.locus.contig=="chrM",keep=False).count()


# In[ ]:


#gnomad_indels_not_in_CADD.count()


# In[ ]:


#gnomad_indels_not_in_CADD.filter(gnomad_indels_not_in_CADD.locus.contig=="chrM",keep=False).anti_join(gnomad_31_CADD_indels).count()


# In[ ]:


#gnomad_31_CADD_indels = gnomad_31_CADD_indels.checkpoint("gs://gnomad-wphu/CADD-indels-gnomad.3.1.ht")


# In[ ]:


#unified = make_unified_CADD()


# In[ ]:


#unified.describe()
#unified.show()


# In[ ]:


#unified.filePathSource["snvs"].show()


# In[ ]:


#unified = unified.checkpoint("gs://gnomad-wphu/complete-CADD-v1.6-annotations.ht")


# In[ ]:


#unified.count()


# In[ ]:


#import requests
#r = requests.get(url = "https://spliceailookup-api.broadinstitute.org/spliceai/?hg=38&variant=chr8-140300616-T-G")
#print(r.json())


# In[ ]:


#test = missing_cadd_scores_no_M.select().head(5).collect()


# In[ ]:


#print(test)


# In[ ]:


#test[0]


# In[ ]:


#primate_ai_info.show()
#primate_ai_info_old.show()


# In[ ]:


#spliceAi_info_skip_invalid.describe()


# In[ ]:


#test = hl.import_vcf("gs://seqr-reference-data/GRCh38/primate_ai/PrimateAI_scores_v0.2.liftover_grch38.vcf.gz", force_bgz=True)


# In[ ]:


#test.describe()


# In[ ]:


#test2 = hl.read('gs://seqr-reference-data/GRCh38/primate_ai/PrimateAI_scores_v0.2.liftover_grch38.vds')


# In[ ]:


#spliceAi_info_skip_invalid.qual.summarize()


# In[7]:


dbNSFP = hl.import_table("gs://gnomad-wphu/dbNSFP/dbNSFP4.1/dbNSFP4.1a_variant.chr*.gz", force=True, missing='.')


# In[22]:


dbNSFP.REVEL_score.summarize()


# In[12]:


revel_information_grch38.show()


# In[18]:


revel_information_grch38.count()


# In[19]:


gnomad_variants.count()


# In[20]:


revel_information_grch38.anti_join(gnomad_variants).count()


# In[8]:


dbNSFP = dbNSFP.repartition(500)
#dbNSFP.n_partitions()
dbNSFP = dbNSFP.checkpoint("gs://gnomad-wphu/dbNSFP.ht")


# In[25]:


82077491-84013093


# In[28]:


revel_information = hl.import_table("gs://gnomad-wphu/revel_grch38_all_chromosomes.csv", delimiter=",", types={'hg19_pos':hl.tint,'grch38_pos':hl.tstr,'REVEL': hl.tfloat64})
revel_information.count()


# In[29]:


82100677-82077491


# In[12]:


dbNSFP.VEST4_score.summarize()


# In[9]:


full_annotations.revel.summarize()


# In[13]:


revel_information_grch38.anti_join(full_annotations).show()


# In[14]:


full_annotations.filter(hl.is_snp(full_annotations.alleles[0], full_annotations.alleles[1])).count()


# In[15]:


revel_information_grch38.anti_join(full_annotations).count()


# In[18]:


revel_information_grch38.aggregate(hl.agg.group_by(revel_information_grch38.locus.contig, hl.agg.count()))


# In[34]:


unique_loci = revel_information_grch38.key_by().key_by("locus").distinct()
a = unique_loci.aggregate(hl.agg.group_by(unique_loci.locus.contig, hl.agg.count()))


# In[35]:


unique_loci_gnomad = gnomad_variants.key_by().key_by("locus").distinct()
b = unique_loci_gnomad.aggregate(hl.agg.group_by(unique_loci_gnomad.locus.contig, hl.agg.count()))


# In[36]:


type(a)
type(b)


# In[7]:


full_annotations.summarize()


# In[8]:


spliceAi_info_skip_invalid.


# In[10]:


full_annotations.filter(hl.is_snp(full_annotations.alleles[0],full_annotations.alleles[1])).splice_ai.summarize()


# In[11]:


full_annotations.filter(~hl.is_snp(full_annotations.alleles[0],full_annotations.alleles[1])).splice_ai.summarize()


# In[12]:


full_annotations.filter(hl.is_indel(full_annotations.alleles[0],full_annotations.alleles[1])).splice_ai.summarize()


# In[13]:


full_annotations.filter(hl.is_indel(full_annotations.alleles[0],full_annotations.alleles[1])).splice_ai.show()


# In[20]:


a = full_annotations.filter(hl.is_indel(full_annotations.alleles[0],full_annotations.alleles[1]))
b = a.filter(hl.is_missing(a.splice_ai))
b.show()
c = a.filter(hl.is_defined(a.splice_ai))
c.show()
b.splice_ai.summarize()
c.splice_ai.summarize()


# In[21]:


print(b.splice_ai.summarize())
print(c.splice_ai.summarize())


# In[22]:


spliceAi_indels_skip_invalid.summarize()


# In[9]:


full_annotations.summarize()


# In[8]:


primate_ai_test= hl.import_table("gs://gnomad-wphu/PrimateAI_scores_v0.2_hg38.tsv.gz",
                                          force_bgz=True, comment="#", skip_blank_lines=True,
                                          types={"pos":hl.tint32, 'primateDL_score':hl.tfloat32, 'ExAC_coverage':hl.tfloat32}
                                         )
primate_ai_test.show()


# In[ ]:




