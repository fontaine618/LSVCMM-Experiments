library(tidyverse)
library(magrittr)

DATA_PATH = "./dmbt1/data/"
GroupGeneTime_PATH = paste0(DATA_PATH, "GroupGeneTime.csv")
SALIVA_OTU_PATH = paste0(DATA_PATH, "SALIVA_OTU.csv")
SALIVA_TAXO_PATH = paste0(DATA_PATH, "SALIVA_TAXO.csv")
DEMO_PATH = paste0(DATA_PATH, "demo.csv")
META_PATH = paste0(DATA_PATH, "meta.csv")



# ==============================================================================
# LOAD OTU TABLE & SAMPLE METADATA
otu = read.table(file=SALIVA_OTU_PATH, header = TRUE, sep=',')
# Format:
# - label (always 0.03)
# - Group (sample descriptor)
# - numOtus (always 1726)
# - OtuXXXX, between 0001 and 1866 (only 1726)
group = otu$Group
subject_id = as.numeric(gsub(".*?([0-9]+).*","\\1", group))
visit_id = as.numeric(gsub("^[^_]*[_]?([0-9]*).*",'\\1', group))
visit_id[is.na(visit_id)] = 0
otu = otu[, -c(1, 2, 3)]  # drop non-otu columns
meta_sample = data.frame(subject_id=subject_id, visit_id=visit_id)
rm(subject_id, visit_id, group)
# ------------------------------------------------------------------------------



# ==============================================================================
# LOAD SUBJECT METADATA
meta_subject = read.table(file=GroupGeneTime_PATH, header = TRUE, sep=',')
# Format:
# - ID
# - Type (WT, KO)
demo = read.table(file=DEMO_PATH, header=T, sep=",")
# Format:
# - ID
# - Gender
# - Group (W, K)
# - Diagnosis (1, 2, 3)
meta_subject %<>% left_join(demo %>% select(ID, Gender), by="ID")
# ------------------------------------------------------------------------------



# ==============================================================================
# MERGE
meta = meta_sample %>% left_join(meta_subject, by=c("subject_id"="ID"))
# Format:
# - subject_id
# - visit_id
# - Type
# - Gender (M/F)
# - Diagnosis (1/2/3) *dropped & replace by second file
# ------------------------------------------------------------------------------



# ==============================================================================
# UPDATED DIAGNOSIS
meta2 = read.table(file=META_PATH, header = TRUE, sep=',')
meta2 %<>% mutate(
  subject_id = as.numeric(gsub(".*?([0-9]+).*","\\1", Group))
)
meta2 %<>% group_by(subject_id) %>% summarize(Diagnosis=head(Diagnosis2, 1))
meta2 %<>% mutate(SCC=ifelse(Diagnosis=="SCC", "SCC", "HP/CIS"))
meta = meta %>% left_join(meta2, by="subject_id")
# ------------------------------------------------------------------------------



# ==============================================================================
# LOAD TAXONOMY
taxonomy = read.table(SALIVA_TAXO_PATH, header=TRUE, sep=",")
# Format:
# - OTU
# - Size
# - Domain, Phylum, Class, Order, Family, Genus
otus = colnames(otu)
rownames(taxonomy) = taxonomy$OTU
taxonomy %<>% filter(OTU %in% otus)
taxonomy = taxonomy[match(otus, taxonomy$OTU), ]
taxonomy %<>% mutate(
  Domain = paste0("d__", Domain),
  Phylum = paste0("p__", Phylum),
  Class = paste0("c__", Class),
  Order = paste0("o__", Order),
  Family = paste0("f__", Family),
  Genus = paste0("g__", Genus)
)
taxonomy %<>% select(-c(Size, OTU))
for(j in 1:ncol(taxonomy)){
  taxonomy[grepl('unclassified', taxonomy[, j], fixed = TRUE), j] = NA
}
# ------------------------------------------------------------------------------



# ==============================================================================
# TO PHYLOSEQ
otu = phyloseq::otu_table(otu, taxa_are_rows=F)
tax = phyloseq::tax_table(as.matrix(taxonomy))
sam = phyloseq::sample_data(meta)
pseq = phyloseq::phyloseq(otu, tax, sam)
# ------------------------------------------------------------------------------



# ==============================================================================
# CLEANUP
rm(demo, meta, meta_sample, meta_subject, meta2, otu, otus, taxonomy, tax,
   sam, DATA_PATH, DEMO_PATH, GroupGeneTime_PATH, META_PATH, SALIVA_OTU_PATH, j, SALIVA_TAXO_PATH)
# only pseq is kept
# ------------------------------------------------------------------------------
