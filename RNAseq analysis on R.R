install.packages("colorspace")
BiocManager::install("rtracklayer")
install.packages("WGCNA", dependencies = T)

BiocManager::install("preprocessCore")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("minet")

library(DESeq2)
library(tidyverse)
library(apeglm)
library(pheatmap)
library(rtracklayer)
library(dplyr)
library(stringr)
library(ggplot2)
library(colorspace)
library(clusterProfiler)
library(BiocParallel)
library(readr)
library(ape)
library(matrixStats)
library(WGCNA)
library(impute)

#
##### BASE SETUP ###############################################################
# Read the count table
setwd("~/Postdoc_files/Experiments_and_projects/3_Results/2_Whole_body_organ_by_organ_RNAseq/2_RNAseq_data_and_analysis/4_Full_resequencing_2nd_RNAseq_round/Secondary_Bioinformatics")

counts <- read.table(
  "featureCounts_gene_counts.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "\t",
  comment.char = "#"
)


# Remove featureCounts annotation columns
counts <- counts[, -(1:5)]


# Export the table for future analysis
#write.csv(counts, "counts.csv", row.names = T)


# Clean sample column names
colnames(counts) <- colnames(counts) |>
  basename() |>
  sub("\\.sorted\\.bam$", "", x = _) |>
  sub(".bam$", "", x = _)

# Create a mapping table
gff <- import("xgBioGlab_hap1_braker-FA_20260209.gff3")
tx <- gff[gff$type %in% c("mRNA", "transcript"), ]

get_mcol <- function(gr, field) {
  m <- mcols(gr)
  if (!(field %in% colnames(m))) {
    return(rep(NA_character_, length(gr)))
  }
  
  x <- m[[field]]
  
  # If it is list-like (AtomicList or list), collapse to one string per row
  if (is(x, "List") || is.list(x)) {
    return(vapply(x, function(v) {
      if (length(v) == 0) NA_character_
      else paste(as.character(v), collapse = ";")
    }, FUN.VALUE = character(1)))
  }
  
  # Otherwise coerce to character safely
  as.character(x)
}


anno_tx <- tibble(
  tx_id    = get_mcol(tx, "makerName"),
  product  = get_mcol(tx, "product"),
  Name     = get_mcol(tx, "Name"),
  ID       = get_mcol(tx, "ID"),
) |>
  # Keep only rows that have something useful
  filter(!is.na(tx_id), tx_id != "") |>
  mutate(
    gene = tx_id |>
      str_replace("^([^;]+);.*$", "\\1") |>
      str_replace("\\.t\\d+$", ""),
    # Clean product to remove transcript-variant clutter
    product_clean = product |>
      replace_na("") |>
      str_replace_all("%2C", ",") |> # in case URL-encoding appears
      str_replace_all("\\s+", " ") |> # normalize spaces
      str_remove(",\\s*transcript variant\\s*X\\d+.*$") |> # remove ", transcript variant X1 ..."
      str_remove("\\s*transcript variant\\s*X\\d+.*$") |> # remove "transcript variant X1 ..."
      str_trim(),
    # A flag that says "this row has a genuinely informative product"
    has_product = product_clean != "" &
      !str_detect(product_clean, regex("^hypothetical protein$", ignore_case = TRUE)) &
      !str_detect(product_clean, regex("^uncharacterized protein$", ignore_case = TRUE)),
    has_name = !is.na(Name) & Name != ""
  ) |>
  group_by(gene) |>
  # Pick the best row per gene:
  # 1) rows with product come first
  # 2) among those, keep the shortest 
  arrange(desc(has_product), desc(has_name), nchar(product_clean)) |>
  slice(1) |>
  ungroup() |>
  transmute(
    gene,
    label = case_when(
      has_product ~ product_clean,
      has_name ~ Name,
      TRUE ~ gene
    ),
    Name
  )

lab_vec <- anno_tx$label
names(lab_vec) <- anno_tx$gene

# To check that lab_vec matches count gene name
matched_labels <- lab_vec[rownames(counts)]
sum(!is.na(matched_labels))
length(matched_labels)
head(matched_labels)



check <- counts
check$gene <- rownames(check)
check <- merge(check, anno_tx, "gene")

check <- check[, c(18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]

Bg_CDA : g6281


#
##### SETUP FOR COMPARING APO VS ALL OTHER ORGANS NON COLLAPSED ################
# Reference table for all samples
coldata_full <- tibble(sample = colnames(counts)) |>
  mutate(
    organ = case_when(
      str_detect(sample, "APO")   ~ "APO",
      str_detect(sample, "AG")    ~ "AG",
      str_detect(sample, "CA")    ~ "CA",
      str_detect(sample, "HEAD")  ~ "HEAD",
      str_detect(sample, "HEART") ~ "HEART",
      str_detect(sample, "HEM")   ~ "HEM",
      str_detect(sample, "HEP")   ~ "HEP",
      str_detect(sample, "HG")    ~ "HG",
      str_detect(sample, "INT")   ~ "INT",
      str_detect(sample, "MAN")   ~ "MAN",
      str_detect(sample, "OVO")   ~ "OVO",
      str_detect(sample, "PU")    ~ "PU",
      str_detect(sample, "SK")    ~ "SK",
      str_detect(sample, "STO")   ~ "STO",
      str_detect(sample, "TK")    ~ "TK",
      str_detect(sample, "(_|^)W(_|$)") ~ "W",
      TRUE ~ NA_character_
    ),
    replicate = case_when(
      str_detect(sample, "_R1") ~ "R1",
      str_detect(sample, "_R2") ~ "R2",
      str_detect(sample, "_R3") ~ "R3",
      TRUE ~ NA_character_
    )
  ) |>
  column_to_rownames("sample")

stopifnot(all(rownames(coldata_full) == colnames(counts)))

coldata_full$replicate  <- factor(coldata_full$replicate)

# Optional: order organs (helps plotting consistency)
coldata_full$organ <- factor(
  coldata_full$organ,
  levels = c("APO","HEM","HEP","HG","INT","STO","CA","HEART","HEAD","AG","MAN","OVO","PU","SK","TK","W")
)


# Building a DESeq2 object and then making a VST matrix for heatmaps
keep <- rowSums(counts >= 10) >= 3
counts_f <- counts[keep, ]

dds_full <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData   = coldata_full,
  design    = ~ organ   # add lane if you want: ~ replicate + organ
)

dds_full <- DESeq(dds_full)

vsd_full <- vst(dds_full, blind = FALSE)
expr_full <- assay(vsd_full)  # genes x samples


# TO CHECK 1 replicate Building a DESeq2 object and then making a VST matrix for heatmaps
keep <- rowSums(counts >= 10) >= 1
counts_f <- counts[keep, ]

dds_full <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData   = coldata_full,
  design    = ~ 1
)

dds_full <- estimateSizeFactors(dds_full)

vsd_full <- vst(dds_full, blind = TRUE)
expr_full <- assay(vsd_full)


#
# PCA with the samples ####
pcaData <- plotPCA(vsd_full, intgroup = "organ", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Organ as a factor
pcaData$organ <- factor(pcaData$organ)

n_organs <- nlevels(pcaData$organ)

# Distinct colors
organ_colors <- setNames(
  qualitative_hcl(n_organs, palette = "Dark 3"),
  levels(pcaData$organ)
)

# Distinct shapes
organ_shapes <- setNames(
  c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4, 5, 6, 9, 10, 11, 12, 13)[seq_len(n_organs)],
  levels(pcaData$organ)
)

# Plot
ggplot(pcaData, aes(PC1, PC2, color = organ, shape = organ)) +
  geom_point(size = 5, stroke = 1.1) +
  scale_color_manual(name = "Organ", values = organ_colors, breaks = levels(pcaData$organ)) +
  scale_shape_manual(name = "Organ", values = organ_shapes, breaks = levels(pcaData$organ)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_classic(base_size = 14)


#
# Visualization of all transcripts differentially expressed ####
pheatmap(expr_full, annotation_col = coldata_full["organ"], scale = "row")

#
# Visualization of top expressed genes in APO (expression-based) #####
apo_samples <- rownames(coldata_full)[coldata_full$organ == "APO"]

# Averaging APO samples
apo_mean <- rowMeans(expr_full[, apo_samples, drop = FALSE])

# Setting the number of top genes to map
topN <- 50
top_apo_genes <- names(sort(apo_mean, decreasing = TRUE ))[1:topN]

mat <- expr_full[top_apo_genes, , drop = FALSE]

# Annotating selected genes
new_names <- lab_vec[rownames(mat)]
new_names[is.na(new_names)] <- rownames(mat)[is.na(new_names)]
rownames(mat) <- new_names

# Heatmap
pheatmap(
  mat,
  annotation_col = coldata_full["organ", drop = FALSE],
  scale = "row",
  show_rownames = TRUE,
  main = paste0("Top ", topN, " genes by VST expression in APO")
)


# TO CHECK 
apo_samples <- rownames(coldata_full)[coldata_full$organ == "APO"]

apo_samples <- intersect(apo_samples, colnames(expr_full))
if (length(apo_samples) == 0) {
  stop("No APO samples found in expr_full. Check sample names and organ assignment.")
}

apo_mean <- rowMeans(expr_full[, apo_samples, drop = FALSE])



#
# Visualization of top APO markers (highest DE genes) #####
dds_bin <- dds_full
coldata_bin <- as.data.frame(coldata_full)
coldata_bin$group <- ifelse(coldata_bin$organ == "HEM", "HEM", "Other")
coldata_bin$group <- factor(coldata_bin$group, levels = c("Other","HEM"))

colData(dds_bin)$group <- coldata_bin$group
design(dds_bin) <- ~ group
dds_bin <- DESeq(dds_bin)

res <- results(dds_bin, contrast = c("group","HEM","Other"))
resultsNames(dds_bin)

resS <- lfcShrink(dds_bin, coef = "group_HEM_vs_Other", res = res, type = "apeglm")
res <- as.data.frame(resS) |> rownames_to_column("gene") |> drop_na(padj)

# choose top APO-up genes
topN <- 50
top_apo_up <- res |>
  arrange(padj) |>
  filter(log2FoldChange > 1) |>
  slice_head(n = topN) |>
  pull(gene)

mat <- expr_full[top_apo_up, , drop = FALSE]

# Annotating selected genes
new_names <- lab_vec[rownames(mat)]
new_names[is.na(new_names)] <- rownames(mat)[is.na(new_names)]
rownames(mat) <- new_names

# Heatmap
pheatmap(
  mat,
  annotation_col = coldata_full["organ", drop = FALSE],  # <-- full organ labels here
  scale = "row",
  show_rownames = TRUE,
  main = paste0("Top ", topN, " APO-up genes (APO vs rest), shown across all organs")
)

#
# My own gene list #####
loc_to_name <- c("LOC106073343"="RXR1", 
                 "LOC106050859"="RXR2", 
                 "LOC106062531"="RXR3",
                 "LOC106056956"="MEIS", 
                 "LOC106061842"="Chromobox1", 
                 "LOC106071927"="Chromobox2",
                 "LOC106053770"="Chromobox3", 
                 "LOC106067727"="Xbox", 
                 "LOC106054679"="THYN1",
                 "LOC106071658"="TLR2", 
                 "LOC106057492"="DDX46", 
                 "LOC106060848"="JIP3",
                 "LOC106064986"="JUN", 
                 "LOC106067432"="JAK2", 
                 "LOC106079295"="JAK2",
                 "LOC106079296"="JAK2",
                 "LOC106076146"="RUVBL2", 
                 "LOC106054037"="RUVBL1",
                 "LOC106063678"="FOSRA2")

# From scalop publication
loc_to_name <- c("LOC106078269"="core-binding factor subunit beta-like", 
                 "LOC106051751"="Elf-1", 
                 "LOC106070023"="Etv2",
                 "LOC106054680"="transcriptional activator Myb-like", 
                 "LOC106060713"="notch 1 1", 
                 "LOC106053973"="notch 1 2",
                 "LOC129926483"="notch 1 3", 
                 "LOC106055368"="notch 1 4", 
                 "LOC106079646"="notch 1 5",
                 "LOC106062514"="notch 1 6", 
                 "LOC106056511"="notch 2", 
                 "LOC106069177"="Jagged-1-like 1",
                 "LOC106073644"="Jagged-1-like 2", 
                 "LOC106075199"="Jagged-1-like 3", 
                 "LOC106056685"="Jagged-1-like 4",
                 "LOC106055572"="hairy/enhancer-of-split protein 1-like",
                 "LOC106067432"="JAK2 1",
                 "LOC106079295"="JAK2 2", 
                 "LOC106079296"="JAK2 3",
                 "LOC106066644"="catenin beta-like", 
                 "LOC106062590"="TCF7-like", 
                 "LOC106050505"="TCF7L2-like",
                 "LOC106065589"="interferon regulatory factor 1-like",
                 "LOC106057526"="interferon regulatory factor 8-like 1", 
                 "LOC106060347"="interferon regulatory factor 8-like 2",
                 "LOC106062041"="receptor-type tyrosine-protein phosphatase alpha-like 1",
                 "LOC106069489"="receptor-type tyrosine-protein phosphatase alpha-like 2",
                 "LOC106063873"="receptor-type tyrosine-protein phosphatase alpha-like 3",
                 "LOC106076987"="receptor-type tyrosine-protein phosphatase alpha-like 4",
                 "LOC106064035"="receptor-type tyrosine-protein phosphatase alpha-like 5",
                 "LOC106064527"="receptor-type tyrosine-protein phosphatase alpha-like 6",
                 "LOC106067448"="receptor-type tyrosine-protein phosphatase alpha-like 7",
                 "LOC106065669"="Signal transducer and activator of transcription 5A",
                 "LOC106057899"="Signal transducer and activator of transcription 5B 1",
                 "LOC106065668"="Signal transducer and activator of transcription 5B 2")

# Genes associated
loc_to_name <- c("LOC106063233"="GATA", #GATA-binding factor 3-like
                 "LOC106069584"="CDC42 1", #homologs nbr5_activated CDC42 kinase 1-like nbr6_CDC42 small effector protein 2-A-like
                 "LOC106063738"="CDC42 2",
                 "LOC106066573"="CDC42 3", 
                 "LOC106072916"="CDC42 4", 
                 "LOC106057399"="CDC42 5",
                 "LOC106071043"="CDC42 6", 
                 "LOC106065739"="HOX 1", #HOXC5 HOX3
                 "LOC106056572"="HOX 2", 
                 "LOC106073122"="CEBPα 1", #CCAAT/enhancer-binding protein alpha-like
                 "LOC106056509"="CEBPα 2",
                 "LOC106054708"="ZBTB 1", #ZBTB17 / 49
                 "LOC106068700"="ZBTB 2", 
                 "LOC106059411"="CD38", #ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase
                 "LOC106080292"="BCL11A", #B-cell lymphoma/leukemia 11A-like
                 "LOC106056054"="MLF1") #myeloid leukemia factor 1-like

# Genes indirectly associated
loc_to_name <- c("LOC106067548"="LIM domain-binding protein 2-like", 
                 "LOC106056619"="NF-kappa-B-activating protein-like", 
                 "LOC106071361"="NF-kappa-B inhibitor-like protein 1",
                 "LOC106076256"="NF-kappa-B inhibitor-interacting Ras-like protein 1", 
                 "LOC106060509"="NF-kappa-B essential modulator-like", 
                 "LOC106055393"="NF-kappa-B inhibitor alpha-like 1",
                 "LOC106050623"="NF-kappa-B inhibitor alpha-like 2", 
                 "LOC106058661"="nuclear factor NF-kappa-B p100 subunit-like")

# APO markers
loc_to_name <- c("LOC106073370"="cAMP", 
                 "LOC106054492"="P450", 
                 "LOC106051983"="mADP",
                 "LOC106073064"="LOC")

# Ancestral
loc_to_name <- c("LOC106067580"="PMS1 protein homolog 1-like", # Cellular organism
                 "LOC106078509"="tRNA wybutosine-synthesizing protein 4-like", 
                 "LOC106063801"="tRNA wybutosine-synthesizing protein 3 homolog",
                 "LOC106065493"="tRNA wybutosine-synthesizing protein 2 homolog", 
                 "LOC106073027"="tRNA wybutosine-synthesizing protein 5-like",
                 "LOC106066063"="Peptidylprolyl_isomerase_domain_and_WD_repeat-containing_protein_1", 
                 "LOC106079168"="Probable_ATP-dependent_RNA_helicase_DDX27", 
                 "LOC106076146"="RuvB-like_2", 
                 "LOC106055973"="GTPase Era, mitochondrial-like",
                 "LOC106068791"="elongation factor Ts, mitochondrial-like", 
                 "LOC106071979"="nucleolar protein 56-like", 
                 "LOC106064670"="DNA replication licensing factor mcm2-like",
                 "LOC106073001"="presequence protease, mitochondrial-like", 
                 "LOC106069890"="DNA replication licensing factor mcm4-A", 
                 "LOC106053831"="polyribonucleotide 5'-hydroxyl-kinase Clp1-like",
                 "LOC129922746"="nucleolar protein 58-like 1",
                 "LOC129922747"="nucleolar protein 58-like 2",
                 "LOC106076945"="nucleolar protein 58-like 3", 
                 "LOC129923744"="nucleolar protein 58-like 4", 
                 "LOC106058876"="DNA-directed RNA polymerase III subunit RPC8-like", 
                 "LOC106076568"="probable ribosome biogenesis protein RLP24", 
                 "LOC106078378"="T-complex protein 1 subunit gamma-like",
                 "LOC106079169"="flap endonuclease 1-like",
                 "LOC106071476"="replication factor C subunit 1-like", 
                 "LOC129926441"="probable replication factor C subunit 1",
                 "LOC106054678"="centrosomal protein of 135 kDa-like", # Eucaryotes
                 "LOC106064884"="general transcription factor IIF subunit 1-like",
                 "LOC106062293"="apoptosis inhibitor 5-like",
                 "LOC106070307"="origin recognition complex subunit 3-like",
                 "LOC106053744"="probable ATP-dependent RNA helicase DDX23",
                 "LOC106061005"="splicing factor 3A subunit 1-like",
                 "LOC106069304"="lupus La protein homolog", #
                 "LOC106064234"="RRP12-like protein", #
                 "LOC106078718"="splicing factor 3B subunit 4-like",
                 "LOC106072343"="U3 small nucleolar RNA-interacting protein 2-like",
                 "LOC106064762"="eukaryotic translation initiation factor 3 subunit B-like", 
                 "LOC106074582"="ESF1 homolog", 
                 "LOC106059859"="eIF-2-alpha kinase activator GCN1-like",
                 "LOC106054037"="ruvB-like 1", 
                 "LOC106068954"="KRR1_small_subunit_processome_component_homolog", 
                 "LOC106057338"="replication factor C subunit 3-like",
                 "LOC106050472"="pre-mRNA-processing-splicing factor 8", 
                 "LOC106066765"="spliceosome RNA helicase Ddx39b", 
                 "LOC106058586"="cell division cycle protein 20 homolog 1",
                 "LOC106072659"="cell division cycle protein 20 homolog 2", 
                 "LOC106061466"="coiled-coil domain-containing protein 134-like", # Filozoan
                 "LOC106068338"="tudor domain-containing protein 5-like", # Metazoans
                 "LOC106077917"="protein virilizer homolog", 
                 "LOC106062761"="heterogeneous nuclear ribonucleoprotein K-like", 
                 "LOC106053627"="acylglycerol kinase, mitochondrial-like", #Opistiokonts
                 "LOC129927271"="39S ribosomal protein L38, mitochondrial-like",
                 "LOC106067135"="serrate RNA effector molecule homolog A-like", 
                 "LOC106078721"="TELO2-interacting protein 1 homolog")

# CDA
loc_to_name <- c("LOC106054955"="cytidine deaminase-like 1",
                 "LOC106065887"="cytidine deaminase-like 2",
                 "LOC106054687" = "True cytidine deaminase")

# GATA
loc_to_name <- c("LOC106070829" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129925975" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC106062953" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129928160" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129928382" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129928412" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129921715" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129923153" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129923437" = "GATA zinc finger domain-containing protein 14-like", 
                 "LOC129926267" = "GATA zinc finger domain-containing protein 15-like", 
                 "LOC129926451" = "GATA zinc finger domain-containing protein 15-like", 
                 "LOC129923160" = "GATA zinc finger domain-containing protein 6-like", 
                 "LOC106071477" = "transcription factor GATA-4-like", 
                 "LOC106063233" = "GATA-binding factor 3-like", 
                 "LOC106051138" = "GATA zinc finger domain-containing protein 1-like")


"LOC106059411"="CD38" #ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase
"LOC106051632"="CD133 (PROM1)" #prominin-1-like
"LOC106053228"="CD13" #aminopeptidase N-like
"LOC106062333"="CD13" #aminopeptidase N-like
"LOC106072368"="CD13" #aminopeptidase N-like
"LOC106074934"="CD13" #aminopeptidase N-like
"LOC106075181"="CD13" #aminopeptidase N-like
"LOC106075904"="CD13" #aminopeptidase N-like
"LOC106079441"="CD13" #aminopeptidase N-like
"LOC106080169"="CD13" #aminopeptidase N-like
"LOC106050756"="CD71" #transferrin receptor protein 1-like
"LOC106050923"="PRDM1" #B lymphocyte-induced maturation protein 1
"LOC106066758"="CD10" #neprilysin-like
"LOC106057406"="NEP1" #neprilysin-1-like
"LOC106058204"="NEP1" #neprilysin-1-like
"LOC106060346"="NEP1" #neprilysin-1-like
"LOC106063401"="NEP1" #neprilysin-1-like
"LOC106077909"="NEP1" #neprilysin-1-like
"LOC106079514"="NEP1" #neprilysin-1-like
"LOC129924726"="NEP1" #neprilysin-1-like
"LOC106064536"="NEP11" #neprilysin-11-like
"LOC106055189"="CD56" #neural cell adhesion molecule 1-like
"LOC106055200"="CD56" #neural cell adhesion molecule 1-like
"LOC106069356"="CD56" #neural cell adhesion molecule 1-like
"LOC106062657"="CD56" #neural cell adhesion molecule 1-like
"LOC106072809"="CD56" #neural cell adhesion molecule 1A-like
"LOC106052465"="CD56" #neural cell adhesion molecule L1-like

"LOC106063233"="GATA3" #GATA-binding factor 3-like
"LOC106062514"="NOTCH1" #neurogenic locus notch homolog protein 1-like
"LOC106053973"="NOTCH1" #neurogenic locus notch homolog protein 1-like
"LOC106060713"="NOTCH1" #neurogenic locus notch homolog protein 1-like
"LOC106079646"="NOTCH1" #neurogenic locus notch homolog protein 1-like
"LOC106055368"="NOTCH1" #neurogenic locus notch homolog protein 1-like
"LOC129926483"="NOTCH1" #neurogenic locus notch homolog protein 1-like
"LOC106056286"="TCF3" #transcription factor E2F3
"LOC106063541"="ikaros" #ikaros family zinc finger protein-like
"LOC106080292"="BCL11A" #B-cell lymphoma/leukemia 11A
"LOC106072317"="Pax-5" #paired box protein Pax-5
"LOC106056509"="CEBPA" #CCAAT/enhancer-binding protein alpha
"LOC106073122"="CEBPA" #CCAAT/enhancer-binding protein alpha
"LOC106061336"="EGR2" #early growth response protein 2
"LOC106059910"="c-Myc" #c-Myc-binding protein
"LOC106051138"="GATA1" #GATA zinc finger domain-containing protein 1
"LOC106078790"="FOXO" #forkhead box protein O



vsd_full <- vst(dds_full, blind = FALSE)
expr_full <- assay(vsd_full)

genes <- names(loc_to_name)  

setdiff(genes, rownames(expr_full))
genes_present <- intersect(genes, rownames(expr_full))
length(genes_present)
genes_present

mat <- expr_full[genes_present, , drop = FALSE]

rownames(mat) <- loc_to_name[rownames(mat)]  

pheatmap(
  mat,
  scale = "row",
  annotation_col = as.data.frame(colData(dds_full)[, "organ", drop = FALSE])
)


#
##### GO enrichment with fgsea/GSEA ############################################################

# Read GO term table found on BioMart
go_raw <- read.delim(
  "bg_go_biomart.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(go_raw)

# Build the gene2go table
gene2go <- go_raw %>%
  filter(
    !is.na(`GO term accession`),
    !is.na(`NCBI Reference Sequences ID`)
  ) %>%
  transmute(
    gene = `NCBI Reference Sequences ID`,  # LOC IDs
    go   = `GO term accession`
  ) %>%
  distinct()

bg <- rownames(dds_full)

gene2go_bg <- gene2go %>%
  filter(gene %in% bg)

TERM2GENE_go <- gene2go_bg %>%
  transmute(term = go, gene = gene) %>%
  distinct()

res_APO_HEM <- as.data.frame(
  results(dds_full, contrast = c("organ", "APO", "HEM"))
) %>%
  rownames_to_column("gene") %>%
  drop_na(log2FoldChange)

gene_ranks <- res_APO_HEM %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  with(setNames(log2FoldChange, gene))

# avoids serialization warnings
register(SerialParam())  

gene_ranks <- sort(gene_ranks, decreasing = TRUE)

gsea_go <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = TERM2GENE_go,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0
)

# Convert to a data frame
gsea_df <- as.data.frame(gsea_go)

# Extract and replace GO number with GO terms
go_id2name <- go_raw %>%
  select(`GO term accession`, `GO term name`, `GO domain`) %>%
  distinct() %>%
  rename(
    GO = `GO term accession`,
    GO_name = `GO term name`,
    GO_domain = `GO domain`
  )

gsea_df <- gsea_df %>%
  left_join(go_id2name, by = c("ID" = "GO"))

# Filter on adjusted p-value
gsea_sig <- gsea_df %>%
  filter(p.adjust < 0.05)

# Split by over or under-representation
gsea_APO_up <- gsea_sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES))

gsea_APO_down <- gsea_sig %>%
  filter(NES < 0) %>%
  arrange(NES)



plot_df <- bind_rows(
  gsea_APO_up  %>% slice_head(n = 15),
  gsea_APO_down %>% slice_head(n = 15)
) %>%
  mutate(Direction = ifelse(NES > 0, "APO-enriched", "HEM-enriched"))

ggplot(plot_df,
       aes(x = NES,
           y = reorder(GO_name, NES),
           size = setSize,
           color = p.adjust)) +
  geom_point() +
  scale_color_viridis_c(trans = "log10", direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Direction, scales = "free_y") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "GO term",
    size = "Gene set size",
    color = "Adjusted p-value",
    title = "GO term enrichment (GSEA): APO vs HEM"
  ) +
  theme_classic(base_size = 12)


# Filter by GO domain
# Biological processes
go_bp <- go_raw %>%
  select(`GO term accession`, `GO domain`) %>%
  distinct()

bp_terms <- go_bp %>%
  filter(`GO domain` == "biological_process") %>%
  pull(`GO term accession`)

gsea_BP <- gsea_sig %>%
  filter(ID %in% bp_terms)

gsea_APO_BP_up <- gsea_BP %>%
  filter(NES > 0) %>%
  arrange(desc(NES))%>%
  distinct(GO_name, .keep_all = T)

gsea_APO_BP_down <- gsea_BP %>%
  filter(NES < 0) %>%
  arrange(NES)%>%
  distinct(GO_name, .keep_all = T)


plot_df_BP <- bind_rows(
  gsea_APO_BP_up  %>% slice_head(n = 100),
  gsea_APO_BP_down %>% slice_head(n = 100)
) %>%
  mutate(Direction = ifelse(NES > 0, "APO-enriched", "HEM-enriched"))

ggplot(plot_df_BP,
       aes(x = NES,
           y = reorder(GO_name, NES),
           size = setSize,
           color = p.adjust)) +
  geom_point() +
  scale_color_viridis_c(trans = "log10", direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Direction, scales = "free_y") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "GO term",
    size = "Gene set size",
    color = "Adjusted p-value",
    title = "Biological Process GO term enrichment (GSEA): APO vs HEM"
  ) +
  theme_classic(base_size = 12)

# Molecular functions
go_mf <- go_raw %>%
  select(`GO term accession`, `GO domain`) %>%
  distinct()

mf_terms <- go_mf %>%
  filter(`GO domain` == "molecular_function") %>%
  pull(`GO term accession`)

gsea_MF <- gsea_sig %>%
  filter(ID %in% mf_terms)

gsea_APO_MF_up <- gsea_MF %>%
  filter(NES > 0) %>%
  arrange(desc(NES))%>%
  distinct(GO_name, .keep_all = T)

gsea_APO_MF_down <- gsea_MF %>%
  filter(NES < 0) %>%
  arrange(NES)%>%
  distinct(GO_name, .keep_all = T)


plot_df_MF <- bind_rows(
  gsea_APO_MF_up  %>% slice_head(n = 15),
  gsea_APO_MF_down %>% slice_head(n = 15)
) %>%
  mutate(Direction = ifelse(NES > 0, "APO-enriched", "HEM-enriched"))

ggplot(plot_df_MF,
       aes(x = NES,
           y = reorder(GO_name, NES),
           size = setSize,
           color = p.adjust)) +
  geom_point() +
  scale_color_viridis_c(trans = "log10", direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Direction, scales = "free_y") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "GO term",
    size = "Gene set size",
    color = "Adjusted p-value",
    title = "Molecular Function GO term enrichment (GSEA): APO vs HEM"
  ) +
  theme_classic(base_size = 12)


# Cellular Component
go_cc <- go_raw %>%
  select(`GO term accession`, `GO domain`) %>%
  distinct()

cc_terms <- go_cc %>%
  filter(`GO domain` == "cellular_component") %>%
  pull(`GO term accession`)

gsea_CC <- gsea_sig %>%
  filter(ID %in% cc_terms)

gsea_APO_CC_up <- gsea_CC %>%
  filter(NES > 0) %>%
  arrange(desc(NES))%>%
  distinct(GO_name, .keep_all = T)

gsea_APO_CC_down <- gsea_CC %>%
  filter(NES < 0) %>%
  arrange(NES)%>%
  distinct(GO_name, .keep_all = T)


plot_df_CC <- bind_rows(
  gsea_APO_CC_up  %>% slice_head(n = 15),
  gsea_APO_CC_down %>% slice_head(n = 15)
) %>%
  mutate(Direction = ifelse(NES > 0, "APO-enriched", "HEM-enriched"))

ggplot(plot_df_CC,
       aes(x = NES,
           y = reorder(GO_name, NES),
           size = setSize,
           color = p.adjust)) +
  geom_point() +
  scale_color_viridis_c(trans = "log10", direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Direction, scales = "free_y") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "GO term",
    size = "Gene set size",
    color = "Adjusted p-value",
    title = "Cellular Component GO term enrichment (GSEA): APO vs HEM"
  ) +
  theme_classic(base_size = 12)




# APO vs all other organs: GO × contrast heatmap

organs <- setdiff(levels(coldata_full$organ), "APO")

gsea_list <- lapply(organs, function(o) {
  
  res <- as.data.frame(
    results(dds_full, contrast = c("organ", "APO", o))
  ) %>%
    rownames_to_column("gene") %>%
    drop_na(log2FoldChange)
  
  ranks <- res %>%
    arrange(desc(log2FoldChange)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    with(setNames(log2FoldChange, gene))
  
  ranks <- sort(ranks, decreasing = TRUE)
  
  gsea <- GSEA(
    geneList = ranks,
    TERM2GENE = TERM2GENE_go,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0
  )
  
  as.data.frame(gsea) %>%
    select(ID, NES, p.adjust) %>%
    mutate(contrast = paste0("APO_vs_", o))
})

gsea_all <- bind_rows(gsea_list)

# Adding GO names and restrict to BP
gsea_all <- gsea_all %>%
  left_join(go_id2name, by = c("ID" = "GO")) %>%
  filter(GO_domain == "cellular_component")

#Building a GO × contrast NES matrix
go_keep <- gsea_all %>%
  filter(p.adjust < 0.05) %>%
  pull(GO_name) %>%
  unique()

heatmap_df <- gsea_all %>%
  filter(GO_name %in% go_keep) %>%
  select(GO_name, contrast, NES)

# Plot
heatmap_mat <- heatmap_df %>%
  pivot_wider(names_from = contrast, values_from = NES) %>%
  column_to_rownames("GO_name") %>%
  as.matrix()

heatmap_mat[!is.finite(heatmap_mat)] <- 0

pheatmap(
  heatmap_mat,
  scale = "none",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  fontsize_row = 8
)

#
##### GO enrichment with GO_MWU (RBGOA) ############################################################

# Read GO term table found on BioMart
go_raw <- read.delim(
  "bg_go_biomart.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(go_raw)

# Build the gene2go table
gene2go <- go_raw %>%
  filter(
    !is.na(`GO term accession`),
    !is.na(`NCBI Reference Sequences ID`)
  ) %>%
  transmute(
    gene = `NCBI Reference Sequences ID`,  # LOC IDs
    go   = `GO term accession`
  ) %>%
  distinct()

bg <- rownames(dds_full) ##########!!!!!!!!########!!!!!!!"########!!!!!!!"#######!!!!!! DDS ET PAS DDS_FULL

gene2go_bg <- gene2go %>%
  filter(gene %in% bg)

# Build GO_MWU GO table
gomwu_go <- gene2go_bg %>%
  group_by(gene) %>%
  summarise(GO = paste(sort(unique(go)), collapse = ";"), .groups = "drop")

# GO_MWU expects TAB-delimited, 2 columns, no header is usually safest
write_tsv(gomwu_go, "bg_LOC2GO_gomwu.tab", col_names = FALSE)


res_APO_HEM <- as.data.frame(
  results(dds_full, contrast = c("organ", "APO", "HEM"))
) %>%
  rownames_to_column("gene") %>%
  drop_na(log2FoldChange)


gomwu_meas <- res_APO_HEM %>%
  select(gene, log2FoldChange) %>%
  filter(gene %in% bg) %>%
  distinct()

write_csv(gomwu_meas, "APO_vs_HEM_logFC_gomwu.csv")

# To continue with GO_MWU, use the GO_MWU script in the GO_MWU folder
# Specify : 
# input="APO_vs_HEM_logFC_gomwu.csv"
# goAnnotations="bg_LOC2GO_gomwu.tab"

##### WGCNA FOR TFs ############################################################
vsd <- vst(dds_full, blind = FALSE)
expr <- assay(vsd)  # genes x samples

gene_var <- rowVars(expr)
expr_f <- expr[gene_var > quantile(gene_var, 0.5), ]  # example: keep top 50% variable

# Sample QC for WGCNA
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

datExpr <- t(expr_f)  # WGCNA expects samples x genes

gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main="Sample clustering", sub="", xlab="")

# Choose soft-threshold power
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers)

# Construct network and modules
net <- blockwiseModules(
  datExpr,
  power = 11,                 # replace with chosen power
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# Module gene membership table
genes <- colnames(datExpr)
mod_df <- data.frame(
  gene = genes,
  module = moduleColors,
  stringsAsFactors = FALSE
)

# Compute kME
kME <- signedKME(datExpr, MEs)  # returns matrix genes × modules
kME_df <- data.frame(gene = rownames(kME), kME, check.names = FALSE)

mod_df <- mod_df %>%
  left_join(kME_df, by = "gene")

# Ensure your GO mapping is in the right format
go_id2name <- go_raw %>%
  select(`GO term accession`, `GO term name`, `GO domain`) %>%
  distinct() %>%
  rename(GO = `GO term accession`,
         GO_name = `GO term name`,
         GO_domain = `GO domain`)

# Recompute module eigengenes using color labels
MEs_color <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs_color <- orderMEs(MEs_color)

# Compute kME against these MEs (now columns are aligned with colors)
kME_color <- signedKME(datExpr, MEs_color)

# Build a clean mod_df with these kME columns
mod_df <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  stringsAsFactors = FALSE
) %>%
  left_join(
    data.frame(gene = rownames(kME_color), kME_color, check.names = FALSE),
    by = "gene"
  )

grep("^kME", colnames(mod_df), value = TRUE)[1:10]








# GO_MWU on module membership (kME)
# Build the GO_MWU inputs for one module
run_gomwu_for_module <- function(module_name, mod_df, gene2go_bg,
                                 outdir = paste0("GO_MWU_", module_name)) {
  
  dir.create(outdir, showWarnings = FALSE)
  
  # GO annotations
  gomwu_go <- gene2go_bg %>%
    group_by(gene) %>%
    summarise(GO = paste(sort(unique(go)), collapse = ";"), .groups = "drop")
  
  write.table(gomwu_go,
              file.path(outdir, "gene2go.tab"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Find the right kME column
  kme_candidates <- c(
    paste0("kME", module_name),
    paste0("kME.ME", module_name),
    paste0("kME_ME", module_name)
  )
  
  kme_col <- kme_candidates[kme_candidates %in% colnames(mod_df)][1]
  if (is.na(kme_col)) {
    stop("Cannot find a kME column for module '", module_name,
         "'. Available kME columns include: ",
         paste(grep("^kME", colnames(mod_df), value = TRUE), collapse = ", "))
  }
  
  meas <- mod_df %>%
    select(gene, all_of(kme_col)) %>%
    rename(measure = all_of(kme_col)) %>%
    tidyr::drop_na(measure)
  
  write.table(meas,
              file.path(outdir, "measure.csv"),
              sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  invisible(list(outdir = outdir, kme_col = kme_col))
}

run_gomwu_for_module("blue", mod_df, gene2go_bg)





# Relate modules to APO trait
meta <- as.data.frame(colData(dds_full))
meta$APO <- as.numeric(meta$organ == "APO")

moduleTraitCor <- cor(MEs, meta$APO, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))




# Functional annotation of modules (hematopoiesis)
mod <- "blue"  # example
modGenes <- names(datExpr)[moduleColors == mod]

# Identify candidate transcription factors within APO modules



































