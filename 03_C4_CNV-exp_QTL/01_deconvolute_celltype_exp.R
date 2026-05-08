################################################################################
# We use Velmeshev's single-cell RNA dataset spanning developmental stages as a reference for deconvolution.
# SCRIPT: Deconvolution of Developmental Brain Bulk RNA-seq
#
# PURPOSE:
# 1. Pre-process bulk RNA-seq data (stage1, stage2, stage3)
#    - Normalize raw counts to log2(TMM-CPM+1); This is different with eQTL matrix we used, because INT is not an expression, but a rank.
#    - Correct for technical batch effects using ComBat
# 2. Merge all 3 processed bulk datasets
# 3. Load and process the scRNA-seq reference profile (Velmeshev et al.)
# 4. Filter for high-confidence marker genes
# 5. Match genes between the bulk (ENSG) and reference (Symbol) datasets
# 6. Run dtangle to estimate cell type proportions
# 7. Run bMIND to estimate cell-type-specific gene expression
#
################################################################################
library(edgeR)
library(limma)
library(sva)     
library(dplyr)
library(ggplot2)
library(pheatmap)
library(dtangle)
library(MIND)
library(gtools)   

################################################################################
#
# Step2: PRE-PROCESS INDIVIDUAL BULK RNA-seq STAGES
#
# We normalize to log2(TMM-CPM+1) and correct for *technical* batches.
# We DO NOT use Quantile Normalization (QN), as that is for eQTL analysis
# and invalidates deconvolution math.
#
################################################################################

# Define base paths
base_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq"
output_path <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/bulk_data"

process_bulk_stage <- function(stage_name) {
  
  message(paste("--- Processing", stage_name, "---"))
  
  # Define file paths
  stage_dir <- file.path(base_path, stage_name)
  count_file <- file.path(stage_dir, "raw_count.txt")
  meta_file <- file.path(stage_dir, "meta.txt")
  samples_file <- file.path(stage_dir, "samples_to_keep.txt") # Here we did expression outlier; and RIN and sex filter.
  genes_file <- file.path(stage_dir, "genes_to_keep.txt")  # keep genes with at least 1 CPM in at least 25% of the individuals; keep only autosomal genes.
  
  # Load data
  count <- read.table(count_file, head=T, row.names=1, check.names=F)
  meta <- read.table(meta_file, head=T, row.names=1)
  samples <- read.table(samples_file, head=F)
  genes <- read.table(genes_file, head=F)
  
  # Filter counts
  count <- count[genes$V1, samples$V1]
  
  # Re-order meta to match count columns
  meta <- meta[colnames(count), ]
  
  # Calculate log2(TMM-CPM+1)
  message("  Calculating log2(TMM-CPM+1)...")
  dge <- DGEList(counts = count)
  dge <- calcNormFactors(dge, method = "TMM")
  count_logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  
  # Correct for technical batch effects using ComBat
  # We assume 'batch' is a column in your 'meta.txt'
  message("  Running ComBat for technical batch correction...")
  # Check if there is a batch variable to correct for
  if ("batch" %in% colnames(meta) && length(unique(meta$batch)) > 1) {
    expr_combat <- ComBat(dat=count_logcpm, 
                          batch=meta$batch, 
                          mod = NULL, 
                          par.prior = TRUE, 
                          prior.plots = FALSE)
  } else {
    message("  Skipping ComBat (no 'batch' variable found or only 1 batch).")
    expr_combat <- count_logcpm
  }

  # Define output file path
  output_file <- file.path(output_path, paste0(stage_name, ".logTMM.ComBat.txt"))
  
  # Write the processed file
  message(paste("  Writing processed file to:", output_file))
  write.table(expr_combat, file=output_file, sep="\t", row.names=T, quote=F, col.names=NA)
  
  return(output_file)
}

# --- Run processing for all 3 stages ---
# This creates the correct input files for deconvolution
process_bulk_stage("stage1")
process_bulk_stage("stage2")
process_bulk_stage("stage3")

message("--- Bulk data pre-processing complete. ---")


################################################################################
#
# Step3: LOAD & MERGE PROCESSED BULK DATA
#
################################################################################

# Load the *corrected* files we just created
d1_path <- file.path(output_path, "stage1.logTMM.ComBat.txt")
d2_path <- file.path(output_path, "stage2.logTMM.ComBat.txt")
d3_path <- file.path(output_path, "stage3.logTMM.ComBat.txt")

d1 <- read.table(d1_path, header=T, sep='\t', row.names=1, check.names=F)
d2 <- read.table(d2_path, header=T, sep='\t', row.names=1, check.names=F)
d3 <- read.table(d3_path, header=T, sep='\t', row.names=1, check.names=F)

# Load the *metadata* files (these are for plotting and analysis)
# NOTE: These are different from the "meta.txt" files used for ComBat
meta1 <- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/covariatesToUse.PCAforQTLfromINT.txt", header=T, sep='\t', row.names=1,check.names=F) # 3genotypePCs+RIN+Sex+AGE(GW)+expressionPCs
meta2 <- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage2/covariatesToUse.PCAforQTLfromINT.txt", header=T, sep='\t', row.names=1,check.names=F)
meta3 <- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage3/covariatesToUse.PCAforQTLfromINT.txt", header=T, sep='\t', row.names=1,check.names=F)

# Find common genes across all 3 datasets
common_genes <- Reduce(intersect, list(rownames(d1), rownames(d2), rownames(d3)))
message(paste("Found", length(common_genes), "common genes across all 3 stages."))

# Filter and merge expression data
d1 <- d1[common_genes, ]
d2 <- d2[common_genes, ]
d3 <- d3[common_genes, ]
expr <- cbind(d1, d2, d3) # This is the final merged bulk matrix

# Create a combined metadata data.frame
# We re-order the metadata to match the expression matrix
meta1 <- t(meta1[,colnames(d1)])
meta2 <- t(meta2[,colnames(d2)])
meta3 <- t(meta3[,colnames(d3)])
meta1<-as.data.frame(meta1)
meta2<-as.data.frame(meta2)
meta3<-as.data.frame(meta3)
# Combine metadata rows
meta_combined <- data.frame(
  sample = c(rownames(meta1), rownames(meta2), rownames(meta3)),
  AGE = c(meta1$GW, meta2$Age, meta3$Age), # Combine GW (gestational week) and Age
  SEX = c(meta1$sex, meta2$sex, meta3$sex),
  RIN = c(meta1$RIN, meta2$RIN, meta3$RIN),
  STAGE = rep(c("stage1", "stage2", "stage3"), c(ncol(d1), ncol(d2), ncol(d3)))
)
rownames(meta_combined) <- meta_combined$sample

# Save the combined objects for future use
save(expr, meta_combined, file=file.path(output_path, 'combined_devbulk_logTMM.RData'))
message("Combined bulk data saved to 'combined_devbulk_logTMM.RData'")


################################################################################
#
# Step4: PROCESS REFERENCE DATA & MARKER GENES
#
################################################################################

# Set working directory to where ref files are
ref_dir <- "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/ref/"
setwd(ref_dir)

# Load reference raw counts
ref_counts <- read.csv("Velmeshev_sum2000_counts.csv", header=T, sep=',', row.names=1)

# ---
# CRITICAL FIX: Correct way to calculate log2(TMM-CPM+1)
# The line `ref<-log2(cpm(ref,method='TMM',log=F)+1)` is incorrect
# as `cpm()` does not use the `method` argument.
# ---
message("Processing reference data to log2(TMM-CPM+1)...")
dge_ref <- DGEList(counts = ref_counts)
dge_ref <- calcNormFactors(dge_ref, method = "TMM")
ref <- cpm(dge_ref, log = TRUE, prior.count = 1) # This is the correct "logTMM"

# Load pre-computed marker genes
mg <- read.csv("velmeshev_mg.csv", header=T, sep=',', row.names=1)

# Filter for strong marker genes
mg1 <- subset(mg, FDR < 0.05 & logFC > 1)
mg1 <- mg1[order(mg1$celltype), ]
mg1$gene <- rownames(mg1)
write.table(mg1,file="/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/resource/cell_maker/cell.mg.velmeshev.csv",col.names=NA,quote=F)
# Get top 50 markers (or all, as in your script)
# Your script has two lines, the second one overrides the first.
# top50_per_celltype <- mg1 %>% ...
# ref2 <- ref[top50_per_celltype$gene,]
ref2 <- ref[mg1$gene, ] # We use all good markers, as in your final script
#ref2 <- ref2[, c(1,2,3,4,5,6,7,9,10,14)] # Filtering to 10 cell types
ref2 <- ref2[!duplicated(rownames(ref2)), ] # Ensure ref marker genes are unique
ref2<-ref[mg1$gene,]
#ref2<-ref2[,c(1,2,3,4,5,6,7,9,10,14)]
ref_pre<- ref[,c(9,10,12,13,14,16)]
#[1] "GLIALPROG"      "prenatal_ExNeu" "prenatal_MG"    "prenatal_OPC"  
#[5] "prenatal_IN"    "prenatal_VASC" 
ref_post<- ref[,c(2,3,4,5,6,7)]
#[1] "postnatal_VASC"  "postnatal_ExNeu" "OL"              "postnatal_IN"   
#[5] "postnatal_MG"    "postnatal_AST" 
message(paste("Filtered reference profile to", nrow(ref2), "marker genes."))


################################################################################
#
# Step5: MATCH BULK (ENSG) & REFERENCE (SYMBOL)
#
################################################################################

# Load annotation file
annot_file <- "/gpfs/hpc/home/chenchao/hanc/database/hg38/gene_symbol.txt"
annot <- read.table(annot_file, header=T, sep='\t')

message("Matching bulk ENSG IDs to reference Gene Symbols...")
bulk <- expr # Use the merged bulk matrix from Stage 3
annot_filt <- annot[annot$geneid %in% rownames(bulk), ]
annot_filt <- annot_filt[!duplicated(annot_filt$geneid), ] # Ensure 1:1 mapping
ensg_ids_in_bulk <- rownames(bulk)
gene_symbols_for_bulk <- annot$gene_name[match(ensg_ids_in_bulk, annot$geneid)]

mapping_df <- data.frame(
  ensg_id = ensg_ids_in_bulk,
  gene_symbol = gene_symbols_for_bulk,
  stringsAsFactors = FALSE
)

# a) Remove rows where Symbol is NA (ENSG not found in annotation)
mapping_df <- mapping_df[!is.na(mapping_df$gene_symbol), ]

# b) Remove rows where Symbol is not in our reference matrix
mapping_df <- mapping_df[mapping_df$gene_symbol %in% rownames(ref2), ]

# c) De-duplicate: Keep only the *first* occurrence of each gene_symbol

message("Resolving duplicates by keeping the first occurrence...")
keep_indices <- !duplicated(mapping_df$gene_symbol)
mapping_final <- mapping_df[keep_indices, ]

message(paste("  Resolved to", nrow(mapping_final), "unique genes."))

# Filter bulk matrix to only the "winning" ENSG IDs
bulk2 <- bulk[mapping_final$ensg_id, ]
rownames(bulk2) <- mapping_final$gene_symbol
# Filter reference matrix to only the "winning" Gene Symbols
ref2 <- ref2[mapping_final$gene_symbol, ]
ref2 <- ref2[rownames(bulk2), ]


# 2. Find common genes (now both are Gene Symbols)
common_genes_final <- intersect(rownames(ref2), rownames(bulk2))
message(paste("Found", length(common_genes_final), "common markers in bulk and ref."))

# 3. Filter both matrices to this final, matched gene set
bulk2 <- bulk2[common_genes_final, ]
ref2 <- ref2[common_genes_final, ]

# 4. Ensure gene order is identical
ref2 <- ref2[rownames(bulk2), ]

message("Bulk and reference matrices are now matched.")


################################################################################
#
# Step6: RUN DECONVOLUTION (dtangle) & VISUALIZE
#
################################################################################

message("Running dtangle to estimate proportions...")
library(dtangle)
prop_estimates <- dtangle(Y = t(as.matrix(bulk2)), references = t(ref2))$estimates
write.csv(prop_estimates, "/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/devbulk_prop.csv")

message("Plotting proportion results...")
library(gridExtra)
library(ggplot2)

prop <- read.csv("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/devbulk_prop.csv", header=T, sep=',', row.names=1)

df_plot <- cbind(meta_combined[rownames(prop), ], prop)

df_plot$AGE <- as.numeric(df_plot$AGE)
setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution")
save.image("devbulk_deconv.RData")


# Generate plots
#pdf("deconvolution_proportion_plots.pdf", width = 5, height = 5)
#ggplot(df_plot, aes(x=AGE, y=postnatal_IN)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=postnatal_ExNeu)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=prenatal_ExNeu)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=prenatal_IN)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=postnatal_AST)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=postnatal_MG)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=postnatal_OPC)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=postnatal_VASC)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=OL)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#ggplot(df_plot, aes(x=AGE, y=GLIALPROG)) + geom_point() + geom_smooth() + facet_wrap(~STAGE, scale='free_x') + theme_bw()
#dev.off()

# Step7: RUN DECONVOLUTION (bMIND) 
#
################################################################################
message("Loading required libraries...")
suppressPackageStartupMessages(library(MIND))
suppressPackageStartupMessages(library(edgeR)) 
suppressPackageStartupMessages(library(limma))

setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution")

message("Loading inputs from  step6 (bulk2, ref2, prop)...")
if (file.exists("bulk_data/combined_devbulk_logTMM.RData")) {
  load("bulk_data/combined_devbulk_logTMM.RData")
} else {
  stop("ERROR: Input file 'bulk_data/combined_devbulk_logTMM.RData' not found. Please create it first.")
}
message("Preparing data for bMIND...")
prop <- read.csv("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/post_devbulk_prop.csv", header=T, sep=',', row.names=1)
ref<-read.csv("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/ref/Velmeshev_sum2000_counts.csv",header=T,sep=',',row.names=1)
annot<-read.table("/gpfs/hpc/home/chenchao/hanc/database/hg38/gene_symbol.txt", header=T, sep='\t')
rownames(expr) <- sub("\\.[0-9]+$", "", rownames(expr))
bulk<- expr[,171:669]
#bulk<- expr[,1:170]
dge_ref <- DGEList(counts = ref)
dge_ref <- calcNormFactors(dge_ref, method = "TMM")
ref <- cpm(dge_ref, log = TRUE, prior.count = 1) 


annot_filt <- annot[annot$geneid %in% rownames(bulk), ]
annot_filt <- annot_filt[!duplicated(annot_filt$geneid), ] # Ensure 1:1 mapping
ensg_ids_in_bulk <- rownames(bulk)
gene_symbols_for_bulk <- annot$gene_name[match(ensg_ids_in_bulk, annot$geneid)]

mapping_df <- data.frame(
  ensg_id = ensg_ids_in_bulk,
  gene_symbol = gene_symbols_for_bulk,
  stringsAsFactors = FALSE
)

mapping_df <- mapping_df[!is.na(mapping_df$gene_symbol), ]
mapping_df <- mapping_df[mapping_df$gene_symbol %in% rownames(ref), ]

message("Resolving duplicates by keeping the first occurrence...")
keep_indices <- !duplicated(mapping_df$gene_symbol)
mapping_final <- mapping_df[keep_indices, ]

message(paste("  Resolved to", nrow(mapping_final), "unique genes."))

bulk2 <- bulk[mapping_final$ensg_id, ]
ref2 <- ref[mapping_final$gene_symbol, ]
rownames(bulk2) <- mapping_final$gene_symbol
ref2 <- ref[rownames(bulk2), ]

common_genes_final <- intersect(rownames(ref2), rownames(bulk2))
message(paste("Found", length(common_genes_final), "common markers in bulk and ref."))

bulk2 <- bulk2[common_genes_final, ]
ref2 <- ref2[common_genes_final, ]
ref2 <- ref2[rownames(bulk2), ]
ref2<- ref2[,colnames(prop)]
message("Bulk and reference matrices are now matched.")

bulk_mind <- bulk2
ref_mind <- ref2
prop_mind <- prop

bulk_mind <- bulk_mind[, rownames(prop_mind)]
sampleid <- rownames(prop_mind)

##For adult,we will skip this non-zero variance filter

ref_sds <- apply(ref_mind, 1, sd)
keep_genes <- ref_sds != 0
bulk_mind <- bulk_mind[keep_genes, ]
ref_mind <- ref_mind[keep_genes, ]
message(paste("Filtered to", nrow(bulk_mind), "non-zero variance genes for bMIND."))
##For adult,we will skip this non-zero variance filter

colnames(bulk_mind) <- rownames(prop_mind) <- paste0('s', 1:nrow(prop_mind))

message("Running bMIND with 8 cores...")
deconv2 = bMIND(as.matrix(bulk_mind), 
                frac = as.matrix(prop_mind), 
                profile = as.matrix(ref_mind), 
                ncore = 8)
message("bMIND run complete.")

message("Processing and saving bMIND results...")
bmind_result <- deconv2$A
for(i in 1:ncol(prop_mind)){
  colnames(bmind_result[, i, ]) <- sampleid 
}

save(bmind_result, file="devbulk_bmind_result.RData")

setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution")
#load("devbulk_bmind_result.RData")
load("pre2_devbulk_bmind_result.RData")
meta<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/all_meta/stage1.filtered.meta.txt",head=T,check.names=F,row.names=1)
setwd("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/04_analysis/05deconvolution/removeversionIDmapping/cts_expr")
if (length(dimnames(bmind_result)[[3]]) != nrow(meta)) {
  stop("Sample size mismatch! Check meta or bmind_result.")
}

dimnames(bmind_result)[[3]] <- meta$sample

meta1<- meta[meta$STAGE=="stage1",]
meta2<- meta[meta$STAGE=="stage2",]
meta3<- meta[meta$STAGE=="stage3",]

cell_types <- dimnames(bmind_result)[[2]]

stage1_expr_list <- list()
stage2_expr_list <- list()
stage3_expr_list <- list()


for (cell in cell_types) {
  
  stage1_expr <- bmind_result[, cell, rownames(meta1)]
  stage2_expr <- bmind_result[, cell, rownames(meta2)]
  stage3_expr <- bmind_result[, cell, rownames(meta3)]
  

  stage1_expr_list[[cell]] <- stage1_expr
  stage2_expr_list[[cell]] <- stage2_expr
  stage3_expr_list[[cell]] <- stage3_expr

write.table(stage1_expr, file = paste0("stage1_", cell, ".txt"),quote=F,col.names=NA)
write.table(stage2_expr,    file = paste0("stage2_", cell, ".txt"),quote=F,col.names=NA)
write.table(stage3_expr,     file = paste0("stage3_", cell, ".txt"),quote=F,col.names=NA)
#saveRDS(stage1_expr, file = paste0("stage1_", cell, ".rds"))
#saveRDS(stage2_expr,    file = paste0("stage2_", cell, ".rds"))
#saveRDS(stage3_expr,     file = paste0("stage3_", cell, ".rds"))
}
message("Done!")

adult_expr_list <- list()

for (cell in cell_types) {
  
  adult_expr <- bmind_result[, cell, rownames(meta)]
  adult_expr_list[[cell]] <- adult_expr

write.table(adult_expr,     file = paste0("adult_", cell, ".txt"),quote=F,col.names=NA)
#saveRDS(stage1_expr, file = paste0("stage1_", cell, ".rds"))
#saveRDS(stage2_expr,    file = paste0("stage2_", cell, ".rds"))
#saveRDS(stage3_expr,     file = paste0("stage3_", cell, ".rds"))
}
