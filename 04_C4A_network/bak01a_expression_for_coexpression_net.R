# raw count: /gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/raw_count.txt
# samples keep : /gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/samples_to_keep.txt
# meta: /gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/meta.filteredsamplessex.txt

##--- Read Prenatal RNA quantification data
#path: /gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq

library(edgeR)
library(reshape)
library(preprocessCore)
library(sva)
library(ggplot2)
library(limma)
library(data.table)

count<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/raw_count.txt",head=T,row.names=1,check.names=F)
meta<-read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/meta.filteredsamplessex.txt",head=T,row.names=1)
samples<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/stage1/samples_to_keep.txt",head=F) #170 samples after RNA preprocessing
annot = fread("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/resource/gencode.v40.annotation.gene.bed")
##--- TMM normalization  
count<- count[,samples$V1]
dge <- DGEList(counts=count)
dge <- calcNormFactors(dge,method="TMM")

# log2 CPM transformation
expr <- cpm(dge,log=TRUE,prior.count=1)

##--- Gene filtering 
# Use raw CPM values for expression filtering
cpm_expr <- cpm(dge)
genes_to_keep <- rowSums(cpm_expr >= 0.1)>= round(0.25*ncol(cpm_expr))
log2cpm.fgene <- expr[genes_to_keep,]
##--- Remove chromosome M and Y genes
annot <- annot[annot$gene_id %in% rownames(log2cpm.fgene),]
# Remove mitochondrial and Y chromosome genes
brainExpressedNoMT<-annot[annot$chr!="chrM"&annot$chr!="chrY"&annot$chr!="chrX",]
log2cpm.fgene <- expr[brainExpressedNoMT$gene_id,]
##---selection Quantile normalization  
expr.qn <- log2cpm.fgene
#expr.qn <- normalize.quantiles(as.matrix(log2cpm.fgene, copy=T))
#rownames(expr.qn) <- rownames(log2cpm.fgene)
#colnames(expr.qn) <- colnames(log2cpm.fgene)

##--- ComBat batch correction 
meta <- meta[colnames(expr.qn),]
# Preserve biological covariates during batch correction
# Prenatal:
# sex, RIN, gestational week
#
# Postnatal:
# sex, RIN, age
expr.qn.comb <- ComBat(dat=expr.qn,batch=meta$batch,mod=NULL,par.prior=TRUE,prior.plots=FALSE)
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2")
write.table(expr.qn.comb,file="prenatal.logTMM.ComBat.txt",sep="\t",row.names=T,quote=F,col.names=NA)

##--- Find covariates
### PCAforQTL R package (Top 3 expression PCs)

suppressPackageStartupMessages({
  library(PCAForQTL)
  library(ggcorrplot)
})

be_perms <- 25
alpha <- 0.05
known_cov_names <- c("sex", "RIN", "GW")
r2_cutoff<- 0.9

expr_t <- t(expr.qn.comb)
prcomp_result <- prcomp(expr_t, center = TRUE, scale. = TRUE)
exp_pcs <- prcomp_result$x
be_result <- PCAForQTL::runBE(expr_t, B = be_perms, alpha = alpha, mc.cores = 1)
k_be <- be_result$numOfPCsChosen
message(paste("BE algorithm chose", k_be, "PCs as the upper limit."))
pdf("prenatal_expressionPC.pdf", width = 8, height = 6)
PCAForQTL::makeScreePlot(prcomp_result, labels = c("BE"), values = c(k_be), titleText = "Scree Plot")
dev.off()
colnames(meta)<- c("sex","RIN","batch","GW")
meta$sex<- as.numeric(as.factor(meta$sex))
known_covs <- meta[, known_cov_names, drop = FALSE]
exp_pcs_top <- exp_pcs[, 1:k_be, drop = FALSE]
common_samples <- intersect(rownames(known_covs), rownames(exp_pcs_top))
known_covs_filtered <- PCAForQTL::filterKnownCovariates(known_covs, exp_pcs_top, unadjustedR2_cutoff = r2_cutoff)
covariates_to_use <- cbind(known_covs_filtered,exp_pcs_top)
message("Final covariate matrix has", ncol(covariates_to_use), "covariates for", nrow(covariates_to_use), "samples.")
write.table(t(covariates_to_use), file = "prenatalcovarites.expressPC_knownforqtl", quote = FALSE, sep = "\t", col.names = NA)

covariates_to_use2<- covariates_to_use[,1:6] #sex  RIN   GW       PC1        PC2        PC3
write.table(t(covariates_to_use2), file = "prenatalcovarites.expressPC_knownfornetwork", quote = FALSE, sep = "\t", col.names = NA)

##--- Read postnatal RNA quantification data
#path: /gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/postnatal

library(edgeR)
library(reshape)
library(preprocessCore)
library(sva)
library(ggplot2)
library(limma)
library(data.table)

count<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/postnatal/raw_count.txt",head=T,row.names=1,check.names=F)
meta<-read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/postnatal/raw_meta.txt",head=F,row.names=1,check.names=F)
colnames(meta)<- c("RIN","sex","disease","Age","batch")
samples<- read.table("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/01_data/RNA-seq/postnatal/samples_to_keep.txt",head=F)
annot = fread("/gpfs/hpc/home/chenchao/hanc/project/2023project_develop_eQTL/resource/gencode.v40.annotation.gene.bed")

##--- TMM normalization  
count<- count[,samples$V1]
dge <- DGEList(counts=count)
dge <- calcNormFactors(dge,method="TMM")

# log2 CPM transformation
expr <- cpm(dge,log=TRUE,prior.count=1)

##--- Gene filtering 
# Use raw CPM values for expression filtering
cpm_expr <- cpm(dge)
# Keep genes expressed at >=1 CPM
# in at least 75% of samples (raw:20%)

genes_to_keep <- rowSums(cpm_expr >=0.1)>= round(0.25*ncol(cpm_expr))
log2cpm.fgene <- expr[genes_to_keep,]

##--- Remove chromosome M and Y genes 
annot <- annot[annot$gene_id %in% rownames(log2cpm.fgene),]
# Remove mitochondrial and Y chromosome genes
brainExpressedNoMT<-annot[annot$chr!="chrM"&annot$chr!="chrY"&annot$chr!="chrX",]
log2cpm.fgene <- expr[brainExpressedNoMT$gene_id,]
##--- Quantile normalization  
#expr.qn <- normalize.quantiles(as.matrix(log2cpm.fgene, copy=T))
#rownames(expr.qn) <- rownames(log2cpm.fgene)
#colnames(expr.qn) <- colnames(log2cpm.fgene)
expr.qn<- log2cpm.fgene
##--- ComBat batch correction 
meta <- meta[colnames(expr.qn),]
# Preserve biological covariates during batch correction
# prenatal:
# sex, RIN, gestational week
#
# Postnatal:
# sex, RIN, age
expr.qn.comb <- ComBat(dat=expr.qn,batch=meta$batch,mod=NULL,par.prior=TRUE,prior.plots=FALSE)
setwd("/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/networkv2")
write.table(expr.qn.comb,file="postnatal.logTMM.ComBat.txt",sep="\t",row.names=T,quote=F,col.names=NA)

##--- Find covariates
### PCAforQTL R package (Top 3 expression PCs)
suppressPackageStartupMessages({
  library(PCAForQTL)
  library(ggcorrplot)
})

be_perms <- 25
alpha <- 0.05
known_cov_names <- c("sex", "RIN", "Age")
r2_cutoff<- 0.9
meta$sex<- as.numeric(as.factor(meta$sex))
expr_t <- t(expr.qn.comb)
prcomp_result <- prcomp(expr_t, center = TRUE, scale. = TRUE)
exp_pcs <- prcomp_result$x
be_result <- PCAForQTL::runBE(expr_t, B = be_perms, alpha = alpha, mc.cores = 1)
k_be <- be_result$numOfPCsChosen
message(paste("BE algorithm chose", k_be, "PCs as the upper limit."))
pdf("postnatal_expressionPC.pdf", width = 8, height = 6)
PCAForQTL::makeScreePlot(prcomp_result, labels = c("BE"), values = c(k_be), titleText = "Scree Plot")
dev.off()
known_covs <- meta[, known_cov_names, drop = FALSE]
exp_pcs_top <- exp_pcs[, 1:k_be, drop = FALSE]
common_samples <- intersect(rownames(known_covs), rownames(exp_pcs_top))
known_covs_filtered <- PCAForQTL::filterKnownCovariates(known_covs, exp_pcs_top, unadjustedR2_cutoff = r2_cutoff)
covariates_to_use <- cbind(known_covs_filtered,exp_pcs_top)
message("Final covariate matrix has", ncol(covariates_to_use), "covariates for", nrow(covariates_to_use), "samples.")
write.table(t(covariates_to_use), file = "postnatalcovarites.expressPC_knownforqtl", quote = FALSE, sep = "\t", col.names = NA)

covariates_to_use2<- covariates_to_use[,1:6] #sex  RIN   Age      PC1        PC2        PC3
write.table(t(covariates_to_use2), file = "postnatalcovarites.expressPC_knownfornetwork", quote = FALSE, sep = "\t", col.names = NA)