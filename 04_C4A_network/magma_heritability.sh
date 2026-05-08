# 1. Annotation step (Mapping SNPs to genes)
# /gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/01_C4_exp/magma/
# export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
magma --annotate window=35,10 \
      --snp-loc pgc_scz_eas_magma.snp.loc \
      --gene-loc NCBI37.3.gene.loc \
      --out EAS_annotation


magma --annotate window=35,10 \
      --snp-loc bip_eas_magma.snp.loc \
      --gene-loc NCBI37.3.gene.loc \
      --out EAS_bip_annotation
      

# 2. Gene Analysis (Calculating gene-level P-values for SCZ)
# PGC3 SCZ GWAS is formatted for MAGMA

magma --bfile g1000_eas \
      --pval pgc_scz_eas_magma.pval use=SNP,P ncol=N \
      --gene-annot EAS_bip_annotation.genes.annot \
      --out BIP_EAS_GeneLevel

# 3. Gene-Set Analysis (networks for heritability enrichment)
/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/01_C4_exp/magma/magma --gene-results gwas/SCZ_EAS_GeneLevel.genes.raw \
      --set-annot C4A_MAGMA_0.3RGeneSets.txt \
      --out /gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/downsampling/magma_results/SCZ_Heritability_Results

/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/01_C4_exp/magma/magma --gene-results gwas/BIP_EAS_GeneLevel.genes.raw \
      --set-annot C4A_MAGMA_0.3RGeneSets.txt \
      --out /gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/02_C4_network/downsampling/magma_results/BIP_Heritability_Results
      