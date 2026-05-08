#!/bin/bash
#SBATCH --job-name=pre_impute
#SBATCH --output=logs/pre_impute_%j.out
#SBATCH --error=logs/pre_impute_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

export SENTIEON_INSTALL_DIR=/gpfs/software/sourcecode/sentieon/sentieon-genomics-202503
export SENTIEON_LICENSE=/gpfs/software/sourcecode/sentieon/lic/Central_South_University_Furong_Lab_cluster.lic
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH

JOINT_DIR=/gpfs/hpc/home/chenchao/share_group_folder/new_chb/adult_gvcf_sentieon_alt_aware/joint_calling
INPUT_VCF="${JOINT_DIR}/MHC.joint.filtered.vcf.gz"
OUTPUT_VCF="${JOINT_DIR}/MHC.upload_ready.vcf.gz"
REFERENCE=/gpfs/hpc/home/chenchao/hanc/share_group_folder/data/gatk_resource_bundle/Homo_sapiens_assembly38.fasta
#echo "chr6 6" > ${JOINT_DIR}/rename_map.txt
export BCFTOOLS_PLUGINS=/gpfs/hpc/home/chenchao/hanc/share_group_folder/tools/bcftools-1.18/plugins

bcftools view -f PASS $INPUT_VCF | \
    bcftools norm -m -any -f $REFERENCE | \
    bcftools view -m2 -M2 -v snps | \
    bcftools +fill-tags - -- -t HWE,F_MISSING | \
    bcftools view -i 'F_MISSING<0.1 && HWE>1e-10 && MAF>0.01' | \
    #bcftools annotate --rename-chrs ${JOINT_DIR}/rename_map.txt | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
    bcftools norm -d all -O z -o $OUTPUT_VCF

bcftools index -t $OUTPUT_VCF
echo "Ready for upload: $OUTPUT_VCF"
