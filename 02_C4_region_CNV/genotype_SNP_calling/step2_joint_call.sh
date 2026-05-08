#!/bin/bash
#SBATCH --job-name=sentieon_joint
#SBATCH --output=logs/joint_%j.out
#SBATCH --error=logs/joint_%j.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=CU

set -euxo pipefail


export SENTIEON_INSTALL_DIR=/gpfs/software/sourcecode/sentieon/sentieon-genomics-202503
#export SENTIEON_LICENSE=/gpfs/software/sourcecode/sentieon/lic/Central_South_University_Furong_Lab_cluster.lic
export SENTIEON_LICENSE=12.12.12.201:8990
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH     
REFERENCE=/gpfs/hpc/home/chenchao/hanc/share_group_folder/data/gatk_resource_bundle/Homo_sapiens_assembly38.fasta
NT=${SLURM_CPUS_PER_TASK:-12}
TARGET_REGION="chr6:24000000-34000000"


RESULT_DIR=/gpfs/hpc/home/chenchao/share_group_folder/new_chb/adult_gvcf_sentieon_alt_aware
JOINT_DIR=${RESULT_DIR}/joint_calling
mkdir -p $JOINT_DIR

RAW_VCF="${JOINT_DIR}/MHC.joint.raw.vcf.gz"
FILTERED_VCF="${JOINT_DIR}/MHC.joint.filtered.vcf.gz"

echo "Listing gVCF files..."

find ${RESULT_DIR} -name "*.gvcf.gz" > ${JOINT_DIR}/gvcf.list

FILE_COUNT=$(wc -l < ${JOINT_DIR}/gvcf.list)
echo "Found ${FILE_COUNT} gVCF files."

# if [ "$FILE_COUNT" -lt 242 ]; then echo "Error: Missing gVCFs!"; exit 1; fi
# === 3. Joint Genotyping ===
# Sentieon GVCFtyper is essentially a combination of GATK CombineGVCFs and GenotypeGVCFs.
# -t $NT

echo "Starting Joint Calling..."
if [ ! -f $RAW_VCF ]; then
sentieon driver -r $REFERENCE -t $NT \
        --interval $TARGET_REGION \
        --algo GVCFtyper \
        --emit_mode variant \
        $RAW_VCF \
        $(cat ${JOINT_DIR}/gvcf.list)
fi

# === 4. Hard Filtering (MHC) ===
echo "Starting Hard Filtering with bcftools..."

if [ ! -f $FILTERED_VCF ]; then  
    bcftools filter \
        -e 'QD < 2.0 || MQ < 30.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
        -s "HardFilter" \
        -m "+" \
        -O z -o $FILTERED_VCF \
        $RAW_VCF
    bcftools index -t $FILTERED_VCF
fi
echo "All Done! Final VCF: $FILTERED_VCF"
