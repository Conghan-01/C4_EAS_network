#!/bin/bash
#SBATCH --job-name=sentieon_array
#SBATCH --output=logs/sentieon_%A_%a.out  
#SBATCH --error=logs/sentieon_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8                
#SBATCH --mem=32G                        
#SBATCH --time=12:00:00
#SBATCH --partition=CU
#SBATCH --array=0-1   
              
set -euxo pipefail

# === setting paths ===
export SENTIEON_INSTALL_DIR=/gpfs/software/sourcecode/sentieon/sentieon-genomics-202503
export SENTIEON_LICENSE=/gpfs/software/sourcecode/sentieon/lic/Central_South_University_Furong_Lab_cluster.lic
export PATH=$SENTIEON_INSTALL_DIR/bin:$PATH

REFERENCE=/gpfs/hpc/home/chenchao/hanc/share_group_folder/data/gatk_resource_bundle/Homo_sapiens_assembly38.fasta
DBSNP=/gpfs/hpc/home/chenchao/share_group_folder/data/genomes/refresh/Homo_sapiens_assembly38.dbsnp138.vcf.gz
KNOWN_1000G_INDELS=/gpfs/hpc/home/chenchao/share_group_folder/data/genomes/refresh/Homo_sapiens_assembly38.known_indels.vcf.gz
KNOWN_MILLS_INDELS=/gpfs/hpc/home/chenchao/share_group_folder/data/genomes/refresh/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

NT=${SLURM_CPUS_PER_TASK:-8}
TARGET_REGION="chr6:24000000-34000000"
SAMPLE_PARENT_DIR=/gpfs/hpc/home/chenchao/share_group_folder/new_chb/adult_fastq
RESULT_DIR=/gpfs/hpc/home/chenchao/share_group_folder/new_chb/adult_gvcf_sentieon_alt_aware


SAMPLES=(${SAMPLE_PARENT_DIR}/PTB522)

SAMPLE_PATH=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
if [ -z "$SAMPLE_PATH" ]; then
    echo "Error: No sample found for ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

SAMPLE_NAME=$(basename $SAMPLE_PATH)
echo "Running sample: $SAMPLE_NAME (Array Job ID: $SLURM_ARRAY_TASK_ID)"

FASTQ_1=$(ls $SAMPLE_PATH/PTB522_1.fq.gz)
FASTQ_2=$(ls $SAMPLE_PATH/PTB522_2.fq.gz)
WORKDIR=$RESULT_DIR/$SAMPLE_NAME
mkdir -p $WORKDIR
cd $WORKDIR

# S1-alignment and sort
if [ ! -f sorted.bam ]; then
    sentieon bwa mem -M -R "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tPL:DNBSEQ" -t $NT -K 10000000 $REFERENCE $FASTQ_1 $FASTQ_2 \
    | sentieon util sort --bam_compression 1 -r $REFERENCE -o sorted.bam -t $NT --sam2bam -i -
fi

# S2-deplicate removal
if [ ! -f deduped.bam ]; then
    sentieon driver -t $NT -i sorted.bam --algo LocusCollector --fun score_info score.txt
    sentieon driver -t $NT -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt --bam_compression 1 deduped.bam
fi

# S3-recalibration
if [ ! -f recal_data.table ]; then
    sentieon driver -r $REFERENCE -t $NT -i deduped.bam \
    --algo QualCal -k $DBSNP -k $KNOWN_MILLS_INDELS -k $KNOWN_1000G_INDELS recal_data.table
fi

# S4-variants call
# You can all --interval $TARGET_REGION \
if [ ! -f ${SAMPLE_NAME}.gvcf.gz.tbi ]; then
sentieon driver -r $REFERENCE -t $NT -i deduped.bam -q recal_data.table \
    --algo Haplotyper \
    -d $DBSNP \
    --emit_conf=10 \
    --call_conf=10 \
    --emit_mode gvcf \
    ${SAMPLE_NAME}.gvcf.gz
fi

echo "Sample $SAMPLE_NAME finished."

