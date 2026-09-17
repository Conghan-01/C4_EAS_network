# We performed localized imputation using the Osprey software suite (https://github.com/broadinstitute/Osprey)
# Genotype imputation was performed using the TOPMed reference panel, and variants with an imputation R2 < 0.3 were excluded to ensure high-quality data.
## Pipeline for C4 Region CNV Imputation using Osprey
## Reference Panel: 1000 Genomes EAS (East Asian)
## Input: Genotype imputed VCFs (TOPMed reference, R2 < 0.3 filtered)

## "Step 1: Filtering variants based on HWE and MAF..."
# Filter fetal cohort
plink2 --vcf chr6.dose.vcf.gz \
      --hwe 1e-6 \
      --maf 0.01 \
      --export vcf bgz\
      --out chr6.fetal.pass2
# Filter adult cohort
plink2 --vcf chr6.adult.vcf.gz \
      --hwe 1e-6 \
      --maf 0.01 \
      --export vcf bgz\
      --out chr6.adult.pass2

# Step 2: Extracting MHC region (chr6:24M-34M) and standardizing chr names(Add "chr": '6' to 'chr6')
bcftools view -r 6:24000000-34000000  chr6.fetal.pass.vcf.gz -Oz -o chr6.region.fetal.vcf.gz 
bcftools view -r 6:24000000-34000000  chr6.adult.pass.vcf.gz -Oz -o chr6.region.adult.vcf.gz 

bcftools annotate --rename-chrs /gpfs/hpc/home/chenchao/share_group_folder/data/CHFB/process/processing/resource/chrname.txt  chr6.region.fetal.vcf.gz -Oz -o chr6.region.fetalchr.vcf.gz 
bcftools annotate --rename-chrs /gpfs/hpc/home/chenchao/share_group_folder/data/CHFB/process/processing/resource/chrname.txt  chr6.region.adult.vcf.gz -Oz -o chr6.region.adultchr.vcf.gz
# Index the extracted VCFs for merging

bcftools index chr6.region.fetalchr.vcf.gz
bcftools index chr6.region.adultchr.vcf.gz

# Step 3: Merge target cohorts with the 1KG EAS Reference Panel
bcftools1.0 merge chr6.region.fetalchr.vcf.gz  /gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/data/1KG_ref/C4_panel_1000G_ext_reference_yuchen_20251019/1KG514EAS.C4.chr6.vcf.gz -Oz -o merged_fetal_raw.vcf.gz
bcftools1.0 merge   chr6.region.adult.vcf.gz  /gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/data/1KG_ref/C4_panel_1000G_ext_reference_yuchen_20251019/1KG514EAS.C4.chr6.vcf.gz -Oz -o merged_adult_raw.vcf.gz

# Remove variants with missing genotypes ('./.') which can interfere with Osprey
# Note: Using bcftools is safer than zcat | grep to preserve headers
bcftools view -e 'GT="mis"' merged_fetal_raw.vcf.gz -Oz -o merged_fetal.vcf.gz
bcftools view -e 'GT="mis"' merged_adult_raw.vcf.gz -Oz -o merged_adult.vcf.gz

#zcat merged_adult.vcf.gz|grep -v './.'|bgzip > merged_adult.vcf.gz2
#zcat merged_fetal.vcf.gz|grep -v './.'|bgzip > merged_fetal.vcf.gz2
#mv  merged_fetal.vcf.gz2  merged_fetal.vcf.gz
#mv  merged_adult.vcf.gz2  merged_adult.vcf.gz

bcftools index merged_adult.vcf.gz
bcftools index merged_fetal.vcf.gz


conda activate osprey_env

# Step 4: VCF Formatting for Osprey (Biallelic SNPs only, strict header)
# Fix missing headers if necessary (missing_header.txt is formatted correctly)
bcftools annotate -h missing_header.txt -Oz -o merged_adult_hfix.vcf.gz merged_adult.vcf.gz
bcftools annotate -h missing_header.txt -Oz -o merged_fetal_hfix.vcf.gz merged_fetal.vcf.gz

bcftools view -m2 -M2 -v snps -Oz -o merged_adult_snps.vcf.gz merged_adult_hfix.vcf.gz
bcftools view -m2 -M2 -v snps -Oz -o merged_fetal_snps.vcf.gz merged_fetal_hfix.vcf.gz

# Strip all FORMAT fields except GT to minimize file size and prevent parsing errors in Osprey

bcftools annotate -x ^FORMAT/GT -Oz -o merged_adult.simple.vcf.gz merged_adult_snps.vcf.gz
bcftools annotate -x ^FORMAT/GT -Oz -o merged_fetal.simple.vcf.gz merged_fetal_snps.vcf.gz

#Or

#bcftools annotate \
#    -x INFO \
#    -Oz -o merged_fetal.vcf.simple.gz \
#    merged_fetal.vcf.gz

# Step 5: Calculate Identity-By-State (IBS) Matrix using ospreyIBS

OSPREY_DIR="/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/Osprey/bin"
REF_SAMPLES="/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/data/1KG_ref/C4_panel_1000G_ext_reference_yuchen_20251019/C4_panel_1000G_514EAS.samples.list"
GENETIC_MAP="/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/data/1KG_ref/phasing_region.shapeit.gmap"
IBS_INTERVAL="chr6:31980001-32046200"

# Calculate IBS for fetal cohort
${OSPREY_DIR}/ospreyIBS \
    --inputFile merged_fetal.simple.vcf.gz \
    --ibsMatrixFile c4.fetal.ibs_matrix.gz \
    --ibsInterval ${IBS_INTERVAL} \
    --geneticMapFile ${GENETIC_MAP} \
    --refSampleFile ${REF_SAMPLES}

# Calculate IBS for adult cohort
${OSPREY_DIR}/ospreyIBS \
    --inputFile merged_adult.simple.vcf.gz \
    --ibsMatrixFile c4.adult.ibs_matrix.gz \
    --ibsInterval ${IBS_INTERVAL} \
    --geneticMapFile ${GENETIC_MAP} \
    --refSampleFile ${REF_SAMPLES}

# Step 6: Impute C4 CNVs using osprey
CNV_REF="/gpfs/hpc/home/chenchao/hanc/project/C4_CNB2025/data/1KG_ref/C4_panel_1000G_ext_reference_yuchen_20251019/C4_panel_1000G_ext_eas_hg38.C4AB_cnvs.vcf.gz"

# Impute adult cohort
${OSPREY_DIR}/osprey \
    --inputFile ${CNV_REF} \
    --outputFile adult_chb_imputed.cnv.vcf.gz \
    --ibsMatrixFile c4.adult.ibs_matrix.gz \
    --threads 3 

# Impute fetal cohort
${OSPREY_DIR}/osprey \
    --inputFile ${CNV_REF} \
    --outputFile fetal_chb_imputed.cnv.vcf.gz \
    --ibsMatrixFile c4.fetal.ibs_matrix.gz \
    --threads 3


# Step 7: Extract target samples from the imputed VCFs (Remove ref samples)

# Generate sample lists from the original pass VCFs
bcftools query -l chr6.adult.pass.vcf.gz > adult.sample.list
bcftools query -l chr6.fetal.pass.vcf.gz > fetal.sample.list 

# Extract only the target samples, excluding the 1KG reference samples
bcftools view -S adult.sample.list adult_chb_imputed.cnv.vcf.gz -Oz -o adult_only_imputed.cnv.vcf.gz 
bcftools view -S fetal.sample.list fetal_chb_imputed.cnv.vcf.gz -Oz -o fetal_only_imputed.cnv.vcf.gz

