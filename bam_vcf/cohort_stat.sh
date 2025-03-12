#PBS -l nodes=1:ppn=4
#PBS -l walltime=08:00:00
#PBS -N stat
#PBS -j oe
#PBS -l mem=8GB

module load java-jdk/1.8.0_92
module load gatk/4.1.3.0
module load gcc/6.2.0
module load bcftools/1.2

odir="/home/t.cri.biowksp13/wk8"

ref="/gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/hg38bundle/Homo_sapiens_assembly38.fasta"

# Run GATK Combine GVCFs to merge per-sample gVCFS
java -Xmx6G -jar $GATK CombineGVCFs \
   -R $ref \
   --variant ${odir}/NA12890.g.vcf.gz \
   --variant ${odir}/NA12889.g.vcf.gz \
   --variant ${odir}/NA12877.g.vcf.gz \
   -L chr20 \
   -O ${odir}/cohort.g.vcf.gz \
   -G StandardAnnotation -G AS_StandardAnnotation 2> ${odir}/combine_gvcf.logs

# Run GATK Joint Genotyping
java -Xmx6G -jar $GATK GenotypeGVCFs \
  -R $ref \
  -O ${odir}/OUTPUT.vcf.gz \
  -G StandardAnnotation -G AS_StandardAnnotation \
  --only-output-calls-starting-in-intervals \
  --use-new-qual-calculator \
  -V ${odir}/cohort.g.vcf.gz \
  -L chr20 \
  --merge-input-intervals


# Run bcftools stats to generate basic statistics of the cohort VCF
bcftools stats -F $ref ${odir}/OUTPUT.vcf.gz > ${odir}/stats_output

