#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -N NA12890-gvcf
#PBS -j oe
#PBS -l mem=8GB

set -e

# Java, Picard, and GATK
module load java-jdk/1.8.0_92 
module load gatk/4.1.3.0

# Path to output directory (yours will be different)
odir="/home/t.cri.biowksp13/wk8"

# Sample name
sample="NA12890"

# Path to sorted bam 
recalBam="${odir}/${sample}.dedup.recal.bam"

# Path to reference files
ref="/gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/hg38bundle/Homo_sapiens_assembly38.fasta"
dbsnp="/gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/hg38bundle/dbsnp_147.hg38.vcf.gz"

gvcf="${odir}/${sample}.g.vcf.gz"

# Make gVCF 
java -Xmx6G -jar $GATK HaplotypeCaller \
-R $ref \
-I $recalBam \
-L chr20 \
-O $gvcf \
-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
-new-qual \
-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
-ERC GVCF > ${odir}/${sample}.gvcf.logs 2>&1
