#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -N NA12877-preprocess
#PBS -j oe
#PBS -l mem=8GB

set -e

# Java, Picard, and GATK
module load java-jdk/1.8.0_92 
module load picard/2.18.29 
module load gatk/4.1.3.0

# Path to output directory (yours will be different)
odir="/home/t.cri.biowksp13/wk8"

# Sample name
sample="NA12877"

# Path to sorted bam 
srtBam="${odir}/${sample}.align.raw.sort.bam"

# Path to reference files
ref="/gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/hg38bundle/Homo_sapiens_assembly38.fasta"
dbsnp="/gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/hg38bundle/dbsnp_147.hg38.vcf.gz"

# Mark duplicates with picard
dupBam="${odir}/${sample}.markdup.bam"
dupMetrics="${odir}/${sample}.markdup_metrics.txt"

java -Xmx6G -jar $PICARD MarkDuplicates TMP_DIR=${odir}/tmp USE_JDK_DEFLATER=true USE_JDK_INFLATER=true I=${srtBam} O=${dupBam} M=${dupMetrics} \
REMOVE_DUPLICATES=false ASSUME_SORTED=true > ${odir}/${sample}.markdup.logs 2>&1

# Index bam
java -Xmx6G -jar $PICARD BuildBamIndex TMP_DIR=${odir}/tmp USE_JDK_DEFLATER=true USE_JDK_INFLATER=true I=${dupBam} > ${odir}/${sample}.markdup.index.logs 2>&1

# Next, we run base recalibration
ogrp="${odir}/${sample}_bqsr.grp"
java -Xmx6G -jar $GATK BaseRecalibrator --use-jdk-inflater --use-jdk-deflater --input ${dupBam} --known-sites ${dbsnp} --reference $ref \
--output $ogrp > ${odir}/${sample}.baserecal.logs 2>&1

# Finally, we apply the recalibration
orecal="${odir}/${sample}.dedup.recal.bam"
java -Xmx6G -jar $GATK ApplyBQSR --use-jdk-inflater --use-jdk-deflater --input ${dupBam} --bqsr-recal-file $ogrp \
--output $orecal > ${odir}/${sample}.apply_recal.logs 2>&1
