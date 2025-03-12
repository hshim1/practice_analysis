#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00
#PBS -N NA12889-align 
#PBS -j oe
#PBS -l mem=8GB

set -e
set -o pipefail

# BWA and samtools
module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.10

# Path to R1 and R2 for the NA12889 sample
fq1="/gpfs/data/mscbmi/hernandez-germline-resources/fastqs/NA12889_1.fq.gz"
fq2="/gpfs/data/mscbmi/hernandez-germline-resources/fastqs/NA12889_2.fq.gz"

# Path to output directory (yours will be differen)
odir="/home/t.cri.biowksp13/wk8"

# Path to reference files for BWA
ref="/gpfs/data/referenceFiles/Homo_sapiens/UCSC/hg38/hg38bundle/Homo_sapiens_assembly38.fasta"

# Align and convert to bam. It is *absolutely* imperative that you make the -R readgroup options correctly.
# We will set ID (identifier) SM (sample) LB (library) and PU (platform unit) to the sample ID for this exampl. The SM is used in downstream
# variant calling to put the sample name in the VCF and the LB is used in the mark duplicates step.
# Ideally you will have access to the necessary metadata to appropriately set these values. 
sample="NA12889"

bwa mem -t 4 -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPU:${sample}\tPL:ILLUMINA" \
$ref $fq1 $fq2 2> ${odir}/align.logs | samtools view -Sb -T $ref -@ 4 - > ${odir}/${sample}.align.raw.bam

# Next, we will sort the bam
samtools sort -@ 4 -T ${odir}/${sample}.align.raw.srt_tmp ${odir}/${sample}.align.raw.bam > ${odir}/${sample}.align.raw.sort.bam

# Finally we index the bam
samtools index -@ 4 -b ${odir}/${sample}.align.raw.sort.bam
