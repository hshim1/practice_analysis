
ch_fastq = Channel.fromFilePairs(params.fastq_pattern)
ch_fastq_star = Channel.fromFilePairs(params.fastq_pattern)
ch_star_index = Channel.fromPath(params.star_index)

ch_ref_fa = Channel.fromPath(params.ref_fa)
ch_ref_rib = Channel.fromPath(params.ribosomal_intervals)
ch_ref_flat = Channel.fromPath(params.ref_flat)
ch_pair_orient = Channel.value(params.pair_orientation)
ch_strand_specificity = Channel.fromPath(params.strand_specificity)

process fastqc {
    tag "fastqc ${sample}"
    cpus 4
    memory '4 GB'
    time { 4.hour * task.attempt }
    module 'java-jdk/1.10.0_1'
    module 'fastqc/0.11.7'

    input:
    tuple val(sample), file(fq) from ch_fastq


    output:
    file "*_fastqc.zip" into ch_fastqc


    script:
    """
        fastqc $fq
    """
}


process star_align {
    tag "STAR Align ${sample}"
    cpus 4
    memory '32 GB'
    time { 4.hour * task.attempt }
    module 'gcc/6.2.0'
    module 'STAR/2.6.1d'
  
    input:
    tuple val(sample), file(fq_files) from ch_fastq_star
    file index_path from ch_star_index.collect()

    output:
    tuple val(sample), file("*Aligned.out.bam") into ch_star_aln_out
    file "*Log.final.out" into ch_star_aln_stats

    script:
	template 'star_aln.sh'
}




process samtools_sort_index {
    tag "samtools sort index ${sample}"
    cpus 4
    memory '4 GB'
    time { 4.hour * task.attempt }
    module 'gcc/6.2.0'
    module 'samtools/1.10'

    input:
    tuple val(sample), file(bam) from ch_star_aln_out

    output:
    tuple val(sample), file ("*.sorted.bam") into ch_sorted_bam
    file "*.sorted.bam.bai" into ch_sorted_bai

    script:
   	 template 'samtools_sort_index.sh'
}

process post_alignment_metrics {
    tag "post_aln_metrics ${sample}"
    publishDir "metrics/"
    cpus 4
    memory '4 GB'
    time { 4.hour * task.attempt }
    module 'java-jdk/1.8.0_92'
    module 'gcc/6.2.0'
    module 'R/3.6.1'
    module 'picard/2.18.29'

    input:
	tuple val(sample), file(bam) from ch_sorted_bam
	val pair_orientation from ch_pair_orient
	file ref_fa from ch_ref_fa.collect()
	file ref_rib from ch_ref_rib.collect()
	file ref_flat from ch_ref_flat.collect()
    
    output:
       
	file "${sample}*" into ch_pval_metrics
    
    script:
    template 'post_alignment_metrics.sh'

}

process multiqc {
    tag "multiqc"
    publishDir "report/"
    cpus 4
    memory '4 GB'
    time { 4.hour * task.attempt }
    module 'gcc/6.2.0'
    module 'python/3.6.0'

    input:
    file "*.metrics" from ch_pval_metrics.collect()
    file "*Log.final.out" from ch_star_aln_stats.collect()
    file "*_fastqc.zip" from ch_fastqc.collect()
    
    output:
    file "multiqc_report.html" into multiqc_report

    script:
    """
        PS1='-'
        . /home/t.cri.biowksp13/wk5/venv/bin/activate
        export LC_ALL=en_US.utf8
        export LANG=en_US.utf8
        multiqc .
     """
}
