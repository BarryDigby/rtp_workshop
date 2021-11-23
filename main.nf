#!/usr/bin/env nextflow

Channel.fromFilePairs("${params.input}", checkIfExists: true)
       .into{ ch_qc_reads; ch_alignment_reads }

ch_fasta = Channel.value(file(params.fasta))
ch_gtf = Channel.value(file(params.gtf))

process FASTQC{
    tag "${base}"
    publishDir params.outdir, mode: 'copy',
        saveAs: { params.save_qc_intermediates ? "fastqc/${it}" : null }

    when:
    params.run_qc

    input:
    tuple val(base), file(reads) from ch_qc_reads

    output:
    tuple val(base), file("*.{html,zip}") into ch_multiqc

    script:
    """
    fastqc -q $reads
    """
}

process MULTIQC{
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    when:
    params.run_qc

    input:
    file(htmls) from ch_multiqc.collect()

    output:
    file("multiqc_report.html") into multiqc_out

    script:
    """
    multiqc .
    """
}

process TX{
    publishDir params.outdir, mode: 'copy',
        saveAs: { params.save_transcriptome ? "reference/transcriptome/${it}" : null }

    when:
    !params.transcriptome && params.fasta

    input:
    file(fasta) from ch_fasta
    file(gtf) from ch_gtf

    output:
    file("${fasta.baseName}.tx.fa") into transcriptome_created

    script:
    """
    gffread -F -w "${fasta.baseName}.tx.fa" -g $fasta $gtf
    """
}

ch_transcriptome = params.transcriptome ? Channel.value(file(params.transcriptome)) : transcriptome_created

process INDEX{
    publishDir params.outdir, mode: 'copy',
        saveAs: { params.save_index ? "reference/index/${it}" : null }

    when:
    !params.kallisto_index

    input:
    file(tx) from ch_transcriptome

    output:
    file("*.idx") into index_created

    script:
    """
    kallisto index -i ${tx.simpleName}.idx $tx
    """
}

ch_index = params.kallisto_index ? Channel.value(file(params.kallisto_index)) : index_created

