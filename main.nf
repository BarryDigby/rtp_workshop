#!/usr/bin/env nextflow

Channel.fromFilePairs("${params.input}", checkIfExists: true)
       .set{ ch_reads }

ch_fasta = Channel.value(file(params.fasta))

process FASTQC{
    tag "${base}"
    publishDir params.outdir, mode: 'copy',
        saveAs: { params.save_qc_intermediates ? "fastqc/${it}" : null }

    when:
    params.run_qc

    input:
    tuple val(base), file(reads) from ch_reads

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

process INDEX{
    publishDir params.outdir, mode: 'copy',
        saveAs: { params.save_index ? "reference/index/${it}" : null }

    when:
    !params.kallisto_index && params.fasta

    input:
    file(fasta) from ch_fasta

    output:
    file("*.idx") into index_created

    script:
    """
    kallisto index -i ${fasta.baseName}.idx $fasta
    """
}

ch_index = params.kallisto_index ? Channel.value(file(params.kallisto_index)) : index_created

ch_index.view()
