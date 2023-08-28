#!/usr/bin/env nextflow

/*
 * pipeline input parameters - these are defaults, we can manually override them in the command line
 * params.reads = location of the reads we want to process
 */

params.reads = "$projectDir/data/*/*_{R1_001,R2_001}.fastq.gz"
params.outdir = 'results/trimmomatic'

// print the reads that will be processed and the output directory

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()


process TRIM {

    params.outdir, params.reads, mode: 'copy' //this is where to put the output files
    container 'quay.io/biocontainers/trimmomatic:0.35--6'

    input:
    tuple genomeName, file(genomeReads) from ch_in_trimmomatic

    output:
    tuple path(fq_1_paired), path(fq_2_paired) into ch_out_trimmomatic

    script:

    fq_1_paired = genomeName + '_R1.p.fastq'
    fq_1_unpaired = genomeName + '_R1.s.fastq'
    fq_2_paired = genomeName + '_R2.p.fastq'
    fq_2_unpaired = genomeName + '_R2.s.fastq'

    """
    trimmomatic \
    PE -phred33 \
    ${genomeReads[0]} \
    ${genomeReads[1]} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}


//workflow

workflow {
        Channel.fromFilePairs(params.reads, checkIfExists: true) //check that the files exist before running the workflow
               .into { ch_in_trimmomatic } //parse the file names into the variable ch_in_trimmomatic??!?
}

// message to tell you if the process worked!
workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}