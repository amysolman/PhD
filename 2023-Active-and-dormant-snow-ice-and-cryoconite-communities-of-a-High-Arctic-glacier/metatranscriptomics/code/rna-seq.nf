/*
 * pipeline input parameters - these are defaults, we can manually override them in the command line
 * params.reads = location of the reads we want to process
 */
params.reads = "$projectDir/data/*/*_{R1_001,R2_001}.fastq.gz"
// params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"

// print the reads that will be processed and the output directory
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()


process FASTQC {
    tag "FASTQC on $sample_id" //print which sample_id is being processed

    input:
    tuple val(sample_id), path(reads) //inputs are the sample_id and path to the reads

    output:
    path "fastqc_${sample_id}_logs" //define output path for the fastqc files

    //make a folder for the fastqc output files
    //then run fastqc and put the output files there

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true) //check that the files exist before running the workflow
        .set { read_pairs_ch } 

    fastqc_ch = FASTQC(read_pairs_ch)
    //MULTIQC(quant_ch.mix(fastqc_ch).collect())
}
