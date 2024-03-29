/*
 * pipeline input parameters - these are defaults, we can manually override them in the command line
 * params.reads = location of the reads we want to process
 * params.multiqc = location of multiqc output
 */
params.reads = "$projectDir/data/*/*_{R1_001,R2_001}.fastq.gz"
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
    publishDir params.outdir, mode:'copy' //this is where to put the output files

    input:
    path '*' // parse in the fastqc output files

    output:
    path 'multiqc_report.html' //output path

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
    MULTIQC(fastqc_ch.collect()) //use collect to process all the files as one task
}

// message to tell you if the process worked!

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}