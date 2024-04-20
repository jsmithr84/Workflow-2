#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "$projectDir/data/*.fastq.gz"
params.assembled = "$projectDir/assembled"
params.qa = "$projectDir/qa"
params.mlst = "$projectDir/mlst"

process AssembleFastA {
    tag "Assembling ${reads.baseName}"

    publishDir "${params.assembled}", mode: 'copy'

    input:
    path reads

    output:
    path "*_assembled.fasta"

    script:
    """
    mkdir -p assembly_output
    spades.py -s ${reads} -o assembly_output --phred-offset 33
    mv assembly_output/contigs.fasta ${reads.baseName}_assembled.fasta
    """
}

process QualityAssessment {
    tag "QA ${assembled_reads.baseName}"

    publishDir "${params.qa}", mode: 'copy'

    input:
    path assembled_reads

    output:
    path "*_qa"

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    quast.py -o ${assembled_reads.baseName}_qa ${assembled_reads}
    """
}

process RunMLST {
    tag "MLST ${assembled_reads.baseName}"

    publishDir "${params.mlst}", mode: 'copy'

    input:
    path assembled_reads

    output:
    path "*_mlst/results.txt"

    script:
    """
    mkdir -p ${assembled_reads.baseName}_mlst
    mlst ${assembled_reads} > ${assembled_reads.baseName}_mlst/results.txt
    """
}

workflow {
    // Create a channel from the reads path
    reads_channel = Channel.fromPath(params.reads)

    // Run assembly on the reads
    AssemblyResults = AssembleFastA(reads_channel)

    // Run QA and MLST in parallel after assembly
    QualityAssessment(AssemblyResults)
    RunMLST(AssemblyResults)
}

