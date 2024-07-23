#! /usr/bin/env nextflow
//vim: syntax=groovy -*- mode: groovy;-*-
// Copyright (C) 2022 IRB Barcelona
// example run: nextflow run ../nextflow_pipeline/Mutation_clustering/main3.nf --serial_genome /g/strcombio/fsupek_cancer1/SV_clusters_project/Test/hg19.fa.p

nextflow.enable.dsl=2

params.input_file = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/GermlineMuts_signatures_indels.tsv"
params.output_folder = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/GLMnb_NoCovs_RawExp"
params.model = "GLMnb"
/*
workflow {
    signatures = Channel.fromPath(params.input_file)
                        .map { file -> file.text }
                        .map { text -> text.readLines().get(0) }
                        .map { header -> header.split('\t')[3..-1] }
                        .flatMap { it.toList() } 
    get_model(signatures)
}
*/

workflow {

    Channel
        .fromPath(params.input_file)
        .first() // Ensures we only process the first line
        .map { file ->
            file.withReader { reader ->
                reader.readLine() // Read only the first line
            }
        }
        .map { header ->
            header.split('\t')[3..-1] // Split the header and take from 4th column
        }
        .flatMap { it.toList() } // Flatten to individual elements
        .set { signatures } // Set to a channel for downstream use

    TestColumns(signatures)

}

process TestColumns {
        input:
        val column

        script:
        """
        echo "Processing column: $column"
        """
}

process get_model {
    tag "${signature}"
    errorStrategy 'retry'
    maxRetries 3
    memory { 15.GB * task.attempt }

    publishDir params.output_folder, mode: 'copy', pattern: '*tsv'

    input:
    val signature

    output:
    path "*.tsv" 

    script:
    """
    Rscript ${baseDir}/get_model.R ${signature} ${params.input_file} ${params.model}
    """
}

