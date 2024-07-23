#! /usr/bin/env nextflow
//vim: syntax=groovy -*- mode: groovy;-*-
// Copyright (C) 2022 IRB Barcelona
// example run: nextflow run ../nextflow_pipeline/Mutation_clustering/main3.nf --serial_genome /g/strcombio/fsupek_cancer1/SV_clusters_project/Test/hg19.fa.p

nextflow.enable.dsl=2

params.input_file = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/test.tsv"
params.input_folder = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/Test"
params.model = "GLMnb"

workflow {
    header_ch = Channel.fromPath(params.input_file)
                        .map { file -> file.text }
                        .map { text -> text.readLines().get(0) }
                        .map { header -> header.split('\t')[3..-1] }
    OperateOnColumns(header_ch)
}

process OperateOnColumns {
    tag "${column}"

    input:
    val column

    script:
    """
    echo "Processing column: ${column}"
    """
}

