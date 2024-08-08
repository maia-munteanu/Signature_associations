#! /usr/bin/env nextflow
//vim: syntax=groovy -*- mode: groovy;-*-
// Copyright (C) 2022 IRB Barcelona

nextflow.enable.dsl=2

params.sig_file = "/g/strcombio/fsupek_cancer3/SV_clusters_project/CNA_genes/Indel_signatures_raw.tsv"
params.cna_file = "/g/strcombio/fsupek_cancer3/SV_clusters_project/CNA_genes/CNA_by_gene_all_samples_restrained.tsv"
params.metadata = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Pipeline_inputs/Hartwig_PCAWG_TCGA_MMRF_CPTAC_OVCARE_MUTes.tsv"
params.output_folder = "/g/strcombio/fsupek_cancer3/SV_clusters_project/CNA_genes/GLMnb_Covs_RawExp"
params.model = "GLMnb"
params.covariates = "TRUE"

workflow {
    signatures = Channel.fromPath(params.sig_file)
        .first()
        .map { file -> 
            file.withReader { reader -> 
                reader.readLine() 
            }
        }
        .map { header -> 
            header.split('\t')[1..-1] 
        }
        .flatMap { it.toList() } 
    get_model(signatures)
}

process get_model {
    tag "${signature}"
    errorStrategy 'retry'
    maxRetries 3
    memory { 10.GB * task.attempt }

    publishDir params.output_folder, mode: 'copy', pattern: "${signature}.tsv"

    input:
    val signature

    output:
    path "*.tsv" 

    shell:
    '''
    Rscript !{baseDir}/get_model_CN.R !{signature} !{params.sig_file} !{params.cna_file} !{params.metadata} !{params.model} !{params.covariates} !{params.output_folder}
    '''
}
