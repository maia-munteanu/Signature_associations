#! /usr/bin/env nextflow
//vim: syntax=groovy -*- mode: groovy;-*-
// Copyright (C) 2022 IRB Barcelona

nextflow.enable.dsl=2

params.input_file = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/GermlineMuts_joint_clustered_signatures_indels_rawExp_wSums.tsv"
params.input_file2 = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/GermlineMuts_joint_clustered_signatures_svs_rawExp_wSums.tsv"
params.output_folder = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/pTweedie_Covs_RawExp_logoSum_pvar_inter"
params.metadata = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Pipeline_inputs/Hartwig_PCAWG_TCGA_MMRF_CPTAC_OVCARE_MUTes.tsv"
params.model = "pTweedie_logoSum_pvar_inter"
params.covariates = "TRUE"

workflow {
    signatures = Channel.fromPath(params.input_file)
        .first()
        .map { file -> 
            file.withReader { reader -> 
                reader.readLine() 
            }
        }
        .map { header -> 
            header.split('\t')[5..-1] 
        }
        .flatMap { it.toList() } 
    get_model(signatures)

    signatures2 = Channel.fromPath(params.input_file2)
        .first()
        .map { file -> 
            file.withReader { reader -> 
                reader.readLine() 
            }
        }
        .map { header -> 
            header.split('\t')[5..-1] 
        }
        .flatMap { it.toList() } 
    get_model(signatures2)
}

process get_model {
    tag "${signature}"
    errorStrategy 'retry'
    maxRetries 3
    memory { 15.GB * task.attempt }

    publishDir params.output_folder, mode: 'copy', pattern: "${signature}.tsv"

    input:
    val signature

    output:
    path "*.tsv" 

    shell:
    '''
    awk -F '\t' -v sig="!{signature}" '
    NR==1 {
        # Store the header line
        header = $0
        for(i=1; i<=NF; i++) {
            if ($i == sig) colnum = i
        }
        print $1, $2, $3, $4, $5, sig
    }
    NR>1 {
        print $1, $2, $3, $4, $5, $colnum
    }
    ' OFS='\t' !{params.input_file} > signature_file.tsv
    Rscript !{baseDir}/get_model.R !{signature} signature_file.tsv !{params.model} !{params.metadata} !{params.covariates}
    '''
}

