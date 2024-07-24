#! /usr/bin/env nextflow
//vim: syntax=groovy -*- mode: groovy;-*-
// Copyright (C) 2022 IRB Barcelona
// example run: nextflow run ../nextflow_pipeline/Mutation_clustering/main3.nf --serial_genome /g/strcombio/fsupek_cancer1/SV_clusters_project/Test/hg19.fa.p

nextflow.enable.dsl=2

params.input_file = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/GermlineMuts_signatures_indels.tsv"
params.output_folder = "/g/strcombio/fsupek_cancer3/SV_clusters_project/Germline/GLMnb_NoCovs_RawExp"
params.model = "GLMnb"

workflow {
    signatures = Channel.fromPath(params.input_file)
        .first()
        .map { file -> 
            file.withReader { reader -> 
                reader.readLine() 
            }
        }
        .map { header -> 
            header.split('\t')[3..-1] 
        }
        .flatMap { it.toList() } 
    get_model(signatures)
}

process test {
        tag "${signature}"
        input:
        val signature

        script:
        """
        echo "Processing column: ${signature}"
        """
}

process get_model {
    tag "${signature}"
    errorStrategy 'retry'
    maxRetries 3
    memory { 5.GB * task.attempt }

    publishDir params.output_folder, mode: 'copy', pattern: "${name}.tsv"

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
        print $1, $2, $3, sig
    }
    NR>1 {
        print $1, $2, $3, $colnum
    }
    ' OFS='\t' !{params.input_file} > signature_file.tsv
    Rscript !{baseDir}/get_model.R !{signature} signature_file.tsv !{params.model}
    '''
}

