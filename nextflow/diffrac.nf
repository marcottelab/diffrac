#!/usr/bin/env nextflow

params.diffrac = '/diffrac/'

params.control_elut = ""
params.experiment_elut = ""
params.annotation_file = ""

params.results_path = "./results"
params.feat_file = 'diffrac.feat'

params.features = "diffrac diffrac_percent pearsonr diffrac_normalized mean_abundance emd zscore sliding_zscore fdr_correct sliding_fdr_correct"

ctrlElutFile = file(params.control_elut)
expElutFile = file(params.experiment_elut)


process runDiffrac {

    publishDir "${params.results_path}", mode: 'copy'

    input:
    file ctrlElut from ctrlElutFile
    file expElut from expElutFile

    output:
    file "${params.feat_file}" into feat_file
    stdout result

    """
    echo 'running Diffrac';

    python ${params.diffrac}/diffrac.py --elution_files $ctrlElut $expElut  --features ${params.features} --output_file ${params.feat_file} --annotated_list ${params.annotation_file}  --use_gmm 
    """
}

result.subscribe {
    println it
}

