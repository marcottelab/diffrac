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


process create_ids_file {

    input:
    file ctrlElut from ctrlElutFile
    file expElut from expElutFile
    
    output:
    file 'ids.txt' into ids_file
    stdout result2

    """
    cat $ctrlElut $expElut |cut -f1|sort|uniq > 'ids.txt'
    """
}

process intermediate_process {

    input:
    file ids_file from ids_file

    output:
    val Channel.from(ids_file.readLines()) into ids

    """
    echo "intermediate_process"
    """
}


process plotSparkline {

    publishDir "${params.results_path}"

    input:
    file ctrlElut from ctrlElutFile
    file expElut from expElutFile
    val id1 from ids
    
    output:
    stdout result3

    """
    echo "plotSparkline: $id1"
    """
}


result2.subscribe {
    println it
}

result3.subscribe {
    println it
}


