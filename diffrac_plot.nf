#!/usr/bin/env nextflow

params.diffrac = '/diffrac/'

params.control_elut = ""
params.experiment_elut = ""
params.annotation_file = ""
params.split_amount = 250

params.results_path = "./results"

ctrlElutFile = file(params.control_elut)
expElutFile = file(params.experiment_elut)
annotationFile = file(params.annotation_file)


process create_ids_file {

    input:
    file ctrlElut from ctrlElutFile
    file expElut from expElutFile
    
    output:
    file 'ids_*' into ids mode flatten
    stdout result2

    """
    cat $ctrlElut $expElut |cut -f1|sort|uniq |grep -v '^\$' > 'ids.txt'
    split -a 4 -l ${params.split_amount} 'ids.txt' ids_
    """
}

process plotSparkline {

    publishDir "${params.results_path}", mode: 'copy'

    input:
    file ctrlElut from ctrlElutFile
    file expElut from expElutFile
    file annotation from annotationFile
    file id1 from ids
    
    output:
    file '*_sparkline.pdf' into sparkline_plots
    stdout result3

    """
    echo "plotSparkline: $id1"
    split -l 1 $id1 id1s_
    for i in id1s_*;
    do
        cat \$i;
        python ${params.diffrac}/evaluation/plots/plot_sparklines.py --filenames $ctrlElut $expElut --proteins `cat \$i` --output_filename `cat \$i`_sparkline.pdf --labels Control Experiment --annotation_file $annotation --parse_fraction_name cell_type col_type condition fraction subindex date
    done;
    """
}


result2.subscribe {
    println it
}

result3.subscribe {
    println it
}


