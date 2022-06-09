#!/usr/bin/env nextflow

params.diffrac = '/diffrac/'

params.complex_file = '/diffrac/data/combined_huMAP_CORUM_complexmerge100percent_everbeke_20180607.txt'

params.control_elut = ""
params.experiment_elut = ""
params.annotation_file = ""
params.split_amount = 50

params.results_path = "./results/complexes"

ctrlElutFile = file(params.control_elut)
expElutFile = file(params.experiment_elut)
annotationFile = file(params.annotation_file)
complexFile = file(params.complex_file)

params.parse_fraction_fields = "cell_type col_type condition fraction subindex date"
params.additional_options = ''


process create_ids_file {

    input:
    file complexFileIn from complexFile
    
    output:
    file 'ids_*' into ids mode flatten
    stdout result2

    """
    split -d -a 4 -l ${params.split_amount} $complexFileIn ids_
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
    bsname=$id1.baseName
    process_id="\${bsname#*_}"
    split -d -a 2 -l 1 $id1 id1s_
    for i in id1s_*;
    do
        #echo \$i
        subprocess_id="\${i#*_}"
        #echo \$subprocess_id
        #kdrew: super narly code, this calculates the index of each entry based on the process_id * split_amount + subprocess_id, the 10# tells bash these are decimals
        final_index=\$((10#\$process_id * ${params.split_amount} + 10#\$subprocess_id))
        #echo \$final_index
        #cat \$i;
        python ${params.diffrac}/evaluation/plots/plot_sparklines.py --filenames $ctrlElut $expElut --proteins `cat \$i` --output_filename complex_\${final_index}_sparkline.pdf --labels Control Experiment --annotation_file $annotation --parse_fraction_name ${params.parse_fraction_fields} ${params.additional_options}
    done;
    """
}


result2.subscribe {
    println it
}

result3.subscribe {
    println it
}


