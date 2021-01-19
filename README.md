# diffrac
Scripts to analyze differential fractionation (DIF-FRAC) experiments. Compares two elut files for significant changes between conditions. 

## Scripts
### diffrac.py
    Primary script to compare experiments
### evaluation/plots/plot_sparklines.py
    Plotting script to visualize changes in elutions across experiments

## Data files:
### ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut 
    Size exclusion chromatography (SEC) fractionation experiment on mouse embryonic stem cells with no treatment (control)

### ./data/MES_SEC_RNAse_20180601.prot_count_mFDRpsm001.elut
    Size exclusion chromatography (SEC) fractionation experiment on mouse embryonic stem cells with RNase A treatment

### ./data/mouse_rna_annotations_go_uniprot_hentze.txt
    Annotation file for identifying known RNA binding proteins, most annotations are from Hentze et al. Nat Rev Mol Cell Biol. 2018 

## Sample Commandline:

```
python diffrac.py --elution_files ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut ./data/MES_SEC_RNAse_20180601.prot_count_mFDRpsm001.elut --features diffrac diffrac_percent pearsonr diffrac_normalized mean_abundance emd zscore sliding_zscore fdr_correct sliding_fdr_correct --output_file MES_SEC_Cntl_RNAse_20180601_diffrac_features.feat --annotated_list ./data/mouse_rna_annotations_go_uniprot_hentze.txt --use_gmm &> MES_SEC_Cntl_RNAse_diffrac_feat.out
```

```
python evaluation/plots/plot_sparklines.py --filenames ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut ./data/MES_SEC_RNAse_20180601.prot_count_mFDRpsm001.elut --proteins P35601 Q9WUK4 Q8R323 Q99J62 Q9D0F6 --output_filename ./RFC_sparklines.pdf --labels Control RNAseA --parse_fraction_name cell_type col_type condition fraction subindex date --annotation_file ./data/uniprot-proteome_mouse_annotations_20180331.csv --annotation_file_sep , --id_column genename
```

#Differential abundance z-score
### diff_abun_zscore.py
    Additional script for comparing experiments.
    
## Commandline format:
```
python zscore_formatter.py [# of replicates] [method for z-score collapse] [control elut file] [treatment elut file]
```

## Sample Commandline:
```
python zscore_formatter.py 1 stouffer ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut ./data/MES_SEC_RNAse_20180601.prot_count_mFDRpsm001.elut
```

## Acknowledgements

Thanks to Ben Liebeskind (@bliebeskind) for elution file reader code. 
