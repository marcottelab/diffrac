# diffrac
Scripts to analyze differential fractionation (DIF-FRAC) experiments. Compares two elut files for significant changes between conditions. 

## Scripts
### diffrac.py
    Primary script to compare experiments

## Data files:
### ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut 
    Size exclusion chromatography (SEC) fractionation experiment on mouse embryonic stem cells with no treatment (control)

### ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut
    Size exclusion chromatography (SEC) fractionation experiment on mouse embryonic stem cells with RNase A treatment

### ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut
    Annotation file for identifying known RNA binding proteins, most annotations are from Hentze et al. Nat Rev Mol Cell Biol. 2018 

## Sample Commandline:

python diffrac.py --elution_files ./data/MES_SEC_Cntl_20180601.prot_count_mFDRpsm001.elut ./data/MES_SEC_RNAse_20180601.prot_count_mFDRpsm001.elut --features diffrac diffrac_percent pearsonr diffrac_normalized mean_abundance emd zscore sliding_zscore fdr_correct sliding_fdr_correct --output_file MES_SEC_Cntl_RNAse_20180601_diffrac_features.feat --annotated_list ./data/mouse_rna_annotations_go_uniprot_hentze.txt --use_gmm &> MES_SEC_Cntl_RNAse_diffrac_feat.out
