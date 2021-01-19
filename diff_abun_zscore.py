import argparse
import csv
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import scipy.stats as st
import math

parser = argparse.ArgumentParser(description='Calculate Z-score for DIFFRAC experiments')
parser.add_argument('replicate_count', metavar = 'N', type=int,
                    help='number of replicates')
parser.add_argument('collapse_method', metavar = 'm', type=str,
                    help='method for collapsing PSM/fraction z-scores to protein z-scores. Either stouffer or max.')
parser.add_argument('control_elut_1', metavar = 'c', type=str,
                    help='control elut file path')
parser.add_argument('treat_elut_1', metavar = 't', type=str,
                    help='treatment elut file path')
parser.add_argument('-c2','--control_elut_2', type=str,
                    help='control elut file path for 2nd replicate',
                    )
parser.add_argument('-t2','--treat_elut_2', type=str,
                    help='treatment elut file path for 2nd replicate',
                    )
parser.add_argument('-c3','--control_elut_3', type=str,
                    help='control elut file path for 3rd replicate',
                    )
parser.add_argument('-t3','--treat_elut_3', type=str,
                    help='treatment elut file path for 3rd replicate',
                    )
args = parser.parse_args()

#make new dataframe with PSM/fraction format for all replicates
ctl_1 = pd.read_csv(args.control_elut_1, sep='\t')
treat_1 = pd.read_csv(args.treat_elut_1, sep='\t')
sample_dict = {'ctl_1_psms':ctl_1, 'treat_1_psms':treat_1}
sample_list  =['ctl_1_psms','treat_1_psms']
ctl_1_df = pd.DataFrame()
treat_1_df = pd.DataFrame()
df_list  =[ctl_1_df,treat_1_df]

if args.replicate_count >= 2:
    ctl_2 = pd.read_csv(args.control_elut_2, sep='\t')
    treat_2 = pd.read_csv(args.treat_elut_2, sep='\t')
    sample_dict['ctl_2_psms'] = ctl_2
    sample_dict['treat_2_psms']= treat_2
    sample_list.append('ctl_2_psms')
    sample_list.append('treat_2_psms')
    ctl_2_df = pd.DataFrame()
    treat_2_df = pd.DataFrame()
    df_list.append(ctl_2_df)
    df_list.append(treat_2_df)

if args.replicate_count == 3:
    ctl_3 = pd.read_csv(args.control_elut_3, sep='\t')
    treat_3 = pd.read_csv(args.treat_elut_3, sep='\t')
    sample_dict['ctl_3_psms'] = ctl_3
    sample_dict['treat_3_psms'] = treat_3
    sample_list.append('ctl_3_psms')
    sample_list.append('treat_3_psms')
    ctl_3_df = pd.DataFrame()
    treat_3_df = pd.DataFrame()
    df_list.append(ctl_3_df)
    df_list.append(treat_3_df)

#next part converts the elut dataframes into psm/fraction dataframes. Could be cleaner but I'm lazy

i=0
for item in sample_list:
    temp_list = []
    try:
        new_df = sample_dict[item].drop(columns='TotalCount')
    except:
        new_df = sample_dict[item]
    for index, row in new_df.iterrows():
        prot=row[0]
        k=0
        for entry in row[1:len(row)]:
            k += 1
            temp_list.append([str(prot)+'$'+str(k), entry])
    df_list[i] = pd.DataFrame(columns = ['protein', str(item)],data = temp_list)
    i+=1

#merge all of the dataframes for each replicate keeping ONLY shared proteins.
'''update 12/14/2020
Fixed the script so that within replicates there should be an union and between replicates an intersection.
This is so that the results better match the DIFFRAC score results.'''

merged_df = df_list[0].merge(df_list[1], on = "protein", how='outer').fillna(0)
if args.replicate_count >= 2:
    #tmp_merged_df = df_list[2].merge(df_list[3], on = "protein", how='outer').fillna(0)
    #merged_df = merged_df.merge(tmp_merged_df, on = "protein", how='inner').fillna(0)
    merged_df = merged_df.merge(df_list[2], on="protein", how='inner').fillna(0)
    merged_df = merged_df.merge(df_list[3], on="protein", how='inner').fillna(0)
if args.replicate_count == 3:
    #tmp_merged_df = df_list[4].merge(df_list[5], on="protein", how='outer').fillna(0)
    #merged_df = merged_df.merge(tmp_merged_df, on="protein", how='inner').fillna(0)
    merged_df = merged_df.merge(df_list[4], on="protein", how='inner').fillna(0)
    merged_df = merged_df.merge(df_list[5], on="protein", how='inner').fillna(0)

merged_df = merged_df.set_index('protein')

#drop any rows where the fractions are all 0s
modified_df = merged_df.loc[~(merged_df==0).all(axis=1)]


#the rest of this script will add a z_score column to the dataframe from everything above
control1_total = np.sum(modified_df[sample_list[0]])
treat1_total = np.sum(modified_df[sample_list[1]])
if args.replicate_count >= 2:
    control2_total = np.sum(modified_df[sample_list[2]])
    treat2_total = np.sum(modified_df[sample_list[3]])
if args.replicate_count == 3:
    control3_total = np.sum(modified_df[sample_list[4]])
    treat3_total = np.sum(modified_df[sample_list[5]])

#estimates the mean PSMs between the control and treatment samples. Repeated for replicates.
dif1 = modified_df.loc[:, sample_list[0]:sample_list[1]]
modified_df['mean_log2_PSMs_1'] = np.log2(dif1.mean(axis=1)+1)

#estimates the fold change between the phosphatase and the control samples. Repeated for replicates.
modified_df['log2_FC_1'] = np.log2(((modified_df[sample_list[1]] +1.0) -
                          (modified_df[sample_list[0]]+1.0))/(modified_df[sample_list[0]]+1.0) +1)

#z-score estimation. See methods in paper for details on estimation.
modified_df['Fc1'] = (modified_df[sample_list[0]] / control1_total)
modified_df['Ft1'] = (modified_df[sample_list[1]] / treat1_total)
modified_df['F1'] = (modified_df[sample_list[1]] +modified_df[sample_list[0]]+1) / (control1_total+treat1_total+1)

modified_df['Z_score_1'] = (modified_df["Ft1"] - modified_df["Fc1"]) / np.sqrt(
(modified_df['F1'] * (1 - modified_df['F1']) / treat1_total) + (modified_df['F1'] *
                                                                (1 - modified_df['F1']) / control1_total))

#if there are replicates z-scores will be estimated for each one. Replicate PSM/fraction Z-scores combined
#using Stouffer's Z-score method.
if args.replicate_count >= 2:
    dif2 = modified_df.loc[:, sample_list[2]:sample_list[3]]
    modified_df['mean_log2_PSMs_2'] = np.log2(dif2.mean(axis=1) + 1)
    modified_df['log2_FC_2'] = np.log2(((modified_df[sample_list[3]] + 1.0) -
                                (modified_df[sample_list[2]] + 1.0)) / (modified_df[sample_list[2]] + 1.0) + 1)
    modified_df['Fc2'] = (modified_df[sample_list[2]] / control2_total)
    modified_df['Ft2'] = (modified_df[sample_list[3]] / treat2_total)
    modified_df['F2'] = (modified_df[sample_list[3]] + modified_df[sample_list[2]] + 1) / \
                        (control2_total + treat2_total + 1)

    modified_df['Z_score_2'] = (modified_df["Ft2"] - modified_df["Fc2"]) / np.sqrt(
        (modified_df['F2'] * (1 - modified_df['F2']) / treat2_total) + (modified_df['F2'] *
                                                                        (1 - modified_df['F2']) / control2_total))

if args.replicate_count == 3:
    dif3 = modified_df.loc[:, sample_list[4]:sample_list[5]]
    modified_df['mean_log2_PSMs_3'] = np.log2(dif3.mean(axis=1) + 1)
    modified_df['log2_FC_3'] = np.log2(((modified_df[sample_list[5]] + 1.0) -
                                (modified_df[sample_list[4]] + 1.0)) / (modified_df[sample_list[4]] + 1.0) + 1)
    modified_df['Fc3'] = (modified_df[sample_list[4]] / control3_total)
    modified_df['Ft3'] = (modified_df[sample_list[5]] / treat3_total)
    modified_df['F3'] = (modified_df[sample_list[5]] + modified_df[sample_list[4]] + 1) / \
                        (control3_total + treat3_total + 1)

    modified_df['Z_score_3'] = (modified_df["Ft3"] - modified_df["Fc3"]) / np.sqrt(
        (modified_df['F3'] * (1 - modified_df['F3']) / treat3_total) + (modified_df['F3'] *
                                                                        (1 - modified_df['F3']) / control3_total))

#collapse replicate Z-scores using Stouffer's Z-score method
if args.replicate_count == 2:
    modified_df['Z_score_S'] = np.abs(modified_df['Z_score_1'] + modified_df['Z_score_2']) / np.sqrt(2)
elif args.replicate_count == 3:
    modified_df['Z_score_S'] = np.abs(modified_df['Z_score_1'] +
                                      modified_df['Z_score_2']+modified_df['Z_score_3']) / np.sqrt(3)

#outfile for psm/fraction Z-scores. Messy but can be informative.
psm_frac_df = modified_df.drop(columns=['Fc1', 'Ft1', 'F1'])
if args.replicate_count == 2:
    psm_frac_df = modified_df.drop(columns=['Fc1', 'Ft1', 'F1', 'Fc2', 'Ft2', 'F2'])
elif args.replicate_count == 3:
    psm_frac_df = modified_df.drop(columns=['Fc1','Ft1','F1','Fc2','Ft2','F2','Fc3','Ft3','F3'])

#write psm/fraction z-score outfile (for those interested in the gritty details).
psm_frac_out = './psm_frac_zscore_out.tab'
psm_frac_df.to_csv(psm_frac_out, sep='\t')

#collapse psm/fraction Z-scores to protein Z-scores using either the max score or Stouffer's Z-score method
protein_list = []
#if list(modified_df.index.values)[0].count('$') >=2:
#    print('You need to remove $ from protein names. Why do you have $s in your protein names?')
[protein_list.append(x[0:x.find('$')]) for x in list(modified_df.index.values)]
set_protein_list = set(protein_list)
final_list = []

#stouffer method for collapse
if args.collapse_method == 'stouffer':
    for prot in set_protein_list:
        temp_list_s = []
        positions = [i for i, x in enumerate(protein_list) if x == prot]
        if args.replicate_count >= 2:
            [temp_list_s.append(modified_df['Z_score_S'][pos]) for pos in positions]
        else:
            [temp_list_s.append(modified_df['Z_score_1'][pos]) for pos in positions]
        s_zscore_s = sum(temp_list_s) / np.sqrt(len(temp_list_s))
        final_list.append([prot, s_zscore_s])

#max method for collapse
elif args.collapse_method == 'max':
    for prot in set_protein_list:
        temp_list_m = []
        positions = [i for i, x in enumerate(protein_list) if x == prot]
        if args.replicate_count >= 2:
            [temp_list_m.append(modified_df['Z_score_S'][pos]) for pos in positions]
        else:
            [temp_list_m.append(modified_df['Z_score_1'][pos]) for pos in positions]
        s_zscore_m = np.max(temp_list_m)
        final_list.append([prot, s_zscore_m])

#write protein z-score outfile
collapsed_df = pd.DataFrame(columns = ['protein', 'zscore'],
                            data = final_list)
collapsed_df = collapsed_df.set_index('protein')
outfile = './protein_zscore_out.tab'
collapsed_df.to_csv(outfile, sep='\t')
