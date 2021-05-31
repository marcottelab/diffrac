
#kdrew: partial original code from ipython notebook: HEK293T_SEC_cntl_RNaseA
#kdrew: earlier version in R from Claire McWhite

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import re

import argparse
import pandas as pd
import itertools as it

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
#rcParams.update({'grid': False})
mpl.rc('pdf', fonttype=42)
import seaborn as sns
sns.set_style("white")

def main():

    parser = argparse.ArgumentParser(description="Plot Sparklines")
    parser.add_argument("--filenames", action="store", dest="filenames", nargs='+', required=True, 
                                    help="Input elution profile")
    parser.add_argument("--output_filename", action="store", dest="output_filename", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--proteins", action="store", dest="proteins", nargs='+', required=True, 
                                    help="Uniprot ACCs of proteins")
    parser.add_argument("--parse_fraction_name", action="store", dest="parse_fraction_name", nargs='+', required=False, default=None,
                                    help="Fields to separate fraction name into, example: cell_type condition col_type fraction date")
    parser.add_argument("--annotation_file", action="store", dest="annotation_file", required=False, default=None,
                                    help="Uniprot mapping file")
    parser.add_argument("--id_column", action="store", dest="id_column", required=False, default='Gene names  (primary )',
                                    help="Column to use for protein identifier in annotation_file, default = 'Gene names  (primary )'")
    parser.add_argument("--colors", action="store", dest="colors", nargs='+', required=False, default=["black","red","green","orange"],
                                    help="Colors to use for plotting, default=(black red green orange)")
    parser.add_argument("--labels", action="store", dest="labels", nargs='+', required=False, default=None,
                                    help="Labels for each input filename")
    parser.add_argument("--fraction_name_sep", action="store", dest="fraction_name_sep", required=False, default='_',
                                    help="Separator to separate fraction name, default = _")
    parser.add_argument("--annotation_file_sep", action="store", dest="annotation_file_sep", required=False, default='\t',
                                    help="Separator to separate annotation file, default = '\t'")
    parser.add_argument("--annotation_file_index", action="store", type=int, dest="annotation_file_index", required=False, default=0,
                                    help="Column id of annotations index (ex. uniprot ACC), default = 0")
    parser.add_argument("--old_elut_format", action="store_true", dest="old_elut_format", required=False, default=False,
                                    help="Input elut files are of the old depricated format, default=False")
    parser.add_argument("--standard_fractions", action="store", dest="standard_fractions", nargs='+', required=False, default=[],
                                    help="Fraction index of standards run on column, example = '18 30 34 37 39 40 45', default = ''")
    parser.add_argument("--standard_mass", action="store", dest="standard_mass", nargs='+', required=False, default=[],
                                    help="Sizes of standards run on column, example = '2000 669 443 200 150 66 29', default = ''")
    parser.add_argument("--plot_y_axis", action="store_true", dest="plot_y_axis", required=False, default=False,
                                    help="Show y axis, default = False")
    parser.add_argument("--window_size", action="store", type=int, dest="window_size", required=False, default=1,
                                    help="Window size if averaging across fractions, default = 1")

    args = parser.parse_args()

    file_dict = dict()
    linspace_max=None
    for i, filename in enumerate(args.filenames):

        df = pd.read_csv(filename, sep='\t', index_col=0)
        if args.old_elut_format:
            df = df[df.columns[1:]]

        #kdrew: preprocess dataframe
        if args.window_size > 1:
            #kdrew: take windowed average across fractions
            df = df.rolling(args.window_size, axis=1, center=True, win_type='gaussian').mean(std=3)


        if args.labels != None and len(args.labels) == len(args.filenames):
            file_dict[args.labels[i]] = df
        elif args.parse_fraction_name != None:
            #kdrew: parse the column name to get the condition as label
            file_dict[df.columns[0].split(args.fraction_name_sep)[args.parse_fraction_name.index('condition')]] = df
        else:
            #kdrew: default to the filename as label
            file_dict[filename] = df


        #kdrew: find max number of fractions
        linspace_max = max(linspace_max, df.shape[1]) 


    if args.annotation_file != None:
        annotations_table = pd.read_csv(args.annotation_file, sep=args.annotation_file_sep, index_col=args.annotation_file_index)
    else:
        annotations_table = None
    #hs_uniprot_table['genename'] = [str(gn).split()[0] for gn in hs_uniprot_table['Gene names'].values]

    #def plot_fractions(input_uids, cntl_df, rnaseA_df):
    uids = []
    for uid in args.proteins:
        for label in file_dict.keys():
            if uid in file_dict[label].index and uid not in uids: 
                uids.append(uid)

    #kdrew: if no ids are present in any dataframe exit
    if len(uids) == 0:
        print "no ids present in dataframes: %s" % ' '.join(args.proteins)
        return 

    std_fractions = [float(x) for x in args.standard_fractions]

    #kdrew: set size of figure so it doesn't get super squished with large number of proteins (an inch per protein seems about right)
    figsize = [6.4, 4.8]
    if len(uids) > 5:
        figsize[1] = len(uids)
    f, axarr = plt.subplots(len(uids), sharex=True, figsize=figsize)

    for i, uid in enumerate(uids):
        x = np.linspace(1,linspace_max,num=linspace_max)
        #print x
        #kdrew: store max value to annotate standards
        max_val = 0
        for label in file_dict.keys():
            try:
                max_val = max(max_val, max(file_dict[label].loc[uid].values)) 
            except KeyError:
                continue

        for j, label in enumerate(file_dict.keys()):
            if args.parse_fraction_name != None:
                #for c in df.columns:
                #    c_split = c.split(args.fraction_name_sep)
                #    print c_split
                #    print c_split[args.parse_fraction_name.index('fraction')]
                #    #[args.parse_fraction_name.index('fraction')]
                x = [c.split(args.fraction_name_sep)[args.parse_fraction_name.index('fraction')] for c in df.columns]
                try:
                    x = [int(xx) for xx in x]
                except  ValueError:
                    x = [int(re.findall(r'\d+', xx)[0]) for xx in x]
                    #re.split('[a-f]+', xx, flags=re.IGNORECASE)


            #kdrew: try to see if id is in data frame
            try:
                #kdrew: try to plot multiple axes
                try:
                    axarr[i].plot(x,file_dict[label].loc[uid].values, color=list(it.islice(it.cycle(args.colors),j+1))[j], label=label)
                    #kdrew: only annotate with standards on the first plot
                    if j == 0:
                        for k, std_frac in enumerate(std_fractions):
                            axarr[i].axvline(std_frac, color="black", ls='dashed', alpha=0.5, linewidth=0.5)
                            if i == 0:
                                axarr[i].annotate(args.standard_mass[k], xy=(std_frac, max_val), xycoords='data', xytext=(-2.5, 0.0), textcoords='offset points', ha="center", rotation=90, fontsize=6,)
                except TypeError:
                    axarr.plot(x,file_dict[label].loc[uid].values, color=list(it.islice(it.cycle(args.colors),j+1))[j], label=label)
                    #kdrew: only annotate with standards on the first plot
                    if j == 0:
                        for k,std_frac  in enumerate(std_fractions):
                            axarr.axvline(std_frac, color="black", ls='dashed', alpha=0.5, linewidth=0.5)
                            axarr.annotate(args.standard_mass[k], xy=(std_frac, max_val+(max_val*0.01)), xycoords='data', xytext=(-2.5, 0.0), textcoords='offset points', ha="center", rotation=90, fontsize=6,)
                #ax = sns.lineplot(x="Time (min)", y="Value (mAU)", data=file_dict[label], color=it.cycle(args.colors)[i], label=label)
            except KeyError:
                continue

        try:
            if annotations_table is not None:
                genename = annotations_table.loc[uid][args.id_column]
            else:
                genename = uid
        except KeyError, AttributeError:
            genename = uid
        try:
            axarr[i].set_ylabel(genename, rotation=0, y=1.08)
            if not args.plot_y_axis:
                axarr[i].get_yaxis().set_ticks([])
        except TypeError:
            axarr.set_ylabel(genename, rotation=0, y=1.08)
            if not args.plot_y_axis:
                axarr.get_yaxis().set_ticks([])

    #kdrew: remove xticks
    try:
        for i in range(len(axarr)):
            #axarr[i].set_xticks([])
            axarr[i].tick_params(axis=u'both', which=u'both',length=0)
    except TypeError:
            #axarr.set_xticks([])
            axarr.tick_params(axis=u'both', which=u'both',length=0)
    
    sns.despine( left=True, bottom=True )

    plt.legend(loc="upper right",fontsize=8)


    #plt.tight_layout()
    #plt.subplots_adjust(hspace = .001)
    plt.savefig(args.output_filename)


if __name__ == "__main__":
    main()


