'''Spectrum Annotator

Reads in a pair of text files (.txt) of annotated fragment ion (MS/MS) peaks. Input must be formatted as three
tab-separated columns: m/z (mass-charge ratio, float), i (intensity in ion current, float), and matches
(comma-separated labels assigned by database search, string).

Generates a figure comparing annotated MS/MS spectra. Annotated peaks
are colored according to fragment ion type, and the most intense matched peak per 100 m/z is annotated with a
text object and arrow.

Supports the following fragment ion types: y, b, a, precursor (M+kH), internal, TMT, immonium.
Supports the following neutral losses (NLs): phosphoric/metaphosphoric acid (Phos), water (H2O), ammonia (NH3).
If desired, user may manually annotate TMT peaks, re-annotate NL of phosphoric/metaphosphoric acid as H3PO4/HPO3,
and manually add/remove annotations for any of the supported types listed above.

Written by: Cameron Flower
Last edited: 4/13/21
Edited: 6/10/21 by Owen Leddy

'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.lines import Line2D
import argparse
from os.path import join
from os import listdir

parser = argparse.ArgumentParser()
parser.add_argument('-b', help = 'path to biological peptide spectra', required = True)
parser.add_argument('-s', help = 'path to synthetic peptide spectra', required = True)
parser.add_argument('-o', help = 'path to output', required = True)

args = parser.parse_args()

########################################################################################################################
# user-defined args




# m/z increment for text annotation (i.e. top peak per this m/z-sized window); small values (<100) risk overcrowding
INCREMENT = 50
# horizontal and vertical dimensions of figure
FIGSIZE = (12,5)

########################################################################################################################

#initialize m/z bounds -> will be updated according to min and max of spectra
min_mz = 1000.
max_mz = 0.

synthetic_files = listdir(args.s)
bio_files = listdir(args.b)

file_pairs = []

for synthetic_file in synthetic_files:
    if synthetic_file != '.DS_Store':
        peptide = synthetic_file.split('_')[0]
        bio_match = None
        for bio_file in bio_files:
            if peptide in bio_file.split('_'):
                bio_match = bio_file
        if bio_match:
            file_pairs.append((synthetic_file, bio_match, peptide))

for synthetic_file, bio_file, peptide in file_pairs:
    TITLE = peptide
    fig = plt.figure(figsize=FIGSIZE)
    print(TITLE)
    for PATH_INPUT, sign in zip([join(args.b, bio_file), join(args.s, synthetic_file)], [1, -1]):
        # read in annotated spectrum file, normalize intensity, reformat
        df_ions = pd.read_csv(PATH_INPUT, sep='\t')
        df_ions.columns = [x.strip() for x in df_ions.columns]
        df_non_background_ions = df_ions.loc[df_ions['matches'].notna()]

        df_ions['i'] = df_ions['i']*100/df_non_background_ions['i'].max()
        df_ions['matches'] = df_ions['matches'].str.replace(' - ', '-', regex=False)
        df_ions['matches'] = df_ions['matches'].str.replace('1+', '+', regex=False)

        # separate ions based on whether or not they matched (are annotated)
        df_unmatched = df_ions[df_ions['matches'].isnull()].reset_index(drop=True)
        df_matched = df_ions[~df_ions['matches'].isnull()].reset_index(drop=True)

        # for each matched ion, re-annotate it with a single label
        label_list = ['']*df_matched.shape[0]
        for i, row in df_matched.iterrows():
            label = row['matches'].split(', ')
            # prioritize TMT
            if 'TMT' in label:
                label_list[i] = 'TMT'
                continue
            # next, prioritize precursor ions
            newlabel = [l for l in label if 'M' in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            # next, prioritize y, b, and a ions (in that order, by expected stability)
            # need to disregard neutral loss ions for now
            newlabel = [l for l in label if 'y ' in l and 'Phos' not in l and 'H2O' not in l and 'NH3' not in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            newlabel = [l for l in label if 'b ' in l and 'Phos' not in l and 'H2O' not in l and 'NH3' not in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            newlabel = [l for l in label if 'a ' in l and 'Phos' not in l and 'H2O' not in l and 'NH3' not in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            # next, prioritize NL phosphate ions
            newlabel = [l for l in label if 'Phos' in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            # next, prioritize other NL ions
            newlabel = [l for l in label if 'H2O' in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            newlabel = [l for l in label if 'NH3' in l]
            if len(newlabel) > 0:
                label_list[i] = newlabel[0]
                continue
            # finally, label internal fragments
            label_list[i] = label[0]
            continue
        df_matched['matches'] = label_list

        # for each top peak per 100 m/z, get label
        new_min = math.floor(df_matched['m/z'].min()/100)*100
        if min_mz > new_min:
            min_mz = new_min
        new_max = math.ceil(df_matched['m/z'].max()/100)*100
        if new_max > max_mz:
            max_mz = new_max
        mz = new_min
        annotations = {}    # key: (x,y) of annotated peak; val: (label,x,y) for text annotation and arrow
        while mz < new_max:
            # don't annotate TMT ions or internal fragments
            df_annotate = df_matched[~df_matched['matches'].str.contains('TMT')]
            df_annotate = df_annotate[~df_annotate['matches'].str.contains('ya')]
            df_annotate = df_annotate[~df_annotate['matches'].str.contains('yb')]
            df_annotate = df_annotate[df_annotate['m/z'].between(mz,mz+INCREMENT)].sort_values('i', ascending=False)
            if not df_annotate.empty:
                df_annotate = df_annotate.reset_index(drop=True).iloc[0]
                x,y,label = df_annotate['m/z'], df_annotate['i']*sign, df_annotate['matches']
                df_allpeaks = df_ions[(df_ions['m/z'] > mz) & (df_ions['m/z'] <= mz+INCREMENT)].sort_values('i', ascending=False)
                df_allpeaks = df_allpeaks.loc[df_allpeaks['matches'].notna()]
                df_maxpeak = df_allpeaks.reset_index(drop=True).iloc[0]
                x_arrow = mz + INCREMENT/2
                y_arrow = df_maxpeak['i']*sign
                annotations[(x, y)] = (label, x_arrow, y_arrow)
            mz += INCREMENT

        # separate ion types for peak coloring
        df_precursor = df_matched[df_matched['matches'].str.contains('[M', regex=False)]
        df_matched = df_matched.drop(df_precursor.index)
        df_tmt = df_matched[df_matched['matches'].str.contains('TMT')]
        df_matched = df_matched.drop(df_tmt.index)
        df_ab = df_matched[df_matched['matches'].str.contains('a ') | df_matched['matches'].str.contains('b ')]
        df_matched = df_matched.drop(df_ab.index)
        df_y = df_matched[df_matched['matches'].str.contains('y ')]
        df_matched = df_matched.drop(df_y.index)
        df_internal = df_matched[df_matched['matches'].str.contains('ya') | df_matched['matches'].str.contains('yb')]
        df_matched = df_matched.drop(df_internal.index)
        df_immonium = df_matched[df_matched['matches'].str.contains('[A-Z][a-z][a-z]', regex = True)]
        df_matched = df_matched.drop(df_immonium.index)

        # plot and color all peaks
        df_list = [df_unmatched, df_internal, df_ab, df_y, df_precursor, df_tmt, df_immonium]
        colors = ['darkgray', 'orange', 'r', 'b', 'g', 'deepskyblue', 'fuchsia']

        for df,c in zip(df_list, colors):
            if not df.empty:
                _, stemlines, _ = plt.stem(df['m/z'], df['i']*sign, c, markerfmt=' ', basefmt=' ')
                plt.setp(stemlines, 'linewidth', 0.8)

        # annotate peaks
        for (x, y), (text, x_arrow, y_arrow) in annotations.items():
            # check if NL is present
            NLs = []
            if 'H2O' in text: NLs.append('H' + '$_{{2}}$' + 'O')
            if 'NH3' in text: NLs.append('NH' + '$_{{3}}$')
            if 'Phos' in text: NLs.append('Phos')               # default annotation from PD for NL of (meta)phosphoric acid
            if 'HPO3' in text: NLs.append('HPO' + '$_{{3}}$')   # user can manutally re-annotate NL of phos to be more specific
            if 'H3PO4' in text: NLs.append('H' + '$_{{3}}$' + 'PO' + '$_{{4}}$')
            NLs = '\n-'.join(NLs)
            text_list = text.split('-')[0].split(' ')
            iontype = text_list[0]
            # check if precursor, which has no position (subscript)
            if 'M' in iontype:
                # charge = text_list[1]
                # text_final = iontype + '$^{{{}}}$'.format(charge[1:-1])
                text_final = iontype
            elif len(text_list) > 1:
                position, charge = text_list[1], text_list[2]
                text_final = iontype + '$_{{{}}}$'.format(position[1:-1]) + '$^{{{}}}$'.format(charge[1:-1])
            else:
                text_final = iontype
            if NLs != '': text_final = text_final + '\n-' + NLs
            if sign == 1:
                offset = 5
            else:
                offset = -24
            plt.annotate(text_final, xy=(x,y), xytext=(x_arrow,y_arrow + offset), ha='center', va='bottom',
                            rotation=0, fontsize=10, arrowprops=dict(facecolor='black', arrowstyle='->'),
                            bbox=dict(pad=1, facecolor='none', edgecolor='none'))

    # format plot

    plot_xmin = min_mz
    plot_xmax = max_mz
    plt.plot([min_mz, max_mz], [0, 0], '-k', linewidth = 1)
    plt.xlim(plot_xmin, plot_xmax)
    plt.ylim(-130, 130)
    plt.xlabel('$\it{m/z}$', size=14)
    plt.ylabel('Intensity (%)', size=14)
    plt.xticks(np.arange(plot_xmin, plot_xmax+100, 100), size=10)
    plt.yticks(np.arange(-100,110,10), size=10)
    ax = plt.gca()
    ax.margins(y = 0.1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    if TITLE == '': TITLE = ' '
    # plt.title(TITLE, y=1.1, size=18)

    # build legend based on which ion types are present
    legend_elements = []
    df_list = [df_ab, df_y, df_precursor, df_internal, df_tmt, df_immonium, df_unmatched]
    colors = ['r', 'b', 'g', 'orange', 'deepskyblue','fuchsia', 'darkgray']
    legend_labels = ['N-terminal\n(a- or b-type)', 'C-terminal\n(y-type)', 'Precursor',
                        'Internal', 'TMT Reporter', 'Immonium', 'Unmatched']
    for df, c, l in zip(df_list, colors, legend_labels):
        if not df.empty: legend_elements.append(Line2D([0], [0], color=c, lw=2, label=l))
    lgd = plt.legend(handles=legend_elements, title=r'$\bf{Ion\ Type}$', loc='center right',
                        ncol=1, frameon=False, prop={'size': 10}, bbox_to_anchor=(1.2,0.5))
    lgd._legend_box.align = 'left'

    # desired path for generated figure
    PATH_OUTPUT = join(args.o, peptide + '.pdf')
    plt.tight_layout()
    # save to output path in png format (adjust resolution as needed)
    plt.savefig(PATH_OUTPUT,
                dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')
