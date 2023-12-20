import pandas as pd
import argparse
from itertools import product
import numpy as np
from os.path import join

parser = argparse.ArgumentParser()
parser.add_argument('-i', help = 'input file (Skyline report)', required = True)
parser.add_argument('-l', help = 'label offsets in Da (comma-separated)', required = True)
parser.add_argument('-c', help = 'charge states (comma-separated)', required = False, default = '2,3')
parser.add_argument('-n', help = 'number of fragment ions to include in pseudospectra', required = False, default = 6)
parser.add_argument('-o', help = 'output path', required = True)

if __name__ == '__main__':

    args = parser.parse_args()
    report = pd.read_csv(args.i)

    labels = args.l.split(',')
    labels = ['[+' + label.strip() +']' for label in labels]

    charges = str(args.c).split(',')
    charges = list(map(int, charges))

    n = int(args.n)

    unique_peptides = list(set(list(report['Modified Sequence'])))

    inclusion_dfs = [[pd.DataFrame() for charge in charges] for label in labels]
    pseudospectrum_dfs = [[pd.DataFrame() for charge in charges] for label in labels]
    
    for peptide in unique_peptides:
        label_index = 0
        for i, label in enumerate(labels):
            if label in peptide:
                label_index = i
        ions = report.loc[report['Modified Sequence'] == peptide]
        unique_ion = ions.loc[ions['Fragment Ion'] == 'precursor']
        charge_index = 0
        inclusion_dfs[label_index][charge_index] = pd.concat([inclusion_dfs[label_index][charge_index], unique_ion])

        ions = ions.loc[ions['Fragment Ion'] != 'precursor']
        ions = ions.sort_values('Area', ascending = False).iloc[0:n]
        pseudospectrum_dfs[label_index][charge_index] = pd.concat([pseudospectrum_dfs[label_index][charge_index], ions])

    #columns needed for inclusion list: Compound, m/z, Intensity Threshold
    for i, label_list in enumerate(inclusion_dfs):
        for j, inclusion_df in enumerate(label_list):
            inclusion_df_out = pd.DataFrame()
            if inclusion_df.shape[0] > 0:
                inclusion_df_out['Compound'] = inclusion_df['Modified Sequence']
                inclusion_df_out['m/z'] = inclusion_df['Precursor Mz']
                inclusion_df_out['Intensity Threshold'] = np.maximum(10000.,inclusion_df['Max Height']*0.01)
                inclusion_df_out.to_csv(join(args.o, labels[i] + '_z=' + str(charges[j]) + '_inclusion_list.csv'))

    #columns needed for pseudospectrum: Compound, m/z, Group ID
    for i, label_list in enumerate(pseudospectrum_dfs):
        for j, pseudospectrum_df in enumerate(label_list):
            pseudospectrum_df_out = pd.DataFrame()
            if pseudospectrum_df.shape[0] > 0:
                pseudospectrum_df_out['Compound'] = pseudospectrum_df['Fragment Ion']
                pseudospectrum_df_out['m/z'] = pseudospectrum_df['Product Mz']
                pseudospectrum_df_out['Group ID'] = pseudospectrum_df['Precursor Mz']
                pseudospectrum_df_out.to_csv(join(args.o, labels[i] + '_z=' + str(charges[j]) + '_pseudospectra.csv'))
