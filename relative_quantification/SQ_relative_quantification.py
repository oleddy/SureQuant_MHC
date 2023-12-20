import pandas as pd
import argparse
import numpy as np
from os.path import join

parser = argparse.ArgumentParser()
parser.add_argument('-d', help = 'path to quantification data', required = True)
parser.add_argument('-t', help = 'path to target list', required = True)
parser.add_argument('-s', help = 'hipMHCs standards table', required = True)
parser.add_argument('-c', help = 'conditions table', required = True)
parser.add_argument('-o', help = 'output path', required = True)
parser.add_argument('-n', help = 'number of fragment ions to use for quantification', required = False, default = 3)

if __name__ == '__main__':

    args = parser.parse_args()

    target_peptides = pd.read_csv(args.t)
    skyline_data = pd.read_csv(args.d)
    conditions_table = pd.read_csv(args.c)

    #hipMHC peptide table has columns light annotated seq, heavy annotated seq, charge
    hipMHC_peptides = pd.read_csv(args.s)

    top_n_fragments = int(args.n)

    #throw out data not pertaining to targets of interest or hipMHCs
    skyline_data = skyline_data.loc[skyline_data['Peptide Modified Sequence'].isin(target_peptides['Peptide']) | skyline_data['Modified Sequence'].isin(hipMHC_peptides['Light Annotated Seq']) | skyline_data['Modified Sequence'].isin(hipMHC_peptides['Heavy Annotated Seq'])]
    sil_data = skyline_data.loc[skyline_data['Isotope Label Type'] == 'heavy']
    targets_data = skyline_data.loc[skyline_data['Isotope Label Type'] == 'light']

    targets_data = targets_data.drop_duplicates(subset = ['Modified Sequence', 'Fragment Ion', 'Precursor Charge', 'Replicate'])
    sil_data = sil_data.drop_duplicates(subset = ['Modified Sequence', 'Fragment Ion', 'Precursor Charge', 'Replicate'])

    conditions_dict = {condition['Filename'] : condition['Condition'] for i, condition in conditions_table.iterrows()}
    charge_dict = {target['Peptide'] : target['Charge'] for i, target in target_peptides.iterrows()}

    reference_condition = conditions_table.loc[conditions_table['Reference'] == True].iloc[0]['Filename']

    top_n_ions_targets = {peptide : [] for peptide in target_peptides['Peptide']}
    top_n_ions_hipMHCs = {peptide[0] : [] for _, peptide in hipMHC_peptides.iterrows()}

    #identify top n ions for targets
    reference_data = targets_data.loc[targets_data['Replicate'] == reference_condition]
    fragment_data = reference_data.loc[reference_data['Fragment Ion'] != 'precursor']
    for i, target in target_peptides.iterrows():
        target_data = fragment_data.loc[fragment_data['Peptide Modified Sequence'] == target['Peptide']]
        target_data = target_data.sort_values(by = 'Area', ascending = False)
        top_n = target_data.iloc[0:top_n_fragments]['Fragment Ion']
        top_n_ions_targets[target['Peptide']] = list(top_n)

    #identify top n ions for hipMHCs
    reference_data = sil_data.loc[sil_data['Replicate'] == reference_condition]
    fragment_data = reference_data.loc[reference_data['Fragment Ion'] != 'precursor']
    for _, (light_seq, heavy_seq, charge) in hipMHC_peptides.iterrows():
        target_data = fragment_data.loc[fragment_data['Modified Sequence'] == heavy_seq]
        target_data = target_data.sort_values(by = 'Area', ascending = False)
        top_n = target_data.iloc[0:top_n_fragments]['Fragment Ion']
        top_n_ions_hipMHCs[light_seq] = list(top_n)


    #get raw values for targets
    raw_target_values = {target : pd.DataFrame(np.zeros((top_n_fragments, len(conditions_table['Filename']))), columns = conditions_table['Filename'], index = top_n_ions_targets[target]) for target in target_peptides['Peptide']} #create dictionary mapping each peptide to a dataframe where fragment ions are rows and conditions are columns

    normalized_target_values = pd.DataFrame(np.zeros((len(target_peptides['Peptide']), len(conditions_table['Filename']))), columns = conditions_table['Filename'], index = target_peptides['Peptide'])
    in_top_n = [row['Fragment Ion'] in top_n_ions_targets[row['Modified Sequence']] for _, row in targets_data.iterrows()] #boolean mask for if the fragment ion is in the top n
    top_n_data = targets_data.loc[in_top_n]
    top_n_data.rename({'Area' : 'Light Area', 'Modified Sequence' : 'Mod Seq'}, inplace = True, axis = 'columns')
    top_n_data = pd.merge(top_n_data, sil_data[['Peptide Modified Sequence', 'Fragment Ion', 'Precursor Charge', 'Replicate', 'Area', 'Modified Sequence']], on = ['Peptide Modified Sequence', 'Fragment Ion', 'Precursor Charge', 'Replicate'])
    top_n_data.rename({'Area' : 'Heavy Area'}, inplace = True, axis = 'columns')
    top_n_data['H/L ratio'] = np.array(top_n_data['Light Area'])/np.array(top_n_data['Heavy Area'])

    for i, row in top_n_data.iterrows():
        if (row['Precursor Charge'] == charge_dict[row['Peptide Modified Sequence']]) and (row['Replicate'] in list(conditions_table['Filename'])):
            raw_target_values[row['Peptide Modified Sequence']][row['Replicate']].loc[row['Fragment Ion']] = row['H/L ratio']
    for target in raw_target_values:
        for ion in top_n_ions_targets[target]:
            raw_target_values[target].loc[ion] = np.array(raw_target_values[target].loc[ion])/float(raw_target_values[target].loc[ion][reference_condition]) #get ratio for each ion
            average_ratios = np.array(raw_target_values[target]).mean(axis = 0)
            normalized_target_values.loc[target] = average_ratios

    #get raw values for hipMHCs
    raw_standard_values = {light_seq : pd.DataFrame(np.zeros((top_n_fragments, len(conditions_table['Filename']))), columns = conditions_table['Filename'], index = top_n_ions_hipMHCs[light_seq]) for _, (light_seq, heavy_seq, charge) in hipMHC_peptides.iterrows()} #create dictionary mapping each peptide to a dataframe where fragment ions are rows and conditions are columns

    for _, (light_seq, heavy_seq, charge) in hipMHC_peptides.iterrows():
        light_data = sil_data.loc[sil_data['Modified Sequence'] == light_seq]
        top_n_light = light_data.loc[light_data['Fragment Ion'].isin(top_n_ions_hipMHCs[light_seq])]
        top_n_light.rename({'Area' : 'Light Area'}, inplace = True, axis = 'columns')

        heavy_data = sil_data.loc[sil_data['Modified Sequence'] == heavy_seq]
        top_n_heavy = heavy_data.loc[heavy_data['Fragment Ion'].isin(top_n_ions_hipMHCs[light_seq])]
        top_n_heavy.rename({'Area' : 'Heavy Area', 'Modified Sequence' : 'Mod Seq'}, inplace = True, axis = 'columns')
        top_n_data = pd.merge(top_n_light, top_n_heavy, on = ['Peptide Modified Sequence', 'Fragment Ion', 'Precursor Charge', 'Replicate'])
        
        top_n_data['H/L ratio'] = np.array(top_n_data['Light Area'], dtype = np.float128)/np.array(top_n_data['Heavy Area'], dtype = np.float128)
        
        for i, row in top_n_data.iterrows():
            if row['Replicate'] in list(conditions_table['Filename']):
                raw_standard_values[light_seq][row['Replicate']][row['Fragment Ion']] = row['H/L ratio']

    hipMHC_light_sequences = [standard[0] for _, standard in hipMHC_peptides.iterrows()]
    #normalize hipMHC intensities to the reference condition
    normalized_standard_values = pd.DataFrame(np.zeros((len(hipMHC_peptides), len(conditions_table['Filename']))), columns = conditions_table['Filename'], index = hipMHC_light_sequences)

    #divide target intensities by average hipMHC relative intensities    
    for _, (light_seq, heavy_seq, charge) in hipMHC_peptides.iterrows():
        for ion in top_n_ions_hipMHCs[light_seq]:
            raw_standard_values[light_seq].loc[ion] = np.array(raw_standard_values[light_seq].loc[ion], dtype = np.float128)/np.float128(raw_standard_values[light_seq].loc[ion][reference_condition]) #get ratio for each ion
            average_ratios = np.array(raw_standard_values[light_seq]).mean(axis = 0)
            normalized_standard_values.loc[light_seq] = average_ratios
    normalized_standard_values.to_csv(join(args.o, 'hipMHC_correction_factors.csv')) #export hipMHC correction factors
    average_norm_factor = np.array(normalized_standard_values).mean(axis = 0)
    normalized_target_values = normalized_target_values / average_norm_factor
    normalized_target_values.columns = conditions_table['Condition']

    #export normalized data
    normalized_target_values.to_csv(join(args.o, 'normalized_target_intensities.csv'))