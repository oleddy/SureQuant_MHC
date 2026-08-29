import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import argparse

def apply_limits(t, x, min_t, max_t):
    min_index = 0
    max_index = 0
    while (t[min_index + 1] <= min_t) and (min_index < (len(t) - 1)):
        min_index += 1
    while (t[max_index] <= max_t) and (max_index < (len(t) - 1)):
        max_index += 1
    return t[min_index:(max_index + 1)], x[min_index:(max_index + 1)]

if __name__ == '__main__':

    buffer = 0.2
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'input chromatogram table', required = True)
    parser.add_argument('-p', help = 'peptide target list', required = True) #each entry must specify time window, mock file, and +Mtb file

    colors = ['red', 'orange', 'gold', 'deepskyblue', 'blue', 'fuchsia']

    #if provided, use a txt file of peptide targets. Otherwise, plot all peptides in the dataset
    args = parser.parse_args()
    data = pd.read_csv(args.i, sep = '\t').drop_duplicates()
    peptides = pd.read_csv(args.p)

    for i, peptide in peptides.iterrows():
        plt.rc('font', size = 18)
        f, ax = plt.subplots(2, 2, figsize = (10, 10))


        intensities = []
        ts = []
        Mtb_data = data.loc[np.logical_and(data['PeptideModifiedSequence'] == peptide['Sequence'], data['FileName'] == peptide['Mtb file'])]
        mock_data = data.loc[np.logical_and(data['PeptideModifiedSequence'] == peptide['Sequence'], data['FileName'] == peptide['mock file'])]

        for j, cond in enumerate(['light', 'heavy']):
            #Get mock data for this condition and peptide
            charge_data = mock_data.loc[mock_data['PrecursorCharge'] == min(mock_data['PrecursorCharge'])]
            condition_rows = charge_data.loc[charge_data['IsotopeLabelType'] == cond]
            fragment_rows = condition_rows.loc[condition_rows['FragmentIon'] != 'precursor']
            mock_fragment_intensity = {(row['FragmentIon'], row['ProductCharge']) : [float(intensity) for intensity in row['Intensities'].split(',')] for k, row in fragment_rows.iterrows()}
            mock_fragment_t = {(row['FragmentIon'], row['ProductCharge']) : [float(t) for t in row['Times'].split(',')] for k, row in fragment_rows.iterrows()}

            min_t = peptide['Min time']
            max_t = peptide['Max time']

            frags = mock_fragment_intensity.keys()
            for k, frag in enumerate(frags):
                t, x = apply_limits(mock_fragment_t[frag], mock_fragment_intensity[frag], min_t - buffer, max_t + buffer)
                ax[j, 0].plot(t, x, color = colors[k])
            ax[j,0].set_xlim(min_t - buffer, max_t + buffer)
            ax[j,0].set_xticks(np.arange(min_t, max_t + 0.5, step = 1))
            ax[j,0].ticklabel_format(axis = 'y', style = 'scientific', scilimits = (0,3), useMathText = True)
            #Get Mtb data for this condition and peptide
            charge_data = Mtb_data.loc[Mtb_data['PrecursorCharge'] == min(Mtb_data['PrecursorCharge'])]
            condition_rows = charge_data.loc[charge_data['IsotopeLabelType'] == cond]
            fragment_rows = condition_rows.loc[condition_rows['FragmentIon'] != 'precursor']
            Mtb_fragment_intensity = {(row['FragmentIon'],row['ProductCharge']): [float(intensity) for intensity in row['Intensities'].split(',')] for k, row in fragment_rows.iterrows()}
            Mtb_fragment_t = {(row['FragmentIon'],row['ProductCharge']) : [float(t) for t in row['Times'].split(',')] for k, row in fragment_rows.iterrows()}

            frags = Mtb_fragment_intensity.keys()
            for k, frag in enumerate(frags):
                t, x = apply_limits(Mtb_fragment_t[frag], Mtb_fragment_intensity[frag], min_t - buffer, max_t + buffer)
                ax[j, 1].plot(t, x, color = colors[k])
            ax[j,1].set_xlim(min_t - buffer, max_t + buffer)
            ax[j,1].set_xticks(np.arange(min_t, max_t + 0.5, step = 1))
            ax[j,1].ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,3), useMathText = True)

        max_heavy_ylim = max([ax[0,0].get_ylim()[1],ax[0,1].get_ylim()[1]])
        max_heavy_ylim = max([max_heavy_ylim, 1e4])
        print(max_heavy_ylim)
        ax[0,0].set_ylim(0,max_heavy_ylim)
        ax[0,1].set_ylim(0,max_heavy_ylim)

        max_light_ylim = max([ax[1,0].get_ylim()[1],ax[1,1].get_ylim()[1]])
        max_light_ylim = max([max_light_ylim, 1e4])
        print(max_light_ylim)
        ax[1,0].set_ylim(0,max_light_ylim)
        ax[1,1].set_ylim(0,max_light_ylim)

        ax[0,0].legend([frag[0] + ' +' + str(frag[1]) for frag in frags])


        shadowaxes = f.add_subplot(111, xticks=[], yticks=[], frame_on=False)
        shadowaxes.set_xlabel('Retention time (min)', labelpad = 36)
        shadowaxes.set_ylabel('Intensity', labelpad = 60)
        f.savefig(peptide['Sequence'] + '_chromatograms.pdf')
