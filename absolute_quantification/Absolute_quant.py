import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy 
import pandas as pd
import math
import seaborn as sns
import os
import argparse

def get_replicate(result_table, replicate_name):
    replicate_table=result_table[result_table['Replicate']=='spleenocyte_pulsed_SQ_20231004']
    return(replicate_table)

def calculation(H1):
    height1=H1['Height'][H1['Peak Rank']==2].values[0]
    height2=H1['Height'][H1['Peak Rank']==3].values[0]
    height3=H1['Height'][H1['Peak Rank']==4].values[0]
    avg_height=(height1+height2+height3)/3
    return(avg_height)

def convert_height(replicate_table,light_seq,H1_seq,H2_seq,H3_seq):
    light=replicate_table[replicate_table['Modified Sequence']==light_seq]
    H1=replicate_table[replicate_table['Modified Sequence']==H1_seq]
    H2=replicate_table[replicate_table['Modified Sequence']==H2_seq]
    H3=replicate_table[replicate_table['Modified Sequence']==H3_seq]
    avg_height_H1=calculation(H1)
    avg_height_H2=calculation(H2)
    avg_height_H3=calculation(H3)
    avg_height_light=calculation(light)
    return([avg_height_H1,avg_height_H2,avg_height_H3,avg_height_light])

def abline(slope, intercept,xmin,xmax):
    """Plot a line from slope and intercept"""
    x_vals = np.array([xmin,xmax])
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

def process_table(result_table,replicate_name,cell_num,filename,light_seq,H1_seq,H2_seq,H3_seq,heavy_input):
    replicate_table=get_replicate(result_table, replicate_name)
    print(replicate_table)
    [avg_height_H1,avg_height_H2,avg_height_H3,avg_height_light]=convert_height(replicate_table,light_seq,H1_seq,H2_seq,H3_seq)
    slope, intercept, r, p, se = scipy.stats.linregress([math.log10(heavy_input/100),math.log10(heavy_input/10),math.log10(heavy_input)],[math.log10(avg_height_H1),math.log10(avg_height_H2),math.log10(avg_height_H3)])
    light_value=10**((math.log10(avg_height_light)-intercept)/slope)
    copies=light_value*10**(-15)*6.02*10**(23)/(cell_num*10**6)
    fig, ax = plt.subplots(figsize=(3,2))
    x_vals=[math.log10(heavy_input/100),math.log10(heavy_input/10),math.log10(heavy_input)]
    yvals=[10**(intercept + slope * i) for i in x_vals]
    ax.loglog([heavy_input/100,heavy_input/10,heavy_input],yvals,'-.',color='#8eb4ba')
    ax.loglog([heavy_input/100,heavy_input/10,heavy_input],[avg_height_H1,avg_height_H2,avg_height_H3],'.',color='#31708e')
    ax.loglog(light_value,avg_height_light,'o',color='#9c3838')
    ax.annotate("r^2 = {:.4f}".format(r), (0.1, 3*10**5))
    ax.set_xlabel('Peptide Abundance(fmol)')
    ax.set_ylabel('Intensity')
    fig.savefig(args.dir+str(filename)+'.png', dpi=300,bbox_inches='tight')
    plt.show()
    print(replicate_name+' contains '+str(copies)+' copies of biological peptide per cell')
    return(copies)

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
# --dir: working directory, --file: Skyline export filename, --rep: Replicate name of the sample to quantify, --output Linear regression plot name, --cell: number of e^6 cells in IP, --light light peptide sequence, --H1-3 heavy peptide with 1-3SIL aa, --input: Highest spike in in fmol
parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, required=True)
parser.add_argument('--file', type=str,required=True)
parser.add_argument('--rep', type=str,required=True)
parser.add_argument('--output',type=str,required=True)
parser.add_argument('--cell',type=float,required=True)
parser.add_argument('--light',type=str,required=True)
parser.add_argument('--H1',type=str,required=True)
parser.add_argument('--H2',type=str,required=True)
parser.add_argument('--H3',type=str,required=True)
parser.add_argument('--input',type=float,required=True)

args = parser.parse_args()
os.chdir(args.dir)
table=pd.read_csv(args.file,sep=',',engine='python')
# print(table)
copies1=process_table(table,args.rep,args.cell,args.output,args.light,args.H1,args.H2,args.H3,args.input)
