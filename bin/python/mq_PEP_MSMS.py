#!/usr/bin/env python3

import pandas as pd
from pandas import Series
import os
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import rfunc
import sys
import yaml

config = yaml.load(open(sys.argv[1]).read())
output = sys.argv[2]

combined = pd.read_csv(output + '/combined.csv')

annotated = combined['Specific annotated peptides - all strains'].dropna().tolist()
annotated = [i.split('\n') for i in annotated ] 
annotated = [val for sublist in annotated for val in sublist if val != '']

novel = combined['Specific novel peptides - all strains'].dropna().tolist()
novel = [i.split('\n') for i in novel ] 
novel = [val for sublist in novel for val in sublist if val != '']

pep_stats = output + '/PEP_MSMS/'
kw        = output + '/PEP_MSMS/PEP_Kruskal-Wallis/'
density   = output + '/PEP_MSMS/PEP_density/'
figures   = output + '/PEP_MSMS/PEP_boxplots/'
tables    = output + '/PEP_MSMS/tables/'

if not os.path.exists(figures):
    os.makedirs(figures)
if not os.path.exists(tables):
    os.makedirs(tables)
if not os.path.exists(kw):
    os.makedirs(kw)
if not os.path.exists(density):
    os.makedirs(density)

def replicates(df):
    for row in df.iterrows():
        reps = 0
        for col in row[1].index:
            if col.startswith('Identification type'):
                if row[1][col] == 'By MS/MS':
                    reps += 1
        df.loc[row[0], 'replicates'] = reps
    return df

def check(value):
    if value in novel:
        return 'n'
    elif value in annotated:
        return 'a'

def pepfig(database, data, names, output): #db name, list of lists, names of lists, output dir, min reps
    #create the figure instance
    #fig = plt.figure(1, figsize=(9, 6))
    fig = plt.figure()
    # Create an axes instance
    a = fig.add_subplot(111)
    a.set_ylabel('Posterior error probability (PEP)')
    # Turn off axis lines and ticks of the big subplot
    a.spines['top'].set_color('none')
    a.spines['bottom'].set_color('none')
    a.spines['left'].set_color('none')
    a.spines['right'].set_color('none')
    a.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax = fig.add_subplot(121)
    # Create the boxplot
    bp = ax.boxplot(data, patch_artist=True)
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=1)
        # change fill color
        box.set( facecolor = '#1b9e77' )
    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=1)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=1)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)
        #median.set(linewidth=2)
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5)

    ## Custom x-axis labels
    ax.set_xticklabels(names)
    
    #ax.set_yticklabels('Posterior Error Probability (PEP)')
    ax.set_title('A. All identified PSMs strain {}'.format(database))
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Save the figure
    fig.savefig(output + '{}_PSM_PEP.png'.format(database), bbox_inches='tight')
    fig.clf()

min_rep = 1

pepdf = pd.DataFrame()

data_dict = {}

count = 0

txt = config['mq_txt'] + '/peptides.txt'

msms_txt = config['mq_txt'] + '/msms.txt'

peptides = pd.read_csv(txt, sep = '\t', engine='python')

global_msms = pd.read_csv(msms_txt, sep='\t',engine='python')

peptides = replicates(peptides)


for strain in config['strains']:
    strain_samples = config['files'][strain]
    msms = global_msms[global_msms['Raw file'].apply(lambda x : x in strain_samples) == True ]
    data_dict[strain] = {}

    data = []
    names= []

    pep = peptides[(peptides['Reverse'].isnull()) & (peptides['Potential contaminant'].isnull())]
    
    msms_pep = msms[(msms['Reverse'].isnull())]

    pep['TEMP'] = pep['Sequence'].apply(check)
    
    msms_pep['TEMP'] = msms['Sequence'].apply(check)

    scont_set  = peptides[(peptides['Potential contaminant'].notnull())&(peptides['replicates'] >=1)]['Sequence'].values.tolist()
    srev_set   = peptides[(peptides['Reverse'].notnull())&(peptides['replicates'] >=1)]['Sequence'].values.tolist()
    snovel_set = pep[(pep['TEMP'] =='n')&(peptides['replicates'] >= 1) ]['Sequence'].values.tolist()
    sannotated_set = pep[(pep['TEMP'] =='a')&(peptides['replicates'] >= 1) ]['Sequence'].values.tolist()

    cont_set =  msms[msms['Sequence'].isin(scont_set)]['PEP'].values.tolist()
    rev_set = msms[msms['Sequence'].isin(srev_set)]['PEP'].values.tolist()
    novel_set = msms_pep[msms_pep['Sequence'].isin(snovel_set)]['PEP'].values.tolist()
    annotated_set = msms_pep[msms_pep['Sequence'].isin(sannotated_set)]['PEP'].values.tolist()
    
    #rev_set = novel_set

    #data_dict[database]['Annotated'] = set([pepseqs[j] for j in annotated_set_]
    
    pepdf.loc['Contaminant PSMs', strain] = len(cont_set)
    pepdf.loc['Contaminant PSMs\nmedian PEP score', strain] = "%E" % np.median(cont_set)
    pepdf.loc['Contaminant PSMs\nPEP score STDEV', strain] = "%E" % np.std(cont_set)

    count += 1
    #data_dict['Annotated hits'] = len(annotated_set)
    
    data.append(annotated_set)
    names.append('Annotated')
    pepdf.loc['Annotated PSMs', strain] = len(annotated_set)
    pepdf.loc['Annotated PSMs\nmedian PEP score', strain] = "%E" % np.median(annotated_set)
    pepdf.loc['Annotated PSMs\nPEP score STDEV', strain] = "%E" % np.std(annotated_set)

    count += 1
    #data_dict['Novel hits'] = len(novel_set)
    data.append(novel_set)
    names.append('Novel')
    pepdf.loc['Novel PSMs', strain] = len(novel_set)
    pepdf.loc['Novel PSMs\nmedian PEP score',strain] = "%E" % np.median(novel_set)
    pepdf.loc['Novel PSMs\nPEP score STDEV', strain] = "%E" % np.std(novel_set)
    count += 1

    #data_dict['Reverse hits'] = len(rev_set)
    data.append(rev_set)
    names.append('Reverse')
    pepdf.loc['Reverse PSMs', strain] = len(rev_set)
    pepdf.loc['Reverse PSMs\nmedian PEP score', strain] = "%E" % np.median(rev_set)
    pepdf.loc['Reverse PSMs\nPEP score STDEV', strain] = "%E" % np.std(rev_set)

    print(names)
    for i in data:
        print(len(i))
    
    outpath = kw + '{}_PEP_score_kw_dunn.txt'.format(strain)
    
    rfunc.list_kw_dunn(names,data,'PEP', 'ANNOTATION',outpath) 
    
    outpath = density
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    #rfunc.list_density( data,
    #                    'All identified peptides\nall PSMs PEP score density curves',
    #                    outpath +'{}_MSMS_PEP_score_density.jpeg'.format('_'.join(strain.split())))

    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    pepfig(strain, data, names, figures)
    pepdf.to_csv(tables +'{}_PSM_PEP_scores.csv'.format(strain))





