#!/home/ptgmat003/ve/bin/python

import config
import pandas as pd
from pandas import Series
import os
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import rfunc

processed = pd.read_csv(config.processed_peptides)

novel     = set(processed[processed['Novel']     == '+']['Sequence'])

annotated = set(processed[processed['Annotated'] == '+']['Sequence'])

pepseqs =Series(processed.reset_index()['index'], index= processed.Sequence).to_dict()

for i in novel:
    assert not i in annotated

pep_stats = config.output + 'PEP_analysis_all_PSMs/'
anova     = config.output + 'PEP_analysis_all_PSMs/PEP_ANOVA/'
kw        = config.output + 'PEP_analysis_all_PSMs/PEP_Kruskal-Wallis/'
density   = config.output + 'PEP_analysis_all_PSMs/PEP_density/'
qqnorm    = config.output + 'PEP_analysis_all_PSMs/PEP_qqnorm/'
figures   = config.output + 'PEP_analysis_all_PSMs/PEP_boxplots/'
tables    = config.output + 'PEP_analysis_all_PSMs/'

if not os.path.exists(figures):
    os.makedirs(figures)
if not os.path.exists(tables):
    os.makedirs(tables)
if not os.path.exists(pep_stats):
    os.makedirs(pep_stats)
if not os.path.exists(anova):
    os.makedirs(anova)
if not os.path.exists(kw):
    os.makedirs(kw)
if not os.path.exists(density):
    os.makedirs(density)
if not os.path.exists(qqnorm):
    os.makedirs(qqnorm)

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

def pepfig(database, data, names, data_,names_, output): #db name, list of lists, names of lists, output dir, min reps
    
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
    ax.set_title('A. All identified PSMs'.format(database))
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Create an axes instance
    ax_ = fig.add_subplot(122)
    # Create the boxplot
    bp_ = ax_.boxplot(data_, patch_artist=True)
    ## change outline color, fill color and linewidth of the boxes
    for box in bp_['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=1)
        # change fill color
        box.set( facecolor = '#1b9e77' )
        

    ## change color and linewidth of the whiskers
    for whisker in bp_['whiskers']:
        whisker.set(color='#7570b3', linewidth=1)

    ## change color and linewidth of the caps
    for cap in bp_['caps']:
        cap.set(color='#7570b3', linewidth=1)

    ## change color and linewidth of the medians
    for median in bp_['medians']:
        median.set(color='#b2df8a', linewidth=2)
        #median.set(linewidth=2)
    ## change the style of fliers and their fill
    for flier in bp_['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5)

    ## Custom x-axis labels
    ax_.set_xticklabels(names_)
    #ax.set_yticklabels('Posterior Error Probability (PEP)')
    
    ax_.set_title('B. At least two replicates')
    ## Remove top axes and right axes ticks
    ax_.get_xaxis().tick_bottom()
    ax_.get_yaxis().tick_left()

    # Save the figure
    fig.savefig(output + '{} PSM PEP.jpg'.format(database), bbox_inches='tight')
    fig.clf()


min_rep = config.peptides['replicates']


pepdf = pd.DataFrame()

data_dict = {}
count = 0

for database in config.searches:
    print database
    data_dict[database] = {}
    txt = config.searches[database]['txt'] +'peptides.txt'
    msms_txt = config.searches[database]['txt'] +'msms.txt'
    data = []
    names= []
    peptides = pd.read_csv(txt, sep = '\t', engine='python')
    msms = pd.read_csv(msms_txt, sep='\t',engine='python')

    peptides = replicates(peptides)

    pep = peptides[(peptides['Reverse'].isnull()) & (peptides['Potential contaminant'].isnull())]
    msms_pep = msms[(msms['Reverse'].isnull())]

    pep['TEMP'] = pep['Sequence'].apply(check)
    msms_pep['TEMP'] = msms['Sequence'].apply(check)


    scont_set       = peptides[(peptides['Potential contaminant'].notnull())&(peptides['replicates'] >=1)]['Sequence'].values.tolist()
    srev_set       = peptides[(peptides['Reverse'].notnull())&(peptides['replicates'] >=1)]['Sequence'].values.tolist()
    snovel_set     = pep[(pep['TEMP'] =='n')&(peptides['replicates'] >= 1) ]['Sequence'].values.tolist()
    sannotated_set = pep[(pep['TEMP'] =='a')&(peptides['replicates'] >= 1) ]['Sequence'].values.tolist()

    scont_set_       = peptides[(peptides['Potential contaminant'].notnull())&(peptides['replicates'] >= min_rep) ]['Sequence'].values.tolist()
    srev_set_       = peptides[(peptides['Reverse'].notnull())&(peptides['replicates'] >= min_rep) ]['Sequence'].values.tolist()
    snovel_set_     = pep[(pep['TEMP'] =='n')&(peptides['replicates'] >= min_rep) ]['Sequence'].values.tolist()
    sannotated_set_ = pep[(pep['TEMP'] =='a')&(peptides['replicates'] >= min_rep) ]['Sequence'].values.tolist()

    cont_set =  msms[msms['Sequence'].isin(scont_set)]['PEP'].values.tolist()
    rev_set = msms[msms['Sequence'].isin(srev_set)]['PEP'].values.tolist()
    novel_set = msms_pep[msms_pep['Sequence'].isin(snovel_set)]['PEP'].values.tolist()
    annotated_set = msms_pep[msms_pep['Sequence'].isin(sannotated_set)]['PEP'].values.tolist()

    cont_set_ =  msms[msms['Sequence'].isin(scont_set_)]['PEP'].values.tolist()
    rev_set_ = msms[msms['Sequence'].isin(srev_set_)]['PEP'].values.tolist()
    novel_set_  = msms_pep[msms_pep['Sequence'].isin(snovel_set_)]['PEP'].values.tolist()
    annotated_set_  = msms_pep[msms_pep['Sequence'].isin(sannotated_set_)]['PEP'].values.tolist()
    print len(rev_set), len(rev_set_)

    count += 1
    #data_dict[database]['Annotated'] = set([pepseqs[j] for j in annotated_set_]
    pepdf.loc['Contaminant PSMs', database] = len(cont_set)
    pepdf.loc['Contaminant PSMs\nmedian PEP score',database] = "%E" % np.median(cont_set)
    pepdf.loc['Contaminant PSMs\nPEP score STDEV',database] = "%E" % np.std(cont_set)
    count += 1

    count += 1
    #data_dict['Annotated hits'] = len(annotated_set)
    data.append(annotated_set)
    names.append('Annotated')
    pepdf.loc['Annotated PSMs', database] = len(annotated_set)
    pepdf.loc['Annotated PSMs\nmedian PEP score', database] = "%E" % np.median(annotated_set)
    pepdf.loc['Annotated PSMs\nPEP score STDEV',database] = "%E" % np.std(annotated_set)

    count += 1
    #data_dict['Novel hits'] = len(novel_set)
    data.append(novel_set)
    names.append('Novel')
    pepdf.loc['Novel PSMs', database] = len(novel_set)
    pepdf.loc['Novel PSMs\nmedian PEP score',database] = "%E" % np.median(novel_set)
    pepdf.loc['Novel PSMs\nPEP score STDEV',database] = "%E" % np.std(novel_set)
    count += 1

    #data_dict['Reverse hits'] = len(rev_set)
    data.append(rev_set)
    names.append('Reverse')
    pepdf.loc['Reverse PSMs',database] = len(rev_set)
    pepdf.loc['Reverse PSMs\nmedian PEP score',database] = "%E" % np.median(rev_set)
    pepdf.loc['Reverse PSMs\nPEP score STDEV',database] = "%E" % np.std(rev_set)

    data_ = []
    names_ = []

    count += 1
    #data_dict[database]['Annotated'] = set([pepseqs[j] for j in annotated_set_]
    pepdf.loc['Contaminant PSMs\n(2 or more replicates)', database] = len(cont_set_)
    pepdf.loc['Contaminant PSMs\nmedian PEP score\n(2 or more replicates)',database] = "%E" % np.median(cont_set_)
    pepdf.loc['Contaminant PSMs\nPEP score STDEV\n(2 or more replicates)',database] = "%E" % np.std(cont_set_)
    count += 1

    #data_dict[database]['Annotated'] = set([pepseqs[j] for j in annotated_set_])
    data_.append(annotated_set_)
    names_.append('Annotated')
    pepdf.loc['Annotated PSMs\n(2 or more replicates)', database] = len(annotated_set_)
    pepdf.loc['Annotated PSMs\nmedian PEP score\n(2 or more replicates)',database] = "%E" % np.median(annotated_set_)
    pepdf.loc['Annotated PSMs\nPEP score STDEV\n(2 or more replicates)',database] = "%E" % np.std(annotated_set_)

    count += 1
    #data_dict[database]['Novel'] = set([pepseqs[j] for j in novel_set_])
    data_.append(novel_set_)
    names_.append('Novel')
    pepdf.loc['Novel PSMs\n(2 or more replicates)', database] = len(novel_set_)
    pepdf.loc['Novel PSMs\nmedian PEP score\n(2 or more replicates)',database] = "%E" % np.median(novel_set_)
    pepdf.loc['Novel PSMs\nPEP score STDEV\n(2 or more replicates)',database] = "%E" % np.std(novel_set_)

    count += 1
    #data_dict[database]['Reverse'] = set([pepseqs[j] for j in rev_set_])
    data_.append(rev_set_)
    names_.append('Reverse')
    pepdf.loc['Reverse PSMs\n(2 or more replicates)',database] = len(rev_set_)
    pepdf.loc['Reverse PSMs\nmedian PEP score\n(2 or more replicates)',database] = "%E" % np.median(rev_set_)
    pepdf.loc['Reverse PSMs\nPEP score STDEV\n(2 or more replicates)',database] = "%E" % np.std(rev_set_)
    count += 1

    #combined = novel_set_ + annotated_set_
    #data_dict[database]['Novel & Annotated'] =set([pepseqs[j] for j in combined]) 

    outpath = anova + '{}_PEP_score_anova_tukey.txt'.format(database)
    rfunc.list_anova(names,data, names_, data_,'PEP', 'ANNOTATION',outpath)
    outpath = kw + '{}_PEP_score_kw_dunn.txt'.format(database)
    rfunc.list_kw_dunn(names,data, names_, data_,'PEP', 'ANNOTATION',outpath)

    #density plots
    outpath = density
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    rfunc.list_density(names,
                       data,
                       'All identified peptides\nall PSMs PEP score density curves',
                       names_,
                       data_,
                       'Only peptides identified in at least two replicates\nall PSMs PEP score density curves',
                       'Posterior Error Probability (PEP)',
                       outpath +'{}_MSMS_PEP_score_density.jpeg'.format('_'.join(database.split())))

    outpath = qqnorm +'_'.join(database.split())
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    rfunc.list_qqnorm(names,
                       data,
                       names_,
                       data_,
                       outpath)


    pepfig(database,data,names, data_, names_, figures)




pepdf.to_csv(tables +'Database PSM PEP scores.csv')







