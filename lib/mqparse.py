#!/usr/bin/env python3
import matplotlib as mpl
mpl.use('agg')
import Bio; from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, translate
import pandas as pd
from pandas import DataFrame,  Series
import json
import collections; from collections import defaultdict
import tempfile
from matplotlib_venn import venn3, venn2, venn3_circles, venn2_circles
from matplotlib import pyplot as plt
import shutil
import numpy as np
import algo
import time
import blast
import re
import os
import rfunc
import subprocess
import yaml
import scikit_posthocs as ph
import scipy
import pickle
import mygene
import rpy2.robjects as ro
import sys
from io import StringIO

# Library to parse MaxQuant txt

def name_dct():
    names = {'sf_novel_m':'Six Frame Novel Only\n(M start)',
            'sf_novel':'Six Frame Novel Only',
            'ref_prot_alt':'Reference Proteome\n(Embl CDS Start AA Alternative)',
            'ref_prot':'Reference Proteome',
            'alt_gms':'GenemarkS Predictions\n(M start)',
            'gms':'GenemarkS Predictions',
            'raw_sf':'Six Frame Translation',
            'alt_sf':'Six Frame Translation\n(M start)'}
    return names

def split_join(value):
    value = value.split(';')
    value = '\n'.join(value)
    return value

def param_table(mapping_dct):
    names = name_dct()
    sum_dct = defaultdict()
    for i in mapping_dct:
        if (i != 'mapping') and (i != 'flatfile'):
            tup =mapping_dct[i]
            fasta = list(SeqIO.parse(tup[0],'fasta'))
            sum = pd.read_csv(tup[1]+'/parameters.txt', sep=None, engine='python')
            sum= sum.set_index('Parameter')
            sum.rename(columns ={'Value':i}, inplace= True)
            sum_dct[i]=sum
            
    summary = pd.concat(sum_dct.values(), axis =1)

    summary = summary[['raw_sf', 
        'alt_sf', 
        'ref_prot', 
        'gms', 
        'alt_gms',
        'sf_novel',
        'sf_novel_m',
        'ref_prot_alt']]
    summary.rename(columns=names, inplace=True)
    return summary

def summary_table(mapping_dct):
    names = name_dct()
    sum_dct = defaultdict()
    for i in mapping_dct:
        if (i != 'mapping') and (i != 'flatfile'):
            tup =mapping_dct[i]
            fasta = list(SeqIO.parse(tup[0],'fasta'))
            sum = pd.read_csv(tup[1]+'/summary.txt', sep=None, engine='python')
            sum['Fasta DB Size'] = len(fasta)
            sum= sum.set_index('Raw file')
            params = sum[['Enzyme', 
                        'Enzyme mode',
                        'Enzyme first search',
                        'Enzyme mode first search',
                        'Use enzyme first search',
                        'Variable modifications',
                        'Variable modifications first search',
                        'Use variable modifications first search',
                        'Multiplicity',
                        'Max. missed cleavages']].fillna(method='pad')
            
            params['Enzyme'] = params['Enzyme'].apply(split_join)
            params['Variable modifications'] = params['Variable modifications'].apply(split_join)
            params = params.loc['Total',:].dropna(how='any')
            series =  sum.loc['Total',:].dropna(how='any')
            params = params.append(series)
            params.name=i
         
            sum_dct[i]=params
    summary = pd.concat(sum_dct.values(), axis =1)
    summary = summary[['raw_sf', 
        'alt_sf', 
        'ref_prot', 
        'gms', 
        'alt_gms',
        'sf_novel',
        'sf_novel_m',
        'ref_prot_alt']]

    summary.rename(columns=names, inplace=True)
    return summary

def verify(df):    #check that the unique peptide occurred in more than onereplicate
    df1 =df.copy()
    columns = []
    for i in df1.columns:
        if i.startswith('Unique peptides '):
            columns.append(i)
    return columns

def list_kw_dunn(names, data, value, group, path):
    kw = scipy.stats.kruskal(*data)
    w = open(path + '/kw.txt', 'w')
    w.write(str(kw))
    w.close()
    post_hoc = pd.DataFrame(ph.posthoc_dunn(data, p_adjust = 'fdr_bh'))
    post_hoc.columns = names
    post_hoc.index = names
    post_hoc.to_csv(path + '/dunn_bh.csv')
    

def parse_protein_ids(val):
    val = val.split(';')
    if len(val) > 1:
        return val[0] +';...'
    else:
        return val[0]
        
def parse_ids(ids):
    ids = ids.split(';')
    new_ids = []
    for i in ids:
        print(i)
        if '|' in i:
            new = i.split('|')[1]
        else:
            new = i
        new_ids.append(new)
    new_ids = ';'.join(new_ids)
    return new_ids

def parse_groups(ids):
    ids = ids.split(';')
    new_ids = []
    for i in ids:
        if '|' in i:
            new = i.split('|')[2].split('_')[1]
        else:
            new = i
        new_ids.append(new)
    new_ids = ';'.join(new_ids)
    return new_ids

def bp(data, names, outfile, title, overlay=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Create the boxplot
    bp = ax.boxplot(data )#, patch_artist=True)
    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=1)
        # change fill color
        #box.set( facecolor = '#1b9e77' )
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
    ax.set_title(title)
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    # Draw the observations
    if overlay == True:
        for i in range(len(data)):
            y = data[i]
            #x = np.random.normal(1 +i, 0.04, size=len(y))
            x = [ i + 1  for _ in y ]
            #print(x)
            #x = np.random.normal(0.5 +i, 0.04, size=len(y))
            plt.plot(x, y, 'r.', alpha=0.2)
    # Save the figure
    fig.savefig(outfile, bbox_inches='tight')
    fig.clf()
    plt.close()

class mq_txt:
    def __init__(self, config, exclude_contaminants=True):
        c = "library('FSA')"; ro.r(c)
        with open(config) as f:
            self.config = yaml.load(f.read())
        print(self.config)
        self.txt_path = self.config['mq_txt']
        self.outdir = self.config['outdir']
        design = self.config['design']
        exclude_contaminants = self.config['exclude_contaminants']
        if not exclude_contaminants == True:
            assert exclude_contaminants == False
        
        
        self.protein_quantification = self.config['protein_quantification']
        assert self.protein_quantification in ['LFQ', 'iBAQ']
        if self.protein_quantification =='LFQ':
            self.protein_quantification ='LFQ intensity'
        self.protein_normalisation = self.config['protein_normalisation']
        self.protein_imputation = self.config['protein_imputation']
        
        self.peptide_quantification = self.config['peptide_quantification']
        assert self.peptide_quantification in ['Intensity']
        self.peptide_normalisation = self.config['peptide_normalisation']
        self.peptide_imputation = self.config['peptide_imputation']


        self.summary = pd.read_csv(self.txt_path +'/summary.txt', sep='\t')
        print(self.summary)
        self.create_folders()
        self.peptides = pd.read_csv(self.txt_path +'/peptides.txt', sep='\t')
        samples = []
        for column in self.peptides.columns:
            if column.startswith('Intensity '):
                name = column.split('Intensity ')[1]
                samples.append(name)

        if not os.path.exists(design):
            self.design = pd.DataFrame()
            self.design['sample'] = pd.Series(samples)
            self.design['rename'] = pd.Series(samples)
            self.design['LEVEL_1'] = 'Group_1'
            self.design.to_csv(design, index=False)
            print("Please edit {}/design.csv, then restart the analysis.".format(self.outdir))
            return
        else:
            self.design=pd.read_csv(design)
            print(self.design)

        rename_columns = {}
        exclude_columns = []
        exclude_samples = []

        group_levels = self.config['group_levels']
        print(design)
        for row in self.design.iterrows():
            rename = row[1]['rename']
            sample = row[1]['sample']
            exclude = str(row[1]['exclude']).strip()
            row_exclude = []
            frm = 'Intensity {}'.format(sample)
            assert frm in self.peptides.columns # Check that the column exists
            to = 'Intensity {}'.format(rename)
            rename_columns[frm] = to
            row_exclude.append(to)
            
            frm = '{} {}'.format(self.protein_quantification, sample)
            to = '{} {}'.format(self.protein_quantification, rename)
            rename_columns[frm] = to
            row_exclude.append(to)
            
            if exclude =='+':
                exclude_columns += row_exclude
                exclude_samples.append(rename)

        for level in group_levels:
            if not level in self.design.columns:
                print('Group level {} is not defined in design {}'.format(level, design))
                return
        
        self.config = self.update_config(self.config, self.design, exclude_samples)
        self.rename_columns = rename_columns
        self.peptides['Identifier'] = self.peptides['Sequence']
        self.proteingroups = self.create_protein_group_identifier(pd.read_csv(self.txt_path +'/proteinGroups.txt', sep='\t'))

        self.proteingroups['Leading Protein'] = self.proteingroups['Protein IDs'].apply(parse_ids).apply(lambda x : x.split(';')[0])
        self.proteingroups['Leading Species'] = self.proteingroups['Protein IDs'].apply(parse_groups).apply(lambda x : x.split(';')[0])
        self.reference_fasta = list(SeqIO.parse(self.config['reference_fasta'], 'fasta'))
        self.search_fasta = list(SeqIO.parse(self.config['search_fasta'],'fasta'))
        self.proteingroups = self.leading_protein_gene(self.proteingroups.copy(), self.search_fasta, self.reference_fasta)
        assert 'Leading.gene' in self.proteingroups.columns
        print(self.proteingroups['Leading.gene'])
        mgfile = self.outdir +'/proteingroups.mygene.csv'
        if not os.path.exists(mgfile):
            self.proteingroups = self.leading_protein_mygene(self.proteingroups)
            self.proteingroups.to_csv(mgfile)
        else:
            self.proteingroups = pd.read_csv(mgfile)
        
        assert len(self.proteingroups['Identifier'].tolist()) == len(set(self.proteingroups['Identifier'].tolist()))
        
        kofile =self.outdir + '/proteingroups.ko.csv'
        if not os.path.exists( kofile):
            self.proteingroups = self.leading_protein_ko(self.proteingroups)
            self.proteingroups.to_csv(kofile)
        else:
            self.proteingroups = pd.read_csv(kofile)
        
        print(list(rename_columns.keys()))
        for sample in samples:
            pepcol = 'Intensity {}'.format(sample)
            protcol = '{} {}'.format(self.protein_quantification, sample)
            if not pepcol in rename_columns:
                del self.proteingroups[protcol]
                del self.peptides[pepcol]
            else:
                r = self.proteingroups[protcol]
                _ = len([i for i in r if not i == 0]) /  len(r)
                print(protcol, _)
                #assert _ * 100 > 20
                #r = self.peptides[pepcol]
                #_ = len([i for i in r if not i == 0]) /  len(r)
                #print(pepcol, _)
                #assert _ * 100 > 20
        
        self.peptides.rename(columns=rename_columns, inplace=True)            
        self.proteingroups.rename(columns=rename_columns, inplace=True)
        self.proteingroups = self.host_proteins(self.proteingroups, self.reference_fasta)
        assert len(self.proteingroups['Identifier'].tolist()) == len(set(self.proteingroups['Identifier'].tolist()))
        self.msms = pd.read_csv(self.txt_path +'/msms.txt', sep='\t')
        self.target_proteingroups = self.exclude_reverse(self.proteingroups)
        self.reverse_msms = self.get_reverse(self.msms)
        self.reverse_peptides = self.get_reverse(self.peptides)
        self.reverse_proteingroups = self.get_reverse(self.proteingroups) 
        self.contaminant_peptides = self.get_contaminants(self.peptides)
        self.contaminant_peptides_list = self.contaminant_peptides['Sequence'].tolist()
        self.contaminant_proteingroups = self.get_contaminants(self.proteingroups)
        self.contaminant_msms = self.msms[self.msms['Sequence'].isin(self.contaminant_peptides_list)]
        self.target_msms = self.exclude_reverse(self.msms)
        self.target_peptides = self.exclude_reverse(self.peptides)
        self.target_proteingroups = self.exclude_reverse(self.proteingroups)
        
        ######
        # QC #
        ######
        self.qc(self.config,self.design,self.qc_dir,self.summary,self.target_peptides,self.target_proteingroups)
        
        if exclude_contaminants == True:
            self.target_peptides = self.exclude_contaminants(self.target_peptides)
            self.target_proteingroups = self.exclude_contaminants(self.target_proteingroups)
        
        self.protein_id_lists(self.target_proteingroups, self.protein_dir +'/protein_ids.txt')
        self.target_peptides_list = self.target_peptides['Sequence'].tolist()
        self.target_msms = self.target_msms[self.target_msms['Sequence'].isin(self.target_peptides_list)]
        self.get_reference_peptides()
        self.reference_peptides=self.target_peptides[self.target_peptides['Sequence'].isin(self.reference_peptides_list)]
        self.non_reference_peptides=self.target_peptides[self.target_peptides['Sequence'].isin(self.non_reference_peptides_list)]
        self.reference_msms=self.target_msms[self.target_msms['Sequence'].isin(self.reference_peptides_list)]
        self.non_reference_msms=self.target_msms[self.target_msms['Sequence'].isin(self.non_reference_peptides_list)]
        self.target_msms_pep = self.target_msms['PEP'].tolist()
        self.target_peptides_pep = self.target_peptides['PEP'].tolist()
        self.reverse_msms_pep = self.reverse_msms['PEP'].tolist()
        self.reverse_peptides_pep = self.reverse_peptides['PEP'].tolist()
        self.reference_msms_pep = self.reference_msms['PEP'].tolist()
        self.non_reference_msms_pep = self.non_reference_msms['PEP'].tolist()
        self.reference_peptides_pep = self.reference_peptides['PEP'].tolist()
        self.non_reference_peptides_pep = self.non_reference_peptides['PEP'].tolist()
        self.contaminant_msms_pep = self.contaminant_msms['PEP'].tolist()
        
        # Get gene names from FASTA
        
        #~######
        # get PEP score statistics
        self.contaminant_pep_median=np.median(self.contaminant_msms_pep)
        self.contaminant_pep_stdev=np.std(self.contaminant_msms_pep)
        self.reverse_pep_median=np.median(self.reverse_msms_pep)
        self.reverse_pep_stdev=np.std(self.reverse_msms_pep)
        self.reference_pep_median=np.median(self.reference_msms_pep)
        self.reference_pep_stdev=np.std(self.reference_msms_pep)
        self.non_reference_pep_median = np.median(self.non_reference_msms_pep)
        self.non_reference_pep_stdev = np.std(self.non_reference_msms_pep)
        self.target_pep_median = np.median(self.target_msms_pep)
        self.target_pep_stdev = np.std(self.target_msms_pep)
        
        
        # Exclude samples according to design.csv
        self.peptide_txt = self.peptide_dir +'target_peptides.txt'
        print('Excluding :')
        print(exclude_columns)
        self.target_peptides = self.target_peptides[[i for i in self.target_peptides.columns if not  i in exclude_columns]]
        self.target_peptides.to_csv(self.peptide_txt, sep='\t')
        self.protein_txt = self.protein_dir + 'target_proteins.txt'

        self.target_proteingroups=self.target_proteingroups[[i for i in self.target_proteingroups.columns if not  i in exclude_columns]]
        self.target_proteingroups.to_csv(self.protein_txt, sep='\t')


        # create summary
        self.create_summary()

        # plot pep scores
        self.pepfig()
    
        # kw analysis
        self.pep_kw()

        # peptide lists
        self.save_peptide_lists()
       
        # fasta export
        fasta_file = self.config['search_fasta']
        prot_ids = self.protein_dir +'/protein_ids.txt'
        outfile = self.fasta_dir +'/leading_proteins.fasta'
        self.export_pg_fasta(fasta_file, prot_ids, outfile)
        
        # interproscan
        if not os.path.exists(self.fasta_dir +'/leading_proteins.fasta.tsv' ):
            infile = self.fasta_dir +'/leading_proteins.fasta'
            outpath  = self.fasta_dir
            self.ips_fasta(infile, outpath) 
    
        # gene sets
        infile = self.fasta_dir +'/leading_proteins.fasta.tsv'
        self.ips_genesets(infile, self.target_proteingroups, self.gsea_dir)
              
        # Create peptide paraameters
        for level in self.config['group_levels']:
            outfile=self.diff_dir + '/peptide_experimental_design_{}.R'.format(level)
            self.create_R_parameters(  config = self.config, 
                                       outfile=outfile, 
                                       quant='Intensity',
                                       group_level=level)
        
        # Normalize peptides
        #norm_method = self.config['normalize']
        #impute_method = self.config['impute']
        outdir = self.diff_dir + 'peptide_normalization'
        infile = self.peptide_txt
        norm_peps=self.diff_dir + '/peptide_normalization/msnbase/normalized.csv'
        if not os.path.exists(norm_peps):
            self.normalize('Intensity.', infile, outdir, self.peptide_normalisation, self.peptide_imputation)
        self.normalized_target_peptides = pd.read_csv(norm_peps)
        self.protein_quant_parameter = '.'.join(self.protein_quantification.split())
        # Create protein paraameters
        for level in self.config['group_levels']:
            outfile=self.diff_dir + '/protein_experimental_design_{}.R'.format(level)
            self.create_R_parameters(  config = self.config, 
                                       outfile=outfile, 
                                       quant=self.protein_quant_parameter,
                                       group_level=level)
        # Normalize proteins
        outdir = self.diff_dir + 'protein_normalization'
        infile = self.gsea_dir +'/ipr_target_proteins.txt'
        norm_prots = self.diff_dir + '/protein_normalization/msnbase/normalized.csv'
        
        if not os.path.exists(norm_prots):
            self.normalize('{}.'.format(self.protein_quant_parameter), infile, outdir, self.protein_normalisation, self.protein_imputation)
        
        self.normalized_target_proteins = pd.read_csv(norm_prots) 
        
        # unipept anslysis
        self.unipept()
        
        # diff analusis
        self.diff_analysis()
        
        # GSEA
        outpath=self.gsea_dir
        for level in self.config['group_levels']:
            gsea_dir = outpath + '/' + level
            if not os.path.exists(gsea_dir):
                os.mkdir(gsea_dir)
            table=self.diff_dir+'/protein_normalization/msnbase/normalized.csv'
            design = self.diff_dir+'/protein_experimental_design_{}.R'.format(level)
            genecol='Mygene.entrez'
            kocol='Leading.Protein.Kegg.Orthology.ID'
            self.ips_gsea(outpath, gsea_dir, design, table, genecol=genecol , kocol=kocol)
            self.summarize_gsea(gsea_dir, gsea_dir + '/summary.txt', self.config, level)
   
    def qc(self, config, design,  outpath, summary, target_peptides, target_proteins):
        print("Starting QC")
        if not os.path.exists(outpath):
            os.mkdir(outpath )
        
        # Protein abundance analysis
        target_proteins = target_proteins.sort_values('Intensity', ascending=False)
        #print(target_proteins.columns.tolist())
        #contaminants = target_proteins[target_proteins['Potential contaminant'] == '+'].head(5)
        target_proteins['Source'] = target_proteins['Protein IDs']
        target_proteins.loc[target_proteins['Potential contaminant'] != '+','Source' ] = 'Non-contaminant'
        target_proteins['Source'] = target_proteins['Source'].apply(lambda x : x.split(';')[0])
        agg_cols = {}
        for col in target_proteins.columns.tolist():
            if col.startswith('Intensity'):
                agg_cols[col] = sum
        agg_cols['Potential contaminant'] = 'first'
        conts = target_proteins.groupby('Source').agg(agg_cols)
        conts = conts.sort_values(['Intensity'], ascending=False)
        print(conts)
        conts = conts.head(5)
        intensity_cols = [col for col in conts.columns if col.startswith('Intensity ') ]
        conts = conts[intensity_cols]
        conts = conts.transpose()
        ax = conts.plot.bar(rot=90, stacked=False, figsize=(10,10))
        ax.set_title('Summed Intensity of all non-contaminant proteins\n compared to top contaminants')
        fig = ax.get_figure()
        fig.savefig( outpath + '/proteingroup_identifications.png')
        fig.clf()
        plt.close()
        #other = other[ibaq_cols].sum(axis=0)
        #print(other)


        # ms vs ms assigned
        summary = summary[summary['Experiment'].notnull()]
        columns = ['MS/MS Submitted','MS/MS Identified']
        ax = summary.set_index('Experiment')[columns].plot.bar(rot=90)
        fig = ax.get_figure()
        fig.savefig( outpath + '/msms_identifications.png')
        fig.clf()
        plt.close()

        ######### 
        # SPLOM #
        #########

        if not os.path.exists(outpath + '/splom'):
            os.mkdir(outpath + '/splom')
        med_intensity = target_proteins[intensity_cols].replace(0,np.nan) #.median(axis=1)
        med_intensity = med_intensity.median(skipna=True,  axis=1)
        for col in intensity_cols:
            newdf = pd.DataFrame()
            newdf['Median'] = med_intensity
            newdf[col] = target_proteins[col]
            newdf = newdf.replace(0, np.nan)
            plt.style.use('ggplot')
            ax = newdf.plot.scatter(x='Median', y=col)
            ax.set_title('Intensity of {}\n compared to median across all samples'.format(col.split()[1]))
            fig = ax.get_figure()
            colname = col.split('Intensity ')[1]
            fig.savefig( outpath + '/splom/{}.png'.format(colname))
            fig.clf()
            plt.close()

        #############
        # Box Plots #
        #############
        title='Number of identified peptides per sample'
        if not os.path.exists(outpath + '/bp'):
            os.mkdir(outpath + '/bp/')
        all_ids = summary['Peptide Sequences Identified'].tolist()
        data = [all_ids]
        names = ['All' ]
        bp(data, names, outpath + '/bp/peptide_identifications_all.png', title,  True)
        for level in config['group_levels']:
            data_dict = defaultdict(list)
            data  = []
            names = []
            mapping=design.set_index('sample')[level].to_dict()
            vals=summary.set_index('Experiment')['Peptide Sequences Identified'].to_dict()
            for key in vals:
                if key in mapping:
                    group = mapping[key]
                    if not str(group) == 'nan':
                        val = vals[key]
                        data_dict[group].append(val)
            for key in data_dict:
                names.append(key)
                data.append(data_dict[key] )
                bp(data, names, outpath + '/bp/peptide_identifications_{}.png'.format(level), title, True)
        
        ##################################
        # % Identified spectra by groups #
        ##################################
        all_ids = summary['MS/MS Identified [%]'].tolist()
        title='Percentage of identified MS/MS per sample'
        data = [all_ids]
        names = ['All' ]
        bp(data, names, outpath + '/bp/msms_percent_identifications_all.png', title, True)
        for level in config['group_levels']:
            data_dict = defaultdict(list)
            data  = []
            names = []
            mapping=design.set_index('sample')[level].to_dict()
            vals=summary.set_index('Experiment')['MS/MS Identified [%]'].to_dict()
            for key in vals:
                if key in mapping:
                    group = mapping[key]
                    if not str(group) == 'nan':
                        val = vals[key]
                        data_dict[group].append(val)
            for key in data_dict:
                names.append(key)
                data.append(data_dict[key] )
                bp(data, names, outpath + '/bp/msms_percent_identifications_{}.png'.format(level), title, True)
        return  

    def update_config(self, config, design, exclude_columns):
        config['samples'] = {}
        for row in design.iterrows():
            rename = row[1]['rename']
            if not rename in exclude_columns:
                print(rename)
                sample = row[1]['sample']
                config['samples'][rename] = {}
                for level in config['group_levels']:
                    assert level in row[1].index
                    config['samples'][rename][level] = row[1][level]
        return config

    def create_protein_group_identifier(self, proteingroups):
        proteingroups['Identifier'] = ' (Protein group ' + proteingroups['id'].apply(str)+')'
        proteingroups['Identifier'] = proteingroups['Protein IDs'].apply(parse_ids).apply(parse_protein_ids) + proteingroups['Identifier']
        return proteingroups

    def create_summary(self):
        w = open(self.outdir +'/summary.txt','w')
        w.write('Target proteingroups: {}\n'.format(len(self.target_proteingroups)))
        w.write('Target peptides: {}\n'.format(len(self.target_peptides)))
        w.write('Target msms: {}\n'.format(len(self.target_msms)))
        d = str(self.summary[-1:].stack())
        w.write(d)
        w.close()
    def host_proteins(self, data, host_fasta):
        ids = []
        for p in host_fasta:
            id = p.id.split('|')[1]
            ids.append(id)
        data['Component'] = np.where(data['Leading Protein'].isin(ids), 'host', 'other')
        return data

    def create_folders(self):
        self.unipept_dir = self.outdir + '/unipept/'
        self.fasta_dir = self.outdir + '/fasta/'
        self.peptide_dir = self.outdir + '/peptides/'
        self.protein_dir = self.outdir + '/proteins/'
        self.pep_dir = self.outdir + '/pep/'
        self.diff_dir = self.outdir + '/diff/'
        self.gsea_dir = self.outdir + '/gsea/'
        self.qc_dir = self.outdir + '/qc/'
        
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir) 
            
        if not os.path.exists(self.pep_dir):
            os.mkdir(self.pep_dir)
        if not os.path.exists(self.fasta_dir):
            os.mkdir(self.fasta_dir)
        if not os.path.exists(self.unipept_dir):
            os.mkdir(self.unipept_dir)
        if not os.path.exists(self.peptide_dir):
            os.mkdir(self.peptide_dir)
        if not os.path.exists(self.protein_dir):
            os.mkdir(self.protein_dir)
        if not os.path.exists(self.diff_dir):
            os.mkdir(self.diff_dir)
        if not os.path.exists(self.qc_dir):
            os.mkdir(self.qc_dir)
        if not os.path.exists(self.gsea_dir):
            os.mkdir(self.gsea_dir)
   
    def get_reference_peptides(self):
        #peptides = self.exclude_contaminants(self.target_peptides)
        seqs = self.target_peptides['Sequence'].tolist()
        refstr = '-'.join([str(s.seq) for s in self.reference_fasta])
        refseqs = []
        nonrefseqs = []
        for i in seqs:
            if i in refstr:
                refseqs.append(i)
            else:
                nonrefseqs.append(i)
        self.reference_peptides_list = refseqs
        self.non_reference_peptides_list = nonrefseqs
        self.target_peptides['ReferenceProteomePeptide']= self.target_peptides['Sequence'].apply(lambda x : str(x in refseqs))
    
    def exclude_contaminants(self, df):
        df = df.copy()
        df = df[df['Potential contaminant'] != '+' ]
        return df
    
    def exclude_reverse(self, df):
        df = df.copy()
        df = df[df['Reverse'] != '+' ]
        return df
    
    def get_contaminants(self, df):
        df =df.copy()
        df = df[df['Potential contaminant'] == '+' ]
        return df
    
    def get_reverse(self, df):
        df = df.copy()
        df = df[df['Reverse'] == '+' ]
        return df

    def save_peptide_lists(self):
        w = open(self.peptide_dir +'/reference_peptides_list.txt', 'w')
        w.write('\n'.join(self.reference_peptides['Sequence'].tolist()))
        w.close()
        w = open(self.peptide_dir +'/non_reference_peptides_list.txt', 'w')
        w.write('\n'.join(self.non_reference_peptides['Sequence'].tolist()))
        w.close()
        w = open(self.peptide_dir +'/reverse_peptides_list.txt', 'w')
        w.write('\n'.join(self.reverse_peptides['Sequence'].tolist()))
        w.close()
        w = open(self.peptide_dir +'/contaminant_peptides_list.txt', 'w')
        w.write('\n'.join(self.contaminant_peptides['Sequence'].tolist()))
        w.close()
        w = open(self.peptide_dir +'/target_peptides_list.txt', 'w')
        w.write('\n'.join(self.target_peptides['Sequence']))
        w.close()

    def pepfig(self):
        data = [self.target_msms_pep, self.contaminant_msms_pep, self.reference_msms_pep, self.non_reference_msms_pep, self.reverse_msms_pep]
        names = ['Target', 'Contaminant', "Reference", "Non-reference", "Reverse"]
        output = self.pep_dir

        fig = plt.figure()
        ax = fig.add_subplot(111)
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
        ax.set_title('PSM PEP Score distributions')
        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        # Save the figure
        fig.savefig(output + '/psm_pep_scores.png', bbox_inches='tight')
        fig.clf()

    def pep_kw(self):
        data = [self.target_msms_pep, self.contaminant_msms_pep, self.reference_msms_pep, self.non_reference_msms_pep, self.reverse_msms_pep]
        names = ['Target', 'Contaminant', "Reference", "Non-reference", "Reverse"]
        outpath = self.pep_dir 
        if not os.path.exists(outpath):
            list_kw_dunn(names, data, 'PEP', 'Category', outpath)
    

    def unipept(self):
        if not (os.path.exists(self.unipept_dir +'unipept.sh')):
            w = open(self.unipept_dir + 'unipept.sh','w')
            w.write('\ncat ../peptides/target_peptides_list.txt | prot2pept | peptfilter | tr I L | sort -u | unipept pept2lca -e -a > pept2lca.txt')
            #w.write('\ncat ../peptides/target_peptides_list.txt | prot2pept | peptfilter | tr I L | sort -u | unipept pept2taxa -e -a > pept2taxa.txt')
            w.close()
            cmd = 'cd {} && chmod 700 unipept.sh && ./unipept.sh > unipept.log 2>&1 && exit'.format(self.unipept_dir)
            process = subprocess.Popen(cmd, shell=True)
            process.wait()
            assert process.returncode == 0

        agg_cols = {'taxon_name':'first', 
                    'taxon_id':'first', 
                    'genus_name':'first', 
                    'genus_id':'first',
                    'species_id':'first',
                    'species_name':'first',
                    'family_id':'first',
                    'family_name':'first',
                    'phylum_name' :'first',
                    'phylum_id':'first',
                    'PeptideCount':sum,
                    'MS.MS.Count': sum }
        
        clean_cols = ['Identifier',
                      'PeptideCount',
                      'MS.MS.Count',
                      'taxon_name', 
                      'taxon_id',
                      'species_name', 
                      'species_id',
                      'genus_name', 
                      'genus_id',
                      'family_name',
                      'family_id',
                      'phylum_name',
                      'phylum_id']

        for col in self.normalized_target_peptides.columns.tolist():
            if col.startswith('Experiment.'):
                agg_cols[col] = np.sum
                clean_cols.append(col)
            elif col.startswith('Intensity.'):
                agg_cols[col] = np.sum
                clean_cols.append(col)
        lca_level_cutoff = 1
        
        # pept2lca
        self.pept2lca = pd.read_csv(self.unipept_dir + '/pept2lca.txt')
        lca_peptides = pd.merge(self.normalized_target_peptides,self.pept2lca,how='left',left_on='Sequence',right_on='peptide')
        lca_peptides['PeptideCount'] = 1
        lca_peptides['taxon_name'].replace(np.nan, "unassigned", inplace=True)

        lca_peptides.to_csv(self.unipept_dir + '/pept2lca_peptides.csv')
        pept2lca_taxon_sc = lca_peptides.groupby(lca_peptides.taxon_name).agg(agg_cols) 
        pept2lca_taxon_sc = pept2lca_taxon_sc[pept2lca_taxon_sc['PeptideCount'] >= lca_level_cutoff ]
        lca_taxon_ids = pept2lca_taxon_sc['taxon_id'].tolist()
        pept2lca_taxon_sc = pept2lca_taxon_sc.sort_values('MS.MS.Count', ascending=False)
        pept2lca_taxon_sc['Identifier'] = pept2lca_taxon_sc['taxon_name']
        newcols = [i for i in clean_cols if i in pept2lca_taxon_sc.columns.tolist()]
        pept2lca_taxon_sc = pept2lca_taxon_sc[newcols]
        pept2lca_taxon_sc.to_csv(self.unipept_dir + '/pept2lca_taxon_sc.csv')

        # pept2lca species
        pept2lca_species_sc = lca_peptides.groupby(lca_peptides.species_name).agg(agg_cols) 
        pept2lca_species_sc = pept2lca_species_sc[pept2lca_species_sc['PeptideCount'] >=  lca_level_cutoff ]
        lca_species_ids = pept2lca_species_sc['species_id'].tolist()
        pept2lca_species_sc['Identifier'] = pept2lca_species_sc['species_name']
        del pept2lca_species_sc['taxon_name']
        del pept2lca_species_sc['taxon_id']
        newcols = [i for i in clean_cols if i in pept2lca_species_sc.columns.tolist()]
        pept2lca_species_sc = pept2lca_species_sc[newcols]
        pept2lca_species_sc.to_csv(self.unipept_dir + '/pept2lca_species_sc.csv')
        
        # pept2lca genus
        pept2lca_genus_sc = lca_peptides.groupby(lca_peptides.genus_name).agg(agg_cols) 
        pept2lca_genus_sc = pept2lca_genus_sc[pept2lca_genus_sc['PeptideCount'] >= lca_level_cutoff ]
        lca_genus_ids = pept2lca_genus_sc['genus_id'].tolist()
        pept2lca_genus_sc['Identifier'] = pept2lca_genus_sc['genus_name']
        del pept2lca_genus_sc['taxon_name']
        del pept2lca_genus_sc['taxon_id']
        del pept2lca_genus_sc['species_name']
        del pept2lca_genus_sc['species_id']
        
        newcols = [i for i in clean_cols if i in pept2lca_genus_sc.columns.tolist()]
        pept2lca_genus_sc = pept2lca_genus_sc[newcols]
        pept2lca_genus_sc.to_csv(self.unipept_dir + '/pept2lca_genus_sc.csv')
        
        # pept2lca family
        pept2lca_family_sc = lca_peptides.groupby(lca_peptides.family_name).agg(agg_cols) 
        pept2lca_family_sc = pept2lca_family_sc[pept2lca_family_sc['PeptideCount'] >= lca_level_cutoff ]
        lca_family_ids = pept2lca_family_sc['family_id'].tolist()
        pept2lca_family_sc['Identifier'] = pept2lca_family_sc['family_name']
        del pept2lca_family_sc['taxon_name']
        del pept2lca_family_sc['taxon_id']
        del pept2lca_family_sc['species_name']
        del pept2lca_family_sc['species_id']
        del pept2lca_family_sc['genus_name']
        del pept2lca_family_sc['genus_id']
        newcols = [i for i in clean_cols if i in pept2lca_family_sc.columns.tolist()]
        pept2lca_famliy_sc = pept2lca_family_sc[newcols]
        pept2lca_family_sc.to_csv(self.unipept_dir + '/pept2lca_family_sc.csv')  
        
        # pept2lca phylum
        pept2lca_phylum_sc = lca_peptides.groupby(lca_peptides.phylum_name).agg(agg_cols) 
        pept2lca_phylum_sc = pept2lca_phylum_sc[pept2lca_phylum_sc['PeptideCount'] >= lca_level_cutoff ]
        pept2lca_phylum_sc['Identifier'] = pept2lca_phylum_sc.index
        del pept2lca_phylum_sc['taxon_name']
        del pept2lca_phylum_sc['taxon_id']
        del pept2lca_phylum_sc['species_name']
        del pept2lca_phylum_sc['species_id']
        del pept2lca_phylum_sc['genus_name']
        del pept2lca_phylum_sc['genus_id']
        del pept2lca_phylum_sc['family_name']
        del pept2lca_phylum_sc['family_id']
        newcols = [i for i in clean_cols if i in pept2lca_phylum_sc.columns.tolist()]
        pept2lca_phylum_sc = pept2lca_phylum_sc[newcols]
        pept2lca_phylum_sc.to_csv(self.unipept_dir + '/pept2lca_phylum_sc.csv')
        
        ## pept2taxa 
        #self.pept2taxa = pd.read_csv(self.unipept_dir + '/pept2taxa.txt')
        #self.pept2taxa['PeptideCount'] = 1
        #pept2taxa = self.pept2taxa[self.pept2taxa['taxon_id'].isin(lca_taxon_ids)]
        #taxa_peptides = pd.merge(self.normalized_target_peptides, pept2taxa, how='inner', left_on='Sequence', right_on='peptide')
        #pept2taxa_sc = taxa_peptides.groupby(taxa_peptides.taxon_name).agg(agg_cols)
        #pept2taxa_sc = pept2taxa_sc.sort_values('MS.MS.Count', ascending=False)
        #pept2taxa_sc['Identifier'] = pept2taxa_sc['taxon_name']
        #newcols = [i for i in clean_cols if i in pept2taxa_sc.columns.tolist()]
        #pept2taxa_sc = pept2taxa_sc[newcols]
        #pept2taxa_sc.to_csv(self.unipept_dir + '/pept2taxa_taxon_sc.csv')
        #
        ## pept2taxa species
        #pept2taxa = self.pept2taxa[self.pept2taxa['species_id'].isin(lca_species_ids)]
        ##self.pept2taxa = self.pept2taxa[['peptide','taxon_id','taxon_name', 'genus_name', 'genus_id']]
        #taxa_peptides = pd.merge(self.normalized_target_peptides, pept2taxa, how='inner', left_on='Sequence', right_on='peptide')
        #pept2taxa_sc = taxa_peptides.groupby(taxa_peptides.species_name).agg(agg_cols)
        ##pept2ltaxa_sc = pept2taxa_sc[['taxon_name','taxon_id','genus_name','genus_id','Intensity','MS/MS Count']]
        #pept2taxa_sc = pept2taxa_sc.sort_values('MS.MS.Count', ascending=False)
        #pept2taxa_sc['Identifier'] = pept2taxa_sc['species_name']
        #del pept2taxa_sc['taxon_name']
        #del pept2taxa_sc['taxon_id']
        #newcols = [i for i in clean_cols if i in pept2taxa_sc.columns.tolist()]
        #pept2taxa_sc = pept2taxa_sc[newcols]
        #pept2taxa_sc.to_csv(self.unipept_dir + '/pept2taxa_species_sc.csv')
        #
        ## pept2taxa genus sc
        #pept2taxa = self.pept2taxa[self.pept2taxa['genus_id'].isin(lca_genus_ids)]
        #taxa_peptides = pd.merge(self.normalized_target_peptides, pept2taxa, how='inner', left_on='Sequence', right_on='peptide')
        #pept2taxa_genus_sc = taxa_peptides.groupby(taxa_peptides.genus_name).agg(agg_cols)
        ##pept2taxa_genus_sc = pept2taxa_genus_sc[['taxon_name','taxon_id', 'genus_name','genus_id','Intensity','MS/MS Count']]
        #del pept2taxa_genus_sc['taxon_name']
        #del pept2taxa_genus_sc['taxon_id']
        #del pept2taxa_genus_sc['species_name']
        #del pept2taxa_genus_sc['species_id']
        #pept2taxa_genus_sc = pept2taxa_genus_sc.sort_values('MS.MS.Count', ascending=False)
        #pept2taxa_genus_sc['Identifier'] = pept2taxa_genus_sc['genus_name']
        #newcols = [i for i in clean_cols if i in pept2taxa_genus_sc.columns.tolist()]
        #pept2taxa_genus_sc = pept2taxa_genus_sc[newcols]
        #pept2taxa_genus_sc.to_csv(self.unipept_dir + '/pept2taxa_genus_sc.csv')

        ## pept2taxa family sc
        #pept2taxa = self.pept2taxa[self.pept2taxa['family_id'].isin(lca_family_ids)]
        #taxa_peptides = pd.merge(self.normalized_target_peptides, pept2taxa, how='inner', left_on='Sequence', right_on='peptide')
        #pept2taxa_sc = taxa_peptides.groupby(taxa_peptides.family_name).agg(agg_cols)
        ##pept2taxa_genus_sc = pept2taxa_genus_sc[['taxon_name','taxon_id', 'genus_name','genus_id','Intensity','MS/MS Count']]
        #del pept2taxa_sc['taxon_name']
        #del pept2taxa_sc['taxon_id']
        #del pept2taxa_sc['species_name']
        #del pept2taxa_sc['species_id']
        #del pept2taxa_sc['genus_name']
        #del pept2taxa_sc['genus_id']
        #pept2taxa_sc = pept2taxa_sc.sort_values('MS.MS.Count', ascending=False)
        #pept2taxa_sc['Identifier'] = pept2taxa_sc['family_name']
        #newcols = [i for i in clean_cols if i in pept2taxa_sc.columns.tolist()]
        #pept2taxa_sc = pept2taxa_sc[newcols]
        #pept2taxa_sc.to_csv(self.unipept_dir + '/pept2taxa_family_sc.csv')
      
    def ips_gsea(self, indir, outpath, design, table, genecol, kocol, keggid='hsa', pval=0.05):
        print(indir, outpath, design, table, genecol, kocol, keggid, pval)
        indir = '/'.join(indir.split('//'))
        outpath = '/'.join(outpath.split('//'))
        
        if not os.path.exists(outpath):
            os.mkdir(outpath)

        table = pd.read_csv(table)
        table.to_csv(indir + '/combined.csv')    
        group_table = indir + '/combined.csv'
        cols = table.columns.tolist()
        
        #print(table['Row names'].tolist())
        cmd = 'gage.R --indir {} --outdir {} --keggid {} --design {} --table {} --genecol {} --kocol {} --pval {}'.format(indir, outpath, keggid, design, group_table, genecol, kocol,  pval)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        assert process.returncode == 0
    
        for name, group in table.groupby('Component'):
            assert len(group) == len(group['Identifier'].dropna())
            #print(group['Row.name'].tolist())
            group_table = indir + '/{}.csv'.format(name)
            group.to_csv(indir + '/{}.csv'.format(name))    
            cmd = 'gage.R --indir {} --outdir {} --keggid {} --design {} --table {} --genecol {} --kocol {} --pval {}'.format(indir, outpath, keggid, design, group_table, genecol, kocol,  pval)
            process = subprocess.Popen(cmd, shell=True)
            process.wait()
            assert process.returncode == 0
    
    def ips_genesets(self, ipr, proteins, outpath, keggid='ko', id_col="Leading Protein"):
        cols = ['ProteinAccession', 
                'MD5', 
                'Length', 
                'Analysis', 
                'SignatureAccession',
                'SignatureDescription',
                'Start',
                'Stop',
                'Score',
                'Status',
                'Date',
                'InterProAnnotationsAccession',
                'InterProAnnotationsDescription',
                'GoAnnotations',
                'PathwaysAnnotations' ]
        data = pd.read_csv(ipr, sep='\t', names = cols,  engine='python')
        lst_col = 'ProteinAccession' 
        x = data.assign(**{lst_col:data[lst_col].str.split('|')})
        data = pd.DataFrame({col:np.repeat(x[col].values, x[lst_col].str.len()) for col in x.columns.difference([lst_col])}).assign(**{lst_col:np.concatenate(x[lst_col].values)})[x.columns.tolist()]
        
        
        # GO Terms
        id2go= defaultdict(set)
        gos = set()
        def go(df):
            go = df['GoAnnotations']
            id = df['ProteinAccession']
            try:
                go=go.split('|')
                id2go[id].update(go)
                for goterm in go:
                    gos.add(goterm)
            except:
                pass
        data.apply(go, axis=1)
        with open( outpath +'/accession2go.p', 'wb') as f:
            pickle.dump( id2go, f)
        go_df = pd.DataFrame()
        go_vals = list(gos)
        go_df['GO_ID'] = pd.Series(go_vals)
        go_df.to_csv(outpath +'/go_terms.csv')

        def gointersect(val):
            vals = val.split(';')
            setlist = []
            for val in vals:
                vset = id2go[val]
                setlist.append(vset)
            union = set.union(*setlist)
            if len(union) > 0:
                return ';'.join(union)
        proteins['_go.term.union']   = proteins[id_col].apply(parse_ids).apply(gointersect)
        
        go2pg = defaultdict(set)
        def go2gene(df):
            go_terms = df['_go.term.union']
            pg = df['Identifier']
            try:
                go_terms = go_terms.split(';')
                for _ in go_terms:
                    go2pg[_].add(pg)
            except:
                pass
        proteins.apply(go2gene,axis=1)

        go_df = pd.DataFrame()
        go_df['GO_ID'] = pd.Series(list(go2pg.keys()))
        go_df['GENES'] = pd.Series(list(go2pg.values())).apply( lambda x  : '|'.join(x))
        go_df.to_csv(outpath +'/go2proteingroups.csv')


        # KEGG
        id2kegg= defaultdict(set)
        keggs = set()
        def kegg(df):
            id = df['ProteinAccession']
            pathways = df['PathwaysAnnotations']
            try:
                kegg = [i for i in pathways.split('|') if i.startswith('KEGG: ')]
                kegg = [i.split('KEGG: ')[1].split('+')[0] for i in kegg]
                id2kegg[id].update(kegg)
                for keggterm in kegg:
                    keggs.add(keggterm)
            except:
                pass
        data.apply(kegg, axis=1)
        with open( outpath +'/accession2kegg.p', 'wb') as f:
            pickle.dump( id2kegg, f)
        kegg_df = pd.DataFrame()
        kegg_vals = list(keggs)
        kegg_df['KEGG_ID'] = pd.Series(kegg_vals)
        kegg_df.to_csv(outpath +'/kegg_terms.csv')
        
        def keggintersect(val):
            vals = val.split(';')
            setlist = []
            for val in vals:
                vset = id2kegg[val]
                setlist.append(vset)
            union = set.union(*setlist)
            if len(union) > 0 :
                return ';'.join(union)
        proteins['_kegg.term.union']   = proteins[id_col].apply(parse_ids).apply(keggintersect)

        kegg2pg = defaultdict(set)
        
        def kegg2gene(df):
            kegg_terms = df['_kegg.term.union']
            pg = df['Identifier']
            try:
                kegg_terms = kegg_terms.split(';')
                for _ in kegg_terms:
                    kegg2pg[_].add(pg)
            except:
                pass
        proteins.apply(kegg2gene,axis=1)
        kegg_df = pd.DataFrame()
        kegg_df['KEGG_ID'] = pd.Series(list(kegg2pg.keys()))
        kegg_df['KEGG_ID'] = kegg_df['KEGG_ID'].apply(str)
        kegg_df['GENES'] = pd.Series(list(kegg2pg.values())).apply( lambda x  : '|'.join(x))
        kegg_df.to_csv(outpath +'/kegg2proteingroups.csv')


        # EC
        id2ec= defaultdict(set)
        ecs = set()
        def ec(df):
            id = df['ProteinAccession']
            pathways = df['PathwaysAnnotations']
            try:
                ec = [i for i in pathways.split('|') if i.startswith('KEGG: ')]
                ec = [i.split('KEGG: ')[1].split('+')[1] for i in ec]
                id2ec[id].update(ec)
                for ecterm in ec:
                    ecs.add(ecterm)
            except:
                pass
        data.apply(ec, axis=1)
        with open( outpath +'/accession2ec.p', 'wb') as f:
            pickle.dump( id2ec, f)
        ec_df = pd.DataFrame()
        ec_vals = list(ecs)
        ec_df['EC_ID'] = pd.Series(ec_vals)
        ec_df.to_csv(outpath +'/ec_terms.csv')
        
        def ecintersect(val):
            vals = val.split(';')
            setlist = []
            for val in vals:
                vset = id2ec[val]
                setlist.append(vset)
            union = set.union(*setlist)
            if len(union) > 0 :
                return ';'.join(union)
        proteins['_ec.term.union']   = proteins[id_col].apply(parse_ids).apply(ecintersect)

        ec2pg = defaultdict(set)
        
        def ec2gene(df):
            ec_terms = df['_ec.term.union']
            pg = df['Identifier']
            try:
                ec_terms = ec_terms.split(';')
                for _ in ec_terms:
                    ec2pg[_].add(pg)
            except:
                pass
        proteins.apply(ec2gene,axis=1)
        ec_df = pd.DataFrame()
        ec_df['EC_ID'] = pd.Series(list(ec2pg.keys()))
        ec_df['EC_ID'] = ec_df['EC_ID'].apply(str)
        ec_df['GENES'] = pd.Series(list(ec2pg.values())).apply( lambda x  : '|'.join(x))
        ec_df.to_csv(outpath +'/ec2proteingroups.csv')

        # REACTOME
        id2reactome= defaultdict(set)
        reactomes = set()
        def reactome(df):
            id = df['ProteinAccession']
            pathways = df['PathwaysAnnotations']
            try:
                reactome = [i for i in pathways.split('|') if i.startswith('Reactome: ')]
                reactome = [i.split('Reactome: ')[1] for i in reactome]
                id2reactome[id].update(reactome)
                for reactometerm in reactome:
                    reactomes.add(reactometerm)
            except:
                pass
        data.apply(reactome, axis=1)
        with open( outpath +'/accession2reactome.p', 'wb') as f:
            pickle.dump( id2reactome, f)
        r_df = pd.DataFrame()
        r_vals = list(reactomes)
        r_df['REACTOME_ID'] = pd.Series(r_vals)
        r_df.to_csv(outpath +'/reactome_terms.csv')
        
        def reactomeintersect(val):
            vals = val.split(';')
            setlist = []
            for val in vals:
                vset = id2reactome[val]
                setlist.append(vset)
            union = set.union(*setlist)
            if len(union) > 0 :
                return ';'.join(union)
        proteins['_reactome.term.union']   = proteins[id_col].apply(parse_ids).apply(reactomeintersect)
        
        reactome2pg = defaultdict(set)
        
        def reactome2gene(df):
            r_terms = df['_reactome.term.union']
            pg = df['Identifier']
            try:
                r_terms = r_terms.split(';')
                for _ in r_terms:
                    reactome2pg[_].add(pg)
            except:
                pass
        proteins.apply(reactome2gene,axis=1)
        
        reactome_df = pd.DataFrame()
        reactome_df['REACTOME_ID'] = pd.Series(list(reactome2pg.keys()))
        reactome_df['REACTOME_ID'] = reactome_df['REACTOME_ID'].apply(str)
        reactome_df['GENES'] = pd.Series(list(reactome2pg.values())).apply( lambda x  : '|'.join(x))
        reactome_df.to_csv(outpath +'/reactome2proteingroups.csv')

        # MetaCyc
        id2metacyc= defaultdict(set)
        metacycs = set()
        def metacyc(df):
            id = df['ProteinAccession']
            pathways = df['PathwaysAnnotations']
            try:
                metacyc = [i for i in pathways.split('|') if i.startswith('MetaCyc: ')]
                metacyc = [i.split('MetaCyc: ')[1] for i in metacyc]
                id2metacyc[id].update(metacyc)
                for metacycterm in metacyc:
                    metacycs.add(metacycterm)
            except:
                pass
        data.apply(metacyc, axis=1)
        with open( outpath +'/accession2metacyc.p', 'wb') as f:
            pickle.dump( id2metacyc, f)

        m_df = pd.DataFrame()
        m_vals = list(metacycs)
        m_df['METACYC_ID'] = pd.Series(m_vals)
        m_df.to_csv(outpath +'/metacyc_terms.csv')
        
        def metacycintersect(val):
            vals = val.split(';')
            setlist = []
            for val in vals:
                vset = id2metacyc[val]
                setlist.append(vset)
            union = set.union(*setlist)
            if len(union) > 0 :
                return ';'.join(union)
        proteins['_metacyc.term.union']   = proteins[id_col].apply(parse_ids).apply(metacycintersect)
        
        metacyc2pg = defaultdict(set)
        
        def metacyc2gene(df):
            m_terms = df['_metacyc.term.union']
            pg = df['Identifier']
            try:
                m_terms = m_terms.split(';')
                for _ in m_terms:
                    metacyc2pg[_].add(pg)
            except:
                pass
        proteins.apply(metacyc2gene,axis=1)
        
        metacyc_df = pd.DataFrame()
        metacyc_df['METACYC_ID'] = pd.Series(list(metacyc2pg.keys()))
        metacyc_df['METACYC_ID'] = metacyc_df['METACYC_ID'].apply(str)
        metacyc_df['GENES'] = pd.Series(list(metacyc2pg.values())).apply( lambda x  : '|'.join(x))
        metacyc_df.to_csv(outpath +'/metacyc2proteingroups.csv')

        # IPR
        id2ipr= defaultdict(set)
        iprs = set()
        def ipr(df):
            id = df['ProteinAccession']
            ipr = str(df['InterProAnnotationsAccession']) +': ' + str(df['InterProAnnotationsDescription'])
            if ipr.startswith('IPR'):
                id2ipr[id].add(ipr)
                iprs.add(ipr)
        data.apply(ipr, axis=1)
        with open( outpath +'/accession2ipr.p', 'wb') as f:
            pickle.dump( id2ipr, f)

        ipr_df = pd.DataFrame()
        ipr_vals = list(iprs)
        ipr_df['IPR_ID'] = pd.Series(ipr_vals)
        ipr_df.to_csv(outpath +'/ipr_terms.csv')

        def iprintersect(val):
            vals = val.split(';')
            setlist = []
            for val in vals:
                vset = id2ipr[val]
                setlist.append(vset)
            union = set.union(*setlist)
            if len(union) > 0 :
                return ';'.join(union)
        proteins['_ipr.term.union']   = proteins[id_col].apply(parse_ids).apply(iprintersect)
        
        ipr2pg = defaultdict(set)
        
        def ipr2gene(df):
            ipr_terms = df['_ipr.term.union']
            pg = df['Identifier']
            try:
                ipr_terms = ipr_terms.split(';')
                for _ in ipr_terms:
                    ipr2pg[_].add(pg)
            except:
                pass
        proteins.apply(ipr2gene,axis=1)
        ipr_df = pd.DataFrame()
        ipr_df['IPR_ID'] = pd.Series(list(ipr2pg.keys()))
        ipr_df['IPR_ID'] = ipr_df['IPR_ID'].apply(str)
        ipr_df['GENES'] = pd.Series(list(ipr2pg.values())).apply( lambda x  : '|'.join(x))
        ipr_df.to_csv(outpath +'/ipr2proteingroups.csv')
        

        cmd = 'mq_genesets.R --outdir {} --keggid {}'.format(outpath, keggid)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        assert process.returncode == 0
        
        proteins.to_csv(outpath +'/ipr_target_proteins.txt',sep='\t')



    def ips_fasta(self, infile, outpath):
        cmd = 'ipr.py {}'.format(infile)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        assert process.returncode == 0
    
    def protein_id_lists(self, proteingroups, outfile):
        lst = proteingroups['Majority protein IDs'].apply(parse_ids).tolist()
        newlst = []
        for l in lst:
            newlst += l.split(';')
        l =pd.Series(newlst)
        l.to_csv(outfile , index = False)

    def export_pg_fasta(self, fasta_file, protein_id_file, outfile):
        fasta = SeqIO.parse(fasta_file,'fasta')
        new_fasta = []
        ids = pd.read_csv(protein_id_file, header = None)[0].tolist()
        for rec in fasta:
            if rec.id.split('|')[1] in ids:
                new_fasta.append(rec)
        SeqIO.write(new_fasta, outfile, 'fasta')

    def create_R_parameters(self, config, outfile, quant, group_level):
        vals=[]
        d = '#!/usr/bin/env R\n\n'
        vals.append(d)
        d = "cols <- c("
        vals.append(d)
        samples = config['samples']
        experiment = defaultdict(list)
        for sample in samples:
            if group_level in samples[sample]:
                group = samples[sample][group_level]
                if isinstance(group, str):
                    experiment[group].append(sample)
        groups = []
        reps = []
        for group in experiment:
            groups.append(group)
            for sample in experiment[group]:
                rep = "'{}.{}'".format(quant, sample)
                reps.append(rep)
        reps = ','.join(reps)
        d = reps
        vals.append(d)
        vals.append(')\n')
        d = 'f <- factor(c('
        vals.append(d)
        reps = []
        for group in groups:
            rep = "rep('{}',{})".format(group, len(experiment[group]))
            reps.append(rep)
        d = ','.join(reps)
        vals.append(d)
        d = "),\n"
        vals.append(d)
        levels = ','.join(["'" +group + "'" for group in groups])
        d ="levels=c({}))\n".format(levels)
        vals.append(d)
        d='design <- model.matrix(~0+f)\n'
        vals.append(d)
        d = 'colnames(design) <- c({})\n'.format(levels)
        vals.append(d)
        d = 'contrast.matrix <- makeContrasts(\n'
        vals.append(d)
        comps = []
        for comparison in config['comparisons'][group_level]:
            comp = '"{}-{}"'.format(comparison[0], comparison[1])
            comps.append(comp)
        d = ','.join(comps)
        vals.append(d)
        d = ',levels=design)\n'
        vals.append(d)
        template = ''.join(vals)
        w = open(outfile, 'w')
        w.write(template)
        w.close()

    def normalize(self, quant, infile, outdir, normalize, impute):
        cmd = 'mq_normalize_intensity.R -q {} -p {} -o {} -n {} -i {} && exit'.format(quant,infile, outdir,normalize, impute)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        assert process.returncode == 0
   

    def create_R_proteingroup_parameters(self):
        config= self.config
        vals=[]
        d = '#!/usr/bin/env R\n\n'
        vals.append(d)
        d = "cols <- c("
        vals.append(d)
        samples = config['samples']
        experiment = defaultdict(list)
        for sample in samples:
            group = samples[sample]['GROUP']
            experiment[group].append(sample)
        groups = []
        reps = []
        for group in experiment:
            groups.append(group)
            for sample in experiment[group]:
                rep = "'{}.{}'".format(self.protein_quant_parameter, sample)
                reps.append(rep)
        reps = ','.join(reps)
        d = reps
        vals.append(d)
        vals.append(')\n')
        d = 'f <- factor(c('
        vals.append(d)
        reps = []
        for group in groups:
            rep = "rep('{}',{})".format(group, len(experiment[group]))
            reps.append(rep)
        d = ','.join(reps)
        vals.append(d)
        d = "),\n"
        vals.append(d)
        levels = ','.join(["'" +group + "'" for group in groups])
        d ="levels=c({}))\n".format(levels)
        vals.append(d)
        d='design <- model.matrix(~0+f)\n'
        vals.append(d)
        d = 'colnames(design) <- c({})\n'.format(levels)
        vals.append(d)
        d = 'contrast.matrix <- makeContrasts(\n'
        vals.append(d)
        comps = []
        for comparison in config['comparisons']:
            comp = '"{}-{}"'.format(comparison[0], comparison[1])
            comps.append(comp)
        d = ','.join(comps)
        vals.append(d)
        d = ',levels=design)\n'
        vals.append(d)
        template = ''.join(vals)
        w = open(self.diff_dir + '/protein_experimental_design.R','w')
        w.write(template)
        w.close()
        peptide_txt = self.peptide_dir +'target_proteins.txt'
        self.target_peptides.to_csv(peptide_txt, sep='\t')
        cmd = 'cd {} && mq_normalize_intensity.R -d protein_experimental_design.R -p {} -o {} && exit'.format(self.diff_dir, peptide_txt,  self.diff_dir + 'peptide_normalization')
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        assert process.returncode == 0
        self.normalized_target_proteins = pd.read_csv(self.diff_dir + '/peptide_normalization/msnbase/normalized_protein_ibaq.csv') 
    
    def diff(self, exp_design, infile, outpath):
        cmd = 'mq_diff.R -d {} -f {} -o {}'.format(exp_design, infile, outpath)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        var = pd.read_csv(outpath +'/group_variance/group_variance.txt',sep='\t')
        var = var.dropna(axis=1, how='all')
        if len(var) > 0:
            names = []
            data = []
            for col in var:
                names.append(col)
                data.append(var[col].tolist())
            list_kw_dunn(names, data, 'Variance', 'Group', outpath + '/group_variance/' )
        

    def summarize_diff(self, inpath, outfile):
        w = open(outfile, 'w')
        
        ##################
        # GROUP VARIANCE #
        ##################
        
        if os.path.exists(inpath + '/group_variance/kw.txt'):
            with open(inpath + '/group_variance/kw.txt') as f:
                kw = f.read()
                kw_pval = float(kw.split('pvalue=')[1].split(')')[0])
            if kw_pval < 0.05:
                w.write('Kruskal-Wallis test showed a significant difference in the variance between groups (p-value {}).\n'.format(np.round(kw_pval,3)))
                
                dunn_bh = pd.read_csv(inpath + '/group_variance/dunn_bh.csv')
                dunn_bh = dunn_bh.set_index('Unnamed: 0')
                comps = []
                for col in dunn_bh.columns:
                    for ind in dunn_bh[col].index:
                        comp = set([ind, col])
                        if not comp in comps:
                            comps.append(comp)
                            pval = float(dunn_bh.loc[ind, col])
                            if  0 <= pval < 0.05:
                                w.write("Dunn's post-hoc test with BH correction showed a significant difference between the variance of {} and {} groups (adjusted p-value {}).\n".format(col, ind, str(np.round(pval,3))))
            else:
                w.write('Kruskal-Wallis test did not show a significant difference in the variance between groups (p-value {}).\n'.format(kw_pval))
            w.write('\n')
        
        #########
        # LIMMA #
        #########
        limma_files = [ inpath + '/limma/' + i for i in os.listdir(inpath + '/limma') if i.endswith('_intensity.txt') ]
        for f in limma_files:
            w.write('Comparison {}:\n'.format(f.split('limma_')[1].split('_intensity.txt')[0]))
            limma = pd.read_csv(f, sep='\t')
            ln = len(limma)
            sig = limma[limma['adj.P.Val'] < 0.05]
            ln_sig = len(sig)
            w.write('Out of {} comparisons, {} were significantly different with qval <  0.05.\n'.format(ln, ln_sig))
            for row in sig.iterrows():
                sig = np.round(row[1]['adj.P.Val'], 3)
                Id = row[1]['Identifier']
                control = row[1]['Control']
                exposed = row[1]['Exposed']
                logfc = np.round(row[1]['logFC'], 3)
                string = '{} was significantly different between {} and {}, with a log2(FC) of {} (qval. {}).\n'.format(Id, control, exposed, logfc, sig) 
                w.write(string)
            w.write('\n')
        w.write('\n')
        w.close()
    
    def summarize_gsea(self, inpath, outfile, config, group_level):
        w = open(outfile, 'w')
        for comp in config['comparisons'][group_level]:
            w.write('Comparison between {} and {}:\n\n'.format(comp[1], comp[0]))
            mcp = '{}_{}'.format(comp[0], comp[1])
            for fl in os.listdir(inpath + '/' + mcp):
                if fl.strip().rstrip().endswith('.csv'):
                    w.write(fl + '\n')
                    data = pd.read_csv(inpath + '/' + mcp + '/' + fl)
                    qdata = data[data['q.val'] < 0.05 ]
                    ln = len(qdata)
                    cmd = '{} terms where differentially enriched with qval < 0.05\n'.format(ln)
                    w.write(cmd)
                    for row in qdata.iterrows():
                        cmd = '{} (pval {}, qval {}, set size {})\n'.format(row[1]['RowName'], np.round(row[1]['p.val'],3),  np.round(row[1]['q.val'],3), row[1]['set.size'])
                        w.write(cmd)
                    w.write('\n')
            w.write('\n')
        w.close()
        
    def diff_analysis(self):
        for level in self.config['group_levels']:
            if not os.path.exists(self.diff_dir + '/' + level): 
                os.mkdir(self.diff_dir + '/' + level)

            # peptides
            diff_dir = self.diff_dir + '/' + level + '/peptide_diff'
            peptide_exp = self.diff_dir + '/peptide_experimental_design_{}.R'.format(level)
            peptides = self.unipept_dir + '/pept2lca_peptides.csv'
            self.diff(peptide_exp, peptides, diff_dir )
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')
            
            
            #pept2lca_phylum_sc
            diff_dir = self.diff_dir + '/' + level + '/pept2lca_phylum'
            peptides=self.unipept_dir + '/pept2lca_phylum_sc.csv'
            self.diff(peptide_exp, peptides, diff_dir )
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')
            
            #pept2lca_species_sc
            diff_dir = self.diff_dir + '/' + level + '/pept2lca_species'
            peptides=self.unipept_dir + '/pept2lca_species_sc.csv'
            self.diff(peptide_exp, peptides, diff_dir )
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')
            
            #pept2lca_genus_sc
            diff_dir = self.diff_dir + '/' + level + '/pept2lca_genus'
            peptides=self.unipept_dir + '/pept2lca_genus_sc.csv'
            self.diff(peptide_exp, peptides, diff_dir )
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')
            
            #pept2lca_genus_sc
            diff_dir = self.diff_dir + '/' + level + '/pept2lca_family'
            peptides=self.unipept_dir + '/pept2lca_family_sc.csv'
            self.diff(peptide_exp, peptides, diff_dir )
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')
            
            #pept2lca_sc
            diff_dir = self.diff_dir + '/' + level + '/pept2lca_taxon'
            peptides=self.unipept_dir + '/pept2lca_taxon_sc.csv'
            self.diff(peptide_exp, peptides, diff_dir ) 
            self.summarize_diff(diff_dir, diff_dir + '/summary.txt')
            
            # proteins
            diff_dir = self.diff_dir + '/' + level + '/protein_diff'
            protein_exp = self.diff_dir + '/protein_experimental_design_{}.R'.format(level)
            table = self.diff_dir + '/protein_normalization/msnbase/normalized.csv'
            self.diff(protein_exp, table,  diff_dir )
            self.summarize_diff(diff_dir ,diff_dir + '/summary.txt')
            
            # genesets diff
            for col in self.normalized_target_proteins.columns.tolist():
                if col.endswith('.term.union'):
                    diff_dir = self.diff_dir + '/' + level + '/' + col
                    self.aggregate_quant(self.normalized_target_proteins, col, '{}.'.format(self.protein_quant_parameter), diff_dir) 
                    table = diff_dir  + '/' + col + '.csv'
                    self.diff(protein_exp, table,  diff_dir )
                    self.summarize_diff(diff_dir , diff_dir + '/summary.txt')
            
            # KEGGG ORTHOLOGY
            col = 'Leading.Protein.Kegg.Orthology'
            diff_dir = self.diff_dir + '/' + level + '/' + col
            self.aggregate_quant(self.normalized_target_proteins, col, '{}.'.format(self.protein_quant_parameter), diff_dir )
            table = diff_dir + '/' + col + '.csv'
            self.diff(protein_exp, table, diff_dir)
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')

            # LEADING SPECIES
            col = 'Leading.Species'
            diff_dir = self.diff_dir + '/' + level + '/' + col
            self.aggregate_quant(self.normalized_target_proteins, col, '{}.'.format(self.protein_quant_parameter), diff_dir )
            table = diff_dir + '/' + col + '.csv'
            self.diff(protein_exp, table, diff_dir)
            self.summarize_diff(diff_dir , diff_dir + '/summary.txt')

    def ec2ko(self, ec):
        ko = rfunc.string2ko('[EC:{}]'.format(ec))
        return ko

    def string2ko(self, ec):
        ko = rfunc.string2ko(ec)
        return ko
    
    def leading_protein_mygene(self, proteingroups):
        genes = proteingroups['Leading Protein'].tolist()
        gene_dict = defaultdict()
        for gene in genes:
            gene_dict[gene] = defaultdict(list)
        mg = mygene.MyGeneInfo()
        df = mg.querymany(genes, as_dataframe=True, scopes='uniprot')
        df = df.reset_index()
        for row in df.iterrows():
            for col in row[1].index:
                gene_dict[row[1]['query']][col].append(row[1][col])
        for row in proteingroups.iterrows():
            entrez = gene_dict[row[1]['Leading Protein']]['entrezgene']
            try:
                entrez = ';'.join([ str(int(i)) for i in entrez ] )
                proteingroups.loc[row[0],'Mygene.entrez'] = entrez
            except:
                pass
        return proteingroups
    
    def leading_protein_gene(self, proteingroups, search_fasta, reference_fasta):
        genes = proteingroups['Leading Protein'].tolist()
        gene_dict = {}
        for p in search_fasta + reference_fasta:
            ID = p.id.split('|')[1]
            if 'GN=' in p.description:
                gn = p.description.split('GN=')[1].split()[0]
                gene_dict[ID] = gn
        proteingroups['Leading.gene'] = proteingroups['Leading Protein'].map(gene_dict)
        return proteingroups

    def leading_protein_ko(self, proteingroups):
        def _(df):
            try:
                ko = rfunc.up2ko(df['Leading Protein'])
                if not ko.ko == '':
                    print(ko.ko, ko.name)
                    df['Leading Protein Kegg Orthology ID'] = ko.ko
                    df['Leading Protein Kegg Orthology Name'] = ko.name
                    df['Leading Protein Kegg Orthology'] = ko.ko +' ' + ko.name
            except:
                pass
            return df
        proteingroups = proteingroups.apply(_, axis = 1)
        return proteingroups

    def up2ko(self, up):
        return rfunc.up2ko(up)
    
    def aggregate_quant(self, df, agg_col, quant_prefix, outpath):
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        data = df[df[agg_col].notnull()]
        x = data.assign(**{agg_col:data[agg_col].str.split(';')})
        data = pd.DataFrame({col:np.repeat(x[col].values, x[agg_col].str.len()) for col in x.columns.difference([agg_col])}).assign(**{agg_col:np.concatenate(x[agg_col].values)})[x.columns.tolist()]
        agg_cols = {}
        for col in data.columns.tolist():
            if col.startswith(quant_prefix):
                agg_cols[col] = np.sum
        agg = data.groupby(agg_col).agg(agg_cols) 
        agg = agg.reset_index()
        agg.rename(columns={agg_col:'Identifier'}, inplace=True)
        agg = agg[agg['Identifier'].notnull()]
        table=outpath +'/' + agg_col + '.csv'
        agg.to_csv(table)
    def parse_ids(self, ids):
        return parse_ids(ids)

def peptides(path):
    pep_dct = {}
    unique_pep_dct = {}
    peptides_path = path +'/peptides.txt'
    peps = pd.read_csv(peptides_path, sep=None, engine='python')
    peps = peps[peps['PEP'] < 0.01]
    for row in peps.iterrows():
        id_count  = 0
        for j in peps.columns:
            if j.startswith('Identification type'):
                if (row[1][j] == 'By MS/MS'):
                    id_count += 1
        peps.loc[row[0], 'Number of Replicates'] = id_count
    peps = peps[peps['Number of Replicates'] > 1 ]
    unique_peps = peps[peps['Unique (Groups)']=='yes']
    for row in peps.iterrows():
        index = row[0]
        sequence = row[1]['Sequence']
        if row[1]['Unique (Proteins)'] =='yes':
            unique_pep_dct[index] = sequence
        else:
            pep_dct[index] = sequence
    return pep_dct, unique_pep_dct

def evidence_table(txt_path, i):
    evp = txt_path +'/evidence.txt'
    ev = pd.read_csv(evp, sep =None, engine='python')
    ev = ev[ev['PEP'] < 0.01]
    sequences = list(set(ev['Sequence'].values))
    seq_dct = {}
    for item in sequences:
        num = len(set(ev[ev['Sequence'] == item]['Experiment']))
        seq_dct[item] = num
    for row in ev.iterrows():
        seq = row[1]['Sequence']
        reps = seq_dct[seq]
        ev.loc[row[0], 'Number of Replicates'] = reps
    ev = ev[ev['Number of Replicates'] > 1]
    return ev

def filteredPG(path, name):   #path to txt, returns filtered
    prot_path = path+'/proteinGroups.txt'
    pg = pd.read_csv(prot_path, sep=None, engine='python')
    potential_contaminants = pg[pg['Potential contaminant']=='+']
    reverse = pg[pg['Reverse'].notnull()]
    if len(reverse)>0:
        reverse = reverse[reverse['Reverse']=='+']
    pg1 = pg[(pg['Reverse'].isnull()) & (pg['Potential contaminant'].isnull())]
    pg1_err = pg1[(pg1['Protein IDs'].str.startswith('REV'))|(pg1['Protein IDs'].str.startswith('CON'))]
    pg1 = pg1[~((pg1['Protein IDs'].str.startswith('REV'))|(pg1['Protein IDs'].str.startswith('CON')))]
    pg1_single = pg1[pg1['Number of proteins'] == 1]
    pg1_multiple = pg1[pg1['Number of proteins']>1]
    
    sum = pd.Series([
        len(pg1),
        len(pg1_single), 
        len(pg1_multiple), 
        len(potential_contaminants),
        len(reverse),
        len(pg)])
    sum.name = name
    sum.index = ['Protein Groups (filtered)',
                'Protein Groups with one member (filtered)',
                'Protein Groups with multiple members (filtered)',
                'Potential Contaminant Protein Groups',
                'Reversed Sequence Hits Protein Groups',
                'Total Protein Groups (unfiltered)']

    return pg1, pg1_multiple, sum

def fasta_dct(path):
    rec_dct = {}
    dct = {}
    ls = list(SeqIO.parse(path, 'fasta'))
    for i in ls:
        id = i.id
        if id not in dct:
            dct[id]=str(i.seq)
            rec_dct[id] = i
    return dct, rec_dct

def add_fasta(df,  fasta_dct, peptide_dct, unique_peps,  db):
    new_df = df.copy()
    fasta_db = fasta_dct[db]
    peptide_db = peptide_dct[db]
    unique_db = unique_peps[db]
    
    old_cols = new_df.columns.tolist()
    prot_cols = []
    peptide_cols = ['Unique_Peptides', 'Non_Unique_Peptides']
    analysis_cols = ['Start_Sites_Found']
    for row in new_df.iterrows():
        index = row[0]
        protein_ids = row[1]['Protein IDs'].split(';')
        peptides = [peptide_db[int(i)] for i in row[1]['Peptide IDs'].split(';') if int(i) in peptide_db]
        unique_peptides=[unique_db[int(i)] for i in row[1]['Peptide IDs'].split(';') if int(i) in unique_db]
        trie = algo.trie_graph(peptides + unique_peptides)
        count = 1
        starts_found = []
        for i in protein_ids:
            colname='Protein_ID_{}'.format(count)
            if colname not in prot_cols:
                prot_cols.append(colname)
            if i in fasta_db:
                rec = fasta_db[i]
                new_id = rec.id
                new_description = rec.description
                new_seq = rec.seq
                positions = algo.trie_matching(trie, new_seq)
                if len(positions) > 0:
                    if positions[0] == 0:
                        starts_found.append(colname)
                new_seq = Seq(algo.trie_upper(trie, str(new_seq)))
                new_rec = SeqRecord(seq = new_seq, id = new_id, description = new_description)
                new_df.loc[index, colname] = new_rec.format('fasta')
            else:
                new_df.loc[index,colname]=i
            count += 1
        new_df.loc[index, 'Unique_Peptides'] = '\n'.join(unique_peptides)
        new_df.loc[index, 'Non_Unique_Peptides'] ='\n'.join(peptides)
        new_df.loc[index, 'Start_Sites_Found']='\n'.join(starts_found)
        new_cols = prot_cols + peptide_cols +  old_cols
        new_df = new_df[new_cols]
    return new_df

def cleavage_site(peptide, fasta_record, genome):
    seq = str(fasta_record.seq)
    description = fasta_record.description
    position = seq.find(peptide) * 3
    coords = [j.split('|') for j in description.split() if (j.startswith('(+)') or j.startswith('(-)'))]
    coords =  [item for sublist in coords for item in sublist] 
    cleavage_pep = True
    if len(coords) != 0:
        cleavage_pep = False
    for coord in coords:
        try:
            strand = coord.split(')')[0][1:]
            datum = coord.split(')')[1].split(':')
        except: 
            strand = coord.split(')')[0][1:]
            datum = coord.split(')')[1].split(':')
        start = int(datum[0])
        end = int(datum[1])
        cleavage = ['R','K']
        if strand == '-':
            pep_pos = end - position
            assert len(Seq(genome[(pep_pos):pep_pos+3]).reverse_complement()) % 3 == 0
            amino = translate(Seq(genome[(pep_pos):pep_pos+3]).reverse_complement(), cds = False, table = 11)
            assert len(genome[pep_pos:pep_pos + 3]) == 3 
        elif strand == '+':
            pep_pos = start + position
            upstream_pos = pep_pos - 3
            amino = translate(Seq(genome[(upstream_pos -1):pep_pos-1]), cds = False, table =11)
            assert len(genome[(upstream_pos -1):pep_pos-1]) == 3
        if amino in cleavage:
            cleavage_pep = True         
    return cleavage_pep

def get_start_coords(seq, coord, genome, starts, stops): # trailing end of a frame (from first found peptide), coordinates of frame. Returns the coordinates with the most downstream start codon for that sequence
    strand = coord.split(')')[0][1:]
    start = int(coord.split(')')[1].split(':')[0])-1
    end = int(coord.split('(')[1].split(':')[1])
    if strand == '+':
        assert len(genome[start:end]) % 3 == 0
        original_seq = translate(Seq(str(genome)[start:end]), table = 11, cds = False)
        position = (original_seq.find(seq[1:])-1) * 3 
        seq_pos = start + position
        #print genome[seq_pos -3: seq_pos + 3], genome[seq_pos:seq_pos + 3], translate(Seq(genome[seq_pos:seq_pos + 9])), seq[:10]
        if str(Seq(genome[seq_pos:seq_pos + 3])) not in starts:
            while (str(genome)[seq_pos - 3:seq_pos] not in starts + stops) and (seq_pos >= 0):
                seq_pos -= 3
                #print genome[seq_pos:seq_pos + 3],
            final_codon = str(Seq(genome[seq_pos -3: seq_pos]))
            #print final_codon
            assert final_codon in starts
            start = seq_pos - 3
        else:
            start = seq_pos
            assert genome[start: start + 3] in starts
    elif strand == '-':
        test_str = str(Seq(genome[end-6:end]).reverse_complement())
        original_seq = translate(Seq(str(genome)[start:end]).reverse_complement(), table = 11, cds = False)
        
        assert len(Seq(str(genome)[start:end]).reverse_complement()) % 3 == 0
        position = (original_seq.find(seq[1:]) -1) *3
        seq_pos = end - position 
        temp_codon = Seq(genome[seq_pos -17:seq_pos]).reverse_complement()
        if str(Seq(genome[seq_pos -3:seq_pos]).reverse_complement()) not in starts:
            while (str(Seq(genome[seq_pos:seq_pos +3]).reverse_complement()) not in starts + stops) and (seq_pos < len(genome)):
                seq_pos += 3
            final_codon= str(Seq(genome[seq_pos:seq_pos + 3]).reverse_complement())
            assert final_codon in starts
            end = seq_pos + 3
        else:
            end = seq_pos
            assert str(Seq(genome[seq_pos - 3 : seq_pos]).reverse_complement()) in starts
    return '({}){}:{}'.format(strand, start + 1, end)

def get_nucs(coord, genome):
    strand = coord.split(')')[0][1:]
    end = int(coord.split(':')[1])
    start = int(coord.split(')')[1].split(':')[0])
    if strand == '+':
        seq = genome[start-1:end]
    elif strand =='-':
        seq = str(Seq(genome[start-1:end]).reverse_complement())
    return seq


def ps(peptides, protein, frame1, frame2, frame3, frame4, frame5, frame6):
    gssps = []
    pssps = []
   
    protein = protein.upper()

    lp = len(protein)
    lp_replace = 'X' * lp

    #print 'lp_replace:', lp_replace
    for m in peptides:
        temp = m
        if protein.find(m) == 0:
            m = m[1:]
        minilocs = 0
        
        locations = [k.start() for k in re.finditer('(?={})'.format(m), frame1)]    
        minilocs += len(locations)
        
        locations = [k.start() for k in re.finditer('(?={})'.format(m), frame2)]
        minilocs += len(locations)
        
        locations = [k.start() for k in re.finditer('(?={})'.format(m), frame3)]
        minilocs += len(locations)
       
        locations = [k.start() for k in re.finditer('(?={})'.format(m), frame4)]
        minilocs += len(locations)
        
        locations = [k.start() for k in re.finditer('(?={})'.format(m), frame5)]
        minilocs += len(locations)
        
        locations = [k.start() for k in re.finditer('(?={})'.format(m), frame6)]
        minilocs += len(locations)
        
        if minilocs == 1:
            gssps.append(temp)

        if minilocs >= 1:
            minilocs = 0
            
            locations = [k.start() for k in re.finditer('(?={})'.format(m), frame1.replace(protein[1:], lp_replace[1:]))]    
            minilocs += len(locations)
            
            locations = [k.start() for k in re.finditer('(?={})'.format(m), frame2.replace(protein[1:], lp_replace[1:]))]
            minilocs += len(locations)
            
            locations = [k.start() for k in re.finditer('(?={})'.format(m), frame3.replace(protein[1:], lp_replace[1:]))]
            minilocs += len(locations)
           
            locations = [k.start() for k in re.finditer('(?={})'.format(m), frame4.replace(protein[1:], lp_replace[1:]))]
            minilocs += len(locations)
            
            locations = [k.start() for k in re.finditer('(?={})'.format(m), frame5.replace(protein[1:], lp_replace[1:]))]
            minilocs += len(locations)
            
            locations = [k.start() for k in re.finditer('(?={})'.format(m), frame6.replace(protein[1:], lp_replace[1:]))]
            minilocs += len(locations)
            
            if minilocs == 0:
                pssps.append(temp)
    return gssps, pssps


def novel_columns(df, genome,cdna):
    features = defaultdict(list)
    #frame1 = str(translate(Seq(genome[:]), cds =False, table = 11))
    #frame2 = str(translate(Seq(genome[1:]), cds = False, table = 11))
    #frame3 = str(translate(Seq(genome[2:]), cds = False, table = 11))
    #frame4 = str(translate(Seq(genome[:]).reverse_complement(), cds = False, table = 11))
    #frame5 = str(translate(Seq(genome[:-1]).reverse_complement(), cds = False, table = 11))
    #frame6 = str(translate(Seq(genome[:-2]).reverse_complement(), cds = False, table = 11))

    for row in df.iterrows():
        evidence = []
        evidence_code = 0
        nblasted = []
        datum = row[1]
        count = row[0]
        if datum['Alternate Start Translations'] == 'True':
            alt = True
        else:
            alt = False
        new_coords = datum['Consensus Coordinates'].split('\n')
        nucs = None
        overlapped = []
        pseudogene = 'False'
        sense_overlapped = []
        for coord in new_coords:
            nucs = get_nucs(coord, genome)
            ndatum = blast.nblast_pretty(nucs, 'MsmegmatisCDNAdb')
            start = int(coord.split(')')[1].split(':')[0])
            end = int(coord.split(':')[1])
            strand = coord.split(')')[0].split('(')[1]
            if strand == '+':
                strand = '1'
            elif strand == '-':
                strand = '-1'
            for f in cdna:
                nstart = int(f.description.split(':Chromosome:')[1].split(':')[0])
                nend = int(f.description.split(':Chromosome:')[1].split(':')[1])
                nstrand = f.description.split(':Chromosome:')[1].split(':')[2].split()[0]
                if ((start >= nstart ) and (start < nend)) or ((end <= nend) and (end > nstart)):
                    overlap = '{} overlaps with: {}'.format(coord, f.description)
                    overlapped.append(overlap)
                    if strand == nstrand:
                        sense_overlapped.append(overlap)
                    if 'pseudogene' in f.description.lower():
                        pseudogene = 'True'
            if len(ndatum) > 0:
                for rec in ndatum:
                    rstrand = rec.split(':Chromosome:')[1].split(':')[2].split()[0]
                    if strand == rstrand:
                        _ = coord+ ' ' + rec
                        nblasted.append(_)
        if len(nblasted) > 0:
            df.loc[count, 'CDNA blastn'] = '\n'.join(nblasted)
            evidence.append('Genomic sense CDNA feature allignment')
        else:
            df.loc[count, 'CDNA blastn'] = np.nan
        assert datum['Annotation_Type'] == 'Novel'
        if len(overlapped) == 0:
            df.loc[count, 'CDNA Overlaps'] = 'None'
        else:
            df.loc[count, 'CDNA Overlaps'] = '\n'.join(overlapped)
        if len(sense_overlapped) == 0:
            df.loc[count, 'CDNA Overlaps (Sense)'] = 'None'
        else:
            df.loc[count, 'CDNA Overlaps (Sense)'] = '\n'.join(sense_overlapped)
            evidence.append('Genomic sense CDNA feature overlap')
        df.loc[count, 'Pseudogene Overlap'] = pseudogene
        df.loc[count, 'nr blastp'] = '_'
        seq = StringIO.StringIO('>' + datum['Consensus Sequence'].split('>')[1])
        seq = SeqIO.read(seq,'fasta')
        longest_seq_str = str(seq.seq).upper() 
        
        #peptides = [''.join(''.join(''.join(p.split('_')).split('(ac)')).split('(ox)')) for p in datum['Non-Redundant Peptide Set'].split('\n')]
        #gss, pss = ps(peptides, longest_seq_str, frame1, frame2, frame3, frame4, frame5, frame6)
        
        gss  = df.loc[count, 'Genome Search Specific Peptides']
        pss  = df.loc[count, 'Protein Search Specific Peptides']
        
        if gss == np.nan:
            gss = ''
        if pss == np.nan:
            pss = ''

        if gss != '':
            gss = gss.split('\n')
        if pss != '':
            pss = pss.split('\n')

        
        num_peps = 'None'
        if len(pss) >= 2:
            num_peps = 'Two or more identified peptides'
            evidence.append(num_peps)
        
        if len(pss) == 1:
            num_peps = 'Only one identified peptide'
            evidence.append(num_peps)

        if len(gss) > 0:
           df.loc[count, 'GSSPs'] = '\n'.join([str(val) for val in gss])
        if len(pss) > 0:
            df.loc[count, 'PSSPs'] = '\n'.join([str(val) for val in pss])
       
        gssm_confirmed = False
        if ((len(gss) > 0) or (len(pss) > 0)):
            df.loc[count, 'GSSM or PSSM Peptides'] = 'True'
            gssm_confirmed = True
        else:
            df.loc[count, 'GSSM or PSSM Peptides'] = 'False'
        
        if gssm_confirmed == True:
            if alt == True:
                res = blast.pblast_pretty(longest_seq_str[1:],'nr')
                df.loc[count, 'nr blastp'] = res
            else:
                res = blast.pblast_pretty(longest_seq_str,'nr')
                df.loc[count, 'nr blastp'] =  res
            #df.loc[count, 'TeST'] = longest_seq_str
            if res != None:
                evidence.append('nr blastp allignment')
                res_accession = res.split('gi|')[1].split('|')[0]
                res_rec = blast.gi2up(res_accession)
                df.loc[count, 'Leading blastp GI id'] = res_accession
                if res_rec != None:
                    for r in res_rec:
                        if not isinstance(res_rec[r], str):
                            df.loc[count,r] = '\n'.join(res_rec[r])
                            for feature in res_rec[r]:
                                features[r].append(feature)
                        else:
                            df.loc[count, r] = res_rec[r]
                            features[r].append(res_rec[r])
        else:
            df.loc[count, 'nr blastp'] ='Sequence excluded (failed GSSP/PSSP search)'
        
        if ('nr blastp allignment' in evidence) and (gssm_confirmed == True) and ('Genomic sense CDNA feature overlap' in evidence) and (num_peps == 'Two or more identified peptides'):
            evidence_code = 1
        
        elif ('nr blastp allignment' in evidence) and (gssm_confirmed == True) and ('Genomic sense CDNA feature allignment' in evidence) and (num_peps == 'Two or more identified peptides'):
            evidence_code = 2

        elif ('nr blastp allignment' in evidence) and (gssm_confirmed == True) and (num_peps == 'Two or more identified peptides'):
            evidence_code = 3

        elif (gssm_confirmed == True) and (num_peps == 'Two or more identified peptides') and ('Genomic sense CDNA feature overlap' in evidence): 
            evidence_code = 4
        
        elif (gssm_confirmed == True) and (num_peps == 'Two or more identified peptides') and ('Genomic sense CDNA feature allignment' in evidence): 
            evidence_code = 5
        
        elif (gssm_confirmed == True) and (num_peps == 'Two or more identified peptides'): 
            evidence_code = 6

        elif ('nr blastp allignment' in evidence) and (gssm_confirmed == True) and ('Genomic sense CDNA feature overlap' in evidence) and (num_peps == 'Only one identified peptide'):
            evidence_code = 7
        
        elif ('nr blastp allignment' in evidence) and (gssm_confirmed == True) and ('Genomic sense CDNA feature allignment' in evidence) and (num_peps =='Only one identified peptide'):
            evidence_code = 8

        elif ('nr blastp allignment' in evidence) and (gssm_confirmed == True) and (num_peps == 'Only one identified peptide'):
            evidence_code = 9
        
        elif (gssm_confirmed == True) and ('Genomic sense CDNA feature overlap' in evidence) and (num_peps == 'Only one identified peptide'):
            evidence_code = 10
        
        elif (gssm_confirmed == True) and ('Genomic sense CDNA feature allignment' in evidence) and (num_peps == 'Only one identified peptide'):
            evidence_code = 11

        elif (gssm_confirmed == True) and (num_peps == 'Only one identified peptide'):
            evidence_code = 12

        elif num_peps == 'None':
            evidence_code = 13
        else:
            assert 4 == 5
        df.loc[count, 'Evidence'] = '\n'.join(evidence)
        df.loc[count, 'Evidence Code (Ranked)'] = evidence_code

    df.sort('Evidence Code (Ranked)', ascending = True, inplace = True)
    return df, features

def PGCombined(mapping_dct, ordered, stage_dct):
    stage1 = stage_dct.keys()[0]
    stage2 = stage_dct.keys()[1]
    stage1_values = stage_dct[stage1]
    stage2_values = stage_dct[stage2]
    txt_dct=defaultdict()
    fasta_dict= {}
    fasta_records = {}
    all_peps = {}
    
    unique_peps = {}
    names= name_dct()
    sum_list = []
    pg_dct = defaultdict()
    pg_mult_dct = defaultdict()
    id_dct = defaultdict(list) 
    ev_dct = {}
    for i in mapping_dct:
        if i == 'mapping':
            f= open(mapping_dct['mapping'])
            mapping = json.load(f)
        elif i == 'flatfile':
            f = mapping_dct['flatfile']
            genome = SeqIO.read(f, 'embl')
            genome = str(genome.seq)
        else:
            tup = mapping_dct[i]
            fasta_path=tup[0]
            seq_dct, rec_dct = fasta_dct(fasta_path)
            fasta_dict[i] = seq_dct
            fasta_records[i] = rec_dct
            txt_path =tup[1]
            pg, pg_mult, sum  = filteredPG(txt_path, i)
            evidence = evidence_table(txt_path, i)
            ev_dct[i] = evidence 
            sum_list.append(sum)
            pg_dct[i] = pg
            pg_mult_dct[i] = pg_mult
            peps, unique = peptides(txt_path)
            all_peps[i]=peps
            unique_peps[i]=unique
    
    frame1 = str(translate(Seq(genome[:]), cds =False, table = 11))
    frame2 = str(translate(Seq(genome[1:]), cds = False, table = 11))
    frame3 = str(translate(Seq(genome[2:]), cds = False, table = 11))
    frame4 = str(translate(Seq(genome[:]).reverse_complement(), cds = False, table = 11))
    frame5 = str(translate(Seq(genome[:-1]).reverse_complement(), cds = False, table = 11))
    frame6 = str(translate(Seq(genome[:-2]).reverse_complement(), cds = False, table = 11))

    pg_summary = pd.concat(sum_list, axis = 1)
    pg_summary.rename(columns=names, inplace =True)
    venn_dct = defaultdict(list)
    for i in pg_dct:
        id_dct[i] = pg_dct[i]['Protein IDs'].tolist()
    combined = DataFrame()
    count = 0
    new_cols = []
    for i in ordered:
        query_name = '_'.join(i.split())
        id_name = 'Leading_Protein_ID\n{}'.format(query_name)
        new_cols.append(id_name)
        unique_peptides = 'Unique_Peptides\n{}'.format(query_name)
        new_cols.append(unique_peptides)
        non_unique_peptides = 'Non_Unique_Peptides\n{}'.format(query_name)
        new_cols.append(non_unique_peptides)
        pep_scores = 'PEP\n{}'.format(query_name)
        new_cols.append(pep_scores)
        group_count = 'Group_Members\n{}'.format(query_name)
        new_cols.append(group_count)

        #start_identified = 'Start_Identified_{}'.format(query_name)
        #new_cols.append(start_identified) 

    new_cols.append('First_Amino_Found')
    new_cols.append('Most_Upstream_ID(s)')
    new_cols.append('Reference IDs Mapped')
    new_cols.append('Consensus Sequence')
    new_cols.append('Consensus Coordinates')
    new_cols.append('Start Codon') 
    new_cols.append('Start Translation') 
    new_cols.append('Identified Start Peptide At Trypsin Cleavage Site')
    new_cols.append('Annotation_Type')
    new_cols.append('Cleavage Site Peptides')
    new_cols.append('Cleavage Isoforms')
    new_cols.append('Acetyl (Protein N-term) Peptides')
    new_cols.append('Oxidation (M) Peptides')
    new_cols.append('Stage Unique: {}'.format(stage1))
    new_cols.append('Stage Unique: {}'.format(stage2))
    new_cols.append('Non-Redundant Peptide Set')
    new_cols.append('Alternate Start Translations')
    new_cols.append('Genome Search Specific Peptides')
    new_cols.append('Protein Search Specific Peptides')
    
    #new_cols.append('Leading blastp gi id')
    #new_cols.append('Annotations')
    #new_cols.append('CDS_Coordinates')
    #new_cols.append('Mulitple Start Sites')
    #new_cols.append('N-term acetylation Peptides (Combined)')
    #new_cols.append('Peptide_IDs_Combined')
    #new_cols.append('Consensus Sequence')
    #new_cols.append('Translation Variants')

    for i in ordered[:]:
        query_ids = id_dct[i]
        df = pg_dct[i]
        ev = ev_dct[i]
        fs = fasta_dict[i]
        fr =fasta_records[i]
        query_name ='_'.join(i.split())
        ordered.remove(i)
        id_name = 'Leading_Protein_ID\n{}'.format(query_name)
        id_col1= 'Unique_Peptides\n{}'.format(query_name)
        id_col2= 'Non_Unique_Peptides\n{}'.format(query_name)
        id_col3='PEP\n{}'.format(query_name)
        id_col4 = 'Group_Members\n{}'.format(query_name)
        for id in query_ids:
            gssps = []
            pssps = []
            non_redundant_peptides = []
            stage1_peptides = []
            stage2_peptides = []
            acetylations = []
            oxidations = []
            covered = []
            covered.append(i)
            new_id = id.split(';')[0]
            ref_seq_present = False
            if i == 'ref_prot':
                ref_seq_present = True
            reference_prot = None
            row_recs = {}
            combined_map = set()
            start_peptides = []
            len_dict = defaultdict(list)    # dict of the most upstream pep found  position per seq 
            starts_found = []
            row = df[df['Protein IDs'] == id]
            assert len(row) == 1
            peptide_ids=row['Peptide IDs'].iget(0).split(';')
            evidence_ids = row['Evidence IDs'].iget(0).split(';')
            evidence_peps = ev[ev['id'].apply(str).isin(evidence_ids)]
            
            combined_peps = []
            
            for pep in list(set(evidence_peps['Modified sequence'])):
                non_redundant_peptides.append(pep)
                if '(ac)' in pep:
                    acetylations.append(pep)
                if '(ox)' in pep:
                    oxidations.append(pep)
                datum = ''.join(pep.split('(ac)'))
                datum = ''.join(datum.split('(ox)'))
                datum = ''.join(datum.split('_'))
                combined_peps.append(datum)

            combined_peps = list(set(combined_peps))
            for pep in list(set(evidence_peps[evidence_peps['Experiment'].isin(stage1_values)]['Modified sequence'])):
                stage1_peptides.append(pep)
            for pep in list(set(evidence_peps[evidence_peps['Experiment'].isin(stage2_values)]['Modified sequence'])):
                stage2_peptides.append(pep)
           
            pep_score = row['PEP'].iget(0)
            group_members = row['Number of proteins'].iget(0)
            unique_peptides = [unique_peps[query_name][int(pep_id)] for pep_id in peptide_ids if int(pep_id) in unique_peps[query_name]]
            non_unique_peps = [all_peps[query_name][int(pep_id)] for pep_id in peptide_ids if int(pep_id) in all_peps[query_name]]
            
            row_peps = combined_peps
            
            if len(combined_peps) > 0:
                pep_trie = algo.trie_graph(combined_peps)
            temp_id = fr[new_id].id
            temp_description = fr[new_id].description
            row_recs[i] = fr[new_id]
            temp_seq = str(fr[new_id].seq)     #the protein sequence
            first_row_seq = temp_seq
            if i == 'ref_prot':
                reference_prot = fr[new_id]

            confirmed_unique = []
            if len(combined_peps) > 0:    
                gss, pss = ps(row_peps, temp_seq, frame1, frame2, frame3, frame4, frame5, frame6) 
                confirmed_unique = pss
                c_u_trie  = algo.trie_graph(confirmed_unique)
                positions = algo.trie_matching(c_u_trie, temp_seq)
                for p in gss:
                    if p not in gssps:
                        gssps.append(p)
                for p in pss:
                    if p not in pssps:
                        pssps.append(p)
            
            else:
                positions = []
            if positions[:1] == 0:
                starts_found.append(query_name)
            verify = temp_seq.lower()

            if len(combined_peps) > 0:
                temp_seq = algo.trie_upper(pep_trie, temp_seq)
            else:
                temp_seq = temp_seq.lower()

            temp_seq_upper = temp_seq.upper()
            #check if peptide starts at a not cleavage site
            for peptide in confirmed_unique:
                if cleavage_site(peptide, fr[new_id], genome) == False:
                    start_peptides.append(peptide)
            
            assert verify == temp_seq.lower()

            if len(positions) > 0:
                upstream = positions[0]
                temp_upstream = verify[upstream:]
            else:
                temp_upstream = ''
            len_upstream = len(temp_upstream)
            len_dict[len_upstream].append(i)
            temp_rec = SeqRecord(seq = Seq(temp_seq), id= temp_id, description = temp_description)
            combined.loc[count, id_name]=temp_rec.format('fasta')
            combined.loc[count,id_col2 ]='\n'.join(non_unique_peps)
            combined.loc[count, id_col1] = '\n'.join(unique_peptides)
            combined.loc[count, id_col3] = pep_score
            combined.loc[count, id_col4] = group_members
            venn_dct[i].append(abs(int(count)))
            found = False
            seq = fs[new_id]
            try:
                map = mapping[seq]
            except:
                try:
                    map = mapping['V' + seq[1:]]
                except:
                    try:
                        map = mapping['L'+seq[1:]]
                    except:
                        map = [10110000000]
            combined_map.update(map)
            for j in ordered[:]:
                found=False
                temp_name = '_'.join(j.split())
                template_ids = id_dct[j]     
                ev1 = ev_dct[j]
                tdf = pg_dct[j]
                tfs = fasta_dict[j]
                for k in template_ids:
                    new_k = k.split(';')[0]
                    temp_seq = tfs[new_k]
                    if found == False:
                        try:
                            temp_map = mapping[temp_seq]
                        except:
                            try: 
                                temp_map = mapping['L'+ temp_seq[1:]]
                            except:
                                try:
                                    temp_map = mapping['V'+temp_seq[1:]]
                                except:
                                    temp_map = [1000000000]
                        intersection = set(map).intersection(temp_map)
                        if len(intersection) >= 1:
                            if j == 'ref_prot':
                                ref_seq_present = True
                                reference_prot = fasta_records[j][new_k]
                            covered.append(j)
                            combined_map.update(temp_map)
                            #combined.loc[count, 'Leading_Protein_ID_{}'.format(temp_name)]=fasta_records[j][k].format('fasta')
                            query_id_name = 'Leading_Protein_ID\n{}'.format(temp_name)
                            query_col1= 'Unique_Peptides\n{}'.format(temp_name)
                            query_col2= 'Non_Unique_Peptides\n{}'.format(temp_name)
                            query_col3 = 'PEP\n{}'.format(temp_name)
                            query_col4 = 'Group_Members\n{}'.format(temp_name)
                            row = tdf[tdf['Protein IDs'] == k]
                            assert len(row) == 1
                            query_peptide_ids=row['Peptide IDs'].iget(0).split(';') 
                            evidence_ids = row['Evidence IDs'].iget(0).split(';')
                            evidence_peps = ev1[ev1['id'].apply(str).isin(evidence_ids)]
                            
                            combined_peps = []

                            for pep in list(set(evidence_peps['Modified sequence'])):
                                non_redundant_peptides.append(pep)
                                if '(ac)' in pep:
                                    acetylations.append(pep)
                                if '(ox)' in pep:
                                    oxidations.append(pep)
                                datum = ''.join(pep.split('(ac)'))
                                datum = ''.join(datum.split('(ox)'))
                                datum = ''.join(datum.split('_'))
                                combined_peps.append(datum)

                            combined_peps = list(set(combined_peps))

                            for pep in list(set(evidence_peps[evidence_peps['Experiment'].isin(stage1_values)]['Modified sequence'])):
                                stage1_peptides.append(pep)
                            for pep in list(set(evidence_peps[evidence_peps['Experiment'].isin(stage2_values)]['Modified sequence'])):
                                stage2_peptides.append(pep)

                            pep_score = row['PEP'].iget(0)
                            group_members = row['Number of proteins'].iget(0)
                            unique_peptides = [unique_peps[temp_name][int(pep_id)] for pep_id in query_peptide_ids if int(pep_id) in unique_peps[temp_name]]
                            non_unique_peps = [all_peps[temp_name][int(pep_id)] for pep_id in query_peptide_ids if int(pep_id) in all_peps[temp_name]] 
                            
                            
                            #combined_peps = list(set(unique_peptides+non_unique_peps))
                            
                            for pep in combined_peps:
                                row_peps.append(pep)
                            
                            if len(combined_peps) > 0:
                                pep_trie = algo.trie_graph(combined_peps)
                            
                            temp_id = fasta_records[j][new_k].id
                            temp_description = fasta_records[j][new_k].description
                            row_recs[j] = fasta_records[j][new_k]
                            temp_seq = str(fasta_records[j][new_k].seq) 
                            confirmed_unique = []

                            if len(combined_peps) > 0:
                                gss, pss = ps(combined_peps, temp_seq, frame1, frame2, frame3, frame4, frame5, frame6) 
                                confirmed_unique = pss
                                c_u_trie  = algo.trie_graph(confirmed_unique)
                                positions = algo.trie_matching(c_u_trie, temp_seq)
                                for p in gss:
                                    if p not in gssps:
                                        gssps.append(p)
                                for p in pss:
                                    if p not in pssps:
                                        pssps.append(p)
                            else:
                                positions = []

                            if positions[:1] == 0:
                                starts_found.append(temp_name)
                            
                            for peptide in confirmed_unique:
                                if cleavage_site(peptide, fasta_records[j][new_k], genome) == False:
                                    start_peptides.append(peptide)
                            
                            if len(positions) > 0:
                                upstream = positions[0]
                                temp_upstream = temp_seq[upstream:]
                            else:
                                temp_upstream = ''
                            len_upstream = len(temp_upstream)

                            len_dict[len_upstream].append(j)
                            
                            if len(combined_peps) > 0:
                                temp_seq = algo.trie_upper(pep_trie, temp_seq)
                            else:
                                temp_seq = temp_seq.lower()

                            temp_rec = SeqRecord(seq = Seq(temp_seq),
                                    description = temp_description,
                                    id = temp_id)
                            combined.loc[count, query_id_name]=temp_rec.format('fasta')
                            combined.loc[count,query_col2 ]='\n'.join(non_unique_peps)
                            combined.loc[count, query_col1] = '\n'.join(unique_peptides)
                            combined.loc[count, query_col3] = pep_score
                            combined.loc[count, query_col4] = group_members
                            venn_dct[j].append(count)
                            id_dct[j].remove(k)
                            found=True

            non_redundant_peptides = list(set(non_redundant_peptides))
            stage1_peptides_new = list(set([pep for pep in stage1_peptides if pep not in stage2_peptides]))
            stage2_peptides_new = list(set([pep for pep in stage2_peptides if pep not in stage1_peptides]))
            acetylations = set(acetylations)
            oxidations = set(oxidations)
            
            combined.loc[count, 'Acetyl (Protein N-term) Peptides'] = '\n'.join(acetylations)
            combined.loc[count, 'Oxidation (M) Peptides'] = '\n'.join(oxidations)
            combined.loc[count, 'Stage Unique: {}'.format(stage1)] = '\n'.join(stage1_peptides_new)
            combined.loc[count, 'Stage Unique: {}'.format(stage2)] = '\n'.join(stage2_peptides_new)
            combined.loc[count, 'Non-Redundant Peptide Set'] = '\n'.join(non_redundant_peptides)
            if len(pssps) > 0:
                combined.loc[count, 'Genome Search Specific Peptides'] = '\n'.join(gssps)
                combined.loc[count, 'Protein Search Specific Peptides'] = '\n'.join(pssps)
            
            covered = ''.join(list(set(covered)))   #checks that the only hit was from a single database

            combined.loc[count,'First_Amino_Found'] = '\n'.join(starts_found)       
            most_upstream = max(len_dict)
            longest_seqs = []
            for item in len_dict[most_upstream]:
                temp = row_recs[item]
                seq = str(temp.seq)[-most_upstream:]
                if seq not in longest_seqs:
                    longest_seqs.append(seq)
            combined.loc[count, 'Alternate Start Translations'] = str(len(longest_seqs) > 1)
            combined.loc[count, 'Most_Upstream_ID(s)'] = '\n'.join(len_dict[most_upstream])
            refs_mapped = []
            for reference_record in fasta_records['ref_prot']:
                rec_seq = str(fasta_records['ref_prot'][reference_record].seq)
                try:
                    intersection = combined_map.intersection(mapping[rec_seq])  
                    if len(intersection) >= 1:
                        refs_mapped.append(fasta_records['ref_prot'][reference_record])
                except:
                    pass
             
            combined_trie = algo.trie_graph(row_peps)
            refs_mapped_upper = [SeqRecord(seq = Seq(algo.trie_upper(combined_trie, str(ref_rec.seq))), id = ref_rec.id, description = ref_rec.description).format('fasta') for ref_rec in refs_mapped]
            #combined.loc[count, 'Reference IDs Mapped'] = '\n'.join([ref_rec.format('fasta') for ref_rec in refs_mapped])
            combined.loc[count, 'Reference IDs Mapped'] = '\n'.join(refs_mapped_upper)

            reference_sequence =list(set([str(ref_rec.seq) for ref_rec in refs_mapped]))
            novel = False        
            
            if (len(reference_sequence) < 1):
                if (len(refs_mapped) == 0) and (ref_seq_present != True):
                    novel = True
                    combined.loc[count, 'Annotation_Type'] = 'Novel'
                elif (len(refs_mapped) == 0) and (covered == 'ref_prot') and (most_upstream == len(first_row_seq)):
                    combined.loc[count, 'Annotation_Type'] = 'Validation'
            
            elif (len(reference_sequence) == 1):
                if most_upstream == len(reference_sequence[0]):
                    combined.loc[count, 'Annotation_Type'] = 'Validation'
                elif most_upstream > len(reference_sequence[0]):
                    combined.loc[count, 'Annotation_Type'] = 'Upstream Start'
            
            elif len(reference_sequence) > 1:
                ref_lens = []
                annotations = []
                longest_reference = 0
                for ref in refs_mapped:
                    ref_lens.append(len(str(ref.seq)))
                    datum = len(str(ref.seq))
                    if datum > longest_reference:
                        longest_reference = datum
                if most_upstream > longest_reference:
                    combined.loc[count, 'Annotation_Type'] = 'Upstream Start'
                elif most_upstream in ref_lens:
                    combined.loc[count, 'Annotation_Type'] = 'Validation'
            
            combined.loc[count, 'Cleavage Site Peptides'] = '\n'.join(list(set(start_peptides)))
            
            longest_seq_str = str(row_recs[len_dict[most_upstream][0]].seq)[-most_upstream:]
            
            cleavage_peps_found = False

            if len(start_peptides) > 0:
                cleavage_peps_found = True
                cleavage_trie = algo.trie_graph(start_peptides)
                clv_lst = algo.trie_matching(cleavage_trie, longest_seq_str) 
                combined.loc[count, 'Cleavage Isoforms'] = len(clv_lst)
            else:
                combined.loc[count, 'Cleavage Isoforms'] = 0
            rec_key = [key for key in row_recs.keys() if ((key != 'ref_prot') and (len(str(row_recs[key].seq)) >= most_upstream) and (longest_seq_str in str(row_recs[key].seq)))]
            
            if len(rec_key) != 0:
                longest_seq_coords =row_recs[rec_key[0]].description.split()
                coords = []
                for item in longest_seq_coords:
                    if (item.startswith('(+)') or item.startswith('(-)')):
                        coords.append(item)
            else:
                coords = []
                rec_key = [key for key in row_recs.keys() if ((key == 'ref_prot') and (len(str(row_recs[key].seq)) >= most_upstream) and (longest_seq_str in str(row_recs[key].seq)))]

                if len(rec_key) != 0:
                    longest_seq_coords =row_recs[rec_key[0]].description.split()
                    coords = []
                    for item in longest_seq_coords:
                        if (item.startswith('(+)') or item.startswith('(-)')):
                            coords.append(item)

            if len(coords) == 1:
                nblasted = [] 
                new_coords = []
                starts = ['ATG', 'GTG', 'TTG']
                stops =  ['TGA', 'TAA', 'TAG']
                coords = coords[0].split('|')
                for coord in coords:
                    new_coord = get_start_coords(longest_seq_str, coord, genome, starts, stops)
                    new_coords.append(new_coord)
                combined.loc[count, 'Consensus Coordinates']  = '\n'.join(new_coords) 
                start_codons = []
                for coord in new_coords:
                    nucs = get_nucs(coord, genome) 
                    start_codon = nucs[:3]
                    start_codons.append(start_codon)
                    #if novel == True:
                    #    combined.loc[count, 'nr blastp'] = '_'
                    #    if len(longest_seqs) > 1:
                    #        res = blast.pblast_pretty(longest_seq_str.upper()[1:],'nr')
                    #        combined.loc[count, 'nr blastp'] = res
                    #    else:
                    #        res = blast.pblast_pretty(longest_seq_str.upper(),'nr')
                    #        combined.loc[count, 'nr blastp'] =  res
                    #    
                    #    if res != None:
                    #        res_accession = res.split('gi|')[1].split('|')[0]
                    #        res_rec = blast.gi2up(res_accession)
                    #        if res_rec != None:
                    #            for r in res_rec:
                    #                combined.loc[count, 'Leading blastp gi id'] = res_accession
                    #                combined.loc[count,r] = '\n'.join(res_rec[r])
                    #    
                    #    datum = blast.nblast_pretty(nucs, 'MsmegmatisCDNAdb')
                    #    if len(datum) > 0:
                    #        for rec in datum:
                    #            _ = coord+ ' ' + rec
                    #            nblasted.append(_)
                #if len(nblasted) > 0:
                    #combined.loc[count, 'CDNA blastn'] = '\n'.join(nblasted)
                assert len(nucs) % 3 == 0
                consensus = translate(Seq(nucs), cds = False, table = 11)
                if consensus.endswith('*'):
                    consensus = consensus[:-1]
                if len(consensus) == len(longest_seq_str):
                    consensuses = []
                    start_amini = []
                    mini_count = 1
                    for seq in longest_seqs:
                        start_amino  = seq[:1] 
                        consensus = algo.trie_upper(combined_trie, seq)
                        consensus_rec = SeqRecord(seq = Seq(consensus), id = 'Consensus_Protein_Sequence_{}.{}'.format(count, mini_count) , description = '|'.join(new_coords))
                        consensuses.append(consensus_rec.format('fasta'))
                        start_amini.append(start_amino)
                        mini_count  += 1
                    start_amini.sort()
                    combined.loc[count, 'Consensus Sequence'] = '\n'.join(consensuses)
                    start_amini.sort()
                    combined.loc[count, 'Start Translation'] = ', '.join(start_amini)
                    if cleavage_peps_found == True:
                        if clv_lst[0] == 0:
                            combined.loc[count,  'Identified Start Peptide At Trypsin Cleavage Site'] = 'False'
                        else:
                            combined.loc[count, 'Identified Start Peptide At Trypsin Cleavage Site'] = 'True'
                    else:
                        combined.loc[count, 'Identified Start Peptide At Trypsin Cleavage Site'] = 'True'

                else:
                    consensus = algo.trie_upper(combined_trie, consensus)
                    consensus_rec = SeqRecord(seq = Seq(consensus), id = 'Consensus_Protein_Sequence_{}'.format(count) , description = '|'.join(new_coords))
                    combined.loc[count, 'Consensus Sequence'] = consensus_rec.format('fasta')
                    combined.loc[count, 'Start Translation'] = 'Not Verified'
                start_codons = list(set(start_codons))
                start_codons.sort()
                combined.loc[count, 'Start Codon'] = ', '.join(start_codons)
            count += 1
            

    mult_df = []
    for i in pg_mult_dct:
        datum = pg_mult_dct[i]
        datum['Database'] = i
        new = ['Database']
        if len(datum) > 0:
            datum = add_fasta(datum, fasta_records, all_peps, unique_peps, i)
            old = filter(lambda a: a != 'Database', datum.columns)
            new += old
            datum = datum[new]
            mult_df.append(datum)
    mult_df = pd.concat(mult_df, axis = 0)
    mult_df = mult_df[new]
    combined = combined[combined['Protein Search Specific Peptides'].notnull()]
    combined.reset_index(inplace=True)
    del combined['index']
    return pg_summary, combined[new_cols], venn_dct, mult_df,genome 

def venn(dct, row, col, worksheet, names):
    temp_dir = tempfile.mkdtemp()
    temp_out = temp_dir+'/out.png'
    sets = []
    new_names = []
    n = name_dct()
    for i in names:
        temp = set(dct[i])
        sets.append(temp)
    names = tuple(names)
    if len(names) == 3:
        v = venn3([i for i in sets], names, alpha=0.4, normalize_to=1.0)
    elif len(names) == 2:
        v = venn2([i for i in sets], names,  alpha=0.4, normalize_to=1.0)
    plt.savefig(temp_out, bbox_inches='tight')
    plt.clf()
    worksheet.insert_image(row, col, temp_out)

class Evidence:
    
    def __init__(self, table):
        self.table  = table
        self.modifications = self.table['Modifications'].tolist()
        self.sequences = self.table['Sequence'].tolist()
        self.experiments = self.table['Experiment'].tolist()
        self.n_term_acetylated = self.table['Modified sequence'].apply(lambda x : '_(ac)' in str(x)) 
        self.methionine_oxidation = self.table['Modified sequence'].apply(lambda x : 'M(ox)' in str(x)) 

        self.modification_map = {}
        self.experiment_modification_map = {}

        for index in range(len(self.sequences)):
            peptide = self.sequences[index]
            experiment = self.experiments[index]
            n_term_acetylated = self.n_term_acetylated[index]
            methionine_oxidation = self.methionine_oxidation[index]

            if not experiment in self.modification_map:
                self.modification_map[experiment] = defaultdict(set)
            if n_term_acetylated == True:
                self.modification_map[experiment]['_(ac)'].add(peptide)
            if methionine_oxidation == True:
                self.modification_map[experiment]['M(ox)'].add(peptide)

            self.modification_map[experiment]['All'].add(peptide)
    
    def export_peptides(self, experiments, modifications):
        export = set()
        for experiment in experiments:
            for modification in modifications:
                export.update(self.modification_map[experiment][modification])
        return export



