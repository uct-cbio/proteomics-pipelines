#!/usr/bin/env python

import pandas as pd
from pandas import DataFrame
import collections
from collections import defaultdict
import numpy as np
import collections
from collections import defaultdict
from collections import Counter

import math
from math import log
from scipy import stats

def spectral_counts(msms):
    spectral = pd.DataFrame()
    count = 0
    raw_cols = list(set(msms['Raw file']))
    for name, group in msms.groupby('Sequence'):
        counts =  Counter(group['Raw file'])
        spectral.loc[count, 'Sequence'] = name
        for key in counts:
            spectral.loc[count, "Spectral count " + key] = counts[key]
        proteing = list(set(group['Protein group IDs']))
        proteing = [str(i) for i in proteing]
        spectral.loc[count, 'Protein group IDs'] = '; '.join(proteing)
        count += 1
    return spectral

#def quantify_sc(spectral_counts, group1, group2):

def spectral_norm(spectral, valcols): # nomalizes on common identified peptides across all replicates
    normind  = spectral[valcols].dropna(how='any').index
    print('Normind: ', len(normind), 'All: ', spectral)
    normdata = spectral.loc[normind].copy()
    for col in valcols:
        colmean = normdata[col].mean()
        normdata[col] = normdata[col]/colmean
    return normdata

def median_normalize(peptides):  # median normalize maxquant peptides or protein raw intensity values
    peps = peptides.copy()
    intensity_cols = []
    all_values = []
    for col in peps.columns:
        if col.startswith('Intensity') and (col != 'Intensity'):
            intensity_cols.append(col)
            peps[col] = peps[col].replace(0,np.nan)
            peps[col] = peps[col].astype('float64')
            all_values += peps[col].dropna().tolist()
            col_med = float(np.median(peps[col].dropna().values))
            peps[col] = peps[col]/col_med
    global_median = np.median(all_values)
    for col in intensity_cols:
        peps[col] = peps[col] * global_median
    return peps

def median_normalize(peptides):  # median normalize maxquant peptides or protein raw intensity values
    peps = peptides.copy()
    intensity_cols = []
    all_values = []
    for col in peps.columns:
        if col.startswith('Intensity') and (col != 'Intensity'):
            intensity_cols.append(col)
            peps[col] = peps[col].replace(0,np.nan)
            peps[col] = peps[col].astype('float64')
            all_values += peps[col].dropna().tolist()
            col_med = float(np.median(peps[col].dropna().values))
            peps[col] = peps[col]/col_med
    global_median = np.median(all_values)
    for col in intensity_cols:
        peps[col] = peps[col] * global_median
    return peps

def quantify(proteingroups, peptides, valcols, level):   #quantify proteins using median normalised peptides
    df = proteingroups.copy()
    new_df = []
    for col in valcols:
        df[col] = df[col].astype(np.float64)
        df[col] = np.nan
    for row in df.iterrows():
        print(row[1].index)
        id = row[1]['id']
        prot_peps = peptides[peptides['Protein group IDs'] == str(id)][valcols].dropna(how='any')
        if len(prot_peps) != 0:
            prot_peps = prot_peps.median()
            for rep in prot_peps.index:
                df.loc[row[0], rep] = prot_peps.loc[rep]
        else:
            df.loc[row[0], 'DELETE'] = 'TRUE'
    df = df[df['DELETE'] != 'TRUE']
    del df['DELETE']
    for name, group in df.groupby(level):
        new = pd.DataFrame()
        med_vals = group[valcols].mean()
        new.loc[0, 'Group'] = level
        new.loc[0, 'Group gene names'] = '\n'.join(group['Gene names'].apply(str).values.tolist())
        new.loc[0, 'Group protein names'] = '\n'.join(group['Protein names'].apply(str).values.tolist())
        new.loc[0, 'Group protein families'] = '\n'.join(group['Protein families'].apply(str).values.tolist())
        new.loc[0, 'Group pathways'] = '\n'.join(group['Pathway'].apply(str).values.tolist())
        new.loc[0, 'Group operon IDs'] = '\n'.join(group['OperonID'].apply(str).values.tolist())
        for rep in med_vals.index:
            new.loc[0, rep] = med_vals.loc[rep]
        new_df.append(new)
    new = pd.concat(new_df)
    new = new.reset_index()
    del new['index']
    return new

def impute(values, num_missing,  width = 0.3, downshift = 1.8):
    std = np.std(values, ddof = 1, dtype=np.float64)       # need to subtract ddof?
    mean =  np.mean(values, dtype=np.float64)
    shifted = mean - (downshift * std)
    generated = np.random.normal(loc=shifted, scale = width * std, size = num_missing)
    return generated

def log2(value):
    try:
        return np.log2(value)
    except:
        return np.nan

def neg_log10(value):
    try:
        return -1 * log(value, 10)
    except:
        return np.nan


def ttest(df_, treated, control):
    df = df_.copy()
    for row in df.iterrows():
        tr = []
        cr = []
        for col in treated:
            datum  = row[1][col]
            if str(datum) != str(np.nan):
                tr.append(datum)
        for col in control:
            datum = row[1][col]
            if str(datum) != str(np.nan):
                cr.append(datum)
        t, p =  stats.ttest_ind(cr, tr,equal_var = True)
        df.loc[row[0],'t.value'] = t
        df.loc[row[0],'p.value'] = p
        fc = np.mean(tr) - np.mean(cr)
        df.loc[row[0], 'logFC'] = fc
    df =df[df['p.value'].notnull()]
    df = df.sort('p.value',ascending = True)
    return df


def welch(df_, treated, control):
    df = df_.copy()
    for row in df.iterrows():
        tr = []
        cr = []
        for col in treated:
            datum  = row[1][col]
            if str(datum) != str(np.nan):
                tr.append(datum)
        for col in control:
            datum = row[1][col]
            if str(datum) != str(np.nan):
                cr.append(datum)
        t, p =  stats.ttest_ind(cr, tr,equal_var = False)
        df.loc[row[0],'t.value'] = t
        df.loc[row[0],'p.value'] = p
        fc = np.mean(tr) - np.mean(cr)
        df.loc[row[0], 'logFC'] = fc
    df = df.sort('p.value',ascending = True)
    return df

def ranksum(df_, treated, control):
    df = df_.copy()
    for row in df.iterrows():
        tr = []
        cr = []
        for col in treated:
            datum  = row[1][col]
            if str(datum) != str(np.nan):
                tr.append(datum)

        for col in control:
            datum = row[1][col]
            if str(datum) != str(np.nan):
                cr.append(datum)
        z, p =  stats.ranksums(cr, tr)
        df.loc[row[0],'z.value'] = z
        df.loc[row[0],'p.value'] = p
        fc = np.mean(tr) - np.mean(cr)
        df.loc[row[0], 'logFC'] = fc
    df = df.sort('p.value',ascending = True)
    return df


def filtervalid(df, groups, minimum):   #dataframe, groups: dict (keys are groupnames, values are groupmembers), cutoff <= minimum fraction of non-np.nan values per group, in each row (0.5 requires more than 50% non-missing values per experimental group)
    groupmin = {}
    for i in groups:
        ln = len(groups[i])
        groupmin[i] = ln*minimum
    df1 = df.copy()
    allocations = {}
    for column in df1.columns:
        seencount = 0
        for group in groups:
            for member in groups[group]:
                if column == member:
                    seencount +=1
                    allocations[column] = group
        assert seencount <= 1
    #print allocations
    for row in df1.iterrows():
        ind = row[0]
        rowcounts = defaultdict(list)
        valid = True
        
        for col in row[1].dropna().index:
            if col in allocations:
                gr = allocations[col]
                rowcounts[gr].append(row[1][col])
        try:
            for key in rowcounts:
                assert len(rowcounts[key]) > groupmin[key]
            for key in groups:
                assert key in rowcounts
        except:
            valid = False
        for key in rowcounts:
            df1.loc[ind, 'Group_Valid_Values_{}'.format(key)] =  len(rowcounts[key])
        if valid != False:
            df1.loc[ind,'Valid'] = 'VALID'
    df1 = df1[df1['Valid'] == 'VALID']
    del df1['Valid'] 
    return df1

#def groupmeans(df, groups):
#    meancols = {}
#    seen = set()
#    for i in df.columns:
#        regcount = 0
#        for j in groupregexlist:
#            if i.startswith(j):
#                regcount += 1
#                seen.add(j)
#                assert i not in groupregexlist
#                meancols[i] = j
#        assert regcount <= 1
#    assert len(seen) == len(groupregexlist)
#    for row in df.iterrows():
#        means = defaultdict(list)
#        for item in row[1].index:
#            try:
#                assert item not in groupregexlist
#                means[meancols[item]].append(row[1][item])
#            except:
#                pass
#        for item in means:
#            mean = np.mean(means[item], dtype = 'float64')
#            df.loc[row[0], 'mean_{}'.format(item)] = mean
#
#    return df

