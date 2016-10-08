#!/usr/bin/env python

import subprocess
import sys
import pandas as pd
import collections
from collections import Counter
import operator
import os

infile = sys.argv[1]
seqcol = sys.argv[2]
outfolder = sys.argv[3]

def isgapped(val):
    if '-' in val:
        return '+'
    else:
        return '-'


def ILequivalence(peptide):
    plist = list(peptide)
    newplist = []
    for p in plist:
        if p == 'I':
            newplist.append('L')
        else:
            newplist.append(p)
     newp = ''.join(newplist)
     return newp

table = pd.read_csv(sys.argv[1], sep=None, engine='python')


table = table[(table['_alignment_rank'] ==1) & (table['_hsp_rank'] == 1) ]
#table = table[:1000]

table['_is.gapped'] = table[seqcol].apply(isgapped)
table['_I2L.equivalent.sequence']  = table[seqcol].apply(ILequivalence)


table2 = table[table['_is.gapped'] != '+']


peps = set(table2[seqcol].tolist())

pepstr = '\n'.join(peps)

w = open(outfolder +'/peptide_set.txt', 'w')
w.write(pepstr)
w.close()


cmd = "cat {}/peptide_set.txt | prot2pept | peptfilter | tr I L | sort -u  > {}/prot2pept.csv".format(outfolder, outfolder)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

cmd = "cat {}/prot2pept.csv | unipept pept2lca --equate > {}/pept2lca.csv".format(outfolder, outfolder)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0


cmd = "cat {}/prot2pept.csv | unipept pept2prot --equate > {}/pept2prot_all.csv".format(outfolder,  outfolder)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

processed = pd.read_csv(outfolder +'/pept2lca.csv').drop_duplicates()
p2 = processed.copy()

pcols = []

for col in processed.columns:
    processed.rename(columns={col : '_' + col},inplace=True)
    pcols.append('_' + col)

tcols = table.columns.tolist()

merged = pd.merge(table, processed, how='left', left_on='_I2L.equivalent.sequence', right_on='_peptide')
newcols = tcols + pcols

assert len(merged) == len(table)

merged = merged[newcols]

merged.to_csv(outfolder + '/merged_blast_unipept.tsv', sep='\t')

counted = Counter(merged['_taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/merged_blast_unipept_taxon_name_counts.tsv', sep='\t')

counted = Counter(processed['_taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2lca_taxon_name_counts.tsv', sep='\t')

processed_non_root = p2[p2['taxon_rank'] != 'no rank']
non_root_peptides = processed_non_root['peptide'].tolist()
w = open(outfolder + '/taxon_name_not_root_pep2lca_peptides.txt','w')
w.write('\n'.join(non_root_peptides))
w.close()

inferred= pd.read_csv(outfolder +'/pept2prot_all.csv').drop_duplicates()

inferred_non_root_peptides = inferred[inferred['peptide'].isin(non_root_peptides)]
non_root_acc_list = inferred_non_root_peptides['uniprot_id'].tolist()
inferred_excl_root_only = inferred[inferred['uniprot_id'].isin(non_root_acc_list)]

counted = Counter(inferred['uniprot_id'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2prot_counts_all.tsv', sep='\t')

#############################
# Excluding root accessions #
#############################


counted = Counter(inferred_excl_root_only['uniprot_id'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
acc_counts = pd.DataFrame(sorted_counted)
acc_counts.to_csv(outfolder + '/pept2prot_counts_excl_root_only_accessions.tsv', sep='\t')
all_acc = acc_counts[0].tolist()

w=open(outfolder + '/pept2prot_excl_root_only_accessions.txt','w')
w.write('\n'.join(all_acc))
w.close()

w=open(outfolder + '/pept2prot_excl_root_only_accessions_pep_export.txt','w')
filt=inferred[inferred['uniprot_id'].isin(all_acc)]
pfilt=list(set(filt['peptide'].tolist()))
w.write('\n'.join(pfilt))
w.close()

pfilt_lca = p2[p2['peptide'].isin(pfilt)]
pfilt_lca.to_csv(outfolder + '/pept2prot_excl_root_only_accessions_filtered_lca.tsv', sep ='\t')

counted = Counter(pfilt_lca['taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2prot_excl_root_only_accessions_filtered_lca_counts.tsv', sep='\t')

############################################################
# Excluding root accession only with more than one peptide #
############################################################

w = open(outfolder + '/pept2prot_excl_root_only_accessions_mult_peptides.txt','w')
one_acc = acc_counts[acc_counts[1] > 1]
one_acc = one_acc[0].tolist()
w.write('\n'.join(one_acc))
w.close()

w=open(outfolder + '/pept2prot_excl_root_only_accessions_mult_peptides_pep_export.txt','w')
filt=inferred[inferred['uniprot_id'].isin(one_acc)]
pfilt=list(set(filt['peptide'].tolist()))
w.write('\n'.join(pfilt))
w.close()
pfilt_lca = p2[p2['peptide'].isin(pfilt)]
pfilt_lca.to_csv(outfolder + '/pept2prot_excl_root_only_accessions_mult_peptides_filtered_lca.tsv', sep ='\t')

counted = Counter(pfilt_lca['taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2prot_excl_root_only_accessions_mult_peptides_filtered_lca_counts.tsv', sep='\t')


########################################
# Analysis excluding all root peptides #
########################################

counted = Counter(inferred_non_root_peptides['uniprot_id'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
acc_counts = pd.DataFrame(sorted_counted)
acc_counts.to_csv(outfolder + '/pept2prot_counts_excl_root_peptides.tsv', sep='\t')
all_acc = acc_counts[0].tolist()

w=open(outfolder + '/pept2prot_excl_root_peptides.txt','w')
w.write('\n'.join(all_acc))
w.close()

w=open(outfolder + '/pept2prot_excl_root_peptides_pep_export.txt','w')
pfilt=list(set(inferred_non_root_peptides['peptide'].tolist()))
w.write('\n'.join(pfilt))
w.close()

pfilt_lca = p2[p2['peptide'].isin(pfilt)]
pfilt_lca.to_csv(outfolder + '/pept2prot_excl_root_peptides_filtered_lca.tsv', sep ='\t')

counted = Counter(pfilt_lca['taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2prot_excl_root_peptides_filtered_lca_counts.tsv', sep='\t')

#############################################################################
# Analysis excluding all root peptides , require multiple non-root peptides #
#############################################################################

w = open(outfolder + '/pept2prot_excl_root_peptides_mult_peptides.txt','w')
one_acc = acc_counts[acc_counts[1] > 1]
one_acc = one_acc[0].tolist()
w.write('\n'.join(one_acc))
w.close()

w=open(outfolder + '/pept2prot_excl_root_peptides_mult_peptides_pep_export.txt','w')
filt=inferred_non_root_peptides[inferred_non_root_peptides['uniprot_id'].isin(one_acc)]
pfilt=list(set(filt['peptide'].tolist()))
w.write('\n'.join(pfilt))
w.close()
pfilt_lca = p2[p2['peptide'].isin(pfilt)]
pfilt_lca.to_csv(outfolder + '/pept2prot_excl_root_peptides_mult_peptides_filtered_lca.tsv', sep ='\t')

counted = Counter(pfilt_lca['taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2prot_excl_root_peptides_mult_peptides_filtered_lca_counts.tsv', sep='\t')

#cmd = "up_acc2fasta.pl {}/pept2prot_excl_root_only_accessions.txt > {}/pept2prot_excl_root_only_accessions.fasta".format( outfolder, outfolder)
#process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
#process.wait()
#assert process.returncode == 0

cmd = "cat {}/pept2prot_excl_root_only_accessions_mult_peptides.txt | sort -u | uniprot --format fasta > {}/pept2prot_excl_root_only_accessions_mult_peptides.fasta".format( outfolder, outfolder)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0


cmd = "cat {}/pept2prot_excl_root_peptides_mult_peptides.txt | sort -u | uniprot --format fasta > {}/pept2prot_excl_root_peptides_mult_peptides.fasta".format( outfolder, outfolder)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

#counted = Counter(inferred_non_root_peptides['uniprot_id'])
#sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
#counts = pd.DataFrame(sorted_counted)
#counts.to_csv(outfolder + '/pept2prot_counts_excl_all_root_peptides.tsv', sep='\t')


w = open(outfolder + '/methods.txt','w')
w.write('Explanantion of unipept analysis and files included.\n\n')
w.write('peptide_set.txt is by default the export of hsp.sbjct sequences from the BLAST results in this pipeline (excluding sequences that are gapped alignemnts).\nWhen using the csv2unipeptlca.py  script in standalone  mode the required field with the peptide sequences needs to be specified .\n\n')
w.write('prot2pept.csv is the unipept prot2pept output of the peptide_set.txt output (in silico tryptic digest).\n\n')
w.write('pept2lca.csv is the unipept pept2lca analysis of the prot2pept.csv file.\n\n')
w.write('pept2prot_all.csv is the unipept pept2prot analysis of prot2pept.csv.\n\n')
w.write('pept2prot_counts_all.csv is the number of peptides per accession from the pep2prot_all.csv results file.\n\n')
w.write('pept2prot_excl_root_only_accessions.txt is the list of accessions that are mapped by at least one peptide that does not have "no rank" in the unipept pep2lca output "taxon_rank" column.\n\n')
w.write('pept2prot_excl_root_only_accessions.fasta is the fasta entries for the above-mentioned file.\n\n')







