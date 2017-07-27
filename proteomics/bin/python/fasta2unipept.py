#!/usr/bin/env python

import subprocess
import sys
import pandas as pd
import collections
from collections import Counter
import operator
import os
import Bio
from Bio import SeqIO

infile = sys.argv[1]
outfolder = sys.argv[2]

fasta = list(SeqIO.parse(infile, 'fasta'))

peps = [ str(i.seq) for i in fasta ]

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


cmd = "cat {}/prot2pept.csv | unipept pept2taxa --all --equate > {}/pept2taxa_all.csv".format(outfolder,  outfolder)
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
assert process.returncode == 0

pept2taxa = pd.read_csv('{}/pept2taxa_all.csv'.format(outfolder)).drop_duplicates()

counted = Counter(pept2taxa['taxon_name'])
sorted_counted = sorted(counted.items(), key=operator.itemgetter(1), reverse=True)
counts = pd.DataFrame(sorted_counted)
counts.to_csv(outfolder + '/pept2taxa_taxon_name_counts.csv')


processed = pd.read_csv(outfolder +'/pept2lca.csv').drop_duplicates()
p2 = processed.copy()

pcols = []

for col in processed.columns:
    processed.rename(columns={col : '_' + col},inplace=True)
    pcols.append('_' + col)

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

inferred_counts = counts
inferred_counts.rename(columns={0:"_Protein.Accession",1:"_Protein.Accession.Tryptic.Peptide.Count"}, inplace=True)

############
## Mergedf pept2taxa, pept2prot, pep2prot_counts

def combined_table(pept2taxa, pept2prot_counts, pept2prot):
     pept2taxa_df = pept2taxa.copy()
     pept2prot_counts_df = pept2prot_counts.copy()
     pept2prot_df = pept2prot.copy()
     del pept2taxa_df['peptide']
     pept2taxa_df = pept2taxa_df.drop_duplicates()
     pept2prot_counts_df.rename(columns={'1':'_Protein.Tryptic.Peptide.Counts' , '0': '_Protein.Accession'}, inplace=True) 
     del pept2prot_df['peptide'] 
     pept2prot_df = pept2prot_df.drop_duplicates() 
     merged_df = pd.merge(pept2prot_df, pept2taxa_df, how = 'left', left_on='taxon_id', right_on='taxon_id') 
     merged_df = pd.merge(merged_df, pept2prot_counts_df, how='left', left_on='uniprot_id', right_on='_Protein.Accession') 
     newcols = [ i for i in merged_df.columns if not i.startswith('Unnamed: ')]
     
     merged_df = merged_df[newcols]

     merged_df.to_csv(outfolder + '/merged_pept2prot_pept2taxa_pept2protCounts.csv')
combined_table(pept2taxa, inferred_counts, inferred)














