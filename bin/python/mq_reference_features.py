#!/usr/bin/env python3
import urllib.parse
import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
import yaml
import blast
from Bio.Blast import NCBIXML
import gff3
import re
config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.Loader)
output = sys.argv[2]

query_fasta = list(SeqIO.parse(output +'/strains/all_mapped_trans_orfs.fasta','fasta'))

strain_samples = defaultdict(list)

design = pd.read_csv(config['design'])

if 'rename' in design.columns:
    scol = 'rename'
else:
    scol = 'sample'
assert scol in design.columns
design = design[design['exclude'] != '+' ] 

rename_dict = sequtils.rename_design(design)

sample2strain = design.set_index(scol)['Strain'].to_dict()
samples = list(sample2strain.keys())

#print(samples)
#print(rename_dict)

#for sample in samples:
#    st = sample2strain[sample]
#    strain_samples[st].append(rename_dict[sample])

#peptide_sequence  = SeqIO.to_dict(list(SeqIO.parse(output + '/blast/peptides2orfs/peptides.fasta','fasta')))
# Per strain analysis
#peptides = pd.read_csv(config['mq_txt'] + '/peptides.txt', sep='\t',engine='python')
#peptides = peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]

#all_peptides = peptides['Sequence'].tolist()
#exp_cols = [i for i in peptides.columns if i.startswith('Experiment')]

# dynamic programming approach to limit the data set for dev


#peptide_annotations =pd.read_csv(output +'/annotations/H37Rv_S5527_peptide_annotations_groups.csv')
#print(peptide_annotations.head(1).stack())
#print(list(query_fasta)[0].format('fasta'))
#quit()

def load_gff3(file):
    # Load GFF3 file, skipping comment lines (lines starting with '#')
    df = pd.read_csv(
        file, 
        sep="\t", 
        comment="#", 
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    )
    return df


for ref in config['reference']:
    assembly_id = config['reference'][ref]['assembly_id']
    for strain in config['strains']:
          
        print("Reference: ", ref, strain)
        print("Loading peptides")
        ref_peptides = load_gff3(output + '/jbrowse/' + strain + '/{}_{}_peptides.gff3'.format(strain, assembly_id))
        #print(ref_peptides.head(1))
        
        print("Loading ORFs")
        ref_orfs = load_gff3(output + '/jbrowse/' + strain + '/{}_{}_orfs.gff3'.format(strain, assembly_id))
        seqid_set = ref_orfs['seqid'].unique()
        assert len(seqid_set) == 1
        seqid = seqid_set[0]
        
        peptide_annotations = pd.read_csv(output +'/annotations/{}_{}_peptide_annotations_groups.csv'.format(ref, strain))
        # get subset of peptide_annotations where VarId is null
        peptide_annotations_null = peptide_annotations[peptide_annotations['VarId'].isnull()]
        assert len(peptide_annotations_null) > 0
        peptide_annotations_dict = peptide_annotations.set_index(['ORF_id', 'PeptideSequence'])['StartType'].to_dict()
        
        peptide_acetylations_dict = peptide_annotations.set_index(['PeptideSequence'])['Strain_Nterm_Acetylated'].to_dict()

        #peptide_annotations = peptide_annotations[peptide_annotations['AnnotationType'] == 'Upstream non-TSS peptide with putative upstream TSS']
        # Step 1: Identify duplicates based on the combination of col1 and col2
        peptide_annotations['is_duplicate'] = peptide_annotations.duplicated(subset=['ORF_id', 'StartPosition'], keep=False)

        # Step 2: Define if var_id is novel (not a duplicate and ref_id is NaN)
        peptide_annotations['is_novel'] = (~peptide_annotations['is_duplicate']) & (peptide_annotations['Ref_id'].isna())

        #peptide_annotations = peptide_annotations[peptide_annotations['VarId'].notnull()]
        #peptide_annotations = peptide_annotations[~peptide_annotations['VarId'].str.contains('.None')]
        #peptide_annotations= peptide_annotations[peptide_annotations['VarId'] == 'S5527_scaffold11_size137936|S5527_scaffold11_size137936_recno_1896.0.P9WMH1.1|(+)40371:41063']
        #print(peptide_annotations['AnnotationType'].unique())
        #print(peptide_annotations.head(1).stack())
        #quit()
        try:
            del peptide_annotations['Unnamed: 0']
        except:
            pass
        ref_peptides['PeptideSequence'] = ref_peptides['attributes'].str.extract('ID=(.*?);')
        ref_peptides= ref_peptides[ref_peptides['PeptideSequence'].notnull()]
        
        ref_orfs['ID'] = ref_orfs['attributes'].str.extract('ID=(.*?);')
        ref_orfs['Parent'] = ref_orfs['attributes'].str.extract('Parent=(.*?);')
        sub_orfs = ref_orfs[ref_orfs['Parent'].notnull()] 
        ref_orfs = ref_orfs[ref_orfs['ID'].notnull()]

        # Regular expression to find occurrences of 'match<anynumber>-string'
        pattern = r"match\d+-"

        # remove the string "match1-" from the start of ORF_id - but it can be any number
        ref_orfs['ORF_id'] = ref_orfs['ID'].apply(lambda x : re.sub(pattern, '', x))
        ref_peptides['PeptideSequence'] = ref_peptides['PeptideSequence'].apply(lambda x : re.sub(pattern, '', x))
        
        rows = ['##gff-version 3']
        feature_rows = ['##gff-version 3']
        seen = []

        for orf, group in peptide_annotations.groupby('ORF_id'):
           #orf = group['ORF_id'].values[0]
           
           #var_start = group['PeptidePosition'].min()
           
           #var_peps = peptide_annotations[peptide_annotations['ORF_id'] == orf]
           #var_pe

           
           blastp_list  = group['BLASTP'].values
           assert len(set(blastp_list)) == 1
           blastp = blastp_list[0]

           orf_id = orf.split('|')[1]
           #print(group) 
           orf_filt = ref_orfs[ref_orfs['ORF_id'] == orf_id]
           if len(orf_filt) == 0:
               continue
           feat_id = orf_filt['ID'].values[0]
           
           sub_orf_filt = sub_orfs[sub_orfs['Parent'] == feat_id]
           #seen=False 
           for mapped_orf in orf_filt.iterrows():
               orf_match_id = mapped_orf[1]['ID']
               orf_start = mapped_orf[1]['start']
               orf_end = mapped_orf[1]['end']
               orf_strand = mapped_orf[1]['strand']

               # upoad orf
               _source='ProteogenomicsPipeline'
               _type='proteoform_orf'
               _start=orf_start
               _end=orf_end
               _score='.'
               _strand=orf_strand
               _phase='.'
               _name=orf_id 

               _attributes1='ID={};Name={};blastp={}'.format(orf_match_id, orf_id, blastp)
               
               gff_row= [seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes1]
               gff_row  = '\t'.join([str(i) for i in gff_row])
               feature_rows.append(gff_row)            
               _attributes1='Parent={};Name={};blastp={}'.format(orf_match_id, orf_id, blastp)
               _type='open_reading_frame'
               gff_row= [seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes1]
               gff_row  = '\t'.join([str(i) for i in gff_row])
               feature_rows.append(gff_row)            
               orf_peptides = [] 
               # iterate through the proteoform peptides
               for row in group.iterrows():
                   var_id = row[1]['VarId']
                     
                   if isinstance(var_id, str):
                       var_no = var_id.split('|')[-2].split('.')[-1]
                   else:
                       assert np.isnan(var_id)
                       var_id= None
                       var_no=None

                   gene_name = row[1]['Gene Name']
                   protein_name = row[1]['Protein Name']
                   gene_name = row[1]['Gene Name']
                   blastp = row[1]['BLASTP']
                   annotation_type = row[1]['AnnotationType']
                   ref_id  = row[1]['Ref_id']
                   
                   try:
                       _ = np.isnan(ref_id)
                       assert _ == True
                       ref_var = False
                   except:
                        ref_var = True

                   start_pos = row[1]['StartPosition']
                   
                   novel = row[1]['is_novel']

                   peptide_pos = row[1]['PeptidePosition'] * 3

                   delta = peptide_pos - start_pos
                   # start peptide must be at max one codon downstream of the start position
                   #assert 3 >= delta >= 0
                   #print('delta', delta)
                   seq = row[1]['PeptideSequence']

                   # make sure the exact peptide occurs in the attributes, not just as a prefix
                   pep_filt = ref_peptides[ref_peptides['PeptideSequence'] == seq]
                   for prow in pep_filt.iterrows():
                       #print(prow)
                       pep_start = prow[1]['start']
                       pep_end = prow[1]['end']
                       pep_strand = prow[1]['strand']
                       #print("ORF: ", orf_start, orf_end, orf_strand)
                       #print("Peptide: ", pep_start, pep_end, pep_strand)
                       #print("Peptide: ", seq)
                       
                       if (pep_start >= orf_start) and (pep_end <= orf_end):
                           if orf_strand == pep_strand:
                               #print("Peptide in ORF")
                               if var_no is None:
                                   # not a variant
                                   proteoform_start = -1
                                   proteoform_end = None
                               elif orf_strand == '+':
                                   proteoform_start =int( pep_start - delta)
                                   proteoform_end =int( orf_end)
                               elif orf_strand == '-':
                                   proteoform_start =int( orf_start)
                                   proteoform_end =int( pep_end + delta)
                               else:
                                   assert False
                               #print("proteoform coordinates: ",proteoform_start, proteoform_end)
                               proteoform_strand = orf_strand
                               
                               for sgroup in sub_orf_filt.iterrows():
                                   sub_start = sgroup[1]['start']
                                   sub_end = sgroup[1]['end']
                                   sub_strand = sgroup[1]['strand']
                                   #print("Sub ORF: ", sub_start, sub_end, sub_strand)
                                   # the end of the matched blast feature, is at or downstream of the start of the peptide, and the start of the feature is downstream or at the calculated TSS site, meaning the TSS site corresponds between reference and proteoform
                                    
                                   if (sub_start <= pep_start) and (sub_end >= pep_start):
                                       # peptide feature
                                       _peptide_annotation = peptide_annotations_dict[(orf, seq)].lower()
                                       _peptide_acetylated = str(peptide_acetylations_dict[seq]).lower()
                                       assert _peptide_acetylated in ['true', 'false']
                                       if _peptide_acetylated == '+':
                                           _peptide_acetylated = 'true'
                                       elif _peptide_acetylated == '-':
                                           _peptide_acetylated = 'false'
                                       
                                       # feature rows
                                       _source='ProteogenomicsPipeline'
                                       _type='polypeptide'
                                       _start=pep_start
                                       _end=pep_end
                                       _score='.'
                                       _strand=pep_strand
                                       _phase='.'
                                       _attributes='Parent={};Name={};start_type={};nterm_acetylated={}'.format(orf_match_id, seq, _peptide_annotation, _peptide_acetylated)
                                       gff_row= [seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes]
                                       gff_row = '\t'.join([str(i) for i in gff_row])
                                       if not seq in orf_peptides:
                                           feature_rows.append(gff_row)      
                                           orf_peptides.append(seq)
                                       
                                       proteoform_valid=False
                                       if (sub_start <= proteoform_start) and (sub_end >= pep_start):
                                           # this peptide linked to a proteoform
                                           proteoform_valid=True
                                           _ID = orf_match_id + '_' + var_id
                                           if not _ID in seen:
                                               _source='ProteogenomicsPipeline'
                                               _type='proteoform_orf'
                                               _start=proteoform_start
                                               _end=proteoform_end
                                               _score='.'
                                               _strand=proteoform_strand
                                               _phase='.'
                                               if ref_var:
                                                   _name=ref_id + '_' + var_no + ' - ' + annotation_type
                                               else:
                                                   _name=orf_id + '_' + var_no + ' - ' + annotation_type

                                               _attributes1='ID={};Name={};Note={};gene_name={};protein_name={};blastp={};annotation_type={}'.format(_ID, _name, annotation_type, gene_name, protein_name, blastp, annotation_type)
                                               _attributes2='Parent={};Name={};Note={};gene_name={};protein_name={};blastp={};annotation_type={}'.format(_ID, _name, annotation_type, gene_name, protein_name, blastp, annotation_type)
                                               
                                               gff_row1= [seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes1]
                                               gff_row1  = '\t'.join([str(i) for i in gff_row1])
                                               _type='proteoform'
                                               gff_row2= [seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes2]
                                               gff_row2  = '\t'.join([str(i) for i in gff_row2])
                                               if ref_var or novel:
                                                   rows.append(gff_row1)            
                                                   rows.append(gff_row2)            
                                               seen.append(_ID)
                                           if ref_var or novel:
                                               # upload the peptide
                                               _source='ProteogenomicsPipeline'
                                               _type='polypeptide'
                                               _start=pep_start
                                               _end=pep_end
                                               _score='.'
                                               _strand=pep_strand
                                               _phase='.'
                                               _name= seq
                                               _attributes='Parent={};Name={};start_type={};nterm_acetylated={}'.format(orf_match_id + '_' + var_id, _name, _peptide_annotation, _peptide_acetylated)
                                               gff_row= [seqid, _source, _type, _start, _end, _score, _strand, _phase, _attributes]
                                               gff_row = '\t'.join([str(i) for i in gff_row])
                                               rows.append(gff_row)      

        gff3 = '\n'.join(rows)
        if not os.path.exists(output + '/jbrowse/' + strain):
            os.makedirs(output + '/jbrowse/' + strain )
        w = open(output + '/jbrowse/' + strain + '/{}_{}_proteoforms.gff3'.format(strain, \
                                                                                  assembly_id), \
                 'w')
        w.write(gff3)
        w.close()
        
        gff3 = '\n'.join(feature_rows)
        if not os.path.exists(output + '/jbrowse/' + strain):
            os.makedirs(output + '/jbrowse/' + strain )
        w = open(output + '/jbrowse/' + strain + '/{}_{}_features.gff3'.format(strain, \
                                                                                  assembly_id), \
                 'w')
        w.write(gff3)
        w.close()
print("Annotation of proteoforms complete")
