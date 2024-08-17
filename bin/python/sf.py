#!/usr/bin/env python3

import sys
import os
import Bio; from Bio import SeqIO
#sys.path.append(os.getcwd())
import sequtils
import yaml

infile  = sys.argv[1]

assembly_name = sys.argv[2]

config = yaml.load(open(sys.argv[3]).read(), Loader=yaml.Loader)

outfile = sys.argv[4]


translation_table = config['translation_table']

start_codons = config['start_codons']

method = config['method']

minimum_translation_length = config['minimum_translation_length']

possible_methods=['stop-to-stop','most upstream start codon']

assert method in possible_methods

if method == 'stop-to-stop':
    method_codons = 'All'
elif method == 'most upstream start codon':
    method_codons = start_codons

seqs = list(SeqIO.parse(infile, 'fasta'))

sf = sequtils.sf_contigs(seqs, assembly_name = assembly_name , table=translation_table, codons=method_codons, peptide_length=minimum_translation_length, translated=False)

if config['alternative_tss_entries'] == True:
    sf = sequtils.alt_starts_recs(sf, starts=start_codons, table= translation_table, translated=False )

nucs = []
prots = []

for rec in sf:
    if config['methionine_start'] == True:
        prot= sequtils.translate_start(rec, table=translation_table, starts=start_codons)
    else:
        prot= sequtils.translate_start(rec, table=translation_table, starts=[])
    
    if len(str(prot.seq)) >= minimum_translation_length:
        prots.append(prot)
        nucs.append(rec)

SeqIO.write(nucs,  outfile + '_sixframe_nuc.fasta' , 'fasta')
SeqIO.write(prots, outfile + '_sixframe_prot.fasta', 'fasta')


