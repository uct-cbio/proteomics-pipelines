#!/usr/bin/env python

import sys
import os
import Bio; from Bio import SeqIO
#sys.path.append(os.getcwd())
import sequtils

if len(sys.argv)==3:
    infile  = sys.argv[1]
    assembly_name = sys.argv[2]
    outfile = sys.stdout


seqs = list(SeqIO.parse(infile, 'fasta'))
sf = sequtils.sf_contigs(seqs, assembly_name = assembly_name, table=11, codons='All', peptide_length=20, translated=False)

SeqIO.write(sf, outfile, 'fasta')
