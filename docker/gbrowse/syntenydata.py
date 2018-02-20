#!/usr/bin/env python

import Bio; from Bio import SeqIO
import sys
from Bio.Seq import Seq

genome = list(SeqIO.parse(sys.argv[1],'fasta'))

new = []

for contig in genome:
    cid = contig.id.split('|')[0]
    seq = str(contig.seq).lower()
    newcid = cid
    contig.id = newcid
    contig.description=''
    contig.seq=Seq(seq)
    new.append(contig)

SeqIO.write(new, sys.argv[2]+'.fasta','fasta')
