#!/usr/bin/python

import Bio; from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

f = list(SeqIO.parse(sys.stdin, 'fasta'))

mlist = []
for i in f:
    desc = i.description.split()[1]
    id = i.id
    seq = i.seq
    mseq = 'M' + str(i.seq)[1:]
    new_seq = SeqRecord(seq = Seq(mseq), id= id, description = desc)
    mlist.append(new_seq)

SeqIO.write(mlist, sys.stdout, 'fasta')

