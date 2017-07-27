#!/usr/bin/env python

import sys; import Bio; from Bio import SeqIO
import sequtils

if len(sys.argv) == 1:
	samfile = sys.stdin.read()
else:
	samfile = sys.argv[1]
	sam = open(samfile)
	samfile = sam.read()
	sam.close()

seqs = sequtils.sam_parse(samfile)
SeqIO.write(seqs,sys.stdout,'fasta')
