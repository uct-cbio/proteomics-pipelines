#!/usr/bin/env python

import sys
import Bio; from Bio import SeqIO
from Bio.Seq import Seq
import algo
import pickle

seqs = [str(i.seq) for i in list(SeqIO.parse(sys.argv[1], 'fasta'))]
Trie = algo.Trie(seqs)

pickle.dump(Trie, open(sys.argv[2], 'wb'))
