#!/usr/bin/env python

import sys
import Bio; from Bio import SeqIO
from Bio.Seq import Seq
import algo
import pickle
import sqlite3

db = sqlite3.connect(sys.argv[1])
cursor = db.cursor()
table = sys.argv[2]
column = sys.argv[3]

cursor.execute("select DISTINCT {} from {}".format(column, table))

seqs= [i[0] for i in cursor.fetchall()]

#seqs = [str(i.seq) for i in list(SeqIO.parse(sys.argv[1], 'fasta'))]
Trie = algo.Trie(seqs)

pickle.dump(Trie, open(sys.argv[4], 'wb'))
