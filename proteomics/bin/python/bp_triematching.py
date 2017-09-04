#!/usr/bin/env python

import collections
from collections import defaultdict
import pandas as pd
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import algo
import re
import tagmatch
import proteomics
import pickle
import numpy as np
import sqlite3
import threading
import time
import multiprocessing as mp
from multiprocessing import Manager
import datrie

NCORE = 8

db = SeqIO.parse(sys.argv[1],'fasta')
print('loaded db')

Trie = pickle.load(open(sys.argv[2],'rb'))
print('loaded trie')



def process(q, iolock):
    global Trie    
    count = 0
    while True:
        rec = q.get()
        if rec is None:
            break
        count += 1
        
        if count % 1000 == 0:
            print(count)
        
        ID = rec.id
        SEQ = str(rec.seq)
        res = Trie.prefixes(SEQ)
        #TrieMatch=algo.TrieMatch(Trie, SEQ)
        #peptides_mapped=TrieMatch.trie_export()

q = mp.Queue(maxsize=NCORE)

iolock = mp.Lock()

pool = mp.Pool(NCORE, initializer=process, initargs=(q, iolock))

for rec in db:
    q.put(rec)  # blocks until q below its max size

for _ in range(NCORE):  # tell workers we're done
    q.put(None)

pool.close()
pool.join()









