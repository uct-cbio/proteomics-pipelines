#!/usr/bin/env python3

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
import os

<<<<<<< HEAD
NCORE = 8
=======
NCORE = mp.cpu_count() 

>>>>>>> e1aa09c7f60c974663baec27314a288147acf6fe
db = SeqIO.parse(sys.argv[1],'fasta')
print('loaded db')

Trie = pickle.load(open(sys.argv[2],'rb'))
print('loaded trie')
database = sys.argv[3]

prectol = float(os.environ["tm_prec_tol"])
specificity=os.environ["tm_specificity"]
enzymes = [i for i in os.environ["tm_enzymes"].split(';') if i != '']
max_missed_cleavages = int(os.environ["tm_max_missed_cleavages"])
gap_tol = float(os.environ["tm_gap_tol"])
fixed_modifications = [i for i in os.environ["tm_fixed_modifications"].split(';') if i != '']
variable_modifications = [i for i in os.environ["tm_variable_modifications"].split(';') if i != '']
sqldb = sqlite3.connect(database)
cursor = sqldb.cursor()
tag_df = pd.read_sql("SELECT DISTINCT scanid, sequence, ngap, cgap from tags", sqldb,columns=['scanid', 'sequence', 'ngap','cgap'])
tag_df = tag_df.drop_duplicates()

seqs = set(tag_df['sequence'].tolist())
tagids = {}

for seq in seqs:
    tagids[seq] = defaultdict(set)
del seqs

def create_dct(table):
    global tagids
    sid = table['scanid']
    seq = table['sequence']
    ngap = table['ngap']
    cgap = table['cgap']
    mw = tagmatch.peptide_mass(seq, fixed_modifications=[], variable_modifications=[], nterm=False, cterm=False)[0]
    tagids[seq][(ngap, mw, cgap)].add(sid)

tag_df.apply(create_dct, axis=1)
del tag_df

print('loaded scans')

def process(q, out, iolock):
    global Trie
    global database
    count = 0
    while True:
        rec = q.get()
        if rec is None:
            out.put(None)
            break
        count += 1
<<<<<<< HEAD
        
        if count % 10000 == 0:
            print(count)
     
        ID = rec.id
        SEQ = str(rec.seq)
        TrieMatch=algo.TrieMatch(Trie, SEQ)
        peptides_mapped=TrieMatch.trie_export()
        
        #TrieMatch=algo.TrieMatch(Trie, SEQ)
        #peptides_mapped=TrieMatch.trie_export()
=======
        ID = rec.id
        SEQ = 'L'.join(str(rec.seq).split('I'))
        TrieMatch=algo.TrieMatch(Trie, SEQ)
        peptides_mapped=list(TrieMatch.trie_export())
        for pepstr in peptides_mapped:
            tseq = pepstr
            prec_ions = tagids[pepstr]
            prec_array = list(prec_ions.keys())
            tm=tagmatch.TagMatch(pepstr, 
                    prec_array, 
                    SEQ, 
                    prec_tol = prectol, 
                    specificity=specificity, 
                    enzymes=enzymes, 
                    max_missed_cleavages=max_missed_cleavages, 
                    gap_tol=gap_tol , 
                    fixed_modifications=fixed_modifications, 
                    variable_modifications=variable_modifications)
        #cursor.execute('SELECT scanid FROM tags WHERE sequence in ({0})'.format(', '.join('?' for _ in peptides_mapped)), peptides_mapped)
        #scanids = cursor.fetchall()
        res = ()
        out.put(res)

def export(out, iolock):
    _ = sqlite3.connect(database)
    cursor = _.cursor()
    count = 0
    close = 0
    while True:
        res = out.get()
        if res is None:
            close += 1
            if close == NCORE:
                break
        count +=1
        if count % 10000 == 0:
            print(count)
>>>>>>> e1aa09c7f60c974663baec27314a288147acf6fe

q = mp.Queue(maxsize=NCORE)
out = mp.Queue(maxsize=NCORE)
iolock = mp.Lock()
pool = mp.Pool(NCORE, initializer=process, initargs=(q, out, iolock))
pool2 = mp.Pool(1, initializer=export, initargs=(out, iolock))
for rec in db:
    q.put(rec)  # blocks until q below its max size
for _ in range(NCORE ):  # tell workers we're done
    q.put(None)
pool.close()
pool2.close()
pool.join()
pool2.join()





