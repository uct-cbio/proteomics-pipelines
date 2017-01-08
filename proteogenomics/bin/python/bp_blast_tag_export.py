#!/usr/bin/env python

import pandas as pd
import sys
import pickle
import tagmatch
from Bio import SeqIO
import sqlite3

table = pd.read_csv(sys.argv[1])

db = sqlite3.connect(sys.argv[2])

cursor = db.cursor()
    
min_identities = int(sys.argv[3])

export=set()

def tag_export(df):
    global export
    subject=df['hsp.sbjct']
    query= df['_query.sequence']
    #ids =df['_query.description'].split()[1].split('scans=')[1].split('|')
    res = cursor.execute("select scanid from combined where seqstr=?",(query,)).fetchall()
    ids =set([i[0] for i in res])
    
    for id in ids:
        res = cursor.execute("select rec from tags where (scanid=? AND source='denovo')",(id,)).fetchall()
        recs = [pickle.loads(i[0]) for i in res]
        tm = tagmatch.blast_tags(subject, recs, min_identities=min_identities)
        for i in tm.newrecords:
            print(i.format('fasta'))
        export.update([i.format('fasta') for i in tm.newrecords])

def fixsubject(val):
    try:
        return ''.join(val.split('-'))
    except:
        return ''


#table = table[table['hsp.match']=='+DNDLR']
#print(table[['hsp.sbjct','hsp.match', 'hsp.query']])

table = table[table['hsp.sbjct'].notnull()]

table['hsp.sbjct'] = table['hsp.sbjct'].apply(fixsubject)

table = table.drop_duplicates(['_query.sequence','hsp.sbjct'])

table.apply(tag_export, axis=1)

w = open(sys.argv[4], 'w')
w.write(''.join(export))
w.close()
