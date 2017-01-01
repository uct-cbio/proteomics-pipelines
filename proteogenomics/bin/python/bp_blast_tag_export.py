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

outfolder = sys.argv[3]


export=[]

def tag_export(df):
    global export
    subject=df['hsp.sbjct']
    ids =df['_query.description'].split()[1].split('scans=')[1].split('|')
    for id in ids:
        print('querying')
        cursor.execute("select rec from tags where scanid=?",(id,))
        recs = [pickle.loads(i[0]) for i in cursor.fetchall()]
        print(recs)
        tm = tagmatch.blast_tags(subject, recs)
        export += tm.newrecords

def fixsubject(val):
    try:
        return ''.join(val.split('-'))
    except:
        return ''


#table = table[table['hsp.match']=='+DNDLR']
#print(table[['hsp.sbjct','hsp.match', 'hsp.query']])
table = table[table['hsp.sbjct'].notnull()]
table['hsp.sbjct'] = table['hsp.sbjct'].apply(fixsubject)
table.apply(tag_export, axis=1)

SeqIO.write(export, sys.argv[3], 'fasta' )
