#!/usr/bin/env python

import sys
import sqlite3
import Bio
from Bio import SeqIO
import pickle
import tagmatch
import pandas as pd
import numpy as np

db = sqlite3.connect(sys.argv[1])
cursor = db.cursor()
prec_tol = float(sys.argv[2])
gap_tol = float(sys.argv[3])
fixed_mods = [i for i in sys.argv[4].split(';') if i!= '']
variable_mods = [i for i in sys.argv[5].split(';') if i!= '']

limit = prec_tol + gap_tol + gap_tol


#cursor.execute("CREATE TABLE IF NOT EXISTS tags (scanid TEXT, tagid TEXT, seqstr TEXT, rec BYTEA, source TEXT, UNIQUE( tagid, seqstr, source ) )")

tag_df = pd.read_sql("SELECT DISTINCT tagid, seqstr from tags", db, columns=['tagid', 'seqstr'])
tag_df = tag_df.drop_duplicates()
print(len(tag_df))
load = []

def validate(table):
    global load
    tid = table['tagid']
    _ = tid.split(';')
    mw = np.round(float(_[2].split('mw=')[1]),4)
    ngap = np.round(float(_[3].split('ngap=')[1]),4)
    cgap = np.round(float(_[4].split('cgap=')[1]),4)
    seq = table['seqstr']
    for mass in tagmatch.peptide_mass(seq, fixed_modifications=fixed_mods, variable_modifications = variable_mods, nterm=True, cterm=True):
        total = mass + ngap + cgap
        diff = np.absolute(total - mw)
        if diff < limit:
            load.append((tid, seq))

tag_df.apply(validate, axis=1)
del tag_df

cursor.execute("CREATE TABLE IF NOT EXISTS valid (tagid TEXT, seqstr TEXT);")

db.commit()   # This creates table, makes sure that tagid, seqstr and source are unique

sql = ("DROP INDEX IF EXISTS index_seqs_valid;")
cursor.execute(sql)
db.commit()

sql = ("DROP INDEX IF EXISTS index_tagids_valid;")
cursor.execute(sql)
db.commit()

sql = ("DROP INDEX IF EXISTS index_tagids_seqs_valid;")
cursor.execute(sql)
db.commit()

sql = "INSERT OR IGNORE INTO valid VALUES(?, ?);"
cursor.executemany(sql, list(set(load)) )
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_seqs_valid ON valid (seqstr);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_tagids_valid ON valid (tagid);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_tagids_seqs_valid ON valid (tagid, seqstr);")
cursor.execute(sql)
db.commit()
print(len(list(set(load))))
