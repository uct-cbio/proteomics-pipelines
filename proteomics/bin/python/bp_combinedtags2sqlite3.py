#!/usr/bin/env python

import sys
import sqlite3
import Bio
from Bio import SeqIO
import pickle

recs = SeqIO.parse(sys.argv[1],'fasta')
source = sys.argv[3]

load = []
for rec in recs:
    scanid = rec.id
    seq = str(rec.seq)
    load.append((scanid, seq, sqlite3.Binary(pickle.dumps(rec,pickle.HIGHEST_PROTOCOL)), source))

db = sqlite3.connect(sys.argv[2])
cursor = db.cursor()

cursor.execute("CREATE TABLE IF NOT EXISTS combined ( scanid TEXT, seqstr TEXT, rec BYTEA, source TEXT );")
#cursor.execute("CREATE TABLE IF NOT EXISTS combined (scanid TEXT, seqstr TEXT, rec BYTEA, source TEXT, UNIQUE( scanid, seqstr, source ) )")

db.commit()   # This creates table, makes sure that scanid, seqstr and source are unique

sql = ("DROP INDEX IF EXISTS index_combined_seqs;")
cursor.execute(sql)
db.commit()

sql = ("DROP INDEX IF EXISTS index_combined_scanids;")
cursor.execute(sql)
db.commit()

sql = "INSERT INTO combined VALUES(?, ?, ?, ?)"
cursor.executemany(sql, load )
db.commit()
    
sql = ("CREATE INDEX IF NOT EXISTS index_combined_seqs ON combined (seqstr);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_combined_scanids ON combined (scanid);")
cursor.execute(sql)
db.commit()
