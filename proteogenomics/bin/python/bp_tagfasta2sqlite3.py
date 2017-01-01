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
    id = rec.id
    scanid = ';'.join(id.split(';')[:2])
    seq = str(rec.seq)
    load.append((scanid, id, seq, sqlite3.Binary(pickle.dumps(rec,pickle.HIGHEST_PROTOCOL)), source))

db = sqlite3.connect(sys.argv[2])
cursor = db.cursor()
cursor.execute("CREATE TABLE IF NOT EXISTS tags (scanid TEXT, tagid TEXT, seqstr TEXT, rec BYTEA, source TEXT, UNIQUE( tagid, seqstr, source ) )")
db.commit()   # This creates table, makes sure that tagid, seqstr and source are unique

sql = ("CREATE INDEX IF NOT EXISTS index_seqs ON tags (seqstr);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_tagids ON tags (tagid);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_scanids ON tags (scanid);")
cursor.execute(sql)
db.commit()

sql = "INSERT OR IGNORE INTO tags VALUES(?, ?, ?, ?, ?)"
cursor.executemany(sql, load )
db.commit()
    

