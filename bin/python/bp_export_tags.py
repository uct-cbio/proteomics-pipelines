#!/usr/bin/env python3

import sqlite3
import sys

output = sys.argv[1]
db = sys.argv[2]

db = sqlite3.connect(sys.argv[2])
cursor = db.cursor()

sql = ("DROP INDEX IF EXISTS index_seqs;")
cursor.execute(sql)
db.commit()

sql = ("DROP INDEX IF EXISTS index_tags;")
cursor.execute(sql)
db.commit()

sql = ("DROP INDEX IF EXISTS index_scanids;")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_seqs ON tags (sequence);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_tags ON tags (tag);")
cursor.execute(sql)
db.commit()

sql = ("CREATE INDEX IF NOT EXISTS index_scanids ON tags (scanid);")
cursor.execute(sql)
db.commit()

cursor.execute("select DISTINCT tag from tags")
tags= [i[0] for i in cursor.fetchall()]

w = open(output, 'w')

w.write('\n'.join(tags))

w.close()

