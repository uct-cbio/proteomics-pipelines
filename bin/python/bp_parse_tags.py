#!/usr/bin/env python3
import os
import sys
import pandas as pd
import re
import sqlite3
import numpy as np

print('Parsing tags')

tagdir = sys.argv[1] 
sample = sys.argv[2]

def process_novor(df):
    global sample 
    sequence = re.sub("[\(\[].*?[\)\]]", "", df[9]).rstrip().strip()
    ngap = 0
    cgap = 0
    scan= str(df[0])
    tag = sequence 
    algorithm = "Novor"
    ID = sample+ '.' + scan
    return (ID, sequence, ngap, cgap, scan, sample, tag, algorithm)

files = os.listdir(tagdir)

combined_tags = []

for f in files:
    if f.endswith('.tags'):
        tags = []
        file = open(tagdir +'/' + f)
        data = file.readlines()
        algorithm = 'DirecTag'
        scans = data[25:]
        scan = '0'
        for row in scans:
            vals = row.split('\n')[0].split('\t')
            if vals[0] == 'S':
                scan = vals[1].split('index=')[1]
            elif vals[0] == 'T':
                sequence = 'O'.join(vals[1].split('0'))
                ngap = str(np.round(float(vals[3]),3))
                cgap = str(np.round(float(vals[4]),3))
                tag = ngap+','+sequence+','+cgap
                if float(ngap) == 0.0:
                    tag = tag.split(',')[1:]
                    tag = ','.join(tag)
                if float(cgap) == 0.0:
                    tag = tag.split(',')[:-1]
                    tag = ','.join(tag)
                ID = sample+'.'+scan
                tags.append((ID, sequence, ngap, cgap, scan, sample, tag, algorithm))
        combined_tags+=tags
    
    elif f.endswith('novor.csv'):
        data = pd.read_csv(tagdir +'/' + f, comment='#', header=None)
        tags = data.apply(process_novor, axis=1).tolist()
        combined_tags+=tags
    
    elif f == sample +'.out':
        tags = []
        regex = re.compile('[^A-Z]')
        file = open(tagdir + '/' + f)
        data = file.read().split('>>')
        data = [d for d in data if not '# No solutions found.' in d]
        data = [d for d in data if not '# too few peaks...' in d]
        data = [d for d in data if not '# Could not process spectrum...' in d]
        data = [d for d in data if not '#Problem reading spectrum...' in d]
        data = [d for d in data if not d == '']
        for res in data:
            res = res.split('\n')
            scan = res[0].split('scan=')[1].split('"')[0]
            for val in res[2:]:
                if val != '':
                    val = val.split('\t')
                    ngap = str(np.round(float(val[3]),3))
                    cgap = str(np.round(float(val[4]),3))
                    sequence = regex.sub('',val[7])
                    tag = ngap +','+sequence +','+cgap
                    if float(ngap) == 0.0:
                        tag = tag.split(',')[1:]
                        tag = ','.join(tag)
                    if float(cgap) == 0.0:
                        tag = tag.split(',')[:-1]
                        tag = ','.join(tag)
                    ID = sample+'.'+scan
                    tags.append((ID, sequence, ngap, cgap, scan, sample, tag, "PepNovo"))
        combined_tags+=tags

db = sqlite3.connect(sys.argv[3], timeout=100)

cursor = db.cursor()
cursor.execute("CREATE TABLE IF NOT EXISTS tags (scanid TEXT, sequence TEXT, ngap REAL, cgap REAL, scan TEXT, sample TEXT, tag TEXT, source TEXT );")
db.commit()

sql = "INSERT INTO tags VALUES(?, ?, ?, ?, ?, ?, ?, ?)"
cursor.executemany(sql, combined_tags )
db.commit()
