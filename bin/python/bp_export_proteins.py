#!/usr/bin/env python3

import sqlite3
import sys
import pandas as pd

print("\n\nSTARTING DATABASE EXPORT...\n\n")


db = sqlite3.connect(sys.argv[1])
cursor1 = db.cursor()
cursor2 = db.cursor()

cursor1.execute("select * from proteins")
names = list(map(lambda x: x[0], cursor1.description))
samples = [i for i in names if i.startswith('SAF ')]

safdict={}

for s in samples:
    cursor1.execute("select SUM(`{}`) from proteins;".format(s))
    for row in cursor1.fetchall():
        safdict[s] = row[0]
        try:
            cmd = "alter table proteins add column '{}' REAL".format('N' + s)
            cursor1.execute(cmd)
            db.commit()
        except:
            pass
try:
    cmd = "alter table proteins add column Summed_NSAF REAL"
    cursor1.execute(cmd)
    db.commit()
except:
    pass


scans = set()

try:
    cursor1.execute('CREATE INDEX accession ON proteins(Accession);')
    db.commit()
except:
    pass

cursor1.execute('SELECT * FROM proteins') 
names = list(map(lambda x: x[0], cursor1.description))
for row in cursor1:
    cols = {}
    for i in range(len(row)):
        col = names[i]
        val = row[i]
        cols[col] = val
    acc = cols['Accession']
    summed_nsaf = 0
    for col in cols:
        val = cols[col]
        if val is not None:
            if col in safdict:
                nsaf = val / safdict[col]
                summed_nsaf += nsaf
                nsaf_col = 'N' + col
                cmd = "UPDATE proteins set `{}`={} where Accession = '{}';".format(nsaf_col, nsaf, acc)
                cursor2.execute(cmd)
            elif col == 'Scans':
                scans.update([ i for i in val.split('\n') if i != ''])
    cmd = "UPDATE proteins set Summed_NSAF={} where Accession = '{}';".format(summed_nsaf, acc)
    cursor2.execute(cmd)
db.commit()


try:
    cmd = "alter table proteins add column Protein_Prob REAL"
    cursor1.execute(cmd)
    db.commit()
except:
    pass

try:
    cmd = "alter table proteins add column Organism_Prob REAL"
    cursor1.execute(cmd)
    db.commit()
except:
    pass

try:
    cmd = "alter table proteins add column MSMS_Percent REAL"
    cursor1.execute(cmd)
    db.commit()
except:
    pass

try:
    cmd = "alter table proteins add column Combined_Prob REAL"
    cursor1.execute(cmd)
    db.commit()
except:
    pass

mapped_scans = len(scans)
scans = set()
cursor1.execute('SELECT * FROM proteins ORDER BY Summed_NSAF DESC') 
names = list(map(lambda x: x[0], cursor1.description))
prev = 0
accs = set()
for row in cursor1:
    cols = {}
    for i in range(len(row)):
        col = names[i]
        val = row[i]
        cols[col] = val
    acc = cols['Accession']
    summed_nsaf = 0
    for col in cols:
        val = cols[col]
        if col == 'Scans':
            scans.update([ i for i in val.split('\n') if i != ''])
    percent = len(scans)/float(mapped_scans) * 100
    if percent != prev:
        prev= percent
        accs.add(acc)
db.commit()

passed_accs = list(accs)
query='SELECT * FROM proteins WHERE Accession in ("{0}")'.format('", "'.join(passed_accs))
df = pd.read_sql_query(sql=query, con=db)

filt = df.groupby(df.Organism).agg({"Summed_NSAF": sum}).sort_values('Summed_NSAF', ascending=False).reset_index()
filt['OrganismProb'] = filt['Summed_NSAF'] / filt['Summed_NSAF'].sum()
filt= filt.set_index('Organism')
org_series = pd.Series(filt.OrganismProb)
org_dict = org_series.to_dict()

cursor1.execute("select SUM(`Summed_NSAF`) from proteins;")
for row in cursor1.fetchall():
    summed_summed_nsaf = row[0]

# Add the organism probability
cursor1.execute('SELECT * FROM proteins ORDER BY Summed_NSAF DESC') 
names = list(map(lambda x: x[0], cursor1.description))
for row in cursor1:
    cols = {}
    for i in range(len(row)):
        col = names[i]
        val = row[i]
        cols[col] = val
    acc = cols['Accession']
    org = cols['Organism']
    prot_prob = cols['Summed_NSAF'] / summed_summed_nsaf
    try:
        org_prob = org_dict[org]
    except:
        org_prob = 0
    combined_prob = prot_prob * org_prob
    cmd = "UPDATE proteins set Organism_Prob={} where Accession = '{}';".format(org_prob, acc)
    cursor2.execute(cmd)
    cmd = "UPDATE proteins set Combined_Prob={} where Accession = '{}';".format(combined_prob, acc)
    cursor2.execute(cmd)
    cmd = "UPDATE proteins set Protein_Prob={} where Accession = '{}';".format(prot_prob, acc)
    cursor2.execute(cmd)
db.commit()

scans = set()
cursor1.execute('SELECT * FROM proteins ORDER BY Combined_Prob DESC, Protein_Prob DESC') 
names = list(map(lambda x: x[0], cursor1.description))
prev = 0
accs = set()
for row in cursor1:
    cols = {}
    for i in range(len(row)):
        col = names[i]
        val = row[i]
        cols[col] = val
    acc = cols['Accession']
    summed_nsaf = 0
    for col in cols:
        val = cols[col]
        if col == 'Scans':
            scans.update([ i for i in val.split('\n') if i != ''])
    percent = len(scans)/float(mapped_scans) * 100
    if percent != prev:
        prev= percent
        accs.add(acc)
passed_accs = list(accs)
query='SELECT * FROM proteins WHERE Accession in ("{0}") ORDER BY Combined_Prob DESC, Protein_Prob DESC'.format('", "'.join(passed_accs))
df = pd.read_sql_query(sql=query, con=db)
output = sys.argv[2]
recs = ''.join(df['Record'].tolist())
w = open( output + '/metanovo.fasta' , 'w' )
w.write(recs)
w.close()
df.to_csv(output + '/metanovo.csv')



