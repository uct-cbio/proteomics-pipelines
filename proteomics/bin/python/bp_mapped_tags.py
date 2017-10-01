#!/usr/bin/env python3

import pandas as pd
import sys
import Bio
from Bio import SeqIO
import tagmatch
import os
from collections import defaultdict
import sqlite3
import gzip

enzymes =os.environ['mn_enzymes'].split(';')
specificity =os.environ['mn_specificity']
max_mc = int(os.environ['mn_max_missed_cleavages'])

f = open(sys.argv[1])
lines = f.readlines()
f.close()

db = sqlite3.connect(sys.argv[2], timeout=100)

cursor = db.cursor()

try:
    query = 'delete from proteins where File=?'
    cursor.execute(query, (sys.argv[1],))
    db.commit()
except:
    pass

seqs = {}
f = sys.argv[1].split('.csv')[0] + '.gz'
with gzip.open(f, "rt") as handle:
        for i in SeqIO.parse(handle, "fasta"):
            acc = i.id.split('|')[1]
            seqs[acc] = i

passed_tags = set()
matched_tags = defaultdict(set)
matched_peptides = defaultdict(set)

for line in lines:
    l = line.split(',')
    tag = ','.join(l[:-3])
    peptide = l[-3]
    acc = l[-2]
    seq = str(seqs[acc].seq)
    pos = int(l[-1])
    if (pos == 2) and (seq[0] == 'M'):
        seq = seq[1:]
        pos = pos - 1
    start = pos-1
    end = start + len(peptide)
    if start > 0:
        amino_acid_before = seq[start-1]
    else:
        amino_acid_before = ""
    first_amino_acid = seq[start]
    
    try:
        last_amino_acid = seq[end-1]
    except:
        continue

    try:
        amino_acid_after = seq[end]
    except:
        amino_acid_after = ""
    
    valid = tagmatch.valid_cleavage(amino_acid_before, first_amino_acid, last_amino_acid, amino_acid_after, enzymes, specificity)
    mc = tagmatch.missed_cleavages(peptide, enzymes)
    if (valid == True) and (mc <= max_mc):
        passed_tags.add(tag)
        matched_tags[acc].add(tag)
        matched_peptides[acc].add(peptide)
passed_tags = list(passed_tags)

ids = defaultdict(set)

from itertools import zip_longest

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

for group in grouper(passed_tags, 500):
    t = [i for i in group if i is not None ]
    cursor.execute('SELECT DISTINCT scanid, tag FROM tags WHERE tag in ({0})'.format(', '.join('?' for _ in t)), t)
    scanids = cursor.fetchall()
    for _ in scanids:
        ids[_[1]].add(_[0])
del scanids

count = 0

def get_organism(val):
    try:
        val = val.split('OS=')[1]
        if '=' in val:
            val = val.split('=')[0]
            val = ' '.join(val.split()[:-1]).strip().rstrip()
        return val
    except:
        return "Unspecified"

_replicates = {}
_accession=[]
_record=[]
_protein_id=[]
_peptides=[]
_peptide_count=[]
_scans=[]
_scan_count=[]
_organism=[]
_seq_length=[]

def create_replicates(df):
    _ = _replicates[df['Accession']]
    for sample in _:
        count = len(_[sample])
        df['Sample {} (msms)'.format(sample)] = count
        df['SAF {}'.format(sample)] = df['Sample {} (msms)'.format(sample)] / df['Length']
    return df

for protein in matched_tags:
    rec = seqs[protein]
    _record.append(rec.format('fasta'))
    organism = get_organism(rec.description)
    tags = matched_tags[protein]
    peptide_set = matched_peptides[protein]
    reps = defaultdict(set)
    scan_set = set()
    for tag in tags:
        scans = ids[tag]
        scan_set.update(scans)
        for scan in scans:
            replicate = scan.split('.mgf')[0]
            reps[replicate].add(scan)
    peptide_count = len(peptide_set)
    peptides = '\n'.join(peptide_set)
    scans='\n'.join(list(scan_set))
    scan_count = len(scan_set)
    accession = protein
    protein_id = rec.id
    seq = str(rec.seq)
    seq_length = len(seq)
    _replicates[protein]=reps
    _accession.append(accession)
    _protein_id.append(protein_id)
    _peptides.append(peptides)
    _peptide_count.append(peptide_count)
    _scans.append(scans)
    _scan_count.append(scan_count)
    _organism.append(organism)
    _seq_length.append(seq_length)
    count += 1

summary = pd.DataFrame()
summary['Accession'] = _accession
summary['Record'] = _record
summary['ID'] = _protein_id
summary['PeptideCount'] = _peptide_count
summary['Peptides'] = _peptides
summary['ScanCount'] = _scan_count
summary['Scans'] = _scans
summary['Organism'] = _organism
summary['Length'] = _seq_length
summary['File'] = sys.argv[1]
summary = summary.apply(create_replicates, axis=1)
sample_cols = [i for i in summary.columns if i.endswith('(msms)') ]
saf_cols = [i for i in summary.columns if i.startswith('SAF ') ]
summary[sample_cols]= summary[sample_cols].fillna(0)
summary[saf_cols]= summary[saf_cols].fillna(0)

for column in summary.columns:
    try:
        cmd = "alter table proteins add column '{}' REAL".format(column)
        cursor.execute(cmd)
        db.commit()
    except:
        pass
summary.to_sql(name="proteins", con=db, if_exists="append")
db.commit()

#print(summary.head())
