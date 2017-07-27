#!/usr/bin/env python

import sys
import sqlite3
import Bio
from Bio import SeqIO
import pickle
import tagmatch
import pandas as pd
import numpy as np
from collections import defaultdict
import shelve

db = sqlite3.connect(sys.argv[1])
cursor = db.cursor()

tag_df = pd.read_sql("SELECT DISTINCT tagid, seqstr from tags", db, columns=['tagid', 'seqstr'])
tag_df = tag_df.drop_duplicates()

load = defaultdict(set)

def load_dict(table):
    global load
    tid = table['tagid']
    seq = table['seqstr']
    load[seq].add(tid)

tag_df.apply(load_dict, axis=1)
del tag_df

pickle.dump(load, open(sys.argv[1] +'.p', 'wb'))
