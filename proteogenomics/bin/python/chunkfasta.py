#!/usr/bin/env python

import Bio
from Bio import SeqIO
from itertools import zip_longest
import pandas as pd
import shutil
import os
import sys
import math

input_fasta=sys.argv[1]
chunknumber=int(sys.argv[2])
output_folder=sys.argv[3]

filename=os.path.basename(input_fasta)

query_folder=output_folder 

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return list(zip_longest(*args, fillvalue=fillvalue))

records =list(SeqIO.parse(input_fasta,'fasta'))

chunksize = math.ceil(len(records)/chunknumber)
chunked = grouper(records, chunksize)


chunks =[]
count=1
for chunk in chunked:
    partname='{}.part.{}.fasta'.format(filename, count)
    chunk = pd.Series(chunk).dropna().tolist()
    count += 1
    SeqIO.write(chunk, query_folder +'/' + partname, 'fasta')



