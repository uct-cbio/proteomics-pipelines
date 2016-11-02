#!/usr/bin/env python

import Bio
from Bio import SeqIO
from itertools import zip_longest
import pandas as pd
import shutil
import os
import sys
import math
import tempfile
import pandas as pd

input_fasta=sys.argv[1]
chunksize=int(sys.argv[2])

output_folder=sys.argv[3]

filename=os.path.basename(input_fasta)

records=SeqIO.parse(input_fasta,'fasta')

count=1
chunkcount=1
recs = []

for rec in records:
    recs.append(rec)
    if chunkcount==chunksize:
        partname='{}.part.{}.fasta'.format(filename, count)
        SeqIO.write(recs, output_folder +'/' + partname, 'fasta')
        recs = []
        count +=1
        chunkcount=1
    else:
        chunkcount+= 1

if len(recs) > 0:
    partname='{}.part.{}.fasta'.format(filename, count)
    SeqIO.write(recs, output_folder + '/' + partname, 'fasta')
    recs = []
    count +=1




