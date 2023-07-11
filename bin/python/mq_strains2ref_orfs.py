#!/usr/bin/env python3

import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
import yaml
import blast
from Bio.Blast import NCBIXML
import gff3

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.Loader)
output = sys.argv[2]

query_fasta = list(SeqIO.parse(output +'/strains/all_mapped_trans_orfs.fasta','fasta'))

for ref in config['reference']:
    assembly_id = config['reference'][ref]['assembly_id']
    gff3.orf_blast_gff3(config, assembly_id, query_fasta , output)
