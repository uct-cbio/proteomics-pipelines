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

config = yaml.load(open(sys.argv[1]).read())
output = sys.argv[2]

ref_genome = list(SeqIO.parse(config['reference_genome'],'fasta'))
SeqIO.write(ref_genome, output + '/jbrowse/{}'.format(config['reference_genome'].split('/')[-1]), 'fasta')

w = open(output +'/jbrowse/jbrowse_upload.sh','w')
w.write('./bin/prepare-refseqs.pl --fasta jbrowse/{}'.format(config['reference_genome'].split('/')[-1]) + '\n')

for strain in config['strains']:
    peptides = "jbrowse/"+ strain +'/{}_reference_peptides.gff3'.format(strain)
    cmd = './bin/flatfile-to-json.pl --gff {} --type polypeptide --getSubs --tracklabel "{}_PEPTIDE_ALIGNMENTS" --key "{}_PEPTIDE_ALIGNMENTS" --cssclass generic_parent --subfeatureClasses '.format(peptides, strain, strain) + "'{ " + '"match_part": "generic_part_a"' + "}'"
    w.write(cmd + '\n')

    orfs = "jbrowse/"+ strain +'/{}_reference_orfs.gff3'.format(strain)
    cmd = './bin/flatfile-to-json.pl --gff {} --type polypeptide --getSubs --tracklabel "{}_ORF_ALIGNMENTS" --key "{}_ORF_ALIGNMENTS" --cssclass transcript --subfeatureClasses '.format(orfs, strain, strain) + "'{ " + '"match_part": "transcript-CDS"' + "}'"
    w.write(cmd + '\n')

w.write('./bin/generate-names.pl -v' + '\n')
w.close()
