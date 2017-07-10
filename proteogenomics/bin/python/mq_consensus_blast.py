#!/usr/bin/env python3

import pandas as pd
import sys
import yaml
import pickle

config = yaml.load(open(sys.argv[1]))
output = sys.argv[2]

for strain in config['strains']:    
    file = output +'/strains/{}/{}_mapped_peptides.p'.format(strain, strain)
    data = pickle.load(open(file, 'rb'))
    print(data.head())
