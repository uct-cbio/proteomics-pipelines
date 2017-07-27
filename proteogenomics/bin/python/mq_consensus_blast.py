#!/usr/bin/env python3

import pandas as pd
import sys
import yaml
import pickle

config = yaml.load(open(sys.argv[1]))
output = sys.argv[2]

for strain in config['strains']:    
    file = output +'/strains/{}/{}_mapped_peptides.csv'.format(strain, strain)
    data = pd.read_csv(file)
    print(data.head())
