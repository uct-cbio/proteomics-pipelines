#!/usr/bin/env python3

import yaml
import json
import pandas as pd
import sys

config = yaml.load(open(sys.argv[1]).read(), Loader=yaml.Loader)
operons = pd.read_csv(config['operons'], sep='\t')

outpath = sys.argv[2]

operons['GI'] = operons['GI'].apply( lambda x : str(int(x)))
operons['OperonID'] = operons['OperonID'].apply( lambda x : str(int(x)))

gi_map =pd.Series(operons.OperonID.values, index=operons.GI).to_dict()

gi_map = json.dumps(gi_map)

w = open(outpath +'/mapping/operons.json', 'w')

w.write(gi_map)

w.close()
