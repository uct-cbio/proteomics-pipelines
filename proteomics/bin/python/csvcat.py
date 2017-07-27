#!/usr/bin/env python

import pandas as pd
import sys
import os

files = []

for file in os.listdir(sys.argv[1]):
    df = pd.read_csv(sys.argv[1] +'/' + file )
    df['File'] = file
    files.append(df)

combined = pd.concat(files)
combined.to_csv(sys.stdout)
