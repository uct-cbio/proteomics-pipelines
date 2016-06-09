#!/usr/bin/env python

import pandas as pd
import sys

frames = []
for key in sys.argv[1:]:
    df = pd.read_csv(key)
    df['File'] = key
    frames.append(df)
combined = pd.concat(frames)
combined.to_csv(sys.stdout)
