#!/usr/bin/python

import pandas as pd
import sys

if len(sys.argv) == 2:
    file = pd.ExcelFile(sys.argv[1])
    outfile = sys.stdout

elif len(sys.argv) == 1:
    file = pd.ExcelFile(sys.stdin)
    outfile = sys.stdout

elif len(sys.argv) == 3:
    file = pd.ExcelFile(sys.argv[1])
    outfile = sys.argv[2]


if len(file.sheet_names) == 1
    df = file.parse('Sheet1')
    df.to_csv(outfile)



