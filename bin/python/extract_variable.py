#!/usr/bin/env python3

import sys
import yaml

infile = sys.argv[1]
var = sys.argv[2]

if infile.endswith('.yml'):
    with open(infile) as f:
        res = yaml.load(f, Loader=yaml.CLoader)

print(res[var])   
