#!/usr/bin/env python3

import itertools
import string
from elasticsearch import Elasticsearch,helpers
import sys
import os
from glob import glob
import pandas as pd
import json

host     = sys.argv[1]
port     = int(sys.argv[2])
alias    = sys.argv[3]

print(host)
print(port)
print(alias)

es = Elasticsearch([{'host': host, 'port': port}])

# create our test index

# Get all csv files in /root/data
files = [y for x in os.walk('/root/data') for y in glob(os.path.join(x[0], '*.csv'))]

count = 0
def clean_field(val):
    val = val.split('.')
    val = [i for i in val if i != '']
    val = '_'.join(val)
    val = val.split()
    val = [i for i in val if i != '']
    val = '_'.join(val)
    val = val.split('/')
    val = [i for i in val if i != '']
    val = '_'.join(val)
    return val

es.indices.delete(index=alias + '*', ignore=[400, 404])

indices = []

for file in files:
    data = pd.read_csv(file, sep=None, engine='python')
    index = alias + '_'.join(file.split('/'))
    index = clean_field(index).lower().split('_csv')[0]
    indices.append(index)
    es.indices.create(index)
    for col in data.columns:
        if col.startswith('Unnamed'):
            del data[col]
        else:
            data.rename(columns= { col : clean_field(col) },inplace=True )
    data = data.reset_index() # Make sure there is no duplicate indexing
    data.rename(columns={'index':'row'},inplace =True)
    data['File'] = file
    data['_id'] = data['File'] + '.{}.'.format(str(count)) + data.reset_index()['index'].apply(str)
    data['_type'] = "document"
    data['_index'] = index
    records = data.to_json(orient='records')
    records = json.loads(records)
    helpers.bulk(es, records, chunk_size=100)
    count += 1
    print(es.count(index=index))

# Create an index table in elasticsearch to locate the files
indices_table = pd.DataFrame()
indices_table['Index'] = pd.Series(indices)
indices_table['File'] = pd.Series(files)
indices_table['Alias'] = alias
indices_table['_id'] = indices_table['Alias'] + '.' + indices_table['File'] 
indices_table['_type'] = "document"
indices_table['_index'] = alias + '_indices'
es.indices.create(alias + '_indices')
records = indices_table.to_json(orient='records')
records = json.loads(records)
helpers.bulk(es, records, chunk_size=100)
print(es.count(index=alias + '_indices'))



