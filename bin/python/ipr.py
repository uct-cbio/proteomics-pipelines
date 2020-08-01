#!/usr/bin/env python3
import time
import sys
import threading
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
import tempfile
import subprocess
import shutil
import pandas as pd
import io

limit = 30
threads = []

jobq=mp.Queue(limit)
assert jobq.empty()
results = []

t1 = time.time()

#names = []
#for i in range(15):
#    names.append(i)
#print(names)

def request(jobq):
    while True:
        rec =jobq.get()
        if rec is None:
            return
        else:
            tempdir = tempfile.mkdtemp()
            tempfasta = tempdir + '/{}.fasta'.format(rec.id.split('|')[-1])
            SeqIO.write(rec, tempfasta, 'fasta')
            t1 = time.time()
            #cmd = 'python3 iprscan5.py --goterms --pathways --appl SMART --appl TMHMM --appl CDD --appl Pfam --appl Phobius --appl ProDom --appl SignalP --appl TIGRFAM --appl COILS --appl Gene3D --appl HAMAP --appl MOBIDB --appl PANTHER --appl PIRSF --appl PRINTS --appl PROSITE --appl SFLD --email=matthys@gmail.com --outfile={} --outformat=tsv --quiet {}'.format(tempfasta, tempfasta)
            cmd = 'iprscan5.py --goterms --pathways --email=matthys@gmail.com --outfile={} --outformat=tsv --quiet {}'.format(tempfasta, tempfasta)
            done=False
            while done == False:
                t1 = time.time()
                try: 
                    p = subprocess.Popen(cmd, shell=True)
                    p.wait()
                    print(p.communicate())
                    assert p.returncode == 0
                    with open(tempfasta+'.tsv.txt') as f:
                        res = f.read()
                        if res != '':
                            print(res)
                            results.append(res)
                        done = True
                    print('Job took {} seconds'.format(str(time.time()-t1)))
                except:
                    time.sleep(60)
                    print('Failed, retrying {}'.format(tempfasta))
            shutil.rmtree(tempdir)

workers = []

for i in range(limit):
    p = threading.Thread(target=request, args=(jobq,))
    p.daemon=True
    p.start()
    workers.append(p)

# Load the recs
recs = SeqIO.parse(sys.argv[1],'fasta')
for i in recs:
    jobq.put(i)


for i in range(limit):
    jobq.put(None)

for w in workers:
    w.join()

assert jobq.empty() # All workers took their None from the queue, ie none were dead

if len(results) > 0:
    results = '\n'.join(results)
    data = io.StringIO(results)
    results = pd.read_csv(data, sep='\t', header=None)
else:
    results = pd.DataFrame()

print(results.head())
results.to_csv(sys.argv[1] + '.tsv', header=False)

print('Results: ', len(results))
print('Time elapsed: ', time.time()-t1)
