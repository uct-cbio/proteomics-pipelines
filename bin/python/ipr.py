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
import os

limit = 30
threads = []

jobq=mp.Queue(limit)
assert jobq.empty()
results = []

t1 = time.time()

names = []
for i in range(15):
    names.append(i)

def request(jobq):
    while True:
        rec =jobq.get()
        if rec is None:
            return
        else:
            tempdir = tempfile.mkdtemp()
            rec_id = rec.id
            if '|' in  rec_id:
                rec_id = rec_id.split('|')[1]
            tempfasta = tempdir + '/{}.fasta'.format(rec_id)
            newrec = Bio.SeqRecord.SeqRecord( id=rec.id, seq =rec.seq)
            #print(newrec.format('fasta'))
            SeqIO.write(newrec, tempfasta, 'fasta')
            t1 = time.time()
            #cmd = 'python3 iprscan5.py --goterms --pathways --appl SMART --appl TMHMM --appl CDD --appl Pfam --appl Phobius --appl ProDom --appl SignalP --appl TIGRFAM --appl COILS --appl Gene3D --appl HAMAP --appl MOBIDB --appl PANTHER --appl PIRSF --appl PRINTS --appl PROSITE --appl SFLD --email=matthys@gmail.com --outfile={} --outformat=tsv --quiet {}'.format(tempfasta, tempfasta)
            cmd = 'iprscan5.py --verbose --goterms --pathways --email=matthys.potgieter@gmail.com --outfile={} --outformat=tsv  {}'.format(tempfasta + '.tsv', tempfasta)
            print(cmd)
            done=False
            #failures = 0
            while done == False:
                t1 = time.time()
                try: 
                    p = subprocess.Popen(cmd, shell=True)
                    #p.wait()
                    print(p.communicate())
                    assert p.returncode == 0
                    #if os.path.exists(tempfasta+'.tsv.txt'):
                    #    infile = tempfasta+'.tsv.txt'
                    #else:
                    infile = tempfasta+'.tsv.tsv.tsv'
                    with open(infile) as f:
                        res = f.read()
                        if res != '':
                            data = io.StringIO(res)
                            df = pd.read_csv(data, names=names, sep='\t', header=None)
                            #print(df)
                            results.append(df)
                        done = True
                    #print('Job took {} seconds'.format(str(time.time()-t1)))
                except Exception as e:
                    #time.sleep(60)
                    print(e)
                    print('Failed interproscan query for  {}'.format(tempfasta))
                #failures += 1
                #if failures > 10:
                #    done=True
                shutil.rmtree(tempdir)

workers = []

for i in range(limit):
    p = threading.Thread(target=request, args=(jobq,))
    p.daemon=True
    p.start()
    workers.append(p)

# Load the recs
count = 0
recs = SeqIO.parse(sys.argv[1],'fasta')
for i in recs:
    #count += 1
    #if count > 10:
    #    continue
    jobq.put(i)


for i in range(limit):
    jobq.put(None)

for w in workers:
    w.join()

assert jobq.empty() # All workers took their None from the queue, ie none were dead

if len(results) == 0:
    res = pd.DataFrame()
else:
    res = pd.concat(results)

if len(res) > 0:
    res.to_csv(sys.argv[1] + '.tsv', sep='\t', index=False, header=False)
print('Results: ', len(results))
print('Time elapsed: ', time.time()-t1)

