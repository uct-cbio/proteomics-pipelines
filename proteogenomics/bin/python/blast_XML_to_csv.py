#!/usr/bin/env python

import Bio
import sys
from Bio.Blast import NCBIXML
import os
import csv
import pandas as pd

combined = pd.DataFrame()
row=0

result_handle = open(sys.argv[1])

blast_records = NCBIXML.parse(result_handle)


output = csv.writer(sys.stdout, delimiter='\t')
headers = [         'blast_record.query',
                    'blast_record.query_id',
                    'blast_record.query_length',
                    'blast_record.query_letters',
                    'alignment.accession',
                     'alignment.hit_def',
                    'alignment.hit_id',
                    'alignment.title',
                    'hsp.align_length',
                    'hsp.bits',
                    'hsp.expect',
                    'hsp.frame',
                    'hsp.gaps',
                    'hsp.identities',
                    'hsp.match',
                    'hsp.num_alignments',
                    'hsp.positives',
                    'hsp.query',
                    'hsp.query_end',
                    'hsp.query_start',
                    'hsp.sbjct',
                    'hsp.sbjct_end',
                    'hsp.sbjct_start',
                    'hsp.score',
                    'hsp.strand',
                    '_visualization',
                    '_alignment_rank',
                    '_hsp_rank',
                    '_alignment.source',
                    '_alignment.entry',
                    '_alignment.entry_name',
                    '_alignment.description',
                    '_blast_record.query.description' ] 

output.writerow(headers)

for blast_record in blast_records:
    aln_count = 1
    top_eval=None
    for alignment in blast_record.alignments:
        hsp_count = 1
        for hsp in alignment.hsps:
            
            if (aln_count == 1) and (hsp_count == 1):
                top_eval = hsp.expect
            else:
                assert hsp.expect >= top_eval

            hsp_match = []
            hsp_match.append(hsp.query)
            hsp_match.append(hsp.match)
            hsp_match.append(hsp.sbjct)
            hsp_match = '\n'.join(hsp_match)
            rd = [  blast_record.query,
                    blast_record.query_id,
                    blast_record.query_length,
                    blast_record.query_letters,
                    alignment.accession,
                    alignment.hit_def,
                    alignment.hit_id,
                    alignment.title,
                    hsp.align_length,
                    hsp.bits,
                    hsp.expect,
                    str(hsp.frame),
                    hsp.gaps,
                    hsp.identities,
                    hsp.match,
                    hsp.num_alignments,
                    hsp.positives,
                    hsp.query,
                    hsp.query_end,
                    hsp.query_start,
                    hsp.sbjct,
                    hsp.sbjct_end,
                    hsp.sbjct_start,
                    hsp.score,
                    str(hsp.strand),
                    hsp_match,
                    aln_count,
                    hsp_count,
                    alignment.hit_def.split('|')[0],
                    alignment.hit_def.split('|')[1],
                    alignment.hit_def.split('|')[2].split()[0],
                    ' '.join(alignment.hit_def.split()[1:]),
                    blast_record.query.split()[1] ] 
            output.writerow(rd)
            hsp_count += 1
            row += 1
        aln_count += 1
