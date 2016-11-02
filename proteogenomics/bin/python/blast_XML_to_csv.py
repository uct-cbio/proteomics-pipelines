#!/usr/bin/env python

import Bio
import sys
from Bio.Blast import NCBIXML
import os
import csv
import pandas as pd
from Bio import SeqIO

combined = pd.DataFrame()
row=0

result_handle = open(sys.argv[1])
fasta_sequence  = SeqIO.to_dict(list(SeqIO.parse(sys.argv[2],'fasta')))

blast_records = NCBIXML.parse(result_handle)


output = open(sys.argv[3],'w')

output = csv.writer(output, delimiter=',')

aln_cutoff = int(sys.argv[4])


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
                    '_blast_record.query.description',
                    '_query.id',
                    '_query.description',
                    '_query.sequence',
                    '_query.sequence.length'] 

output.writerow(headers)

for blast_record in blast_records:
    top_eval=None
    
    aln_count = 1
        
    qdesc = fasta_sequence[blast_record.query.split()[0]].description
    qseq = str(fasta_sequence[blast_record.query.split()[0]].seq)
    qlen = len(qseq)
    blast_record.alignments.sort(key=lambda align: max(hsp.score for hsp in align.hsps), reverse=True)
    for alignment in blast_record.alignments:
        hsp_count = 1
        
        alignment.hsps.sort(key=lambda hsp :  hsp.score , reverse=True)

        for hsp in alignment.hsps:
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
                    ' '.join(blast_record.query.split()[1:]),
                    blast_record.query.split()[0],
                    qdesc,
                    qseq,
                    qlen ]
            assert len(headers) == len(rd)
            output.writerow(rd)
            hsp_count += 1
            row += 1
        if aln_count == aln_cutoff:
            break
        aln_count += 1
        
