#!/usr/bin/env python

import Bio
import sys
from Bio.Blast import NCBIXML
import os
import pandas as pd

combined = pd.DataFrame()
row=0

result_handle = open(sys.argv[1])

blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    aln_count = 1
    top_eval=None
    
    for alignment in blast_record.alignments:
        hsp_count = 1
        for hsp in alignment.hsps:
            
            if aln_count == 1:
                top_eval = hsp.expect
                assert hsp_count == 1
            else:
                assert hsp.expect >= top_eval

            hsp_match = []
            hsp_match.append(hsp.query)
            hsp_match.append(hsp.match)
            hsp_match.append(hsp.sbjct)
            hsp_match = '\n'.join(hsp_match)
            
            combined.loc[row, 'blast_record.query'] = blast_record.query
            combined.loc[row, 'blast_record.query_id'] = blast_record.query_id
            combined.loc[row, 'blast_record.query_length'] = blast_record.query_length
            combined.loc[row, 'blast_record.query_letters'] = blast_record.query_letters
            
            combined.loc[row, 'alignment.accession'] = alignment.accession
            combined.loc[row, 'alignment.hit_def'] = alignment.hit_def
            combined.loc[row, 'alignment.hit_id'] = alignment.hit_id
            combined.loc[row, 'alignment.title'] = alignment.title
         
            combined.loc[row, 'hsp.align_length'] = hsp.align_length
            combined.loc[row, 'hsp.bits'] = hsp.bits
            combined.loc[row, 'hsp.expect'] = hsp.expect
            combined.loc[row, 'hsp.frame'] = str(hsp.frame)
            combined.loc[row, 'hsp.gaps'] = hsp.gaps
            combined.loc[row, 'hsp.identities'] = hsp.identities
            combined.loc[row, 'hsp.match'] = hsp.match
            combined.loc[row, 'hsp.num_alignments'] = hsp.num_alignments
            combined.loc[row, 'hsp.positives'] = hsp.positives
            combined.loc[row, 'hsp.query' ] = hsp.query
            combined.loc[row, 'hsp.query_end' ] = hsp.query_end
            combined.loc[row, 'hsp.query_start' ] = hsp.query_start
            combined.loc[row, 'hsp.sbjct' ] = hsp.sbjct
            combined.loc[row, 'hsp.sbjct_end' ] = hsp.sbjct_end
            combined.loc[row, 'hsp.sbjct_start' ] = hsp.sbjct_start
            combined.loc[row, 'hsp.score' ] = hsp.score
            combined.loc[row, 'hsp.strand' ] = str(hsp.strand)
            combined.loc[row, '_visualization' ] = hsp_match
            combined.loc[row, '_alignment.rank' ] = aln_count
            combined.loc[row, '_hsp.rank' ] = hsp_count
            combined.loc[row, '_alignment.source'] = alignment.hit_def.split('|')[0]
            combined.loc[row, '_alignment.entry'] = alignment.hit_def.split('|')[1]
            combined.loc[row, '_alignment.entry_name'] = alignment.hit_def.split('|')[2].split()[0]
            combined.loc[row, '_alignment.description'] = ' '.join(alignment.hit_def.split()[1:])
            combined.loc[row, '_blast_record.query.description'] = blast_record.query.split()[1]

            hsp_count += 1
            row += 1
        aln_count += 1

combined.to_csv(sys.argv[2])
