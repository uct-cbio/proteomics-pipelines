#!/usr/bin/python

import sys
import os
import Bio; from Bio import SeqIO
import sequtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


if len(sys.argv)==3:
    infile  = sys.argv[1]
    outfile = sys.argv[2]

elif len(sys.argv) == 1:
    infile = sys.stdin
    outfile = sys.stdout

elif len(sys.argv)==2:
    infile = sys.argv[1]
    outfile = sys.stdout

def alt_starts_recs(records, starts = ['ATG','GTG','TTG']):
    all= []
    for seq in records:
        alt_seqs = alt_tss(seq, starts = starts)
        all += alt_seqs
    return all

def alt_tss(seq_record, starts = ['ATG', 'GTG', 'TTG']):
    alt_seqs = []
    #print seq_record.format('fasta')
    alt_seq_count = 1
    desc_ = seq_record.description
    id_   = seq_record.id
    seq_  = str(seq_record.seq)
    num_  =  desc_.split('_')[-1].split('.')[0]
    base_desc = '_'.join(desc_.split('_')[:-1])
    coords  =id_.split('|(')[-1]
    #print id_
    base_id = '_'.join(id_.split('|(')[:-1])
    strand = coords[0]

    coords = coords.split(')')[1].split(':')
    start  = int(coords[0])
    end    = int(coords[1])
    orig_rec = SeqRecord(seq = Seq(seq_) , id = id_, description = desc_)
    alt_seqs.append(orig_rec)
    count = 3
    while count <= len(seq_)-3:
        codon = seq_[count:count +3]
        if codon in starts:
            alt_seq = Seq(seq_[count:])
            assert len(alt_seq) % 3 == 0
            if strand == '+':
                new_start = start + count
                new_end   = end
            elif strand == '-':
                new_start = start
                new_end = end - count
            new_rec = SeqRecord(id = base_id + '|({}){}:{}'.format(strand, new_start, new_end),
                                description = base_desc +'_{}.{}'.format(str(num_), str(alt_seq_count)),
                                seq = alt_seq)
            alt_seqs.append(new_rec)
            alt_seq_count += 1
        count += 3
    return alt_seqs

seqs = list(SeqIO.parse(infile, 'fasta'))
alternative = alt_starts_recs(seqs)
SeqIO.write(alternative, outfile, 'fasta')