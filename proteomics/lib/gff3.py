#!/usr/bin/env python3
import pandas as pd
import sys
import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
from io import StringIO
import yaml
import blast
from Bio.Blast import NCBIXML
import gff3
 

def gff3_contig_export(recs):
    contig_columns = []
    for rec in recs:
        id = rec.id
        reference_sequence = id.split('|')[0]
        size = len(str(rec.seq))
        start = '1'
        end = str(len(str(rec.seq)))
        strand = '.'
        source = 'Assembly'
        type='contig'
        score ='.'
        strand='.'
        phase='.'
        attributes='Name={}'.format(reference_sequence)
        row = reference_sequence + '\t' + source + '\t' + type + '\t' + start + '\t' + end     + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes
        contig_columns.append(row)
    return contig_columns

def orfblastp2gff3(blast_records, query_fasta, orf_fasta,  eval_cutoff, subject_cutoff):
    mapdict = {}
    subject_count = {}

    recdict = defaultdict(list)
    for blast_record in blast_records:
        top_eval=None
        aln_count = 1
        qdesc = query_fasta[blast_record.query.split()[0]].description
        qid = query_fasta[blast_record.query.split()[0]].id
        qseq = str(query_fasta[blast_record.query.split()[0]].seq)
        qlen = len(qseq)
        blast_record.alignments.sort(key=lambda align: max(hsp.score for hsp in align.hsps), reverse=True)
        rec_count = 1
        for alignment in blast_record.alignments:
            alignment.hsps.sort(key=lambda hsp :  hsp.score , reverse=True)
            if not qseq in mapdict:
                mapdict[qseq] = None
                subject_count[qseq] = 0
            for hsp in alignment.hsps:
                if mapdict[qseq] != hsp.sbjct:
                    mapdict[qseq] = hsp.sbjct
                    subject_count[qseq] += 1
                if (subject_count[qseq] <= subject_cutoff) and (hsp.expect <= eval_cutoff):
                    orf_rec = orf_fasta[alignment.hit_def.split()[0]]
                    orf_seq = str(orf_rec.seq)
                    orf_coords = orf_rec.id.split('|')[-1]
                    orf_start = int(orf_coords.split(':')[0].split(')')[1])
                    orf_end = int(orf_coords.split(':')[1])
                    orf_strand = orf_coords.split('(')[1].split(')')[0]
                    hitseq = hsp.sbjct.replace('-','')
                    amino_offset = orf_seq.find(hitseq)  
                    b2g = blast.blast2genome(fasta_id = qid, query_seq = hsp.query, match_seq=hsp.match, hit_seq = orf_seq, hit_from=orf_start, hit_to=orf_end, hit_strand=orf_strand, amino_offset=amino_offset)  
                    
                    contig = orf_rec.id.split('_')[-3]
                    source = "BLASTP"
                    type="polypeptide"
                    score = hsp.expect
                    phase= '.'
                    attributes="ID=match{}-{};Name={}".format(rec_count, qseq, qseq)
                    start = b2g.alignment_start
                    end = b2g.alignment_end
                    strand = orf_strand
                    row = contig +'\t' + source +'\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + strand + '\t' + phase + '\t' + attributes
                    
                    recdict[qseq].append(row)
                    
                    for rec in b2g.records:
                        start = rec[1]
                        end = rec[2]
                        strand = rec[3]
                        score = hsp.expect
                        phase= '.'
                        attributes="Parent=match{}-{};Name={}".format(rec_count, qseq, qseq)
                        row = contig +'\t' + source +'\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + strand + '\t' + phase + '\t' + attributes
                        recdict[qseq].append(row)
                    rec_count += 1
    return recdict

def tblastn2gff3(blast_records, query_fasta, eval_cutoff, subject_cutoff):
    mapdict = {}
     
    subject_count = {}

    recdict = defaultdict(list)
    
    query_dict = SeqIO.to_dict(query_fasta)

    for blast_record in blast_records:
        top_eval=None
        aln_count = 1
        qid = blast_record.query.split()[0]
        newids = blast_record.query.split()[1].split(';')
        newids = [i.split('|')[1] for i in newids]
        qseq = str(query_dict[blast_record.query.split()[0]].seq)
        qlen = len(qseq)
        blast_record.alignments.sort(key=lambda align: max(hsp.score for hsp in align.hsps), reverse=True)
        
        rec_count = 1
        
        for alignment in blast_record.alignments:
            alignment.hsps.sort(key=lambda hsp :  hsp.score , reverse=True)
            if not qseq in mapdict:
                mapdict[qseq] = None
                subject_count[qseq] = 0

            for hsp in alignment.hsps:
                
                if mapdict[qseq] != hsp.sbjct:
                    mapdict[qseq] = hsp.sbjct
                    subject_count[qseq] += 1

                if (subject_count[qseq] <= subject_cutoff) and (hsp.expect <= eval_cutoff):
                    start = hsp.sbjct_start
                    end = hsp.sbjct_end
                    frame  = hsp.frame[1]
                    if frame < 0:
                        strand='-'
                    else:
                        strand='+'
                    hitseq = hsp.sbjct.replace('-','')
                    b2g = blast.blast2genome(fasta_id = qid, query_seq = hsp.query, match_seq=hsp.match, hit_seq = hitseq, hit_from=start, hit_to=end, hit_strand=strand, amino_offset=0)  
                    contig = alignment.hit_def
                    source = "BLASTP"
                    type="polypeptide"
                    score = hsp.expect
                    phase= '.'
                    start = b2g.alignment_start
                    end = b2g.alignment_end
                    for newid in newids:
                        attributes="ID=match{}-{};Name={}".format(rec_count, newid, newid)
                        row = contig +'\t' + source +'\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + strand + '\t' + phase + '\t' + attributes
                        recdict[newid].append(row)
                    for rec in b2g.records:
                        contig = alignment.hit_def
                        source = "BLASTP"
                        type="polypeptide"
                        start = rec[1]
                        end = rec[2]
                        strand = rec[3]
                        score = hsp.expect
                        phase= '.'
                        for newid in newids:
                            attributes="Parent=match{}-{};Name={}".format(rec_count, newid, newid)
                            row = contig +'\t' + source +'\t' + type + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + strand + '\t' + phase + '\t' + attributes
                            recdict[newid].append(row)
                    rec_count += 1
    return recdict
