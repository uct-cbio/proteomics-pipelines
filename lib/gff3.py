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
        reference_sequence = str(rec.id).split(' ')[0]
        
        #reference_sequence = id.split('|')[0]
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
                    contig = orf_rec.id.split('|')[0]
                    if '_ENA_' in contig:
                        contig = 'ENA_' + contig.split('_ENA_')[1]
                        #contig = '|'.join(contig.split('_'))
                        contig = contig.split('_')[1]
                    else:
                        contig = contig.split('_')[-1]
                    
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
        try:
            newids = blast_record.query.split()[1].split(';')
            newids = [i.split('|')[1] for i in newids]
        except:
            newids = [ blast_record.query.split()[0]]
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
                    contig = alignment.hit_def.split()[0]
                    if 'ENA' in contig:
                        contig = contig.split('ENA')[1]
                        contig = contig.split('|')[1]
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
                        contig = alignment.hit_def.split()[0]
                        if 'ENA' in contig:
                            contig = contig.split('ENA')[1]
                            contig = contig.split('|')[1]
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

def peptide_blast_gff3(blast_name, output, peptide_sequence, config, strain_samples, peptides):
    orf_sequence  = SeqIO.to_dict(list(SeqIO.parse(output + '/blast/peptides2orfs/{}_reforfs.fasta'.format(blast_name),'fasta')))
    blast_results = output + '/blast/peptides2orfs/{}.xml'.format(blast_name)
    result_handle = open(blast_results)
    blast_records = NCBIXML.parse(result_handle)
    print("PARSED BLAST RESULTS FOR REFERENCE GENOMES")
    mapdict = {}
    recdict = defaultdict(list)
    recdict = gff3.orfblastp2gff3(blast_records=blast_records, query_fasta=peptide_sequence, orf_fasta=orf_sequence,  eval_cutoff=200000, subject_cutoff=1)

    ref_genome =list(SeqIO.parse(output + '/ena/{}/{}.fasta'.format(blast_name, blast_name),'fasta'))
    contigs = gff3.gff3_contig_export(ref_genome)


    for strain in config['strains']:
        samples = strain_samples[strain]
        sample_columns = []
        
        for sample in samples:
            sample_columns.append('Experiment {}'.format(sample))
        
        strain_filt = peptides[peptides[sample_columns].sum(axis=1) >=1]
        strain_peptides = list(set(strain_filt['Sequence'].tolist()))
        strain_rows = ['##gff-version 3\n##Index-subfeatures 1\n']
        #strain_rows += contigs
        for peptide in strain_peptides:
            if peptide in recdict:
                results = recdict[peptide]
                strain_rows += results
        gff_str = '\n'.join(strain_rows)

        if not os.path.exists(output + '/jbrowse/' + strain):
            os.mkdir(output + '/jbrowse/' + strain)     

        # Write the peptide features 
        w =open(output+ '/jbrowse/' + strain + '/{}_{}_peptides.gff3'.format(strain,blast_name), 'w')
        w.write(gff_str)
        w.close()


def orf_blast_gff3(config, blast_name, query_fasta, output):
    blast_results = output + '/blast/orfs2genome/{}.xml'.format(blast_name)

    result_handle = open(blast_results)

    blast_records = NCBIXML.parse(result_handle)

    recdict = gff3.tblastn2gff3(blast_records=blast_records, query_fasta=query_fasta, eval_cutoff=200000, subject_cutoff=1)

    strains = defaultdict(list)

    for key in recdict:
        s = key.split('_')[0]
        strains[s] += recdict[key]

    quant = pd.read_csv(output+ '/diff/protein_normalization/msnbase/normalized.csv')

    id_dict = {}

    for row in quant.iterrows():
        try: 
            orfs = row[1]['ORF.ids...all.strains'].split('\n')
        except:
            continue
        orfs = [ i.split('|')[1] for i in orfs ] 
        for orf in orfs:
            if not orf in id_dict:
                id_dict[orf] = defaultdict()
        for col in row[1].index:
            if col.startswith('iBAQ.'):
                s = col.split('.')[1]
                intensity = row[1][col]
                for orf in orfs:
                    if not s in id_dict[orf]:
                        id_dict[orf][s] = intensity         
                    elif not intensity > id_dict[orf][s]:
                        id_dict[orf][s] = intensity         
    orf_coords = {}

    for strain in config['strains']:
        
        orf_coords[strain] = defaultdict(list)

        strain_rows = ['##gff-version 3\n##Index-subfeatures 1\n']
        strain_rows += strains[strain]

        gff_str = '\n'.join(strain_rows)
        
        if not os.path.exists(output + '/jbrowse/' + strain):
            os.mkdir(output + '/jbrowse/' + strain)     
        # Write the peptide features 
        
        w =open(output+ '/jbrowse/' + strain + '/{}_{}_orfs.gff3'.format(strain,blast_name), 'w')
        w.write(gff_str)
        w.close()
        
        for row in strains[strain]:
            _ = row.split('\t')
            if _[8].startswith('ID='):
                start = int(_[3])
                end = int(_[4])
                coords = [i for i in range(start, end + 1)]
                Id = _[8].split(';Name=')[1]
                orf_coords[strain][Id] += coords

        samps = defaultdict(list)
        for Id in orf_coords[strain]:
            coords = orf_coords[strain][Id]
            if Id in id_dict: 
                for samp in id_dict[Id]:
                    intensity = id_dict[Id][samp]
                    for coord in coords:
                        tup = (coord, intensity)
                        samps[samp].append(tup)
        
        for samp in samps:
            tups = samps[samp]
            ints = {}
            for tup in tups:
                if tup[0] in ints:
                    ints[tup[0]] += tup[1]
                else:
                    ints[tup[0]] = tup[1]
            vals = [(k, v) for k, v in ints.items()]
            vals = sorted(vals, key=lambda x: x[0])
            vals=  [ ' '.join([str(v[0]), str(v[1])]) for v in vals ] 
            header = 'variableStep chrom=Chromosome\n' + '\n'.join(vals)
            w =open(output+ '/jbrowse/' + strain + '/{}_{}_reference_orfs.wig'.format(strain,samp), 'w')
            w.write(header)
            w.close()



