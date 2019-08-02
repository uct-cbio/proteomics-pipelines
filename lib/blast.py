#!/usr/bin/env python3

import Bio
import numpy as np
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import tempfile
import shutil
import urllib3
import re
import subprocess
from collections import defaultdict
from natsort import natsorted

def pblast(seqstr, db, threshold = 0.04, path = '/home/thys/biotools/blast/', outfmt = 5):
    dbpath = path + db
    rec = SeqRecord(seq = Seq(seqstr), id = 'query')
    temp_dir = tempfile.mkdtemp()
    temp_fasta = temp_dir +'/temp.fasta'
    temp_out = temp_dir + '/temp.xml'
    SeqIO.write(rec, temp_fasta, 'fasta')
    blastp_cline = NcbiblastpCommandline(query = temp_fasta,
                                        db = dbpath,
                                        evalue = threshold,
                                        outfmt = outfmt,
                                        out = temp_out)
    stdout, stderr = blastp_cline()
    result_handle = open(temp_out)
    blast_record = NCBIXML.read(result_handle)
    shutil.rmtree(temp_dir)
    return blast_record

def pblast_pretty(seqstr, db, threshold = 0.04):
    rec = pblast(seqstr, db)
    output_list = []
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < threshold:
                output = [
                'Alignment E value: {}'.format(hsp.expect), 
                'length: {}'.format(alignment.length),
                'Feature: {}'.format(alignment.title)]
                output = ', '.join(output)
                output_list.append(output) 
    for i in output_list:
        return i
        break

def nblast(seqstr, db, threshold = 0.04, path = '/home/thys/biotools/cdna_blast/', outfmt = 5):
    dbpath = path + db
    rec = SeqRecord(seq = Seq(seqstr), id = 'query')
    temp_dir = tempfile.mkdtemp()
    temp_fasta = temp_dir +'/temp.fasta'
    temp_out = temp_dir + '/temp.xml'
    SeqIO.write(rec, temp_fasta, 'fasta')
    blastn_cline = NcbiblastnCommandline(query = temp_fasta,
                                        db = dbpath,
                                        evalue = threshold,
                                        outfmt = outfmt,
                                        out = temp_out)
    stdout, stderr = blastn_cline()
    result_handle = open(temp_out)
    blast_record = NCBIXML.read(result_handle)
    shutil.rmtree(temp_dir)
    return blast_record


def tblastn(seqstr, db, threshold = 0.04, path = '/home/thys/biotools/tblastn/', outfmt = 5):
    dbpath = path + db
    rec = SeqRecord(seq = Seq(seqstr), id = 'query')
    temp_dir = tempfile.mkdtemp()
    temp_fasta = temp_dir +'/temp.fasta'
    temp_out = temp_dir + '/temp.xml'
    SeqIO.write(rec, temp_fasta, 'fasta')
    blastn_cline = NcbiblastnCommandline(query = temp_fasta,
                                        db = dbpath,
                                        evalue = threshold,
                                        outfmt = outfmt,
                                        out = temp_out,
                                        strand = plus)
    stdout, stderr = blastn_cline()
    result_handle = open(temp_out)
    blast_record = NCBIXML.read(result_handle)
    shutil.rmtree(temp_dir)
    return blast_record

def nblast_pretty(seqstr, db, threshold = 0.04):
    rec = nblast(seqstr, db)
    output_list = []
    count =1
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < threshold:
                output = [
                'Alignment E value: {}'.format(hsp.expect), 
                'length: {}'.format(alignment.length),
                'Feature: {}'.format(alignment.title)]
                output = ', '.join(output)
                output_list.append(output)
                count +=1
    return output_list

def gi2up(gi_accession_number):
    try:
        url = 'http://www.uniprot.org/mapping/?from=P_GI&to=ACC&query={}'.format(gi_accession_number)
        p = urllib2.urlopen(url).read()
        splitter = 'xml:lang="en"><head><title>'
        ls = p.split(splitter)[1].split(' ')[0].split(':')[1]
        url2 = 'http://www.uniprot.org/uniprot/?query=yourlist:{}&sort=yourlist:{}&columns=yourlist%28{}%29,id%2Centry%20name%2Creviewed%2Cprotein%20names%2Cgenes%2Corganism%2Clength%2Cexistence%2Ccomment%28PATHWAY%29%2Cgo%2Cgo%28biological%20process%29%2Cgo%28molecular%20function%29%2Cgo%28cellular%20component%29%2Cgo-id'.format(ls,ls,ls)
        x = urllib2.urlopen(url2).read()
        datum = x.split('class="addRemoveColumn mid"')[0].split('</script></td></tr></thead><tbody><tr ')[1]
        datum = datum.split('class=')
        biorec = {}
        biorec['Entry (Uniprot)'] = datum[5].split('"entryID"><a href="/uniprot/')[1].split('"')[0]
        biorec['Entry Name (Uniprot)'] = datum[6].split('>')[1].split('<')[0]
        biorec['Protein Names'] = datum[11].split('title="')[1].split('"')[0]
        biorec['Gene Names'] = ''.join(''.join('>'.join(datum[14].split('>')[1:]).split('<strong>')).split('</strong>')).split('</div>')[0].strip() 
        biorec['Organism'] = ''.join(datum[15].split('">')[2:]).split('</a><')[0]
        biorec['Length'] = datum[16].split('>')[1].split('<')[0]
        biorec['Protein Existence'] = datum[17].split('>')[1].split('<')[0] 
        biorec['Pathway'] = datum[18].split('td>')[1].split('<td')[0] 
        go = datum[19].split('<td style=""')[0].split('</a>')
        biorec['Gene Ontology (GO)'] = [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        go = datum[20].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology (biological process)'] = [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        #return entryID1, entryID2, protein_name, gene_names, organisms
        go = datum[21].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology (molecular function)'] =  [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        go = datum[22].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology (cellular component)'] = [i.split('>')[-1] for i in go if len(i.split('>')[-1]) > 0]
        go = datum[23].split('<td style=""')[0].split('</a>')
        biorec['Gene ontology IDs'] = [i.split('>')[-1] for i in go if ((len(i.split('>')[-1]) > 0) and (i.split('>')[-1] != '<td ')) ] 
        return biorec
    except:
        return None

def bparse(rec, threshold = 0.04):
    output_list = []
    for alignment in rec.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < threshold:
                output = ['****Alignment****\n',
                    'Sequence: {}'.format(alignment.title),
                    'Length: {}'.format(alignment.length),
                    'E value: {}'.format(hsp.expect),
                    'HSP Query:   {}'.format(hsp.query),
                    'HSP Match:   {}'.format(hsp.match),
                    'HSP Subject: {}'.format(hsp.sbjct)]
                output = '\n'.join(output)
                output_list.append(output)
    return output_list

class blast2genome:
    def __init__(self, fasta_id, query_seq, match_seq, hit_seq, hit_from, hit_to, hit_strand, amino_offset=0):
        #amino offset is the position in hit_seq from which to start searching for match segments (ie peptide alignment to ORF protein. all seqs must be amino acids, hit_from and hit_to are genomics coords, and query_from and query_to are amino acid coords (1 based).
        self.fasta_id = fasta_id
        self.query_seq = query_seq
        self.match_seq = match_seq
        self.hit_seq = hit_seq
        self.hit_to = hit_to
        self.hit_from = hit_from
        self.segments = self.match_seq.replace('+',' ').split()
        target = self.hit_seq.replace('-','')
        self.records = []
        amino_offset = 0
        self.hit_strand= hit_strand
        self.alignment_start = np.inf
        self.alignment_end = 0

        for segment in self.segments: 
            newstart, newend, newoffset = self.evaluate_start(segment, self.hit_seq, self.hit_from, self.hit_to, self.hit_strand, amino_offset)
            start_position = newstart
            end_position = newend
            
            if start_position < self.alignment_start:
                self.alignment_start = start_position
            if end_position > self.alignment_end:
                self.alignment_end = end_position

            strand=self.hit_strand
            record = (segment, start_position, end_position, strand)
            self.records.append(record)
            amino_offset = newoffset

    def evaluate_start(self, segment, hit_seq, hit_start, hit_end, hit_strand, hit_amino_offset):
        starts = [m.start() for m in re.finditer('(?={})'.format(segment), hit_seq)]
        assert hit_amino_offset % 1 == 0
        starts =  [ i for i in starts if i >= hit_amino_offset ]
        start = starts[0]
        end = start + len(segment)
        new_hit_amino_offset = end
        if hit_strand=='+':
            newstart = ((start ) * 3) + hit_start
            newend   = (((end-start )* 3) + newstart )-1
            return newstart, newend, new_hit_amino_offset
        elif hit_strand=='-':
            newend = hit_end - ((start ) * 3) 
            newstart   = newend - (((end-start )* 3) ) + 1
            return newstart, newend, new_hit_amino_offset

class pblast_mult:
    def __init__(self, query_recs, target_recs, outdir, evalue=0.0000000001, num_threads=5, max_target_seqs=500, max_hsps=1, matrix='BLOSUM62', gapopen=11, word_size=3, gapextend=1, comp_based_stats=2, window_size=15, threshold=11, seg='no', pairwise=True, perc_identity_cutoff=90):
        
        self.query_recs = query_recs
        self.target_recs = target_recs
        self.outdir = outdir
        self.num_threads = num_threads
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.max_hsps = max_hsps
        self.matrix = matrix
        self.gapopen = gapopen
        self.word_size =  word_size
        self.gapextend = gapextend
        self.comp_based_stats = comp_based_stats
        self.window_size = window_size
        self.threshold = threshold
        self.seg = seg
        
        blastdir = outdir + '/db'
        self.outfile = outdir + '/results.xml'
        os.mkdir(blastdir)
        SeqIO.write( target_recs, blastdir + '/target.fasta','fasta')
        cmd="cd {} && makeblastdb -in target.fasta -dbtype 'prot' -out 'target'".format(blastdir)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        assert process.returncode == 0
        cmd ="cd {} && blastp -query target.fasta -db target -out={} -outfmt 5 -max_target_seqs {} -max_hsps {} -num_threads {} -evalue {} -matrix {} -gapopen {} -word_size {} -gapextend {} -comp_based_stats {} -window_size {} -threshold {} -seg {}".format(blastdir, self.outfile, self.max_target_seqs, self.max_hsps, self.num_threads, self.evalue, self.matrix, self.gapopen, self.word_size, self.gapextend, self.comp_based_stats, self.window_size, self.threshold, self.seg)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) 
        process.wait()
        assert process.returncode == 0
        result_handle = open(self.outfile)
        self.blast_records = list(NCBIXML.parse(result_handle))
        result_handle.close() 
        self.aln_dict = defaultdict(set)

        for rec in self.blast_records:
            for aln in rec.alignments:
                for hsp in aln.hsps:
                    perc = hsp.identities/float(hsp.align_length) * 100
                    if perc >= perc_identity_cutoff:
                        self.aln_dict[rec.query.split()[0]].add(aln.hit_def.split()[0])

        if pairwise == True:
            new_dict = defaultdict(set)
            for rec in self.aln_dict:
                for match in self.aln_dict[rec]:
                    if rec in self.aln_dict[match]:
                        new_dict[rec].add(match)
            self.aln_dict = new_dict

        self.sorted_map = {}
        for rec in self.query_recs:
            recs = natsorted(list(self.aln_dict[rec.id]))
            self.aln_dict[rec.id] = recs
            recs = ';'.join(recs)
            self.sorted_map[rec.id] = recs







