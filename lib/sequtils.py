#/usr/bin/env python
import rpy2.robjects as ro
import collections
from collections import Counter
import multiprocessing
import pandas as pd
import numpy as np
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, translate
from Bio.Data import CodonTable
from numpy import NaN
import Bio; from Bio import SeqIO
import pandas as pd
import os
import natsort
import subprocess
import tempfile
import collections; from collections import defaultdict
from io import StringIO
import re
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import Phylo
import contextlib
import concurrent.futures
import sys
import algo
import blast

def alt_starts_recs(records, starts = ['ATG','GTG','TTG'], peptide_length=20):
    all= []
    for seq in records:
        alt_seqs = alt_tss(seq, starts = starts)
        all += alt_seqs
    
    newall = []

    for rec in all:
        temp = str(rec.seq)
        if len(temp) >= peptide_length:
            newall.append(rec)
    
    return newall

def Three_Frames(sequence): # All three frames on the forward strand (returns three sequences and includes possible stop codons)
    sequence_length = len(sequence)
    Frame1 = sequence[0:int(sequence_length/3)*3]
    Frame2 = sequence[1:int((sequence_length+2)/3)*3-2]
    Frame3 = sequence[2:int((sequence_length+1)/3)*3-1]
    Frames = (Frame1, Frame2, Frame3)
    return Frames

def Six_Frames_and_Reverse_Complement(sequence):# Do three frame translation on forward and reverse-complement of forward strand
    seq = Seq(sequence)
    reversed_seq = seq.reverse_complement()
    rev_seq = str(reversed_seq).rstrip()
    return Three_Frames(sequence) + Three_Frames(rev_seq), rev_seq

def Stop_List(table): # Return possible stop codons given translation table
    Codon_Table = table
    codon_table = CodonTable.unambiguous_dna_by_id[Codon_Table]
    Stop_Codons = codon_table.stop_codons
    return Stop_Codons

def Start_List(table): # Possible starts
    Codon_Table = table
    codon_table = CodonTable.unambiguous_dna_by_id[Codon_Table]
    Start_Codons = codon_table.start_codons
    return Start_Codons

def codon_list(table):
    codons = CodonTable.unambiguous_dna_by_id[table].forward_table.keys()
    return codons

def six_frame(Genome, assembly_name, contig_name, table, peptide_length, codons = ['ATG','GTG','TTG'], translated = True, contig_id='sixframe', methionine_start=True, translate_stop=False): 
    Minimum_Peptide_Length = peptide_length
    Codon_Table = table
    Start_Codons = codons

    codon_table = CodonTable.unambiguous_dna_by_id[Codon_Table]
    
    Stop_Codons = codon_table.stop_codons
    
    if (translated==True) and (methionine_start==True):
        for i in Start_Codons:
            assert i in codon_table.start_codons                       # make sure the selected start codons are found in the selected translation table

    Six_Frame, Reverse_Complement_Genome = Six_Frames_and_Reverse_Complement(Genome)

    ORFs_by_Frame = {}

    for i in Six_Frame:
        count = 1
        Frame = Six_Frame.index(i)+1
        [a] = [Frame - 3 if Frame > 3 else Frame] # Adjust for frameshift to find index.
        b = a-1         # Shift the value of y to == index of last nucleotide (-1), and adjust for frame shift with 'a'.

        ORFs_by_Frame[Frame] = []
        starts = False
        index = 0
        ORF = []
        Start = NaN

        if Frame <= 3:
            String = 1

        elif Frame >3:
            String = -1

        while index < (len(i))-2:
            x = index
            y = index + 3
            if starts == True:
                if i[x:y] in Stop_Codons:
                    ORF.append(i[x:y])
                    temp = ''.join(ORF)
                    if Frame <= 3:
                        location_start = Start + a
                        location_end = y + b
                    elif Frame > 3:
                        location_start = len(Genome)-(y + b)+1 # Index by direct strand
                        location_end = len(Genome)-(Start+a)+1  # Index by direct strand

                    ORFs_by_Frame[Frame].append((count,Frame, String, location_start, location_end, temp))
                    count += 1
                    ORF = []
                    starts = False
                elif x == (len(i) - 3):
                    ORF.append(i[x:y])
                    if (len(Genome)-(y+b))<=0:
                        assert len(Genome)-(y+b) == 0
                    temp = ''.join(ORF)
                    if Frame <= 3:
                        location_start = Start + a
                        location_end = y + b
                    elif Frame > 3:
                        location_start = len(Genome)-(y + b)+1 # Index by direct strand
                        location_end = len(Genome)-(Start+a)+1  # Index by direct strand
                    ORFs_by_Frame[Frame].append((count,Frame, String,location_start, location_end, temp))
                    count += 1
                    ORF = []
                    starts = False
                else:
                    ORF.append(i[x:y])
            if starts == False:
                if i[x:y] in Start_Codons:
                    starts = True
                    Start = x
                    ORF.append(i[x:y])
            index += 3
    
    Six_Frame = []
    rec_count = 1
    for i in ORFs_by_Frame:
        for j in ORFs_by_Frame[i]:                   #check by frame
            if j[5][-3:] in Stop_Codons:             #check if strand ends in a stop codon
                Seq_len = (len(j[5])-3)/3                  #subtract 3 if above true, and divide by three to get length
            elif j[5][-3:] not in Stop_Codons:       #if last codon is not a stop codon, then just divide by 3 to get petide length
                Seq_len = len(j[5])/3
            if Seq_len >= Minimum_Peptide_Length:    #compare to minimum peptide length
                seq = str(Seq(j[5]))
                if j[1] < 4:
                    frame_direction = '+'
                elif j[1] > 3:
                    frame_direction = '-'
                
                seq = Seq(seq)
                if translated == False:
                    record = SeqRecord(seq = seq,id = "{}_{}|{}_{}_recno_{}|({}){}:{}".format(assembly_name, contig_name, assembly_name, contig_name,  float(rec_count), frame_direction,j[3],j[4]), description='Six_Frame_ORF')

                elif translated == True: 
                    prot_seq = str(translate(seq, table=table))
                    if (prot_seq[-1] == '*') and (translate_stop == False):            
                        prot_seq = prot_seq[:-1]  

                    if (methionine_start == True):
                        seqstart =seq[:3]
                        if seqstart in Start_Codons:
                            prot_seq = 'M' + prot_seq[1:]
                    
                    prot_seq = Seq(prot_seq)
                    record = SeqRecord(prot_seq,id = "{}_{}|{}_{}_recno_{}|({}){}:{}".format(assembly_name, contig_name, assembly_name, contig_name, float(rec_count), frame_direction,j[3],j[4]),description='Six_Frame_Translated_ORF')
                
                if frame_direction == '+':
                    preceding = Genome[j[3]-4:j[3]-1]
                if frame_direction == '-':
                    preceding = str(Seq(Genome[j[4]: j[4]+3]).reverse_complement())
                
                if len(preceding) == 0:
                    preceding="None"

                record.description = record.description + ' Preceding_codon={}'.format(preceding)

                Six_Frame.append(record)
                rec_count += 1
    return Six_Frame

def get_frame(seq, genome_seq, stops): #'''BUG'''translated seqrecord with coords, genomic seq record, possible stops eg ['TAA', 'TAG', 'TGA']
    '''In put a translated seq_record with description with format eg. (+)4454545:4454548 and a single genomic seq record to output the stop-to-stop nucleic acid sequence that contains that seq_record. NB input a list of possible stop codons, eg. stops = ['TAA', 'TAG', 'TGA'] (for bacterial table 11)'''
    frames = []
    coords = seq.description.split()[1].split('|')
    for i in coords:
        coord = i[3:].split(':')
        start = int(coord[0])
        end = int(coord[1])
        strand = i[:3]
        if strand == '(-)':
            while str(genome_seq.seq[end:end+3].reverse_complement()) not in stops:
                end = end + 3
            end = end + 3
            seq_ = genome_seq.seq[start-1:end].reverse_complement()
            description = '(-){}:{}'.format(start,end)
        elif strand == '(+)':
            while str(genome_seq.seq[start-4:start-1]) not in stops:
                start = start - 3
            start = start - 3
            seq_ = genome_seq.seq[start-1:end]
            description = '(+){}:{}'.format(start, end)
        record = SeqRecord( seq_, id = seq.id+' Genomic sequence option (stop to stop) at {}'.format(description))
        frames.append(record)
    return frames

def gmhmmp(model, path, sequence):   # --mode (native, combined, heuristic), --org (path), seq_record object
    temp_dir = tempfile.mkdtemp()
    temp_out = temp_dir + '/out'
    temp_faa = temp_dir +'/faa'
    temp_fnn = temp_dir + '/fnn'
    f = temp_dir+'/fasta'                 #temporary file to write the fasta record to
    # SeqIO.write(sequence, f, 'fasta')
    command = 'gmhmmp -m {}{} {} -r -a -s d -n -g 11 -A {}gmhmmp_proteins.fasta'.format(path, model, f, temp_dir)
    process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    process.wait()
    fasta = list(SeqIO.parse(temp_dir +'gmhmmp_proteins.fasta', 'fasta'))
    return fasta                                                      # protein sequence of predicted genes

def gmsn(genome, outpath = '/tmp'):  #prokaryotes only, requires a seq-record object for query sequence
    if outpath[-1] =='/':
        outpath = outpath[:-1]
    temp_dir = tempfile.mkdtemp()
    f = temp_dir+'/fasta'                 #temporary file to write the fasta record to
    SeqIO.write(genome, f, 'fasta')
    command = 'cd {} && gmsn.pl --prok --faa {} && rm -rf {}/GeneMarkS_output_files && mv {} {}/GeneMarkS_output_files'.format(temp_dir, f, outpath,temp_dir, outpath)
    process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
    process.wait()
    fasta = list(SeqIO.parse('{}/{}/fasta.faa'.format(outpath, 'GeneMarkS_output_files'), 'fasta'))
    return fasta                                                      # protein sequence of predicted genes

def m_start_list(proteins):     #inport a list of amino acid sequences, and return all starting with M
    m = []
    for protein in proteins:
        seq = 'M'+str(protein.seq)[1:]
        record = SeqRecord(seq = Seq(seq), id = protein.id, description = protein.description)
        m.append(record)
    return m

def sf_contigs(contigs, assembly_name, table=11, peptide_length = 20, codons = 'All', translated = False, methionine_start=False): # contigs fasta file, output genomic six frame sequences
    
    if (codons == 'All'):
        codons = codon_list(table)
        if translated ==True:
            assert methionine_start==False # can only translate known start codons with methionine
    six_frame_seqs = []
    for i in contigs:
        cid = '_'.join(i.id.split('|'))  #get the contig id
        #print(cid)
        six = six_frame(str(i.seq).upper(),assembly_name= assembly_name, contig_name = cid, table = table, peptide_length = peptide_length, codons = codons, translated=translated, methionine_start = methionine_start)
        
        for j in six:
            #print(j)
            six_frame_seqs.append(j)
    return six_frame_seqs

def translate_start(rec, table=11, starts = ['ATG','GTG','TTG']): # list of Met-initiatiator codons    
    seq = rec.seq
    id = rec.id
    desc = rec.description

    if seq[:3] in starts:
        translated = 'M' + str(translate(seq, table = table))[1:]
    else:
        translated=str(translate(seq,table=table, cds=False))

    if translated.endswith('*'):
        translated = translated[:-1]

    newrec = SeqRecord(id = id, description = desc, seq = Seq(translated))
    return newrec


def alt_starts_recs(records, starts = ['ATG','GTG','TTG'], table=11, translated=False):
    all= []
    for seq in records:
        alt_seqs = alt_tss(seq, starts = starts, table=table, translated=translated)
        all += alt_seqs
    return all

# returns a lists of records of all alternative start codons for an ORF (including the ORF given)
def alt_tss(seq_record, table=11, starts = ['ATG', 'GTG', 'TTG'], translated=True):
    alt_seqs = []
    alt_seq_count = 1
    id_   = seq_record.id
    seq_  = str(seq_record.seq)
    desc_ = seq_record.description.split()[1]
    num_  =  int(id_.split('|')[1].split('recno_')[1].split('.')[0])

    coords  =id_.split('|')[2].split('(')[1]
    base_id = id_.split('recno_')[0]
    strand = coords[0]
    
    coords = coords.split(')')[1].split(':')

    start  = int(coords[0])
    end    = int(coords[1])
    orig_rec = SeqRecord(seq = Seq(seq_) , id = id_, description = 'Six_Frame_ORF ' + desc_)
    alt_seqs.append(orig_rec)
    count = 3
    while count <= len(seq_)-3:
        codon = seq_[count:count +3]
        if codon in starts:
            alt_seq = Seq(seq_[count:])
            prev_codon =seq_[count-3:count ]
            newdesc = 'Alternative_TSS_Sequence ' + desc_.split('=')[0] +'=' + prev_codon
            assert len(alt_seq) % 3.0 == 0
            if strand == '+':
                new_start = start + count
                new_end   = end
            elif strand == '-':
                new_start = start
                new_end = end - count
            new_rec = SeqRecord(id = base_id + 'recno_{}.{}|({}){}:{}'.format(str(num_), str(alt_seq_count), strand, new_start, new_end),
                                description = newdesc,
                                seq = alt_seq)
            alt_seqs.append(new_rec)
            alt_seq_count += 1
        count += 3

    if translated == True:
        translated_recs = []
        for rec in alt_seqs:
            rec = translate_start(rec, table=table, starts=starts)
            translated_recs.append(rec)
        alt_seqs = translated_recs
    return alt_seqs

class gff3:
    '''Basic class to parse GFF3'''
    gffcols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    def __init__(self, GFF3):
        self.table = pd.read_csv(os.path.abspath(GFF3),sep='\t',comment='#')
        self.table = self.table.iloc[:,:9]
        self.table.columns = self.gffcols
        #self.attributes_columns()

    def attributes_columns(self):
        pos_cols = ['ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Derives_from','Note', 'Dbxref', 'Ontology_term', 'Is_circular']
        for pos in pos_cols:
            self.table['_attributes.{}'.format(pos)] = self.table['attributes'].apply(lambda x : x.split('{}='.format(pos))[1].split(';')[0] if '{}='.format(pos) in x else None)
        self.table = self.table.dropna(how='all',axis=1)

    def split_attributes(self):
        for row in self.table.iterrows():
            attr = row[1]['attributes']
            fields = attr.split(';')
            for f in fields:
                f = f.split('=')
                name = f[0]
                field = f[1]
                if (name =='Dbxref'):
                    subfields = field.split(',')
                    for sub in subfields:
                        sub = sub.split(':')
                        sname = sub[0]
                        sfield = sub[1]
                        self.table.loc[row[0], '_attributes.{}.{}'.format(name, sname)] = sfield
                else:
                    self.table.loc[row[0], '_attributes.{}'.format(name)] = field

    def entryCount(self):
        return len(self.table)
    
    def expand_table(self, mapping):
        new_df = pd.DataFrame()
        new_df['_mapped.id'] = pd.Series(list(mapping.keys()))
        new_df['seqid'] = new_df['_mapped.id'].apply(lambda x : mapping[x])
        
        merged = pd.merge(self.table, new_df)
        merged['seqid'] = merged['_mapped.id'].apply(lambda x : x.split('|')[1])
        self.table=merged
    


class vcf:
    '''Basic class for vcf data'''
    def __init__(self, VCF):
        path = os.path.abspath(VCF)
        w = open(path).readlines()
        headers = [line.strip() for line in w if line.startswith('#CHROM')][0].split('\t')
        self.table = pd.read_csv(path,sep='\t',comment='#', header=None)
        self.table.columns = headers

    def entryCount(self):
        return len(self.table)

class genome:
    '''Basic class of genome object'''
    def __init__(self, FASTA):
        self.reflist = list(SeqIO.parse(os.path.abspath(FASTA),'fasta'))    
        self.refdict = {}
        for seq in self.reflist:
            self.refdict[seq.id] = seq
        
class variant_genome:
    def __init__(self, VCF, FASTA, GFF3 = None, name=None):  # vcf file, name, genome,gff3
        self.ref_genome = genome(FASTA)
        self.vcf = vcf(VCF)
        self.name = name
        self.var_genome = {} 
        self.var_mapping = {}
        self.ref_genome = SeqIO.to_dict(SeqIO.parse(os.path.abspath(FASTA),'fasta'))

        for chrom in self.ref_genome:      
            start_pos = 0
            step = 0
            seen_pos = []
            rec = self.ref_genome[chrom]  
            rseq = str(rec.seq)

            _ = self.vcf.table
        
            new_seq_list = []
            new_seq_map  = {}
            
            subset = _[_['#CHROM']==chrom]  # get only vcf rows relevant to the fasta id
            for row in subset.iterrows():
                pos = row[1]['POS'] - 1
                assert pos not in seen_pos
                assert pos >=start_pos
                seen_pos.append(pos)
                ref = row[1]['REF']
                alt = row[1]['ALT']
                alt = alt.split(',')
                assert not isinstance(alt, str)
                alt = alt[0]
                assert ',' not in ref
                assert ',' not in alt

                for p in range(start_pos, pos):
                    new_seq_map[p] = p + step
                    new_seq_list.append(rseq[p])

                new_seq_map[pos] = pos + step
                new_seq_list.append(alt)
                step += len(alt) - len(ref) 
                start_pos = pos + len(ref) 
            
            for p in range(start_pos, len(rseq)):
                new_seq_map[p] = p + step
                new_seq_list.append(rseq[p])

            vseq = ''.join(new_seq_list)
            vid  = 'variant_' + rec.id
            vdescription  = 'variant_' + rec.description
            vsrecord = SeqRecord(seq = Seq(vseq), id = vid, description = vdescription)

            self.var_genome[chrom] = vsrecord
            self.var_mapping[chrom]= new_seq_map


    def var_feature_table(self, GFF3):
        gff = gff3(os.path.abspath(GFF3))
        cds_table = gff.table[gff.table['type']=='CDS']
        start_affected = 0
        end_affected = 0

        for row in cds_table.iterrows():
            seqid = row[1]['seqid']
            start = row[1]['start'] -1
            end = row[1]['end']
            strand = row[1]['strand']

        
            ref_genome = str(self.ref_genome[seqid].seq)
            var_genome = str(self.var_genome[seqid].seq)
            mapping = self.var_mapping[seqid]
            
            ref_cds = ref_genome[start : end]
            
            var_start = None
            var_end = None
            
            if start in mapping:
                var_start = mapping[start]
            if end in mapping:
                var_end   = mapping[end]

            if (var_start != None) and (var_end != None):
                var_region = var_genome[var_start:var_end]
                if strand == '+':
                    var_cds = Seq(var_region)
                elif strand =='-':
                    var_cds = Seq(var_region).reverse_complement()
                
                var_cds_len = len(str(var_cds))
                if var_cds_len % 3 != 0:
                    print(row[1]['_attributes.Name'])
            
            if var_start == None:
                start_affected += 1
            if var_end == None:
                end_affected += 1
        print('Start affected: {}'.format(start_affected))
        print('End affected: {}'.format(end_affected))

def get_end(cigar,start): # subtract 1 from start
    flags={'D':1,'I':0,'M':1,'S':0}
    end = start 
    num =[]
    for char in cigar:
        if char not in flags:
            num.append(char)
        elif char in flags:
            new_num = int(''.join(num))*flags[char]
            num=[]
            end += new_num
    return end

def sam_parse(sam):
    count = 0
    sequences = []
    orfs = defaultdict(list)
    seqdict = defaultdict(list)
    sam = sam.split('\n')
    if sam[0].startswith('@'):
        sam = sam[1:]
        #reference_title = sam[0].split('\tSN:')[1].split('\tLN:')[0].split('.')[0]+'_aligned.fasta'
        #length = sam[0].split('\tLN:')[1]

    for aln in sam:
        if len(aln.split('\t')) > 1:
            datum =  aln.split('\t')
            q_name = datum[0]
            sam_flag = datum[1]
            ref_name = datum[2]
            start = datum[3]
            
            mapq = datum[4]
            cigar = datum[5]
            
            end = get_end(cigar, int(start)-1)

            seq = Seq(datum[9])
            
            if sam_flag =='16':
                seq = seq.reverse_complement()
                strand = '-'
            else:
                strand = '+'
            
            header = '({}){}:{}'.format(strand, start, end)
            mapping = ref_name+'_'+start+'_'+sam_flag+'_'+cigar+'_'+mapq+'_'+header
            seqdict[q_name].append(str(seq))

            if sam_flag == '4':
                record = SeqRecord(seq=seq,id=q_name,description='Unmapped')
                sequences.append(record)
            else:
                record = SeqRecord(seq=seq,id=q_name,description=mapping)
                sequences.append(record)

                #orfs[q_name].append(record)
    #new_seqs = []
    #for i in orfs:
    #    description = []
    #    id = i
    #    assert len(list(set(seqdict[i])))==1
    #    seq = Seq(seqdict[i][0])
    #    for j in orfs[i]:
    #        description.append(j.description)
    #        new_seqs.append(j)   
    #    record = SeqRecord(seq=seq, description=' '.join(description),id=id)
    #   sequences.append(record)
    return sequences


class peptides2genome:
    '''Basic class for peptides to genome mapping - list of BioPython seqrecords for genome contigs'''
    def __init__(self, genome, assembly_name, translation_table, peptides_list, outdir, threads=1, start_codons=['ATG','GTG','TTG']):
        self.threads=threads
        self.translation_table = translation_table
        self.stop_codons = Stop_List(translation_table)
        self.start_codons = start_codons
        self.outdir = outdir
        
        os.mkdir(self.outdir +'/sixframe_blast')
        #os.mkdir(self.outdir +'/sixframe_blast')

        mstart_peps=[i[1:] for i in peptides_list if i.startswith('M')]
        non_mstart_peps=[i for i in peptides_list if not i.startswith('M')]

        self.nonMetTrie = algo.Trie(non_mstart_peps)
        self.MetTrie = algo.Trie(mstart_peps) #excluding N-term M
        self.peptide2orf = defaultdict(list)

        orfs = sf_contigs(genome, assembly_name = assembly_name, table=translation_table, codons='All', peptide_length=1, translated=False )
        
        self.orfs = orfs

        orf_counts = defaultdict(list)  # This baby is to hold the values of the id's of ORFs that occur one or more times in the genome
        preceding_codon = defaultdict() 
        
        orf_seq_list= [str(rec.seq) for rec in orfs]
        #orf_trans_list = [str(translate(rec.seq, table=self.translation_table, cds=False)) for rec in orfs]

        self.trans_orfs = []
        orf_trans_list = []

        for orf in self.orfs:
            newrec = SeqRecord(id=orf.id, seq = translate(orf.seq, table = self.translation_table, cds=False))
            self.trans_orfs.append(newrec)
            orf_trans_list.append(str(newrec.seq))

        orf_id_list = [rec.id for rec in orfs]

        preceding_codons = [rec.description.split('Preceding_codon=')[1].split('.')[0] for rec in orfs]
        
        most_upstream_peptide_pos=[]
        trie_peptide_dict = []

        orf_series = pd.Series(orf_id_list, name='ORF_id')
        orf_series.index=orf_seq_list
        orf_series.index.name = 'ORF_sequence'
        
        orf_df = orf_series.reset_index()
        #orf_df['ORF_genome_count'] = orf_df['ORF_ids'].apply(lambda x : len(x.split(';')))
        orf_df['ORF_translation'] = pd.Series(orf_trans_list) # DEV change the way in-sequence stops are translated (eg. use selenocysteine)
        orf_df['ORF_preceding_codons'] = pd.Series(preceding_codons)
        print("Calculating cterm fragments")
        orf_df['ORF_Cterm_fragment'] = orf_df['ORF_sequence'].apply(lambda x : self.cterm_fragment(x))
        print("Calculating neterm fragments")
        orf_df['ORF_Nterm_fragment'] = orf_df['ORF_preceding_codons'].apply(lambda x : self.nterm_fragment(x))
        orf_df['Assembly'] = assembly_name
        print("Calculating most upstream start")
        orf_df['Most_Upstream_Inferred_Start'] = orf_df.apply(self.most_upstream_start, axis=1) 
        orf_df = orf_df[orf_df['Most_Upstream_Inferred_Start'].notnull()]
        print("Calculating most upstream inferred codon")
        orf_df['Most_Upstream_Inferred_Codon'] = orf_df.apply(self.most_upstream_codon, axis=1)
        orf_df['Most_Upstream_Inferred_Translation'] = orf_df.apply(self.most_upstream_translation, axis=1)
        
        self.mapped_orfs = [orf for orf in self.orfs if orf.id in list(orf_df['ORF_id']) ]
        self.mapped_trans_orfs =  [ orf  for orf in self.trans_orfs if orf.id in list(orf_df['ORF_id']) ]
        
        blasted_orfs = blast.pblast_mult(self.mapped_trans_orfs, self.trans_orfs, evalue=0.0000000001, outdir=self.outdir + '/sixframe_blast')
        orf_df['Orthologous_ORFs'] = orf_df['ORF_id'].apply(lambda x : blasted_orfs.sorted_map[x])
        orf_df['Orthologous_ORFs_Cluster_Size'] = orf_df['ORF_id'].apply(lambda x : len(blasted_orfs.aln_dict[x]))

        self.orf_trans_list = orf_df['ORF_translation'].tolist()
        self.orf_df = orf_df
        print("ready to analyze the peptides")
        self.peptides = self.process_peptides(peptides_list)
    
    def most_upstream_start(self, df):
        
        pep_positions = []
        non_met_position_aminos = defaultdict(list)

        orf_trans = df['ORF_translation']
        orf_nucs = df['ORF_sequence']
        orf_id = df['ORF_id']

        nonMetTM = algo.TrieMatch(self.nonMetTrie, orf_trans)
        pep_positions += nonMetTM.trie_matching()
        for pep in nonMetTM.trie_export():
            self.peptide2orf[pep].append(orf_id)
            non_met_positions= [m.start() for m in re.finditer('(?={})'.format(pep), orf_trans)]
            for pos in non_met_positions:
                non_met_position_aminos[pos].append(pep[0])
        
        MetTM = algo.TrieMatch(self.MetTrie, orf_trans)
        met_peptides = MetTM.trie_export()
        for metpep in met_peptides:
            met_positions= [m.start() for m in re.finditer('(?={})'.format(metpep), orf_trans)]
            for metpos in met_positions:
                if metpos != 0:
                    first_pos = metpos-1
                    nuc_pos= first_pos * 3
                    first_codon = orf_nucs[nuc_pos:nuc_pos + 3]
                    first_amino = orf_trans[first_pos]
                    if (first_amino == 'M') or (first_codon in self.start_codons):
                        pep_positions.append(first_pos )
                        #position_aminos[first_pos].append('M')
                        self.peptide2orf["M" + metpep].append(orf_id)

        if len(pep_positions)> 0:
            most_upstream_codon_pos  = min(pep_positions) * 3
            most_upstream_codon = orf_nucs[most_upstream_codon_pos:most_upstream_codon_pos + 3]
            while ((most_upstream_codon not in self.start_codons) and (most_upstream_codon_pos != 0)) or ((most_upstream_codon_pos/3 in non_met_position_aminos) and (most_upstream_codon_pos !=0)): # Keep searching upstream if a non-Met amino was mapped to the most upstream position.
                most_upstream_codon_pos -= 3
                most_upstream_codon = orf_nucs[most_upstream_codon_pos:most_upstream_codon_pos + 3]
            #print(most_upstream_codon)
            #most_upstream_aminos = []
            #if most_upstream_codon_pos/3 in non_met_position_aminos:
            #    most_upstream_aminos += non_met_position_aminos[most_upstream_codon_pos/3]
            return most_upstream_codon_pos + 1

    def most_upstream_codon(self, df):
        upstream_start = df['Most_Upstream_Inferred_Start']
        nuc_pos = int(upstream_start) - 1
        orf_sequence = df['ORF_sequence']
        upstream_codon = orf_sequence[nuc_pos : nuc_pos + 3]
        return upstream_codon

    def most_upstream_translation(self, df):
        upstream_start = int(df['Most_Upstream_Inferred_Start'])
        aa_start  = int((upstream_start-1)/3)
        pstream_codon = df['Most_Upstream_Inferred_Codon']
        orf_trans = df['ORF_translation']
        most_upstream = orf_trans[aa_start:]
        return most_upstream        

    def ORF_set_count(self):
        print(self.orf_df[self.orf_df['ORF_Cterm_fragment']=='True'].head())
        return len(self.orf_df)
    
    def cterm_fragment(self,orf_nucs):
        last_codon = orf_nucs[-3:]
        if last_codon in self.stop_codons:
            return 'False'
        else:
            return 'True'

    def nterm_fragment(self, val): #supply a list of preceding codons/fragments of codons.
        for v in val.split(';'):
            if v in self.stop_codons:
                return 'False'
            else:
                return 'True'
    
    def locate_peptide(self, prot, peptide):
        shift = 0
        if peptide.startswith('M'):
            temp_pep = peptide[1:]   # because of non-standard start codons, disregard 'M'
            shift -= 1 # position adjustment for disregarding 'M'
        else:
            temp_pep = peptide
        starts  = [k.start() for k in re.finditer('(?={})'.format(temp_pep), prot)]
        starts = [i + 1 + shift for i in starts] # 1 based numbering
        starts = ';'.join([str(_) for _ in starts])
        return starts
    
    def tryptic_nterm_peptide(self, df):
        starts = df['Peptide_starts']
        starts = [int(i) for i in starts.split(';') if i != '']
        prot  = df['ORF_translation']

        cleavages = ['R', 'K']
        tryptic = []
        starts = [i-1 for i in starts ]
    
        for s in starts:
            before = prot[s-1:s]
            if before in cleavages:
                tryptic.append('True')
            else:
                tryptic.append('False')
        tryptic = ';'.join(tryptic)
        return tryptic

    def tryptic_cterm_peptide(self, df):
        starts = df['Peptide_starts']
        starts = [int(i) for i in starts.split(';') if i != '']
        prot  = df['ORF_translation']
        peptide = df['Peptide_sequence']
        
        cleavages = ['R', 'K']
        tryptic = []
        starts = [i-1 + len(peptide) for i in starts]
        for s in starts:
            before = prot[s-1:s]
            if before in cleavages:
                tryptic.append('True')
            else:
                tryptic.append('False')
        tryptic = ';'.join(tryptic)
        return tryptic
    
    def first_codon_peptide(self, df):
        starts = df['Peptide_starts']
        starts = [int(i) for i in starts.split(';') if i != '']
        orf = df['ORF_sequence']
        starts = [(i-1)*3 for i in starts]
        codons=[]
        for s in starts:
            codon = orf[s:s+3]
            codons.append(codon)
        codons = ';'.join(codons)
        return codons

    def previous_codon_peptide(self, df):
        starts = df['Peptide_starts']
        starts = [int(i) for i in starts.split(';') if i != '']
        orf = df['ORF_sequence']
        starts = [(i-1)*3 for i in starts]
        codons=[]
        for s in starts:
            codon = orf[s-3:s]
            codons.append(codon)
        codons = ';'.join(codons)
        return codons

    def strfind(self, translated, peptide):
        if peptide in translated:
            return True
        else:
            return False

    def amino_after(self, df):
        peptide = df['Peptide_sequence']
        starts = df['Peptide_starts']
        starts = [int(i) -1 + len(peptide) for i in starts.split(';') if i != '']
        translated = df['ORF_translation']
        aminos=[]
        for s in starts:
            amino = translated[s:s+1]
            if amino == '':
                amino = 'None'
            aminos.append(amino)
        aminos = ';'.join(aminos)
        return aminos

    def amino_before(self, df):
        peptide = df['Peptide_sequence']
        starts = df['Peptide_starts']
        starts = [int(i) -1 for i in starts.split(';') if i != '']
        translated = df['ORF_translation']
        aminos=[]
        for s in starts:
            amino = translated[s-1:s]
            if amino == '':
                amino = 'None'
            aminos.append(amino)
        aminos = ';'.join(aminos)
        return aminos

    def process_peptides(self, peptides_list):       
        pool = multiprocessing.Pool(self.threads)
        peptides = pd.concat(pool.map(self.peptide_df, peptides_list))
        peptides = peptides.reset_index()
        del peptides['index']
        return peptides

    def peptide_df(self, peptide):
        #print(peptide)
        ids = self.peptide2orf[peptide]
        df = self.orf_df.copy()
        df = df[df['ORF_id'].isin(ids)]
        if len(df) != 0:
            df['Peptide_sequence'] = peptide
            df['Peptide_starts'] = df['ORF_translation'].apply(lambda x : self.locate_peptide(x, peptide))
            df['Peptide_tryptic_nterm'] = df.apply(self.tryptic_nterm_peptide, axis=1)
            df['Peptide_tryptic_cterm'] = df.apply(self.tryptic_cterm_peptide, axis=1)
            df['Peptide_previous_codon'] = df.apply(self.previous_codon_peptide,axis=1)
            df['Peptide_first_codon'] = df.apply(self.first_codon_peptide,axis=1)
            df['Peptide_amino_acid_before'] =  df.apply(self.amino_before,axis=1)
            df['Peptide_amino_acid_first'] =  df['Peptide_sequence'].apply(lambda x : x[:1])
            df['Peptide_amino_acid_last'] =  df['Peptide_sequence'].apply(lambda x : x[-1:])
            df['Peptide_amino_acid_after'] =  df.apply(self.amino_after,axis=1)
            
            distinct_translated_ORF_count = len(df.drop_duplicates(['Peptide_first_codon', 'Peptide_previous_codon', 'Most_Upstream_Inferred_Translation'], keep='first'))
            df['Peptide_inferred_translated_sequence_count'] = distinct_translated_ORF_count
            mapped_orfs = len(list(set(df['ORF_id'].tolist())))
            if distinct_translated_ORF_count  == 1:
                df['Peptide_inferred_translated_sequence_specific']='+'
            elif distinct_translated_ORF_count > 1:
                df['Peptide_inferred_translated_sequence_specific']='-'
            df['Peptide_genome_ORF_count'] = mapped_orfs 
            
            # paralogs
            ortholog_df = df[df['Orthologous_ORFs_Cluster_Size'] > 0 ]

            distinct_translated_ORF_cluster_count = len(ortholog_df.drop_duplicates(['Orthologous_ORFs'], keep='first'))

            df['Peptide_translated_ORF_cluster_count'] = distinct_translated_ORF_cluster_count
            
            if (distinct_translated_ORF_cluster_count  == 1) or (mapped_orfs == 1):
                df['Peptide_translated_ORF_cluster_specific']='+'
            else:
                df['Peptide_translated_ORF_cluster_specific']='-'
            
            return df

class peptides2proteome:
    '''Basic class for peptide to proteome mapping - list of BioPython seqrecords for genome contigs'''
    def __init__(self, proteome, peptides_list, threads=1):
        self.threads=threads
        self.peptides_list=peptides_list  # list of peptides
        self.proteome=proteome            # Bio.SeqIO list of fasta records
        self.pepdict = self.process_peptides()           # map peptides
    
    def process_peptides(self):       
        peptides_list =list(set(self.peptides_list))
        pool = multiprocessing.Pool(self.threads)
        temps = pool.map(self.peptide_search, peptides_list)
        new_dict = {}
        for t in temps:
            for key in t:
                new_dict[key] = t[key]
        return new_dict
        
    def peptide_search(self, peptide):
        print('Searching {} in fasta'.format(peptide))
        temp =defaultdict(list)
        for rec in self.proteome:
            if peptide in str(rec.seq):
                temp[peptide].append(rec.id)
        return temp

class mapping2peptides:
    '''Basic class to export classes of peptides from peptide2genome mapping table'''
    # consider that when calculating the unique tanslated sequences, * at the end should be excluded.  
    
    def __init__(self, mapping, translation_table):
        self.mapping=mapping
        self.translation_table = translation_table
        self.paralogous_peptides = self.paralogous()
        self.orf_mapping = self.orf2peptides()

    def non_specific(self):
        map = self.mapping 
        ns =map[map['Peptide_translated_ORF_cluster_specific'] == '-']['Peptide_sequence'].tolist()
        ns = set(ns)
        return ns
    
    def specific(self):
        map = self.mapping
        s =map[map['Peptide_translated_ORF_cluster_specific'] == '+']['Peptide_sequence'].tolist()
        s = set(s)
        return s

    def paralogous(self):
        map = self.mapping
        s =map[ (map['Peptide_translated_ORF_cluster_specific'] =='+') & (map['Peptide_genome_ORF_count'] > 1) ]['Peptide_sequence'].tolist()
        s = set(s)
        return s

    #def specific(self):
    #    map = self.mapping
    #    s =map[ (map['Peptide_genome_ORF_count'] == 1) ]['Peptide_sequence'].tolist()
    #    s = set(s)
    #    return s

    #def non_specific(self):
    #    map = self.mapping
    #    s =map[ (map['Peptide_genome_ORF_count'] > 1) ]['Peptide_sequence'].tolist()
    #    s = set(s)
    #    return s

    def check_for_met(self, codons):
        met = False
        for codon in codons.split(';'):
            t = translate(Seq(codon), table = self.translation_table, cds=False)
            if t == 'M':
                met = True
        return met

    def non_atg_m(self):
        map = self.mapping
        mstarts = map[map['Peptide_amino_acid_first'] =='M']
        mstarts['MetCodon'] = mstarts['Peptide_first_codon'].apply(self.check_for_met)
        nonatgm = mstarts[mstarts['MetCodon'] == False]
        return set(nonatgm['Peptide_sequence'].tolist())
    
    def orf2peptides(self):
        orf_mapping = defaultdict(list)
        #specific_map=map[map['Peptide_sequence'].isin(peptides)].copy()
        ids = self.mapping['ORF_id'].tolist()
        peps = self.mapping['Peptide_sequence'].tolist()
        for i in range(len(self.mapping)):
            _id  = ids[i]
            _pep = peps[i]
            orf_mapping[_id].append(_pep)
        return orf_mapping 
    
    def get_peptides(self, orfs):
        peps = []
        for orf in orfs:
            peps += self.orf_mapping[orf]
        return peps

    def return_fasta(self, peptides):
        map = self.mapping
        orfs_fasta=[]
        prots_fasta=[]
        identical_translated= False
        specific_map=map[map['Peptide_sequence'].isin(peptides)].copy()
        specific_map = specific_map.drop_duplicates(['ORF_id'])
        orf_id_list =specific_map['ORF_id'].tolist()
        orf_nucs_list =specific_map['ORF_sequence'].tolist()
        orf_trans_list =specific_map['ORF_translation'].tolist()
        
        #if len(set(orf_trans_list)) == 1:
        #    identical_translated=True
        
        for item in range(len(orf_id_list)):
            orf_id = orf_id_list[item]
            orf_nucs = orf_nucs_list[item]
            orf_trans = orf_trans_list[item]
            
            nucs=SeqRecord(id=orf_id, seq=Seq(orf_nucs))
            trans=SeqRecord(id=orf_id, seq=Seq(orf_trans))
            
            #nucs=SeqRecord(id='|'.join(orf_id.split('|')[1:]),seq=Seq(orf_nucs))
            #trans=SeqRecord(id='|'.join(orf_id.split('|')[1:]),seq=Seq(orf_trans))
            
            orfs_fasta.append(nucs)
            prots_fasta.append(trans)
        
        #if identical_translated == True: # In case the translated ORF sequence is identical
        #    prot_ids = []
        #    for rec in prots_fasta:
        #        prot_ids.append(rec.id)
        #    prots_fasta = [SeqRecord(id = ';'.join(prot_ids), seq=rec.seq)]
        return orfs_fasta, prots_fasta


class pairwise_blast:
    def __init__(self, query, target, temp_folder, max_evalue=0.0001, matrix='BLOSUM62', word_size=3 ):
        assert not isinstance(query, list)
        assert not isinstance(target, list)
        
        SeqIO.write(query,  temp_folder +"/query.fasta", "fasta")
        SeqIO.write(target, temp_folder +"/target.fasta", "fasta")
        
        output = NcbiblastpCommandline(query=temp_folder +"/query.fasta", subject=temp_folder +"/target.fasta", outfmt=5, evalue = max_evalue, matrix=matrix, word_size=word_size)()[0] 
        blast_result_record = NCBIXML.read(StringIO(output)) 
        
        results = [] 
        variants = []
        hsps = []
        positions = []
        self.query = query
        self.target = target

        self.querylength=len(str(query.seq))
        self.targetlength = len(str(target.seq))
        
        var_dict={}

        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < max_evalue:
                    hsps.append(hsp)
                    h = []
                    h.append('****Alignment****')
                    h.append('sequence:         ' + str(alignment.title))
                    h.append('length:           ' + str(alignment.length))
                    h.append('e value:          ' + str(hsp.expect))
                    h.append('align_length:     ' + str(hsp.align_length))
                    h.append('frame:            ' + str(hsp.frame))
                    h.append('gaps:             ' + str(hsp.gaps))
                    h.append('identities:       ' + str(hsp.identities))
                    h.append('num alignments:   ' + str(hsp.num_alignments))
                    h.append('positives:        ' + str(hsp.positives))
                    h.append('hsp query start:  ' + str(hsp.query_start))
                    h.append('hsp query end:    ' + str(hsp.query_end))
                    h.append('hsp query:        ' + str(hsp.query))
                    h.append('hsp match:        ' + str(hsp.match))
                    h.append('hsp subject:      ' + str(hsp.sbjct))
                    h.append('hsp align_length: ' + str(hsp.align_length))
                    h.append('hsp gaps:         ' + str(hsp.gaps))
                    h.append('hsp sbjct_start:  ' + str(hsp.sbjct_start))
                    h.append('hsp sbjct_end:    ' + str(hsp.sbjct_end))
                    h.append('hsp score:        ' + str(hsp.score))
                    h.append('hsp strand:       ' + str(hsp.strand))
                    
                    h = '\n'.join(h)
                    match = str(hsp.match)
                    start = int(hsp.sbjct_start)
                    
                    qstart = int(hsp.query_start)
                    tstart = int(hsp.sbjct_start)
                    
                    diff = tstart-qstart
                    
                    target_len = len(str(target.seq))

                    for i, letter in enumerate(str(query.seq)):
                        location = i + 1 + diff
                        matchpos= location-tstart
                        if (matchpos >= 0) and (matchpos < len(match)) and (location <= target_len):
                            if match[matchpos] == '+':
                                positions.append(location)
                                var='{}->{}'.format(str(target.seq)[i+diff], str(query.seq)[i])
                                var_dict[location]=var
                        elif (location >= 1) and (location <= len(str(target.seq))):
                            if str(target.seq)[i+diff] != str(query.seq)[i]:
                                positions.append(location)
                                var='{}->{}'.format(str(target.seq)[i+diff], str(query.seq)[i])
                                var_dict[location]=var
                    results.append(h)
                break

        os.remove(temp_folder +'/query.fasta')
        os.remove(temp_folder +'/target.fasta')
        
        self.results = '\n'.join(results)
        self.hsps = hsps
        self.differences = list(set(positions))
        self.variants = var_dict

    def feature_overlap(self, gff3):  # sequtils gff3 object (must be from UniProt! Target must be a UniProt fasta record with id sp|P9WQP5|P9WQP5_MYCTU etc..
        if len(self.differences) > 0:
            features = gff3.table[gff3.table['seqid'] == self.target.id.split('|')[1]]
            overlap  = features[(features['start'] < max(self.differences))  & (features['end'] > min(self.differences))]
            
            res = []
            for row in overlap.iterrows():
                fields = []
                for ind in row[1].index:
                    if not ind.startswith('_'):
                        datum = '{}: {}'.format(ind, row[1][ind])
                        fields.append(datum)
                fields = '; '.join(fields)
                res.append(fields)
            res = '\n\n'.join(res)
            return res
        else:
            return None

class frameshift_peptides:
    def __init__(self, recs, peptides, tempfolder):
        mapped_orfs = {}
        self.reclist= recs.copy()
        frameshift_peptides=[]
        icds_map=defaultdict(list)
        peptide_stt_map = defaultdict(list)
        peptide_end_map = defaultdict(list)
        frameshifts=defaultdict(list)
        #blast_res=defaultdict()
        
        frameshift_results=defaultdict(list)
        frameshift_report = []

        for rec in self.reclist[:]:
            self.reclist.remove(rec)
            
            for trec in self.reclist:
                lst = [rec, trec]
                assert len(lst) ==2 
                if len(icds_blast(lst, tempfolder)) > 0:
                    
                    icds_map[rec.id].append(trec.id)
                    icds_map[trec.id].append(rec.id)
                    
                    start_rec=None
                    end_rec=None

                    for peptide in peptides:
                        
                        pep_res={}

                        start_found=False
                        end_found=False

                        fs_start =None
                        
                        fake_pep = 'X'*len(peptide)
                        for protrec in lst:

                            protseq = str(protrec.seq).upper()
                            
                            pstart = peptide[:4]
                            pend = peptide[-4:]
                            
                            id = protrec.id.split('|')
                            
                            nucs_coords = False
                            pg = None

                            if len(id) == 3:
                                try: 
                                    strand = id[2].split('(')[1].split(')')[0]
                                    assert ((strand=='+') or (strand =='-'))
                                    coords = id[2].split(')')[1].split(':')
                                    start = int(coords[0])
                                    end   = int(coords[1])
                                    contig = id[0]
                                    nucs_coords = True
                                except:
                                    pass

                            if (pstart in protseq) & (pend not in protseq):
                                start_found=True
                                fake_pep = pstart + fake_pep[len(pstart):]
                                start_rec=protrec.id
                                for i, letter in enumerate(peptide[4:]):
                                    if pstart + letter in protseq:
                                        pstart = pstart + letter
                                        fake_pep =  pstart + fake_pep[len(pstart):] 
                                    else:
                                        break
                                
                                if nucs_coords == True:
                                    pg = pep2genome(pstart, protseq, strand, start, end, contig)
                                    pep_res['start_rec_nuc_start'] = pg.genomic_start
                                    pep_res['start_rec_nuc_end'] = pg.genomic_end
                                    pep_res['start_rec_contig'] = pg.contig
                                    pep_res['start_rec_strand'] = strand

                                pep_res['start_rec_prot_start'] = protseq.find(pstart) + 1 
                                pep_res['start_rec_prot_end'] = protseq.find(pstart) + len(pstart) 
                                assert protseq[pep_res['start_rec_prot_start']-1:pep_res['start_rec_prot_end']] == pstart
                                pep_res['start_rec'] = protrec.id
                                pep_res['start_pep'] = pstart

                            elif (pend in protseq) & (pstart not in protseq):
                                end_found=True
                                fake_pep = fake_pep[:-len(pend)] + pend 
                                end_rec = protrec.id
                                
                                for i, letter in enumerate(peptide[:-4][::-1]):
                                    if letter + pend in protseq:
                                        pend = letter + pend
                                        fake_pep = fake_pep[:-len(pend)] + pend
                                    else:
                                        break
                                if nucs_coords == True:
                                    pg = pep2genome(pend, protseq, strand, start, end, contig)
                                    pep_res['end_rec_nuc_start'] = pg.genomic_start
                                    pep_res['end_rec_nuc_end'] = pg.genomic_end
                                    pep_res['end_rec_contig'] = pg.contig
                                    pep_res['end_rec_strand'] = strand
                                pep_res['end_rec_prot_start'] = protseq.find(pend) + 1
                                pep_res['end_rec_prot_end'] = protseq.find(pend) + len(pend) 
                                pep_res['end_rec'] = protrec.id
                                pep_res['end_pep'] = pend
                                assert protseq[pep_res['end_rec_prot_start']-1:pep_res['end_rec_prot_end']] == pend
                        if (end_found&start_found) == True:
                            frameshift_peptides.append(peptide)
                            fs='{} {}'.format(start_rec, end_rec)
                            frameshifts[fs].append(peptide)
                            
                            pep_res['covered_sequence'] = fake_pep
                            pep_res['peptide_length'] = len(peptide)
                            pep_res['covered_length'] = len(''.join(fake_pep.split('X')))
                            pep_res['peptide'] = peptide
                            
                            results = []
                            results.append('***Frameshift validation peptide {}***'.format(peptide))
                            desc_start = '{} matched at aa postion {}-{} in {}'.format(pstart, pep_res['start_rec_prot_start'], pep_res['start_rec_prot_end'], pep_res['start_rec'])
                            if 'start_rec_strand' in pep_res:
                                desc_start = desc_start + ' (genomic coordinates: ({}){}:{} in contig {})'.format(pep_res['start_rec_strand'], pep_res['start_rec_nuc_start'], pep_res['start_rec_nuc_end'], pep_res['start_rec_contig'])
                            results.append(desc_start)

                            desc_end = '{} matched at aa postion {}-{} in {}'.format(pend, pep_res['end_rec_prot_start'], pep_res['end_rec_prot_end'], pep_res['end_rec'])
                            if 'end_rec_strand' in pep_res:
                                desc_end = desc_end + ' (genomic coordinates: ({}){}:{} in contig {})'.format(pep_res['end_rec_strand'], pep_res['end_rec_nuc_start'], pep_res['end_rec_nuc_end'], pep_res['end_rec_contig'])
                            results.append(desc_end)
                            results = '\n'.join(results)
                            frameshift_results[fs].append(pep_res)
                            frameshift_report.append(results)

        self.frameshift_results = frameshift_results
        self.frameshift_peptides = frameshift_peptides
        self.peptide_stt_map = json.loads(json.dumps(peptide_stt_map))
        self.peptide_end_map = json.loads(json.dumps(peptide_end_map))
        self.frameshifts=json.loads(json.dumps(frameshifts))
        self.frameshift_report = '\n'.join(frameshift_report)

class pep2genome:
    def __init__(self, peptide, protein, strand, start, end, contig):
        if peptide in protein:
            if strand == '+':
                pstart = protein.find(peptide) 
                nstart = (pstart * 3) + start 
                nend = (pstart + len(peptide)) * 3 + start -1
                
                self.genomic_start = nstart
                self.genomic_end = nend

            elif strand == '-':
                pstart = protein.find(peptide) 
                
                nstart = end - (pstart * 3) 

                nend = end - (pstart + len(peptide)) * 3  + 1 
                
                self.genomic_start = nend
                self.genomic_end = nstart

        self.contig = contig
        self.peptide =peptide

class proteogenomics:
    
    def __init__(self, peptides, orf, reference, n_term_acetylated=[], table=11, start_codons=['ATG','GTG','TTG'], cleavage_aminos=['R','K']): #SeqRecord xx|xx|(+)123:1231 orf, ref entry, as BioPython SeqRecord
        self.peptides = peptides

        for pep in n_term_acetylated:
            assert pep in peptides
        
        self.table=table
        self.orf = orf
        self.orf_sequence = str(orf.seq)
        self.translated_orf_sequence = str(translate(orf.seq, table=self.table, cds=False))
        self.orf_id = self.orf.id
        self.translated_orf_rec = SeqRecord(id = self.orf_id, seq = Seq(self.translated_orf_sequence))
        
        #print(self.translated_orf_rec.format('fasta'))

        self.orf_contig = self.orf_id.split('|')[0]
        self.orf_number = self.orf_id.split('|')[1]
        self.orf_coords = self.orf_id.split('|')[2]

        self.reference = reference
        self.reference_sequence = ''.join(str(reference.seq).split('*'))
        
        
        temp_dir = tempfile.mkdtemp()
        
        self.pairwise_blast = pairwise_blast(self.translated_orf_rec, self.reference, temp_dir)
        
        os.removedirs(temp_dir)

        self.start_codons = start_codons
        
        self.cleavage_aminos = cleavage_aminos

        self.annotated_peptides = [pep for pep in peptides if pep in self.reference_sequence]
        
        self.novel_peptides = [pep for pep in peptides if not pep in self.annotated_peptides]

        self.acetylated = n_term_acetylated
        
        self.orf_strand, self.orf_start, self.orf_end = self.get_coords()

        pep_starts = []
        
        positions = defaultdict(set)        
        
        tss_validation_type = defaultdict(list)
        
        self.annotation_type = []
        
        self.tss_peptides = []
        
        self.validated_tss_sequences = []
        
        self.met_ap_peptides = []
        
        self.identified_tss_sites = set()
        
        self.variant_sequences = []
        
        self.other_peptidase_sites = []
        
        self.other_peptidase_type = defaultdict(list)
        
        self.other_peptidase = []
        
        self.other_peptidase_sequences = []
        
        self.upstream_inferred_tss_sequences = []
        
        self.variant_sequences_trie = []
        
        self.mapped_peptides = []
        
        self.unmapped_peptides = []

        for peptide in self.peptides:
            if peptide.startswith("M"):
                if peptide[1:] in self.translated_orf_sequence:
                    self.mapped_peptides.append(peptide)
                else:
                    self.unmapped_peptides.append(peptide)
            else:
                if peptide in self.translated_orf_sequence:
                    self.mapped_peptides.append(peptide)
                else:
                    self.unmapped_peptides.append(peptide)

        if len(self.pairwise_blast.hsps) == 0: # end the analysis if translated orf and ref dont align
            print('ORF and Reference not alligned')
            return
        
        if len(self.mapped_peptides) == 0: # end the analysis if translated orf and ref dont align
            print('No peptides for annotation')
            return
        
        orf_ref_pos = self.pairwise_blast.hsps[0].query_start - (self.pairwise_blast.hsps[0].sbjct_start - 1) -1 # 0 based position of annotated start site in translated orf sequence

        for peptide in self.mapped_peptides:
            enzymatic     = True
            met_init      = False   #  non-atg start codon translated with Met
            non_enzymatic = False   #  peptide not after anyzmatic cleavage site
            start_site    = False   #7  peptide starts at known start codon
            met_cleaved   = False   #  peptide evidence for initiator methionine cleavage
            nterm_acetyl  = False 

            if peptide.startswith('M'):
                pepstart = (self.translated_orf_sequence.find(peptide[1:]) - 1) # 0 based
            else:
                pepstart = self.translated_orf_sequence.find(peptide) 
            if pepstart == orf_ref_pos:
                relative = 'at annotated TSS site'
            elif pepstart > orf_ref_pos:
                relative = 'downstream of annotated TSS site'
            elif pepstart < orf_ref_pos:
                relative = 'upstream of annotated TSS site'
            first_codon = self.orf_sequence[pepstart * 3: pepstart * 3 + 3]
            first_amino = str(translate(Seq(first_codon),table=11, cds= False))
            previous_codon = self.orf_sequence[pepstart * 3-3: pepstart * 3 ]
            previous_amino = str(translate(Seq(previous_codon),table=11, cds=False))
            
            if previous_amino not in self.cleavage_aminos:
                enzymatic = False
            if peptide[0] =='M' and (first_codon in self.start_codons) and (first_codon != 'ATG'):
                met_init = True
            if first_codon in self.start_codons:
                start_site = True
            if (first_amino != 'M') and (previous_codon in self.start_codons):
                met_cleaved = True
            if peptide in self.acetylated:
                nterm_acetyl = True
            
            # Check if peptide at TSS

            if (nterm_acetyl == True) and (start_site == True) and (peptide[0] == 'M'):
                self.tss_peptides.append(peptide)
                self.identified_tss_sites.add(pepstart)
                tss_validation_type[pepstart].append('N-terminal acetylated start site peptide {} beginning with Methionine {}'.format(peptide, relative))

            elif (start_site == True) and (enzymatic == False) and (peptide[0] == 'M'): 
                self.tss_peptides.append(peptide)
                self.identified_tss_sites.add(pepstart)
                tss_validation_type[pepstart].append('Start site peptide {} beginning with Methionine with a non-enzymatic N-terminus {}'.format(peptide,relative))

            elif (met_init == True):
                self.tss_peptides.append(peptide)
                self.identified_tss_sites.add(pepstart )
                tss_validation_type[pepstart].append('Non-ATG start site ({}) peptide {} beginning with Methionine {}'.format(first_codon, peptide, relative))

            elif (met_cleaved == True):
                self.met_ap_peptides.append(peptide)
                self.tss_peptides.append(peptide)
                self.identified_tss_sites.add(pepstart - 1 )
                tss_validation_type[pepstart-1].append('Non-tryptic N-terminus peptide {} following an {} start-codon (initiator Methionine cleavage)'.format(peptide, previous_codon))

            elif (nterm_acetyl == True) or (enzymatic == False):
                other_peptidase_message = []
                if nterm_acetyl == True:
                    self.other_peptidase_type[pepstart].append('N-term acetylated peptide {} {}'.format(peptide, relative))
                if enzymatic == False:
                    self.other_peptidase_type[pepstart].append('Non-enzymatic N-terminus peptide {} {}'.format(peptide, relative))
                self.other_peptidase.append(peptide)
                self.other_peptidase_sites.append(pepstart)

            positions[pepstart].update(peptide[0])
            pep_starts.append(pepstart)
        
        most_upstream    = min(pep_starts)
        upstream_start = True
        if len(self.identified_tss_sites) > 0:
            if most_upstream > min(self.identified_tss_sites):
                upstream_start = False
            elif most_upstream == min(self.identified_tss_sites):
                _ = list(positions[most_upstream])
                non_M = [i for i in _ if i != 'M']
                if len(non_M) == 0:
                    upstream_start = False

        if most_upstream > orf_ref_pos:  # downstream of annotated start site
            upstream_start = False

        elif most_upstream == orf_ref_pos:
            _ = list(positions[most_upstream])
            non_M = [i for i in _ if i != 'M']
            if len(non_M) == 0:
                upstream_start = False

        variant_count = 1
        if upstream_start == True:
            current_pos = most_upstream 
            current_orf_start = current_pos * 3
            #current_codon = self.orf_sequence[current_pos * 3: current_pos * 3 + 3 ]
            current_codon = None 
            previous_codon = self.orf_sequence[(current_pos -1)* 3 : (current_pos-1) * 3 + 3]
            while (current_codon not in start_codons) and (previous_codon != ''):
                current_pos -= 1
                current_orf_start = current_pos * 3
                current_codon = self.orf_sequence[current_orf_start: current_orf_start + 3 ]
                previous_codon = self.orf_sequence[(current_pos -1)* 3 : (current_pos-1) * 3 + 3]
            if current_orf_start > 0:
                new_rec = self.get_var_rec(current_orf_start, variant_count, cds = True, variant_description="(Next upstream TSS of peptide identified upstream of mapped annotated sequence TSS)")
                self.annotation_type.append('Upstream non-TSS peptide with putative upstream TSS')
            else:
                new_rec = self.get_var_rec(current_orf_start, variant_count, cds = False, variant_description= "(No upstream ORF TSS of peptide identified upstream of mapped annotated sequence TSS)")
                self.annotation_type.append('Upstream non-TSS peptide with no putative upstream TSS')
            self.variant_sequences.append(new_rec)
            self.upstream_inferred_tss_sequences.append(new_rec)
            variant_count += 1
        
        for tss in self.identified_tss_sites:
            desc = "Identified TSS by {}.".format('; '.join(tss_validation_type[tss]))
            annotation_type = None
            if tss == orf_ref_pos:
                annotation_type = 'Annotated TSS validated'
            elif tss < orf_ref_pos:
                annotation_type = 'Upstream TSS identified'
            elif tss > orf_ref_pos:
                annotation_type = 'Downstream TSS identified'
            desc = desc + ' {}.'.format(annotation_type)
            new_rec = self.get_var_rec(tss * 3, variant_count, cds = True, variant_description=desc)
            self.variant_sequences.append(new_rec)
            self.validated_tss_sequences.append(new_rec)
            self.annotation_type.append(annotation_type) 
            variant_count += 1
        
        for other in self.other_peptidase_sites:
            new_rec = self.get_var_rec(other * 3, variant_count, cds = False, variant_description="(Other peptidase site identified by {})".format(' (AND) '.join(self.other_peptidase_type[other])))
            self.variant_sequences.append(new_rec)
            self.other_peptidase_sequences.append(new_rec)
            if other == orf_ref_pos:
                self.annotation_type.append('Annotated start site non-enzymative cleavage site')
            elif other < orf_ref_pos:
                self.annotation_type.append('Upstream non-enzymative cleavage site')
            elif other > orf_ref_pos:
                self.annotation_type.append('Downstream non-enzymative cleavage site')
            variant_count += 1
        self.variant_sequences_trie += list_trie_upper(self.variant_sequences, self.peptides)

    def get_var_rec(self, orf_base_position, variant_number, cds = True, variant_description='Sequence Isoform'):
        variant_orf = self.orf_sequence[orf_base_position:]
        try:
            translated_variant_orf = str(translate(Seq(self.orf_sequence[orf_base_position:]) , table = self.table, cds= cds))
        except:
            translated_variant_orf = str(translate(Seq(self.orf_sequence[orf_base_position:]) , table = self.table, cds= False))
        translated_variant_orf = ''.join(translated_variant_orf.split('*'))
        us_start, us_end = self.nuc_base_position_to_coords(orf_base_position)
        new_number = self.orf_number +'.{}'.format(variant_number)
        new_id = self.orf_contig + '|' + new_number + '|({}){}:{}'.format(self.orf_strand, us_start, us_end)
        variant_record = SeqRecord(seq = Seq(translated_variant_orf), id = new_id, description= variant_description)
        return variant_record

    def nuc_base_position_to_coords(self, orf_base_position):
        new_start = self.orf_start
        new_end = self.orf_end
        if self.orf_strand == '+':
            new_start = self.orf_start + orf_base_position
        elif self.orf_strand == '-':
            new_end = self.orf_end - orf_base_position
        return new_start, new_end
    
    def get_coords(self):
        coords =self.orf.id.split('|')[2]
        strand = coords.split(')')[0].split('(')[1]
        start = int(coords.split(')')[1].split(':')[0])
        end = int(coords.split(')')[1].split(':')[1])
        return strand, start, end
                    
def icds_blast(fasta, temp_folder, max_evalue=0.0001):
    
    tempfasta = fasta.copy()
    non_alligned = []
    blasted = []

    for qrec in tempfasta[:]:
        qid = qrec.id
        tempfasta.remove(qrec) # dont blast it against itself or again after the loop

        for trec in tempfasta:
            tid = trec.id
            datum = pairwise_blast(qrec, trec, temp_folder, max_evalue=max_evalue)
            if len(datum.hsps) == 0:
                non_alligned.append('{}, {} (evalue cutoff {})'.format(qid, tid, str(max_evalue)))
    return non_alligned

def reference_mapping_blast(fasta_orfs, fasta_refs, temp_folder, max_evalue=0.0001):
    
    non_alligned = []
    min_hsp_e = np.inf
    blast_pair = None

    for qrec in fasta_orfs:
        qid = qrec.id
        for trec in fasta_refs:
            tid = trec.id
    
            datum = pairwise_blast(qrec, trec,  temp_folder, max_evalue=max_evalue)
            for hsp in datum.hsps:
                if hsp.expect < min_hsp_e:
                    min_hsp_e = hsp.expect

                    blast_pair = '\n'.join(["Query id: {}".format(datum.query.id), "Target id: {}".format(datum.target.id), datum.results])

    return min_hsp_e, blast_pair

class annotate:
    def __init__(self, reference_protein, orf_sequence, specific_peptides, translation_table):
        self.reference_protein = reference_protein
        self.orf_sequence = orf
        self.peptides = specific_peptides
        self.reference_length = len(str(reference_protein.seq))
        self.translation_table = translation_table
        self.translated = translate(orf.seq, cds=False, table=self.translation_table)
    
    def annotation_type(self): 
        pass


def remove_asterisk(fasta):
    new_fasta = []
    for rec in fasta:
        rs = str(rec.seq)
        if rs[-1] == '*':
            rs = rs[:-1]
            rec.seq = Seq(rs)
        new_fasta.append(rec)
    return new_fasta



@contextlib.contextmanager
def stdout_redirect(where):
    sys.stdout = where
    try:
         yield where
    finally:
        sys.stdout = sys.__stdout__
        
def clustalw(output_file, fasta):
     new_fasta = fasta.copy()

     #for rec in new_fasta:
     #    new_ids = []
     #    ids = rec.id.split(';')
     ##    for id in ids:
     #        print(id)
     #        new_ids.append(id.split('|')[1])
     #    new_ids = ';'.join(new_ids)
     #    rec.id = new_ids

     SeqIO.write(new_fasta, output_file, 'fasta')
     cline = ClustalwCommandline("clustalw2", infile=output_file)
     stdout, stderr = cline()
     aln = output_file.split('.fasta')[0]

     try:
         align = AlignIO.read(aln +'.aln', "clustal")
         tree = Phylo.read(aln +".dnd", "newick")
         with stdout_redirect(StringIO()) as new_stdout:
             Phylo.draw_ascii(tree)
         new_stdout.seek(0)
         tree = new_stdout.read()
         return align, tree
     except:
         return None, None


def muscle(fasta):
    new_fasta = fasta.copy()
    #for rec in new_fasta:
    #    new_ids = []
    #    ids = rec.id.split(';')
    #    for id in ids:
    #        print(id)
    #        new_ids.append(id.split('|')[1])
    #    new_ids = ';'.join(new_ids)
    #    rec.id = new_ids    
    cline = MuscleCommandline(clwstrict=True)
    child = subprocess.Popen(str(cline),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True,
    shell=(sys.platform!="win32"))
    SeqIO.write(new_fasta, child.stdin, "fasta")
    child.stdin.close()
    lines=[]
    try:
        align = AlignIO.read(child.stdout, "fasta")
        return align
    except:
        return None

def list_trie_upper(fasta, peptide_list):
     new_fasta = fasta.copy()
     peptide_list = list(peptide_list)
     upper = []
     Trie = algo.Trie(peptide_list)
     
     for rec in new_fasta:
         id = rec.id
         temp_seq = str(rec.seq)
         description= rec.description
         TrieMatch=algo.TrieMatch(Trie, temp_seq)
         temp_seq = TrieMatch.trie_upper()
         assert(len(temp_seq)) == len(str(rec.seq))
         newrec = SeqRecord(seq = Seq(temp_seq), id = id, description = description)
         upper.append(newrec.format('fasta'))
     return upper




'''
def alignment_code(query, query_start, alignment, ref_prot):
    query_len = len(query[query_start:]) # query len is the length of the translated region
    query_translated = str(translate(Seq(query[query_start:])))
    assert ref_prot.endswith('*')#to account for stop codon in algorithm, ensure that the function recieves the translated sequence from BioPython
    ref_prot_len = len(ref_prot)*3 #nucleotide length
    if query == alignment:
        if query_len == ref_prot_len:
            code=5        # when whole sequence match and translation matches reference proteome (after first amino acid)
        else:
            code=1     #whole sequence match but translation does not excactly match reference proteome   
    elif query[-query_len:]==alignment[-query_len:]:
        if query_len == ref_prot_len:
            code=6 #translated sequence nucleotide sequence match to genomeand matches reference protein
        else:
            code=2    # translated sequence nucleotide sequence match  but not to reference protein
    elif query[-ref_prot_len:]==alignment[-ref_prot_len:]:
        code=3   # only regions matching are length of ref prot nuc sequence
    else:
        code=4 # mapps to a region coding for a known protein but not exact match
    return code
'''
'''
def sam_cigar(description):
    new_desc = description.split('_')
    ref = new_desc[0]
    start = int(new_desc[1])
    flag = new_desc[2]
    if flag == '0':
        strand = '+'
    elif flag == '16':
        strand ='-'
    cigar=new_desc[3]
    datum=(ref,strand,start,cigar)
    return datum
'''
'''
def sf_map(six_frame, genome,ref_prot, starts=['ATG','GTG','TTG'], table=11):
    stops=Stop_List(table)
    count=0
    sf_count=0
    y=set()
    for i in six_frame:
        query = str(i.seq)
        if query[-3:] in stops:
            x = 0
            found = False
            prot_query = ''
            while x<= len(query) and found == False:
                if query[x:x+3] in starts:
                    prot_query = translate(Seq(query[x:]),
                            table=11,
                            cds=False)         
                    found = True
                else:
                    x+=3
            seq_id = i.id
            seq_description = i.description.split()[1:] 
            for j in seq_description:
                prots_found = []
                datum = sam_cigar(j)
                strand=datum[1]
                start =int(datum[2])-1
                cigar = datum[3]
                end = get_end(cigar,start)
                if strand =='+':
                    aligned = str(genome[start:end].seq)
                elif strand =='-':
                    aligned = str(genome[start:end].seq.reverse_complement())
                template =  Seq(aligned).translate(table=11, cds=False)
                if aligned[-3:] in stops:
                    prot_found = False
                    for k in ref_prot:
                        if prot_found == False:
                            ref_record = k
                            prot = str(k.seq)+'*'
                            query_prot = prot[1:]
                            if query_prot in template:
                                prot_start = template.find(query_prot)
                                code = alignment_code(query, 
                                        x ,  
                                        aligned,
                                        prot)
                                prot_found = True
                                if strand == '+':
                                    orf_start = start+x
                                    orf_end = end
                                    prot_id = k.id
                                    prot_start=template.find(str(k.seq)[1:])*3-3    
                                    prot_header = '{}|{}:{}'.format(prot_id,
                                            start+ prot_start +1, 
                                            orf_end)
                                    count += 1
                                elif strand=='-':
                                        pass
                    if prot_found == False:
                        code = '*'
                        prot_header = None
                else:
                    pass
                if prot_found == False:
                    header = '({}){}:{}|{}|{}|*'.format(strand, start+1, end, cigar,code)
                else:
                    header = '({}){}:{}|{}|{}|{}'.format(strand, start+1, end, cigar,code,prot_header)
                if code == 4 and strand == '+':
                    print 
                    print header
                    print 'Query:\n', query
                    print 'Alignment:\n',aligned
                    print str(aligned)==str(query)
                    print len(query), len(aligned)
                    print [i for i in xrange(len(query)) if query[i] != aligned[i]]
                    #print query[494], aligned[494]

                    print '\nTranslated Query:\n', prot_query
                    print 
                    print ref_record.format('fasta')
                    print '########################################'
        elif query[-3:] not in stops:
             pass
        sf_count+=1 
      
'''


'''
def annotate_gm(gm, sf):#not used, but can be used to annotate a gm prediction
    annotated = []
    count = 0
    found = []
    for i in gm:
        lst = i.description.split('\t>')[0].split('|')
        done = False
        strand = lst[3]
        start = lst[4]
        end = lst[5]
        description = '({}){}:{}'.format(strand, start, end)
        annos = []
        count += 1
        while done == False:
            for j in sf:
                if done == False and j not in found:
                    desc = j.description.split()[1].split('|')
                    for coord in desc:
                        j_strand = j.description.split()[1].split(')')[0][1:]
                        j_start  = j.description.split()[1].split(')')[1].split(':')[0]
                        j_end = j.description.split()[1].split(')')[1].split(':')[1]
                        if strand == j_strand:
                            if strand == '-':
                                if start == j_start:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
                            elif strand == '+':
                                if end == j_end:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
            done = True
        annotations = description+' '+'|'.join(annos)
        record = SeqRecord(i.seq, id = i.id, description = annotations)
        annotated.append(record)
return annotated
'''

'''
def m_start_dict(proteins):
	m_dict = defaultdict()
	for prot in proteins:
		assert not isinstance(proteins[prot], (list, tuple))
		record =proteins[prot]
		if prot[:1] !='M':
			id = record.id
			seq = 'M' + str(record.seq)[1:]
			description = record.description
			s_record = SeqRecord(id = id, seq = Seq(seq), description = description)
			if seq not in proteins:
				m_dict[seq] = s_record   #discard elements already known to exist
		elif prot[:1]=='M':
			s_record = SeqRecord(id = record.id, seq = Seq(str(record.seq)), description = record.description)
			m_dict[prot] =s_record     #add known reference proteins starting with 'M'
        for i in m_dict:
		rec = m_dict[i]
		rec.id = rec.id.split('|')[0]+ '|' + 'M_START_'+ rec.id.split('|')[1]
	return m_dict
'''
'''
def alt_ref(ref, seq):     # refs annotated by cds location, as dict, and a SeqRecord object (fasta, embl etc). Returns alternative start amino acid
    fasta = seq.format('fasta')
    genome = SeqIO.read(StringIO.StringIO(fasta),'fasta')
    M_refs = defaultdict(list)
    M_refs_id = defaultdict(list)
    alt_ref = []
    for i in ref:
        for j in ref[i]:
            if len(j.description.split(' Coordinates: ')) > 1:
                coords =  ref[i][0].description.split(' Coordinates: ')[1]
                for coord in coords.split('|'):
                    try:
                        start = int(coord.split(')')[1].split(':')[0])-1
                        end = int(coord.split(')')[1].split(':')[1])
                        strand = coord.split(')')[0][1:]
                        nucs = genome.seq[start:end]

                        if strand == '+':
                            prot =  str(translate(nucs, table = 11, cds = False))
                            assert prot[-1] == '*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)

                        elif strand == '-':
                            prot = str(translate(nucs.reverse_complement(), table= 11, cds = False))
                            assert prot[-1] =='*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)
                    except:
                        pass
    for i in M_refs_id:
        for j in M_refs_id[i]:
            j.id = j.id +' Alternative Start'
            j.description = j.description.split(' Coordinates: ')[0] +' Non-ATG start at coordinates: '+'|'.join(M_refs[i])
            j.seq = Seq(i)
            alt_ref.append(j)
            break
    return alt_ref
'''

'''
#5. GET NOVEL SEQUENCES
def get_novel(proteins):    # input a list of annotated and non annotated sequences, and output only the non-annotated ones
    novel = []
    for i in proteins:
        if len(i.description.split()) == 1:
            novel.append(i)
    return novel

def ldict(seqs):
	seqs_dict = defaultdict(list)
	for i in seqs:
		seqs_dict[str(i.seq)]
'''
'''
def annotate_gm(gm, sf):#not used, but can be used to annotate a gm prediction
    annotated = []
    count = 0
    found = []
    for i in gm:
        lst = i.description.split('\t>')[0].split('|')
        done = False
        strand = lst[3]
        start = lst[4]
        end = lst[5]
        description = '({}){}:{}'.format(strand, start, end)
        annos = []
        count += 1
        while done == False:
            for j in sf:
                if done == False and j not in found:
                    desc = j.description.split()[1].split('|')
                    for coord in desc:
                        j_strand = j.description.split()[1].split(')')[0][1:]
                        j_start  = j.description.split()[1].split(')')[1].split(':')[0]
                        j_end = j.description.split()[1].split(')')[1].split(':')[1]
                        if strand == j_strand:
                            if strand == '-':
                                if start == j_start:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
                            elif strand == '+':
                                if end == j_end:
                                    done = True
                                    anno = j.id
                                    annos.append(anno)
                                    found.append(j)
            done = True
        annotations = description+' '+'|'.join(annos)
        record = SeqRecord(i.seq, id = i.id, description = annotations)
        annotated.append(record)
return annotated
'''

'''
def m_start_dict(proteins):
	m_dict = defaultdict()
	for prot in proteins:
		assert not isinstance(proteins[prot], (list, tuple))
		record =proteins[prot]
		if prot[:1] !='M':
			id = record.id
			seq = 'M' + str(record.seq)[1:]
			description = record.description
			s_record = SeqRecord(id = id, seq = Seq(seq), description = description)
			if seq not in proteins:
				m_dict[seq] = s_record   #discard elements already known to exist
		elif prot[:1]=='M':
			s_record = SeqRecord(id = record.id, seq = Seq(str(record.seq)), description = record.description)
			m_dict[prot] =s_record     #add known reference proteins starting with 'M'
        for i in m_dict:
		rec = m_dict[i]
		rec.id = rec.id.split('|')[0]+ '|' + 'M_START_'+ rec.id.split('|')[1]
	return m_dict
'''
'''
def alt_ref(ref, seq):     # refs annotated by cds location, as dict, and a SeqRecord object (fasta, embl etc). Returns alternative start amino acid
    fasta = seq.format('fasta')
    genome = SeqIO.read(StringIO.StringIO(fasta),'fasta')
    M_refs = defaultdict(list)
    M_refs_id = defaultdict(list)
    alt_ref = []
    for i in ref:
        for j in ref[i]:
            if len(j.description.split(' Coordinates: ')) > 1:
                coords =  ref[i][0].description.split(' Coordinates: ')[1]
                for coord in coords.split('|'):
                    try:
                        start = int(coord.split(')')[1].split(':')[0])-1
                        end = int(coord.split(')')[1].split(':')[1])
                        strand = coord.split(')')[0][1:]
                        nucs = genome.seq[start:end]

                        if strand == '+':
                            prot =  str(translate(nucs, table = 11, cds = False))
                            assert prot[-1] == '*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)

                        elif strand == '-':
                            prot = str(translate(nucs.reverse_complement(), table= 11, cds = False))
                            assert prot[-1] =='*'
                            prot = prot[:-1]
                            if prot[:1] != 'M':
                                M_refs[prot].append(coord)
                                M_refs_id[prot].append(j)
                    except:
                        pass
    for i in M_refs_id:
        for j in M_refs_id[i]:
            j.id = j.id +' Alternative Start'
            j.description = j.description.split(' Coordinates: ')[0] +' Non-ATG start at coordinates: '+'|'.join(M_refs[i])
            j.seq = Seq(i)
            alt_ref.append(j)
            break
    return alt_ref
'''

'''
#5. GET NOVEL SEQUENCES
def get_novel(proteins):    # input a list of annotated and non annotated sequences, and output only the non-annotated ones
    novel = []
    for i in proteins:
        if len(i.description.split()) == 1:
            novel.append(i)
    return novel

def ldict(seqs):
	seqs_dict = defaultdict(list)
	for i in seqs:
		seqs_dict[str(i.seq)]
'''
def proteome2string(fasta):
    proteome_strings = []
    for rec in fasta:
        proteome_strings.append(str(rec.seq))
    combined_string = 'X'.join(proteome_strings)
    return combined_string

def genome2aa(fasta):
    genome_list = []
    for contig in fasta:
        g = str(contig.seq)
        frame1 = str(translate(Seq(g[:]),  cds = False, table = 11))
        frame2 = str(translate(Seq(g[1:]), cds = False, table = 11))
        frame3 = str(translate(Seq(g[2:]), cds = False, table = 11))
        frame4 = str(translate(Seq(g[:]).reverse_complement(),  cds = False, table = 11))
        frame5 = str(translate(Seq(g[:-1]).reverse_complement(),cds = False, table = 11))
        frame6 = str(translate(Seq(g[:-2]).reverse_complement(),cds = False, table = 11))
        frame_list = [frame1, frame2, frame3, frame4, frame5, frame6]
        genome_list.append('X'.join(frame_list))
    genome_str = 'X'.join(genome_list)
    return genome_str


