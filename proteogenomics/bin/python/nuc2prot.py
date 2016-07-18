#!/usr/bin/env python

import sequtils
import sys
import Bio; from Bio import SeqIO; from Bio.Seq import translate
from Bio.Seq import translate, Seq
from Bio.SeqRecord import SeqRecord
import time

genome = None
proteome = None

if len(sys.argv) == 3:
    genome = sys.argv[1]
    proteome = sys.argv[2]
    genome = SeqIO.read(genome, 'fasta')
    proteome  = list(SeqIO.parse(proteome, 'fasta'))


starts = ['ATG','GTG','TTG']
stops  = ['TAG','TGA', 'TAA']

prots = []

fasta =  list(SeqIO.parse(sys.stdin,'fasta'))


for i in fasta:
    
    count = 0
    maps = ''

    start = int(i.description.split()[1].split('_')[-1].split(')')[1].split(':')[0])
    
    end = int(i.description.split()[1].split('_')[-1].split(')')[1].split(':')[1])
    strand  = i.description.split()[1].split('_')[-1][1]

    cigar = i.description.split()[1].split('_')[-3]
    
    additions = []

    full_match = False
    
    letters = []
    for letter in cigar:
        if letter in 'DSMI':
            letters.append(letter)
    
    if proteome != None:
        if strand == '+':
            if (str(i.seq)[-3:] not in stops) and (letters[-1:][0] == 'M'):
                while (str(genome.seq[end-3:end]) not in stops) and (end <= len(str(genome.seq)) - 3):
                    end += 3
                    add = str(genome.seq[end-3:end])
                    assert len(add) == 3                
                    additions.append(add)
            
            if letters[-1:][0] == 'M':
                full_match = True
                start = end-6 
                while (str(genome.seq[start-3:start]) not in stops) and (start >= 3):
                    start -= 3
                nucs = genome.seq[start:end]
                assert len(str(nucs)) % 3 == 0
        
        elif strand == '-': 
            if (str(i.seq)[-3:] not in stops) and (letters[:1][0] == 'M'):
                
                while (str(genome.seq[start-1:end].reverse_complement())[-3:] not in stops) and (start >= 3):
                    start -= 3
                    add =str(genome.seq[start-1:end].reverse_complement())[-3:] 
                    assert len(add) == 3
                    additions.append(add)
            
            
            if letters[:1][0] == 'M':
                full_match = True
                end = start + 5
                while (str(genome.seq[end:end+3].reverse_complement()) not in stops) and (end <= len(str(genome.seq))-3):
                    end += 3
                nucs = genome.seq[start-1:end].reverse_complement()
                start -= 1

        if full_match == True:
            nuc_trans = str(translate(nucs, table=11, cds=False))
            maps = []
            for j in proteome:     
                protein = j.seq
                pr_id = j.id
                if str(protein)[1:] in nuc_trans:
                    maps.append(pr_id)
            maps = ';'.join(maps)
    
    if len(maps) == 0:
        maps = '*'
    
    seq = str(i.seq)
    while (seq[count:count + 3] not in starts) and (count < len(seq)):
        count += 3

    extend = len(additions) * 3
    maps = "({}){}:{}_{}_3PrimeExtended:{}".format(strand, start + 1, end, maps, extend)
    
    seq = seq[count:] + ''.join(additions)
    

    seq = str(translate(Seq(seq), cds =False, table = 11)) 
    if seq.endswith('*'):
        seq = seq[:-1]
   
    if len(seq) >= 20:
        rec = SeqRecord(seq = Seq(seq), description = i.description+'_'+maps, id = i.id)
        prots.append(rec)
    
SeqIO.write(prots, sys.stdout, 'fasta')

