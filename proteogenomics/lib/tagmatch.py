#!/usr/bin/env python

import pandas as pd
import collections
from collections import defaultdict
import re
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np

def peptide_mass(peptide, fixed_modifications=["Carbamidomethylation of C"], variable_modifications=[], nterm=True, cterm=True):
    if 'X' in peptide:
        return [np.inf]
    if 'Z' in peptide:
        return [np.inf]
    if 'B' in peptide:
        return [np.inf]

    
    masses = []

    mass_h20=18.01056
    nterm_mass =    1.007825035
    cterm_mass =    1.007825035 + 15.99491463


    monoisotopic = {'A':   71.037114,
                    'C':   103.009185,
                    'D':   115.026943,
                    'E':   129.042593,
                    'F':   147.068414,
                    'G':   57.021464,
                    'H':   137.058912,
                    'U':   150.95363,
                    'O':   237.147735,
                    'I':   113.084064,
                    'K':   128.09496,
                    'L':   113.084064, 
                    'J':   113.084064,
                    'M':   131.04049,
                    'N':   114.042927,
                    'P':   97.05276,
                    'Q':   128.05858,
                    'R':   156.101111,
                    'S':   87.03203,
                    'T':   101.04768,
                    'V':   99.06841,
                    'W':   186.07931,
                    'Y':   163.06333  }
    ##############
    # Fixed mods #
    ##############
    options = ["Carbamidomethylation of C"]
    for i in fixed_modifications:
        assert i in options
    if "Carbamidomethylation of C" in fixed_modifications:
        monoisotopic['C'] = monoisotopic['C'] + 57.021464
    
    for aa in peptide:
        masses.append(monoisotopic[aa])

    # < insert directed graph here > #
    masslist = [sum(masses)]
    # < insert directed graph here > #
    
    if nterm == True:
        masslist = [i + nterm_mass for i in masslist]
    if cterm == True:
        masslist = [i + cterm_mass for i in masslist]

    return masslist

def mz2mw(mz, charge):
    proton_mass=1.007276
    mw = (mz * charge) - (proton_mass * charge) 
    return mw


class TagMatch:
    def __init__(self, query, tag_mass_list, target, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], cleavage_rule='tryptic', specificity='specific', prec_tol=0.02, gap_tol=0.5):
        self.target = target
        self.target_length = len(target)

        self.tag_mass_list = tag_mass_list

        self.query=query
        self.query_length = len(query)
        self.fixed_modifications=fixed_modifications
        self.variable_modifications=variable_modifications
        self.query_positions = self.query_positions()
        self.prec_tol = prec_tol
        self.gap_tol  = gap_tol

        self.validated=False
        self.specificity = specificity

        if cleavage_rule=='tryptic':
            self.tryptic_positions = self.tryptic_match()
            self.tryptic_peptides  = self.tryptic_peptides()
            self.tryptic_masses = self.tryptic_masses()
            
            self.amino_gaps = self.amino_gaps()
            self.validate_masses(self.tryptic_masses)
            
    def query_positions(self):
        positions = [(m.start(), m.start()+self.query_length ) for m in re.finditer('(?={})'.format(self.query), self.target)]
        return positions
    
    def validate_masses(self, mass_dict):
        validated_precursor = []
        precursor_peptides = defaultdict(list)
        validated_peptides = set()
        
        for position_key in mass_dict:
            target_masses = mass_dict[position_key]
            target_peptide = self.tryptic_peptides[position_key]
            
            for target_mass in target_masses:
                for tag_mass_set in self.tag_mass_list:
                    assert len(tag_mass_set) == 3
                    query_mass = tag_mass_set[1]
                    diff = abs(target_mass - query_mass)
                    if diff < self.prec_tol:
                        nterm_mass_gap = tag_mass_set[0]
                        cterm_mass_gap = tag_mass_set[2]
                        for amino_gaps in self.amino_gaps[target_peptide]:
                            nterm_gap_match=False
                            cterm_gap_match=False

                            nterm_amino_gaps = peptide_mass(amino_gaps[0], fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=False, cterm=False)
                            for nterm_amino_gap in nterm_amino_gaps:
                                diff = abs(nterm_mass_gap - nterm_amino_gap)
                                if diff < self.gap_tol:
                                    nterm_gap_match=True

                            cterm_amino_gaps = peptide_mass(amino_gaps[1], fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=False, cterm=False)
                            for cterm_amino_gap in cterm_amino_gaps:
                                diff = abs(cterm_mass_gap - cterm_amino_gap)
                                if diff < self.gap_tol:
                                    cterm_gap_match=True
                            
                            if (cterm_gap_match==True) and (nterm_gap_match==True):
                                validated_precursor.append(tag_mass_set)
                                precursor_peptides[tag_mass_set].append(self.tryptic_peptides[position_key])
                                validated_peptides.add(self.tryptic_peptides[position_key])
        self.validated_precursor = validated_precursor
        self.precursor_peptides  = precursor_peptides
        self.validated_peptides = list(validated_peptides)
        if len(validated_precursor) > 0:
            self.validated = True

    def tryptic_masses(self):
        tryptic_masses = {}
        for key in self.tryptic_peptides:
            peptide = self.tryptic_peptides[key]
            masslist = peptide_mass(peptide, fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=True, cterm=True)
            tryptic_masses[key] = masslist
        return tryptic_masses

    def tryptic_peptides(self):
        tryptic_peptides = {}
        for coord in self.tryptic_positions:
            p = self.target[coord[0]:coord[1]]
            
            key='{}:{}'.format(str(coord[0]), str(coord[1]))
            tryptic_peptides[key] = p
        return tryptic_peptides

    def tryptic_match(self):
        tryptic_coords = []
        for coord in self.query_positions:
            nterms = self.tryptic_nterm(coord[0])
            cterms = self.tryptic_cterm(coord[1]-1)
            for nterm in nterms:
                for cterm in cterms:
                    tryptic_coords.append((nterm, cterm))
        return tryptic_coords


    def tryptic_nterm(self, pos):
        nterm = pos
        cleavage_aminos=['R','K']
        nterms=[]
        while True:
            if nterm == 0:
                nterms.append(nterm)
                return nterms
            elif self.target[nterm-1] in cleavage_aminos:
                nterms.append(nterm)
                return nterms
            else: 
                if self.specificity == 'semi-specific':
                    nterms.append(nterm)
                nterm -= 1

    def tryptic_cterm(self, pos):
        cterm = pos
        cleavage_aminos=['R','K']
        cterms=[]
        while True:
            if cterm == self.target_length-1:
                cterms.append(cterm + 1)
                return cterms
            elif self.target[cterm] in cleavage_aminos:
                cterms.append(cterm + 1)
                return cterms
            else: 
                if self.specificity == 'semi-specific':
                    cterms.append(cterm + 1)
                cterm += 1
    
    def amino_gaps(self):
        amino_gap_dict = defaultdict(list)
        for key in self.tryptic_peptides:
            peptide = self.tryptic_peptides[key]
            positions = [(m.start(), m.start()+self.query_length ) for m in re.finditer('(?={})'.format(self.query), peptide)]
            for coord in positions:
                tag_start = coord[0]
                tag_end = tag_start + self.query_length

                nterm_gap = peptide[:tag_start]
                cterm_gap = peptide[tag_end:]
                
                amino_gap_dict[peptide].append( (nterm_gap, cterm_gap) )

        return amino_gap_dict


def character_strip(sequence, character):
    started=False
    start_pos=0
    while (started==False) and (start_pos < len(sequence)):
        if sequence[start_pos]==character:
            start_pos += 1
        else:
            started = True

    end_pos=len(sequence)-1
    ended=False
    while (ended==False) and (end_pos > 0):
        if sequence[end_pos]==character:
            end_pos -= 1
        else:
            ended = True
    return sequence[start_pos: end_pos + 1]

class gap_sequence:
    def __init__(self, gap, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02, gap_tol=0.5):
        self.gap = gap
        self.prec_tol = prec_tol
        self.gap_tol  = gap_tol
        self.fixed_modifications = fixed_modifications
        self.variable_modifications = variable_modifications
        self.gap_sequences=[]
        
        self.possible = set([0])
        self.valid = False

        self.aminos     =           ['A', 
                                     'C', 
                                     'D', 
                                     'E', 
                                     'F' ,
                                     'G', 
                                     'H', 
                                     'U', 
                                     'O', 
                                     'I', 
                                     'K', 
                                     'L', 
                                     'J',
                                     'M', 
                                     'N', 
                                     'P' , 
                                     'Q', 
                                     'R', 
                                     'S', 
                                     'T' , 
                                     'V' , 
                                     'W', 
                                     'Y']

    def sequences(self, mass=0, sequence=''):
        diff = np.absolute(mass - self.gap)

        if diff < self.prec_tol:
            self.gap_sequences.append(sequence)
        
        elif (mass < self.gap):
            for amino in self.aminos:
                masses=peptide_mass(amino, variable_modifications = self.variable_modifications, fixed_modifications=self.fixed_modifications, nterm=False, cterm=False)
                for m in masses:
                    newmass=mass + m
                    newsequence = sequence + amino
                    self.sequences(newmass, newsequence)
    
    def validate(self):
        newpossible=set()
        for mass in self.possible:
            diff = np.absolute(mass - self.gap)

            if diff < self.prec_tol:
                self.valid = True
        
            elif (mass < self.gap) and (self.valid==False):
                for amino in self.aminos:
                    masses=peptide_mass(amino, variable_modifications = self.variable_modifications, fixed_modifications=self.fixed_modifications, nterm=False, cterm=False)
                    for m in masses:
                        newmass = mass + m
                        newpossible.add(np.round(newmass,6))
        if (len(newpossible) > 0) and (self.valid ==False):
            self.possible=newpossible
            self.validate()

class blast_tags:
    def __init__(self, subject, tags, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02):

        self.subject=subject
        self.tags=tags 
        self.fixed_modifications = fixed_modifications
        self.variable_modifications = variable_modifications
        self.prec_tol  = prec_tol

        self.mweights = set()
        self.nmasses  = set()
        self.cmasses  = set()
        
        self.samples  = set()
        self.scans    = set()

        for qtag in self.tags:
            qseq = str(qtag.seq)
            qsubst = self.longest_substring(qseq, self.subject)
            qnmass_offsets=self.nmass_gap(self.subject,self.substring_coordinates(self.subject, qsubst))
            qid = qtag.id
            self.samples.add(qid.split(';')[0])
            self.scans.add(qid.split(';')[1])


            qid_mw = float(qid.split(';')[2].split('mw=')[1])
            self.mweights.add(qid_mw)

            qid_ng = float(qid.split(';')[3].split('ngap=')[1])
            qid_cg = float(qid.split(';')[4].split('cgap=')[1])
            
            # Calculate mass gaps
            qcoords = self.substring_coordinates(qseq, qsubst)
            qnmassgaps = self.nmass_gap(qseq, qcoords)
            for qngap in qnmassgaps:
                for offset in qnmass_offsets:
                    self.nmasses.add(qid_ng + (qngap - offset))
            
            qcmassgaps = self.cmass_gap(qseq, qcoords)
            qcmass_offsets=self.cmass_gap(self.subject,self.substring_coordinates(self.subject, qsubst)) 
            for qcgap in qcmassgaps:
                for offset in qcmass_offsets:               
                    self.cmasses.add( qid_cg +  (qcgap - offset) )
            
            #for ttag in self.tags:
            #    tseq = str(ttag.seq)
            #    tsubst = self.longest_substring(tseq, self.subject)
            #    tsubst = self.longest_substring(tseq, self.subject)
            #    tid = ttag.id
            #    
                
            #    tid_mw = float(tid.split(';')[2].split('mw=')[1])
            #    tid_ng = float(tid.split(';')[3].split('nmg=')[1])
            #    tid_cg = float(tid.split(';')[4].split('cmg=')[1])
            #    
            #    # Calculate mass gaps
            #    tcoords = self.substring_coordinates(tseq, tsubst)
        
        self.cmasses = set([np.round(i,6) for i in self.cmasses if i >= 0])
        self.nmasses = set([np.round(i,6) for i in self.nmasses if i >= 0])

        self.newtags = set()
        
        assert len(self.samples) ==1
        assert len(self.scans) == 1

        for mw in self.mweights:
            for nmass in self.nmasses:
                for cmass in self.cmasses:
                    newtag = (self.subject, mw, nmass, cmass )
                    if self.valid_tag(newtag) == True:
                        self.newtags.add(newtag)
        
        self.newrecords = self.new_records()

    def valid_tag(self, tag):
        sequence_masses = peptide_mass(tag[0], fixed_modifications=self.fixed_modifications, variable_modifications = self.variable_modifications, cterm=True, nterm=True)
        for mass in sequence_masses:
            total = mass + tag[2] + tag[3]
            diff = np.absolute(total - tag[1])
            if diff < self.prec_tol:
                return True

    def longest_substring(self, reference, match_sequence):
        substrings=[]
        for charpos in range(len(match_sequence)):
            mismatch=False
            pos= charpos
            while (mismatch == False) & (pos < len(match_sequence)):
                if match_sequence[charpos:pos +1]  in reference:
                    pos += 1
                else:
                    mismatch=True
                subs = match_sequence[charpos:pos ]
            substrings.append(subs)
        return max(substrings, key=len)
    
    def substring_coordinates(self, reference, substr):
        coordinates=[(m.start(), m.start() + len(substr)) for m in re.finditer('(?={})'.format(substr), reference)]
        return coordinates

    def cmass_gap(self, reference, coordinates):
        mass_gaps=[]
        for coordinate in coordinates:
            amino_gap = reference[coordinate[1]:]
            amino_masses = peptide_mass(amino_gap, fixed_modifications = self.fixed_modifications, variable_modifications = self.variable_modifications, nterm=False, cterm=False)
            mass_gaps += amino_masses
        return mass_gaps

    def nmass_gap(self, reference, coordinates):
        mass_gaps=[]
        for coordinate in coordinates:
            amino_gap = reference[:coordinate[0]]
            amino_masses = peptide_mass(amino_gap, fixed_modifications = self.fixed_modifications, variable_modifications = self.variable_modifications, nterm=False, cterm=False)
            mass_gaps += amino_masses
        return mass_gaps
    
    def new_records(self):
        newrecords = []
        for newtag in self.newtags:
            rec = SeqRecord(seq=Seq(newtag[0]), id = '{};{};mw={};ngap={};cgap={}'.format(list(self.samples)[0], list(self.scans)[0], str(newtag[1]), str(newtag[2]), str(newtag[3])))
            newrecords.append(rec)
        return newrecords
