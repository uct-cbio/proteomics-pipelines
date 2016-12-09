#!/usr/bin/env python

import pandas as pd
import collections
from collections import defaultdict
import re

def peptide_mass(peptide, fixed_modifications=["Carbamidomethylation of C"], variable_modifications=[], nterm=True, cterm=True):
    if 'X' in peptide:
        return 0
    if 'Z' in peptide:
        return 0
    if 'B' in peptide:
        return 0

    
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
    def __init__(self, query, tag_mass_list, target, fixed_modifications=[], variable_modifications=[], cleavage_rule='Tryptic', prec_tol=0.02):
        self.target = target
        self.target_length = len(target)

        self.tag_mass_list = tag_mass_list

        self.query=query
        self.query_length = len(query)
        self.fixed_modifications=fixed_modifications
        self.variable_modifications=variable_modifications
        self.query_positions = self.query_positions()
        self.tolerance = prec_tol
        self.validated=False
        
        if cleavage_rule=='Tryptic':
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
                    if diff < self.tolerance:
                        nterm_mass_gap = tag_mass_set[0]
                        cterm_mass_gap = tag_mass_set[2]
                        for amino_gaps in self.amino_gaps[target_peptide]:
                            nterm_gap_match=False
                            cterm_gap_match=False

                            nterm_amino_gaps = peptide_mass(amino_gaps[0], fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=True, cterm=False)
                            for nterm_amino_gap in nterm_amino_gaps:
                                if nterm_amino_gap > (nterm_mass_gap - self.tolerance):
                                    nterm_gap_match=True

                            cterm_amino_gaps = peptide_mass(amino_gaps[1], fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=False, cterm=True)
                            for cterm_amino_gap in cterm_amino_gaps:
                                if cterm_amino_gap > (cterm_mass_gap - self.tolerance):
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
            nterm = self.tryptic_nterm(coord[0])
            cterm = self.tryptic_cterm(coord[1]-1)
            tryptic_coords.append((nterm, cterm))
        return tryptic_coords


    def tryptic_nterm(self, pos):
        nterm = pos
        cleavage_aminos=['R','K']
        while True:
            if nterm == 0:
                return nterm
            elif self.target[nterm-1] in cleavage_aminos:
                return nterm
            else: 
                nterm -= 1

    def tryptic_cterm(self, pos):
        cterm = pos
        cleavage_aminos=['R','K']
        while True:
            if cterm == self.target_length-1:
                return cterm + 1
            elif self.target[cterm] in cleavage_aminos:
                return cterm + 1
            else: 
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

