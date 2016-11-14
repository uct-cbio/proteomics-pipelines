#!/usr/bin/env python

import pandas as pd
import collections
from collections import defaultdict
import re

def peptide_mass(peptide):
    masses = []
    mass_h20=18.01056

    monoisotopic = {'A':   71.03711,
                    'C':   103.00919,
                    'D':   115.02694,
                    'E':   129.04259,
                    'F':   147.06841,
                    'G':   57.02146,
                    'H':   137.05891,
                    'U':   150.95363,
                    'O':   255.158295,
                    'I':   113.08406,
                    'K':   128.09496,
                    'L':   113.08406,
                    'M':   131.04049,
                    'N':   114.04293,
                    'P':   97.05276,
                    'Q':   128.05858,
                    'R':   156.10111,
                    'S':   87.03203,
                    'T':   101.04768,
                    'V':   99.06841,
                    'W':   186.07931,
                    'Y':   163.06333  }
    for aa in peptide:
        masses.append(monoisotopic[aa])
    mass = sum(masses) + mass_h20
    return mass

def mz2mw(mz, charge):
    proton_mass=1.007276
    mw = (mz * charge) - (proton_mass * charge) 
    return mw


class TagMatch:
    def __init__(self, query,  precursor_mass_list, target, modifications=[], cleavage_rule='Tryptic', tol=0.02):
        self.target = target
        self.target_length = len(target)

        self.precursor_mass_list = precursor_mass_list
        self.query=query
        self.query_length = len(query)
        self.modifications=modifications
        self.query_positions = self.query_positions()
        self.tolerance = tol
        self.validated=False
        
        if cleavage_rule=='Tryptic':
            self.tryptic_positions = self.tryptic_match()
            self.tryptic_peptides  = self.tryptic_peptides()
            self.tryptic_masses = self.tryptic_masses()
        
            self.validate_masses(self.tryptic_masses)

    def query_positions(self):
        positions = [(m.start(), m.start()+self.query_length ) for m in re.finditer('(?={})'.format(self.query), self.target)]
        return positions
    
    def validate_masses(self, mass_dict):
        validated_precursor = []
        precursor_peptides = defaultdict(list)
        validated_peptides = set()
        for position_key in mass_dict:
            target_mass = mass_dict[position_key]
            for query_mass in self.precursor_mass_list:
                diff = abs(target_mass - query_mass)
                if diff < self.tolerance:
                    validated_precursor.append(query_mass)
                    precursor_peptides[query_mass].append(self.tryptic_peptides[position_key])
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
            mass = peptide_mass(peptide)
            tryptic_masses[key] = mass
        return tryptic_masses

    def tryptic_peptides(self):
        tryptic_peptides = {}
        for coord in self.tryptic_positions:
            p = self.target[coord[0]:coord[1]]
            tryptic_peptides['{}:{}'.format(str(coord[0]), str(coord[1]))] = p
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

            
