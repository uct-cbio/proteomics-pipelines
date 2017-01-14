#!/usr/bin/env python
import time
import pandas as pd
import collections
from collections import defaultdict
import re
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import re

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
    return np.round(mw, 4)

class TagMatch:
    def __init__(self, query, tag_mass_list, target, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], enzymes=['Trypsin'], specificity='specific', prec_tol=0.02, gap_tol=0.5, max_missed_cleavages=3):
        self.possible_enzymes = ['Trypsin', 'Trypsin, no P rule','Whole protein']
        self.enzymes = enzymes
        self.max_missed_cleavages = max_missed_cleavages
        for enzyme in self.enzymes:
            assert enzyme in self.possible_enzymes

        self.possible_specificity = ['specific', 'semi-specific', 'unspecific']
        self.specificity = specificity
        assert self.specificity in self.possible_specificity
        
        self.target = target
        self.target_length = len(target)

        self.tag_mass_list = tag_mass_list
        
        self.gap_mass_dict = defaultdict(set)

        if not ((self.specificity=='specific') and (self.enzymes==['Trypsin'])):
            self.gaps = self.gaps()
            self.max_gaps = (max(self.gaps[0]), max(self.gaps[1]))
            self.maxngap = self.max_gaps[0]
            self.maxcgap = self.max_gaps[1]

        self.query=query
        self.query_length = len(query)
        self.fixed_modifications=fixed_modifications
        self.variable_modifications=variable_modifications
        self.query_positions = self.query_positions()
        self.prec_tol = prec_tol
        self.gap_tol  = gap_tol

        self.validated=False
        self.specificity = specificity

        self.positions = self.match()
        self.peptides  = self.peptides()
        self.masses = self.masses()
        
        self.amino_gaps = self.amino_gaps()
        self.validate_masses(self.masses)
            
    def query_positions(self):
        positions = [(m.start(), m.start()+self.query_length ) for m in re.finditer('(?={})'.format(self.query), self.target)]
        return positions
    
    def validate_masses(self, mass_dict):
        validated_precursor = []
        precursor_peptides = defaultdict(list)
        validated_peptides = set()
        
        for position_key in mass_dict:
            target_masses = mass_dict[position_key]
            target_peptide = self.peptides[position_key]
            
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
                            if not amino_gaps[0] in self.gap_mass_dict:
                                nterm_amino_gaps = peptide_mass(amino_gaps[0], nterm=False, cterm=False, fixed_modifications=self.fixed_modifications, variable_modifications = self.variable_modifications)
                            else:
                                nterm_amino_gaps = self.gap_mass_dict[amino_gaps[0]]
                            for nterm_amino_gap in nterm_amino_gaps:
                                if abs(nterm_amino_gap - nterm_mass_gap) < self.gap_tol:
                                    nterm_gap_match=True
                            if not amino_gaps[1] in self.gap_mass_dict:
                                cterm_amino_gaps = peptide_mass(amino_gaps[1], nterm=False, cterm=False, fixed_modifications = self.fixed_modifications, variable_modifications=self.variable_modifications)
                            else:
                                cterm_amino_gaps = self.gap_mass_dict[amino_gaps[1]]
                            for cterm_amino_gap in cterm_amino_gaps:
                                if abs(cterm_amino_gap - cterm_mass_gap) < self.gap_tol:
                                    cterm_gap_match=True
                            
                            if (cterm_gap_match==True) and (nterm_gap_match==True):
                                validated_precursor.append(tag_mass_set)
                                precursor_peptides[tag_mass_set].append(self.peptides[position_key])
                                validated_peptides.add(self.peptides[position_key])
        self.validated_precursor = validated_precursor
        self.precursor_peptides  = precursor_peptides
        self.validated_peptides = list(validated_peptides)
        if len(validated_precursor) > 0:
            self.validated = True

    def masses(self):
        masses = {}
        for key in self.peptides:
            peptide = self.peptides[key]
            masslist = peptide_mass(peptide, fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=True, cterm=True)
            masses[key] = masslist
        return masses

    def peptides(self):
        peptides = {}
        for coord in self.positions:
            p = self.target[coord[0]:coord[1]]
            key='{}:{}'.format(str(coord[0]), str(coord[1]))
            peptides[key] = p
        return peptides

    def match(self):
        coords = []
        for coord in self.query_positions:
            nterms = self.nterms(coord[0])
            cterms = self.cterms(coord[1]-1)

            for nterm in nterms:
                for cterm in cterms:
                    amino_acid_before = self.target[nterm -1: nterm]
                    first_amino_acid  = self.target[nterm]
                    last_amino_acid   = self.target[cterm -1]
                    amino_acid_after   = self.target[cterm : cterm + 1]
                    valid = self.valid_cleavage( amino_acid_before, first_amino_acid, last_amino_acid, amino_acid_after)
                    if valid == True:
                        mc = self.missed_cleavages(self.target[nterm:cterm])
                        if (mc <= self.max_missed_cleavages):
                            coords.append((nterm, cterm))
        return coords


    def nterms(self, pos):
        nterm = pos
        nterms=[]
        
        if ((self.specificity=='specific') and (self.enzymes==['Trypsin'])):
            valid_gaps=True
            missed=0
            while valid_gaps==True:
                if nterm == 0:
                    nterms.append(nterm)
                    return nterms
                else: 
                    current_amino = self.target[nterm]
                    previous_amino = self.target[nterm-1]
                    if (current_amino != 'P') and (previous_amino in {'R', 'K'}):
                        nterms.append(nterm)
                        missed += 1
                        if missed > self.max_missed_cleavages:
                            return nterms
                    nterm -= 1
                    assert nterm >= 0
        else:
            valid_gaps=True
            limit = self.maxngap + self.gap_tol
            while valid_gaps==True:
                amino_gap = self.target[nterm:pos]
                validated=False
                gaps = peptide_mass(amino_gap, fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=False, cterm=False)
                self.gap_mass_dict[amino_gap].update(gaps)
                for gap in gaps:
                    diff_gaps = set([(abs(i-gap) < self.gap_tol) for i in self.gaps[0]])
                    if True in diff_gaps:
                        validated=True

                if nterm == 0:
                    if validated == True:
                        nterms.append(nterm)
                    return nterms
                else: 
                    if validated == True:
                        nterms.append(nterm)
                    nterm -= 1
                    assert nterm >= 0
                if min(gaps) > limit:
                    return nterms

    def cterms(self, pos):
        cterm = pos
        cterms=[]
        
        if ((self.specificity=='specific') and (self.enzymes==['Trypsin'])):
            valid_gaps =True
            missed=0
            while valid_gaps ==True:
                amino_gap = self.target[pos +1:cterm+1]

                if cterm == self.target_length-1:
                    cterms.append(cterm + 1)
                    return cterms
                else: 
                    current_amino = self.target[cterm]
                    next_amino = self.target[cterm +1]
                    if (next_amino != 'P') and (current_amino in {'R', 'K'}):
                        cterms.append(cterm + 1)
                        missed += 1
                        if missed > self.max_missed_cleavages:
                            return cterms
                    cterm += 1
                    assert cterm <= self.target_length -1

        else:
            valid_gaps =True
            limit = self.maxcgap + self.gap_tol
            while valid_gaps ==True:
                amino_gap = self.target[pos +1:cterm+1]
                gaps = peptide_mass(amino_gap, fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications, nterm=False, cterm=False)
                self.gap_mass_dict[amino_gap].update(gaps)
                validated=False
                
                for gap in gaps:
                    diff_gaps = set([(abs(i-gap) < self.gap_tol) for i in self.gaps[1]])
                    if True in diff_gaps:
                        validated=True

                if cterm == self.target_length-1:
                    if validated == True:
                        cterms.append(cterm + 1)
                    return cterms
                else: 
                    if validated == True:
                        cterms.append(cterm + 1)
                    cterm += 1
                    assert cterm <= self.target_length -1
                if min(gaps) > limit:
                    return cterms

    def gaps(self):
        ngaps = []
        cgaps = []
        for tag in self.tag_mass_list:
            n = tag[0]
            c = tag[2]
            ngaps.append(n)
            cgaps.append(c)
        return (list(set(ngaps)), list(set(cgaps)))

    def amino_gaps(self):
        amino_gap_dict = defaultdict(list)
        for key in self.peptides:
            peptide = self.peptides[key]
            positions = [(m.start(), m.start()+self.query_length ) for m in re.finditer('(?={})'.format(self.query), peptide)]
            for coord in positions:
                tag_start = coord[0]
                tag_end = tag_start + self.query_length
                nterm_gap = peptide[:tag_start]
                cterm_gap = peptide[tag_end:]
                amino_gap_dict[peptide].append( (nterm_gap, cterm_gap) )
        return amino_gap_dict

    def valid_cleavage(self, amino_acid_before, first_amino_acid, last_amino_acid, amino_acid_after):
        valid = False
        nterm_valid=False
        cterm_valid=False
        for enzyme in self.enzymes:
            if enzyme=='Trypsin, no P rule':
                cleavage_aminos=['K','R']
                if amino_acid_before == '':
                    nterm_valid=True
                elif amino_acid_before in cleavage_aminos:
                    nterm_valid=True
                if amino_acid_after =='':
                    cterm_valid=True
                elif last_amino_acid in cleavage_aminos:
                    cterm_valid=True

            elif enzyme == 'Trypsin':
                cleavage_aminos=['K','R']
                
                if amino_acid_before == '':
                    nterm_valid=True
                elif (amino_acid_before in cleavage_aminos) and (first_amino_acid != 'P'):
                    nterm_valid=True
                
                if amino_acid_after =='':
                    cterm_valid=True
                elif last_amino_acid in cleavage_aminos:
                    cterm_valid=True
            
            elif enzyme == 'Whole protein':
                
                if amino_acid_before == '':
                    nterm_valid=True
                
                if amino_acid_after =='':
                    cterm_valid=True
        if ((nterm_valid and cterm_valid) == True) and (self.specificity == 'specific'):
            valid = True
        elif ((nterm_valid or cterm_valid) == True) and (self.specificity == 'semi-specific'):
            valid = True
        elif self.specificity == 'unspecific':
            valid=True
        return valid
    
    def missed_cleavages(self, peptide):
        
        new_peptides = [peptide]

        for enzyme in self.enzymes:
            holder=[]
            for datum in new_peptides:
                if enzyme=='Trypsin, no P rule':
                    peptides = re.sub(r'(?<=[RK])','\n', datum).split('\n')
                    holder += peptides

                elif enzyme == 'Trypsin':
                    
                    peptides = re.sub(r'(?<=[RK])(?=[^P])','\n', datum).split('\n')
                    holder += peptides

                elif enzyme == 'Whole protein':
                    peptides = [datum]
                    holder += peptides
            new_peptides = holder
        return len(new_peptides) - 1

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
                        newpossible.add(np.round(newmass,4))
        if (len(newpossible) > 0) and (self.valid ==False):
            self.possible=newpossible
            self.validate()

class blast_tags:
    def __init__(self, subject, tags, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02, gap_tol=0.5, min_identities = 4, IL_equivalence=True):

        self.subject=subject
        self.tags=tags 
        self.fixed_modifications = fixed_modifications
        self.variable_modifications = variable_modifications
        self.prec_tol  = prec_tol
        self.gap_tol = gap_tol

        self.min_identities = min_identities
        self.IL_equivalence = IL_equivalence

        self.mweights = set()
        self.nmasses  = set()
        self.cmasses  = set()
        
        self.samples  = set()
        self.scans    = set()

        record_holder = set()
        self.newrecords = []


        for qtag in self.tags:
            qseq = str(qtag.seq)
            qsubstrings = self.substrings(qseq, self.subject)
            
            modified_tags = self.segment_matching(qseq, self.subject, min_identities=self.min_identities, IL_equivalence=self.IL_equivalence)
            idlist = qtag.id.split(';')
            file = idlist[0]
            scan = idlist[1]
            mw = float(idlist[2].split('mw=')[1])
            ngap = float(idlist[3].split('ngap=')[1])
            cgap = float(idlist[4].split('cgap=')[1])

            for i in modified_tags:
                newseq=i[0]
                newngap = np.round(ngap + i[1],6)
                newcgap = np.round(cgap + i[2],6)
                newseq, newngap, newcgap = self.fix_negative_gaps(newseq, newngap, newcgap)
                valid_tag=True
                if newngap <= -(self.gap_tol):
                    valid_tag=False
                elif newcgap <= -(self.gap_tol):
                    valid_tag=False
                if valid_tag == True:
                    newid ='{};{};mw={};ngap={};cgap={}'.format(file, scan, mw, newngap, newcgap)
                    newrec=SeqRecord(seq=Seq(newseq),id=newid,description='Segment correction of {}-{}-{}'.format(ngap,qseq,cgap))
                    temp = newrec.format('fasta')
                    if not temp in record_holder:
                        record_holder.add(temp)
                        self.newrecords.append(newrec)

    def valid_tag(self, tag):
        sequence_masses = peptide_mass(tag[0], fixed_modifications=self.fixed_modifications, variable_modifications = self.variable_modifications, cterm=True, nterm=True)
        for mass in sequence_masses:
            total = mass + tag[2] + tag[3]
            diff = np.absolute(total - tag[1])
            if diff < self.prec_tol:
                return True

    def substrings(self, reference, match_sequence):
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
            if subs != '':
                substrings.append(subs)
        return substrings

    def longest_substring(self, substrings):    
        maxlen = len(max(substrings, key=len))
        maxes = [i for i in substrings if len(i) == maxlen]
        return maxes
    
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
        test = []
        for newtag in self.newtags:
            rec = SeqRecord(seq=Seq(newtag[0]), id = '{};{};mw={};ngap={};cgap={}'.format(list(self.samples)[0], list(self.scans)[0], str(newtag[1]), str(newtag[2]), str(newtag[3])))
            trec = rec.format('fasta')
            if not trec in test:
                test.append(trec)
                newrecords.append(rec)
        return newrecords
    
    def match_segments(self, tagsequence, blastsequence, substrings, min_identities, matchsequence='', matchngap=None, matchcgap=None, identities=[], matches = set()):
        if (len(substrings) == 0):
            if len(identities) > 0:
                if (len(''.join(identities)) >= min_identities):
                    for cg in peptide_mass(blastsequence, cterm=False, nterm=False, fixed_modifications = self.fixed_modifications, variable_modifications = self.variable_modifications):
                        for cg2 in peptide_mass(tagsequence, cterm=False, nterm=False, fixed_modifications = self.fixed_modifications, variable_modifications = self.variable_modifications ) :
                            matches.add((matchsequence + blastsequence, matchngap, cg2-cg, '-'.join(identities)))
            return matches

        subst = substrings[0]
        newsubstrings = substrings[1:]
        refpos = [(m.start(), m.start()+len(subst)) for m in re.finditer('(?={})'.format(subst), blastsequence)]
        tagpos = [(m.start(), m.start()+len(subst)) for m in re.finditer('(?={})'.format(subst), tagsequence)]
        
        newtags = self.match_segments(tagsequence, blastsequence, newsubstrings, min_identities, matchsequence=matchsequence, matchngap=matchngap, matchcgap=matchcgap, identities=identities, matches=matches)
        matches.update(newtags)
        
        for r in refpos:
            for t in tagpos:
                refgap = blastsequence[:r[0]]
                matchgap = tagsequence[:t[0]]
                newblastsequence = blastsequence[r[1] :]
                newtagsequence = tagsequence[t[1] :]

                refgapmasses=peptide_mass(refgap, nterm=False, cterm=False, fixed_modifications=self.fixed_modifications, variable_modifications=self.variable_modifications)
                matchgapmasses = peptide_mass(matchgap, nterm=False, cterm=False, fixed_modifications=self.fixed_modifications, variable_modifications = self.variable_modifications)

                for rm in refgapmasses:
                    for tm in matchgapmasses:
                        if matchngap==None:
                            newmatchngap = tm - rm
                            newmatchngap = np.round(newmatchngap, 6)
                            newmatchsequence = refgap + subst
                            newtags = self.match_segments(newtagsequence, newblastsequence, newsubstrings, min_identities, matchsequence=newmatchsequence, matchngap=newmatchngap, matchcgap = matchcgap, identities=identities + [subst] , matches=matches)
                            matches.update(newtags)
                        elif rm==tm==0:
                            return matches
                        else:
                            diff = abs(rm - tm)
                            if diff < self.gap_tol:
                                newmatchsequence = matchsequence + refgap + subst
                                newtags = self.match_segments(newtagsequence, newblastsequence, newsubstrings, min_identities, matchsequence=newmatchsequence, matchngap=matchngap, matchcgap=matchcgap, identities=identities + [subst], matches=matches)
                                matches.update(newtags)
        return matches

    def segment_matching(self, tagsequence, blastsequence, min_identities=4, IL_equivalence=True):
        if IL_equivalence==True:
            newtagsequence = 'L'.join(tagsequence.split('I'))
            newblastsequence = 'L'.join(blastsequence.split('I'))
        else:
            newtagsequence = tagsequence
            newblastsequence = blastsequence

        substrings = self.substrings(newblastsequence, newtagsequence)
        newtags = self.match_segments(newtagsequence, newblastsequence, substrings, min_identities=min_identities, matchsequence='', matchngap=None, matchcgap=None, identities=[], matches=set())
        besttags=[]
        bestlen=0
        
        for tag in newtags:
            identities = ''.join(tag[3].split('-'))
            if len(identities) > bestlen:
                bestlen = len(identities)
        
        for tag in newtags:
            identities = ''.join(tag[3].split('-'))
            if len(identities) == bestlen:
                assert tag[0] == newblastsequence
                tag = (blastsequence, tag[1], tag[2], tag[3])
                besttags.append(tag)
        return set(besttags)
        
    def fix_negative_gaps(self, sequence, ngap, cgap):
        
        if ngap < 0:
            ngap_fixed=False
            for position in range(1, len(sequence)):
                newsequence = sequence[position:]
                newgapsequence = sequence[:position]
                ng = peptide_mass(newgapsequence, nterm=False, cterm=False,  fixed_modifications=self.fixed_modifications, variable_modifications= self.variable_modifications)
                mng = min(ng)
                
                if mng + ngap >= self.gap_tol:
                    break    
                
                for _ in ng:
                    newg = ngap + _
                    if np.absolute(newg) < self.gap_tol:
                        ngap = newg
                        sequence = newsequence
                        ngap_fixed=True
                        break
                if ngap_fixed==True:
                    break

        if cgap < 0:
            cgap_fixed=False
            for position in range(1, len(sequence)):
                newsequence = sequence[:-position]
                newgapsequence = sequence[-position:]
                cg = peptide_mass(newgapsequence, nterm=False, cterm=False,  fixed_modifications=self.fixed_modifications, variable_modifications= self.variable_modifications)
                
                mcg = min(cg)
                
                if mcg + cgap >= self.gap_tol:
                    break    
    #            
                for _ in cg:
                    newg = cgap + _
                    if np.absolute(newg) < self.gap_tol:
                        cgap = newg
                        sequence = newsequence
                        cgap_fixed=True
                        break
                if cgap_fixed==True:
                    break

        return sequence, np.round(ngap,4), np.round(cgap,4)

