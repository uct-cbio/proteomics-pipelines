#!/usr/bin/env python

import pandas as pd
import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
from collections import defaultdict

class peptidesmapped:
    def __init__(self, peptides, reference_lineage, taxon_peptides, taxon_data):
        self.peptides = peptides
        self.reference_lineage = reference_lineage
        self.taxon_peptides = taxon_peptides
        self.taxon_data = taxon_data
        self.mapped = self.mapped_entries()
        self.counted = Counter(self.mapped)
        self.best, self.best_count, self.best_fasta, self.peptides_counted = self.best_mapped()

    def mapped_entries(self):
        mapped = []
        for peptide in self.peptides:
            if peptide in self.taxon_peptides:
                entries = self.taxon_peptides[peptide]
                assert len(entries) == len(list(set(entries)))
                mapped += entries
        return mapped

    def best_mapped(self):
        if len(self.mapped) > 0:
            mapped_fastas = {}
            counted = self.counted
            biggest = max(counted.values())
            overlap_max = defaultdict(list)
            for key in counted:
                if counted[key] == biggest:
                    dict = self.taxon_data[key]
                    val = '{} - {} Taxid: {}'.format(dict['Entry'] ,dict['Organism'],dict['Organism ID'])
                    seq = Seq(dict['Sequence'])
                    rec = SeqRecord(seq = seq, id = dict['Entry'], description = dict['Organism'    ])    
                    mapped_fastas[dict['Entry']] = rec
                    lineage = dict['Taxonomic lineage IDs (all)'].split(', ')
                    lineage  = [int(i) for i in lineage]
                    overlap  = set(lineage) & set(self.reference_lineage)
                    overlap_max[len(overlap)].append(val)
            new_biggest = max(overlap_max.keys())
            newmapped = '\n'.join(overlap_max[new_biggest])
            spec_fastas = []
            for datum in overlap_max[new_biggest]:
                d = datum.split()[0]
                fs = mapped_fastas[d]
                spec_fastas.append(fs)
            return newmapped, biggest, spec_fastas, counted
        else:  
            return None, None, None, None

    def count_best_mapped(self, reference_list): # Return highest peptide count
        max_ref = []
        for ref in reference_list:
            refc = self.peptides_counted[ref]
            max_ref.append(refc)
        if len(max_ref) > 0:
            max_ref = max(max_ref)
        else:
            max_ref = 0
        return max_ref

    def pep2proteome(self, proteome_id):
        mapped_proteome = []
        
        if len(self.mapped) > 0:
            for entry in self.mapped:
                try:
                    dict = self.taxon_data[entry]
                    p = dict['Proteomes'].split(':')[0]
                    protein = dict['Entry']
                    genes = dict['Gene names']
                    count = self.counted[entry]
                    map = (protein, genes, count)
                    if p == proteome_id:
                        mapped_proteome.append(map)
                except:
                    pass

        return mapped_proteome           
            

def map2fasta(ids, mapping):
    recs = []
    taxonomy = None
    for id in ids:
        header = id
        id = id.split('|')[1]
        data = mapping[id]
        seq = data['Sequence']
        organism = data['Organism']
        name = data['Protein names']
        description = "{} {}".format(organism, name)
        rec = SeqRecord(seq = Seq(seq), id = header, description = description)
        recs.append(rec)
    return recs


