#!/usr/bin/env python3

import shutil
import sequtils
import Bio; from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO
import unittest
import os
from Bio.Seq import Seq, translate
import collections
import tempfile
import time
import numpy as np
import blast


class blast2genome(unittest.TestCase):

    def setUp(self):
        self.genome = SeqIO.read("testdata/H37Rv.fasta","fasta")

    def tearDown(self):
        pass

    def test_blast2genome(self):
        fasta_id = "orf_sequence_1"
        query='PATGAVTGRWGSERAVSVPRVGVLALQGDTREHLAALRECGAEPTTVRRRDELDAVDALVIPGGESTTMSHXXXXXXXXXXXXXXXXXXXPAYGSCAGMILLASEILDAGAAGRQALPLRAMNMTVRRNAFGSQVDSFEGDIEFAGLDDPVRAVFIRAPWVERVGDGVQVLARAAGHIVAVRQGAVLATAFHPEMTGDRRIHQLFVDIVTS'
        hit='PATGAVTGRWGSERAVSVPRVGVLALQGDTREHLAALRECGAEPMTVRRRDELDAVDALVIPGGESTTMSHLLLDLDLLGPLRARLADGLPAYGSCAGMILLASEILDAGAAGRQALPLRAMNMTVRRNAFGSQVDSFEGDIEFAGLDDPVRAVFIRAPWVERVGDGVQVLARAAGHIVAVRQGAVLATAFHPEMTGDRRIHQLFVDIVTS'
        match='PATGAVTGRWGSERAVSVPRVGVLALQGDTREHLAALRECGAEP TVRRRDELDAVDALVIPGGESTTMSHLLLDLDLLGPLRARLADGLPAYGSCAGMILLASEILDAGAAGRQALPLRAMNMTVRRNAFGSQVDSFEGDIEFAGLDDPVRAVFIRAPWVERVGDGVQVLARAAGHIVAVRQGAVLATAFHPEMTGDRRIHQLFVDIVTS'
        
        hit_from=2931702
        hit_to=2932334
        hit_strand='-'
        
        b2g = blast.blast2genome(fasta_id = fasta_id, query_seq = query,match_seq=match, hit_seq = hit, hit_from=hit_from, hit_to=hit_to, hit_strand=hit_strand) 
        self.assertEqual(len(b2g.records), 2)
    
    def test_evaluate_start(self):
        fasta_id = "orf_sequence_1"
        nucs = 'gacttgctcttgcccacgatcgcaccgcagcgggcg'
        self.assertEqual(len(str(nucs)) % 3 , 0 )
        prot = translate(Seq(nucs[3:-3]))
        self.assertEqual(str(prot), 'LLLPTIAPQR')
        hit_from=3
        hit_to=33
        hit_strand='+'
        query='LLLPTIAPQR'
        hit='LLLPTIAPQR'
        match='LLLPTIA QR'
        
        # Forward strand
        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None, segment="LLLPTIA", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=0)
        self.assertEqual(newstart, 4)
        self.assertEqual(newend, 24)
        self.assertEqual(newoffset, 7)
        self.assertEqual(len(str(nucs[newstart-1:newend])) % 3 , 0 )
        self.assertEqual(translate(nucs[newstart-1:newend]), "LLLPTIA")
        
        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None, segment="QR", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=7)
        self.assertEqual(newstart,  28)
        self.assertEqual(newend, 33)
        self.assertEqual(len(str(nucs[newstart-1:newend])) % 3 , 0 )
        self.assertEqual(translate(nucs[newstart-1:newend]), "QR")
        self.assertEqual(newoffset, 10) 
        
        # Reverse strand
        hit_strand='-'
        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None,segment="LLLPTIA", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=0)
        self.assertEqual(newoffset, 7)
        self.assertEqual(newend, 33)
        self.assertEqual(newstart, 13)
        self.assertEqual(len(str(nucs[newstart-1:newend])) % 3 , 0 )
        
        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None, segment="QR", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=7)
        self.assertEqual(newoffset, 10) 
        self.assertEqual(newstart, 4)
        self.assertEqual(newend, 9)
        self.assertEqual(len(str(nucs[newstart-1:newend])) % 3 , 0 )

    def test_evaluate_start_2(self):
        fasta_id = "orf_sequence_1"
        hit_from=3
        hit_to=33
        hit_strand='+'
        query='TTTTTTTTTT'
        hit=  'TTTTVVTTTT'
        match='TTTT  TTTT'
        
        # Forward strand
        newstart, newend, newoffset= blast.blast2genome.evaluate_start(None, segment="TTTT", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=0)
        self.assertEqual(newstart, 4)
        self.assertEqual(newend, 15)
        self.assertEqual(newoffset, 4)
        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None, segment="TTTT", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=4)
        self.assertEqual(newstart,  22)
        self.assertEqual(newend, 33)
        
        # Reverse strand
        hit_strand='-'
        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None, segment="TTTT", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=0)
        self.assertEqual(newend, 33)
        self.assertEqual(newstart, 22)
        self.assertEqual(newoffset, 4)

        newstart, newend, newoffset = blast.blast2genome.evaluate_start(None, segment="TTTT", hit_seq=hit, hit_start=4, hit_end=33, hit_strand=hit_strand, hit_amino_offset=4)
        self.assertEqual(newstart, 4)
        self.assertEqual(newend, 15)
    
    def test_peptide_to_orf_alignment_positive_strand(self):

        orf_seq = 'DEVHNLLHRAKPSSQCSLRHVIALLAPGQGSQTEGMLSPWLQLPGAADQIAAWSKAADLDLARLGTTASTEEITDTAVAQPLIVAATLLAHQELARRCVLAGKDVIVAGHSVGEIAAYAIAGVIAADDAVALAATRGAEMAKACATEPTGMSAVLGGDETEVLSRLEQLDLVPANRNAAGQIVAAGRLTALEKLAEDPPAKARVRALGVAGAFHTEFMAPALDGFAAAAANIATADPTATLLSNRDGKPVTSAAAAMDTLVSQLTQPVRWDLCTATLREHTVTAIVEFPPAGTLSGIAKRELRGVPARAVKSPADLDELANL'
        hit_strand = '+'
        hit_start=2516727
        hit_end=2517695 

        fasta_id='PEPTIDE_1'
        query_seq='IALLAPGQGSQTEGMLSPWLQLPGAADQIPPWSK'
        match_seq='IALLAPGQGSQTEGMLSPWLQLPGAADQI++WSK'
        hit_seq='IALLAPGQGSQTEGMLSPWLQLPGAADQIAAWSK'
        amino_offset = orf_seq.find(hit_seq)
        b2g = blast.blast2genome(fasta_id = fasta_id, query_seq = query_seq, match_seq=match_seq, hit_seq = orf_seq, hit_from=hit_start, hit_to=hit_end, hit_strand=hit_strand, amino_offset=amino_offset)
        self.assertEqual(len(b2g.records), 2)
        for rec in b2g.records:
            segment = rec[0]
            strand = rec[3]
            start = rec[1]
            end = rec[2]
            nucs = str(self.genome.seq)[start-1:end]
            test_aminos = translate(nucs, table=11, cds=False)
            #test_aminos = translate(nucs.reverse_complement(), table=11, cds=False)
            self.assertEqual(segment, test_aminos)

    def test_peptide_to_orf_alignment_negative_strand(self):
        orf_seq='GARAHQGPQDPYPQVQEQDWLPQTAGTPSAADGPEGHRHRITEATDMAHKKGASSSRNGRDSAAQRLGVKRYGGQVVKAGEILVRQRGTKFHPGVNVGRGGDDTLFAKTAGAVEFGIKRGRKTVSIVGSTTA'
        hit_strand='-'
        hit_start=2739772
        hit_end=2740170

        fasta_id="PEPTIDE_1"

        query_seq   ="YGGPVVK"
        hit_seq     ="YGGQVVK"
        match_seq   ="YGG+VVK"
        
        amino_offset = orf_seq.find(hit_seq)
        b2g = blast.blast2genome(fasta_id = fasta_id, query_seq = query_seq, match_seq=match_seq, hit_seq = orf_seq, hit_from=hit_start, hit_to=hit_end, hit_strand=hit_strand, amino_offset=amino_offset)
        self.assertEqual(len(b2g.records), 2)
        for rec in b2g.records:
            segment = rec[0]
            strand = rec[3]
            start = rec[1]
            end = rec[2]
            nucs = str(self.genome.seq)[start-1:end]
            test_aminos = translate(Seq(nucs).reverse_complement(), table=11, cds=False)
            self.assertEqual(segment, test_aminos)

        
class pblast_mult(unittest.TestCase):

    def setUp(self):
        self.genome = list(SeqIO.parse("testdata/H37Rv.fasta","fasta"))
        self.orfs = sequtils.sf_contigs(self.genome, assembly_name = 'H37Rv', table=11 , codons='All', peptide_length=1, translated=True )
        self.orfs = [ i for i in self.orfs if len(str(i.seq)) >= 100 ]
        self.orfs = self.orfs[:100]
        self.outdir = tempfile.mkdtemp()

    def tearDown(self):     
        shutil.rmtree(self.outdir)

    def test_pblast_mult(self):
        res = blast.pblast_mult(self.orfs, self.orfs, self.outdir)
        self.assertEqual(len(res.blast_records), 100)
        self.assertEqual(len(list(res.aln_dict.keys())), 100)
        test_id = self.orfs[0].id
        self.assertEqual(res.aln_dict[test_id], ['H37Rv_Chromosome|H37Rv_Chromosome_recno_1.0|(+)1:1524'] )
        self.assertEqual(res.sorted_map[test_id], 'H37Rv_Chromosome|H37Rv_Chromosome_recno_1.0|(+)1:1524' )


if __name__ == '__main__': unittest.main()
