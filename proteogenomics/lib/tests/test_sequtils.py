#!/usr/bin/env python

import sequtils
import Bio; from Bio import SeqIO
from io import StringIO
import unittest
import os
from Bio.Seq import Seq
import collections
import tempfile

testData = os.getcwd() + '/testdata'

contig  = SeqIO.read(testData + '/genome_contig.fasta','fasta')
contigs = list(SeqIO.parse(testData + '/S507_comb_assmbly_18_03_16.fasta','fasta'))

class Three_Frame(unittest.TestCase):


    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_three_frame_translate_with_all_codons(self):

        # Sanity check of the logic for the three frame breakdown of a nucleotide sequence
        # We do not want to throw away any bases by mistake :)

        gene_str1 = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAATG'  #64 bases, 21.333333 codons,
        # checks the case with one base extra in first frame
        s1f1      = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAAT'
        s1f2      = 'CGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAATG'
        s1f3      = 'GGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAA'

        gene_str2 = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAAT'   #63 bases, 21 codons
        # checks the case with no bases extra in the first frame
        s2f1      = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAAT'
        s2f2      = 'CGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTA'
        s2f3      = 'GGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAA'

        gene_str3 = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAA'    #62 bases, 20.66666 codons
        # checks the case with two bases extra in the first frame
        s3f1      = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCT'
        s3f2      = 'CGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTA'
        s3f3      = 'GGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAA'

        three_frame_template = {gene_str1: (s1f1, s1f2, s1f3),
                               gene_str2: (s2f1, s2f2, s2f3),
                               gene_str3: (s3f1, s3f2, s3f3),
                               }
        for raw_string in three_frame_template:
            three_frame_results = sequtils.Three_Frames(raw_string)
            template_results    = three_frame_template[raw_string]
            set_three_frame_results = set(three_frame_results)
            set_template_results = set(template_results)
            intersect_set = set_template_results.intersection(set_three_frame_results)

            for frame in three_frame_results:
                self.assertEqual(len(frame) % 3, 0)  # makes sure each each frame is divisible by three
            self.assertEqual(len(three_frame_results), len(template_results))
            self.assertEqual(len(three_frame_results), 3, 'Check that three frames are returned')
            self.assertEqual(len(intersect_set), 3, 'Makes sure that the template results match to the Three_Frame function')

class Six_Frames_and_Reverse_Complement(unittest.TestCase):

    def setUp(self):
            pass

    def tearDown(self):
        pass

    def test_Six_Frames_and_Reverse_Complement(self):
        gene_str1 = 'CCGGGTACTCTCGGGAAGTCCATGGAGGCGTGCCAGAGCGGCCGAATGGGGCTCACTGCTAATG'
        template_fwd = sequtils.Three_Frames(gene_str1)
        template_rev_seq = str(Seq(gene_str1).reverse_complement())
        template_rev = sequtils.Three_Frames(template_rev_seq)

        out_combined, out_rev_seq = sequtils.Six_Frames_and_Reverse_Complement(gene_str1)

        self.assertEqual(template_fwd + template_rev, out_combined)
        self.assertEqual(out_rev_seq, template_rev_seq)

class six_frames(unittest.TestCase):
    '''This code tests the six frame database against EMBOSS sixpack and sms2 with a few different options'''
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def assert_two_sequence_lists_are_identical(self, database1, database2):
        self.assertEqual(sorted(database1) == sorted(database2), True, '{} and {} are not the same.'.format(database1,database2))

    def get_strings(self, database):
        database = [str(rec.seq) for rec in database]
        return database

    def strip_trailing_char(self, list, char):
        list = [st[:-1] if st.endswith(char) else st for st in list ]
        return list

    def test_six_frame_database_any_size_protein_all_codons_prok(self):
        emboss = list(SeqIO.parse(testData +'/emboss_sixpack_min1_table11_translated.fasta', 'fasta'))
        genome = str(contig.seq)
        table=11
        codons = sequtils.codon_list(table) # all possible codons
        translated = True
        emboss_seqs = self.get_strings(emboss)
        emboss_seqs = self.strip_trailing_char(emboss_seqs, 'X')

        test = sequtils.six_frame(genome, peptide_length=1, table=table, codons=codons, translated=translated)
        test_seqs = self.get_strings(test)

        self.assert_two_sequence_lists_are_identical(test_seqs, emboss_seqs)

    def test_six_frame_database_min20aa_protein_all_codons_prok(self):
        emboss = list(SeqIO.parse(testData +'/genome_contig_emboss_sixpack_prok_allCodons_min20aa.fasta', 'fasta'))
        genome = str(contig.seq)
        table=11
        peptide_length=20
        codons = sequtils.codon_list(table) # all possible codons
        emboss_seqs = self.get_strings(emboss)
        emboss_seqs = self.strip_trailing_char(emboss_seqs, 'X')
        emboss_seqs = [rec for rec in emboss_seqs if len(rec) >= peptide_length]

        test = sequtils.six_frame(genome, peptide_length=peptide_length, table=table, codons=codons, translated=True)
        test_seqs = self.get_strings(test)

        self.assert_two_sequence_lists_are_identical(test_seqs, emboss_seqs)

    def test_six_frame_database_min20aa_protein_ATG_prok(self):
        genome = str(SeqIO.read(testData + '/sms2genome.fasta', 'fasta').seq)
        codons = ['ATG', 'GTG', 'CTG', 'TTG']
        table = 11
        peptide_length = 20
        stops = ['TAA', 'TGA', 'TAG' ]
        fwd = list(SeqIO.parse(testData + '/sms2genome_prok_20aa_ctg_gtg_ttg_atg_fwdStrand.fasta','fasta'))
        rev = list(SeqIO.parse(testData +'/sms2genome_prok_20aa_ctg_gtg_ttg_atg_revStrand.fasta','fasta'))
        combined_sms = rev + fwd
        sms_ORF  = [rec for rec in combined_sms if rec.id.startswith('ORF')]
        sms_ORF_strings = self.get_strings(sms_ORF)
        sms_PROT = [rec for rec in combined_sms if rec.id.startswith('Translation')]
        sms_PROT_strings = self.get_strings(sms_PROT)
        sms_PROT_strings = self.strip_trailing_char(sms_PROT_strings, '*')
        sms_PROT_strings = [str for str in sms_PROT_strings if len(str) >= 20]
        new = []
        for seq in sms_ORF_strings:
            if seq[-3:] in stops:
                if len(seq) >= 63: #this makes sure that 30 codons are coding (= 20 aa)
                    new.append(seq)
            else:
                if len(seq) >= 60: # no stop codon at the end so equivalent to 20 aa
                    new.append(seq)
        sms_ORF_strings = new

        test_raw = sequtils.six_frame(genome, peptide_length=peptide_length, table=table, codons=codons, translated=False)
        test_translated = sequtils.six_frame(genome, peptide_length=peptide_length, table=table, codons=codons, translated=True)
        test_raw_strings = self.get_strings(test_raw)
        test_translated_strings = self.get_strings(test_translated)

        self.assert_two_sequence_lists_are_identical(test_translated_strings, sms_PROT_strings)
        self.assert_two_sequence_lists_are_identical(test_raw_strings, sms_ORF_strings)


def main():
    unittest.main()

if __name__ == '__main__':
    main()

