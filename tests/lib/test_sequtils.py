#!/usr/bin/env python

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

testData = os.getcwd() + '/testdata'

contig  = SeqIO.read(testData + '/genome_contig.fasta','fasta')
contig_peptides = open(testData +'/genome_contig_peptides.txt').read().split('\n')

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
        self.assertEqual(sorted(database1) == sorted(database2), True) # '{} and {} are not the same.'.format(database1,database2))

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

        test = sequtils.six_frame(genome, assembly_name = 'testassembly', contig_name='test_contig', peptide_length=1, table=table, codons=codons, translated=translated, methionine_start=False)
        test_seqs = self.get_strings(test)

        self.assert_two_sequence_lists_are_identical(emboss_seqs, test_seqs)

    def test_six_frame_database_min20aa_protein_all_codons_prok(self):
        emboss = list(SeqIO.parse(testData +'/genome_contig_emboss_sixpack_prok_allCodons_min20aa.fasta', 'fasta'))
        genome = str(contig.seq)
        table=11
        peptide_length=20
        codons = sequtils.codon_list(table) # all possible codons
        emboss_seqs = self.get_strings(emboss)
        emboss_seqs = self.strip_trailing_char(emboss_seqs, 'X')
        emboss_seqs = [rec for rec in emboss_seqs if len(rec) >= peptide_length]

        test = sequtils.six_frame(genome, assembly_name='assembly_name', contig_name='contig_name', peptide_length=peptide_length, table=table, codons=codons, translated=True, methionine_start=False)
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

        test_raw = sequtils.six_frame(genome,assembly_name= 'assembly_name',contig_name='contig_name', peptide_length=peptide_length, table=table, codons=codons, translated=False)
        test_translated = sequtils.six_frame(genome, assembly_name='assembly_name',contig_name='contig_name', peptide_length=peptide_length, table=table, codons=codons, translated=True, methionine_start=False)
        test_raw_strings = self.get_strings(test_raw)
        test_translated_strings = self.get_strings(test_translated)

        self.assert_two_sequence_lists_are_identical(test_translated_strings, sms_PROT_strings)
        self.assert_two_sequence_lists_are_identical(test_raw_strings, sms_ORF_strings)
    
    def test_six_frame_coordinates_correct(self):

        codons = ['ATG', 'GTG', 'TTG']
        table = 11
        peptide_length = 20
        stops = ['TAA', 'TGA', 'TAG' ]
        genome = str(SeqIO.read(testData + '/genome_contig.fasta', 'fasta').seq)
        test_translated = sequtils.six_frame(genome, assembly_name='testassembly',contig_name='testcontig', peptide_length=peptide_length, table=table, codons=codons, translated=True, methionine_start=True)
        
        for seq in test_translated:
            prot = str(seq.seq)
            strand = seq.id.split(')')[0].split('(')[1]
            coords = seq.id.split(')')[1].split(':')
            start = int(coords[0]) - 1
            end = int(coords[1])

            genome_seq = Seq(genome[start:end])
            
            if strand == '+':
                start_codon = genome_seq[:3]
                prot_seq = ''.join(str(translate(genome_seq, cds = False, table=11)).split('*'))
            elif strand =='-':
                start_codon = genome_seq.reverse_complement()[:3]
                prot_seq = ''.join(str(translate(genome_seq.reverse_complement(),cds=False,table=11)).split('*'))
            if start_codon in codons:
                prot_seq = 'M' + prot_seq[1:]
            self.assertEqual(prot_seq, prot, 'Coordinates do not match sequence')
    
    def test_translated_and_untranslated_six_frame_correspond(self):
        table = 11
        peptide_length = 20
        codons =['ATG','GTG','TTG']

        contig = str(list(SeqIO.parse(testData + '/genome_contig.fasta', 'fasta'))[0].seq)
        six_frame_translated = sequtils.six_frame(contig,assembly_name='testassembly', contig_name='testcontig',codons=codons,peptide_length=peptide_length, table=11, translated=True, methionine_start=True)
        six_frame_untranslated = sequtils.six_frame(contig,assembly_name='assembly_name', contig_name='contig_name',  codons=codons, peptide_length=peptide_length, table=11, translated=False)

        translist = []
        untranslist = []

        for rec in six_frame_translated:
            translist.append(str(rec.seq))
        for rec in six_frame_untranslated:
            prot = str(translate(rec.seq, table=11, cds=False))
            if prot.endswith('*'):
                prot = prot[:-1]
            self.assertTrue(rec.seq[:3] in codons, 'Start codon not in given start list')
            prot = 'M' + prot[1:]
            untranslist.append(prot)
        self.assert_two_sequence_lists_are_identical(translist, untranslist)

    def test_alt_tss_gets_right_records(self):
        table = 11
        peptide_length = 20
        assembly_name='assembly_name'
        contigs = list(SeqIO.parse(testData + '/genome_contig.fasta', 'fasta'))
        six_frame = sequtils.sf_contigs(contigs, assembly_name=assembly_name, table=11)

        self.assertTrue(len(six_frame) > 0)
        
        contig_dict = {}
        for rec in contigs:
            id = assembly_name +'_' + '_'.join(rec.id.split('|'))    
            contig_dict[id] = str(rec.seq)

        starts = ['ATG', 'GTG', 'TTG']
        for sfrec in six_frame:
            raw = str(sfrec.seq)
            alt_recs = sequtils.alt_tss(sfrec, starts = starts, translated=False)
            start = 3
            alt_starts = [raw]
            test_alt_starts = []
            while start <= len(raw)-3:
                if raw[start:start+3] in starts:
                    alt_starts.append(raw[start:])
                start += 3

            for rec in alt_recs:
                #print(rec.format('fasta'))
                
                self.assertTrue(str(rec.seq) in alt_starts)
                test_alt_starts.append(str(rec.seq))

                strand = rec.id.split(')')[0].split('(')[1]
                coords = rec.id.split(')')[1].split(':')
                start = int(coords[0]) - 1
                end = int(coords[1])
                contigmap = rec.id.split('|')[0]
                genome = contig_dict[contigmap]
                genome_seq = Seq(genome[start:end])
                
                if strand == '+':
                    seq_ = genome_seq
                elif strand =='-':
                    seq_ = genome_seq.reverse_complement()
                self.assertTrue(len(str(seq_)) % 3.0 == 0)
                self.assertEqual(str(seq_), str(rec.seq), 'Coordinates do not match sequence')
            self.assertEqual(set(alt_starts), set(test_alt_starts),'Incorrect alternative TSS identified')
        
    def test_translate_starts(self):

        table = 11
        peptide_length = 20
        contigs = list(SeqIO.parse(testData + '/genome_contig.fasta', 'fasta'))
        six_frame = sequtils.sf_contigs(contigs = contigs, assembly_name = 'testgenome',table=table)
            
        starts = ['ATG', 'GTG', 'TTG']
        for sfrec in six_frame:
            seq = sfrec.seq

            if seq[:3] in starts:
                tseq = 'M' + str(translate(seq, table = table))[1:]
            else:
                tseq = str(translate(seq, table = table))
            if tseq.endswith('*'):
                tseq = tseq[:-1]
            
            test_seq =str(sequtils.translate_start(sfrec, table=table, starts=starts).seq)
            self.assertEqual(tseq,test_seq)
            self.assertTrue(len(sfrec.seq) % 3.0 == 0)
        self.assertTrue(len(six_frame) > 0)

    def test_sf_contigs_and_six_frame_correspond(self):

        table = 11
        peptide_length = 20
        contigs = list(SeqIO.parse(testData + '/genome_contig.fasta', 'fasta'))
        sixframe = sequtils.sf_contigs(contigs = contigs, assembly_name = 'testgenome',  table=11)
        
        for contig in contigs:
            pass

class pairwise_blast(unittest.TestCase):
    
    def setUp(self):
        self.tempdir=tempfile.mkdtemp() 
    
    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_pairwise_blast_non_matched_first_amino_query_shorter_query(self):
        query_str='PGGAVATSYAQRVSIFTK'
        target_str='NLGGAVATSYAQRVSIFTKIINRELPGRFVYEDDDVVAFLTIEPMTQGHTLVVPRAEIDHWQNVDPALFGRVMSVSQLIGKAVCRAFSTQRAGMIIAGLEVPHLHIHVFPTRSLSDFGFANVDRNPSPGSLDEAQAKIRAALAQLA'
        res=sequtils.pairwise_blast(SeqRecord(id='Query', seq=Seq(query_str)), SeqRecord(id='Target', seq=Seq(target_str)),temp_folder=self.tempdir)
        
        self.assertEqual(res.differences, [2])
        self.assertEqual(res.variants, {2: 'L->P'})

    def test_pairwise_blast_non_matched_inner_amino_query_shorter_query(self):
        query_str='YGTATTFSQAR'
        target_str='RAGVFGAIPPGSRRRPARCRVPVRPVRGADGRGGPLVPRRTDLHHVLVIGSGPIVIGQACEFDYSGTQACRVLRAEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVERVIAQQAERGNKIDALLATLGGQTALNTAVALYESGVLEKYGVELIGADFDAIQRGEDRQRFKDIVAKAGGESARSRVCFTMAEVRETVAELGLPVVVRPSFTMGGLGSGIAYSTDEVDRMAGAGLAASPSANVLIEESIYGWKEFELELMRDGHDNVVVVCSIENVDPMGVHTGDSVTVAPAMTLTDREYQRMRDLGIAILREVGVDTGGCNIQFAVNPRDGRLIVIEMNPRVSRSSALASKATGFPIAKIAAKLAIGYTLDEIVNDITGETPACFEPTLDYVVVKAPRFAFEKFPGADPTLTTTMKSVGEAMSLGRNFVEALGKVMRSLETTRAGFWTAPDPDGGIEEALTRLRTPAEGRLYDIELALRLGATVERVAEASGVDPWFIAQINELVNLRNELVAAPVLNAELLRRAKHSGLSDHQIASLRPELAGEAGVRSLRVRLGIHPVYKTVDTCAAEFEAQTPYHYSSYELDPAAETEVAPQTERPKVLILGSGPNRIGQGIEFDYSCVHAATTLSQAGFETVMVNCNPETVSTDYDTADRLYFEPLTFEDVLEVYHAEMESGSGGPGVAGVIVQLGGQTPLGLAHRLADAGVPIVGTPPEAIDLAEDRGAFGDLLSAAGLPAPKYGTATTFAQARRIAEEIGYPVLVRPSYVLGGRGMEIVYDEETLQGYITRATQLSPEHPVLVDRFLEDAVEIDVDALCDGAEVYIGGIMEHIEEAGIHSGDSACALPPVTLGRSDIAKVRKATEAIAHGIGVVGLLNVQYALKDDVLYVLEANPRASRTVPFVSKATAVPLAKACARIMLGATIAQLRAEGLLAVTGDGAHAARNAPIAVKEAVLPFHRFRRADGAAIDSLLGPEMKSTGEVMGIDRDFGSAFAKSQTAAYGSLPAQGTVFVSVANRDKRSLVFPVKRLADLGFRVLATEGTAEMLRRNGIPCDDVRKHFEPAQPGRPTMSAVDAIRAGEVNMVINTPYGNSGPRIDGYEIRSAAVAGNIPCITTVQGASAAVQGIEAGIRGDIGVRSLQELHRVIGGVER'
        res=sequtils.pairwise_blast(SeqRecord(id='Query', seq=Seq(query_str)), SeqRecord(id='Target', seq=Seq(target_str)),temp_folder=self.tempdir)
        
        self.assertEqual(res.differences, [749])
        self.assertEqual(res.variants, {749: 'A->S'})

    def test_pairwise_blast(self):
        query_str='TYLHGIVGEIVEPERIDAYLDRGPEMLSFVLKHTPLKMCWVPGYSDYYPEAPGGRPGGRSIEPKPFNARKLGADMAGLEPAYGKVPLNVVVMQQDYVRLNQLKRHPRGVLRSMKVGARTMWAKATGKNLVGMGRALIGPLRIGLQRAGVPVELNTAFTDLFVENGVVSGVYVRDSHEAESAEPQLIRARRGVILACGGFEHNEQMRIKYQRAPITTEWTVGASANTGDGILAAEKLGAALDLMDDAWWGPTVPLVGKPWFALSERNSPGSIIVNMSGKRFMNESMPYVEACHHMYGGEHGQGPGPGENIPAWLVFDQRYRDRYIFAGLQPGQRIPSSWLDSGVIVQADTLAELAGKAGLPADELTATVQRFNAFARSGVDEDYHRGESAYDRYYGDPSNKPNPNLGEVGHPPYYGAKMVPGDLGTKGGIRTDVNGRALRDDGSIIDGLYAAGNVSAPVMGHTYPGPGGTIGPAMTFGYLAALHIADQAGKRP'
        target_str='MTVQEFDVVVVGSGAAGMVAALVAAHRGLSTVVVEKAPHYGGSTARSGGGVWIPNNEVLKRRGVRDTPEAARTYLHGIVGEIVEPERIDAYLDRGPEMLSFVLKHTPLKMCWVPGYSDYYPEAPGGRPGGRSIEPKPFNARKLGADMAGLEPAYGKVPLNVVVMQQDYVRLNQLKRHPRGVLRSMKVGARTMWAKATGKNLVGMGRALIGPLRIGLQRAGVPVELNTAFTDLFVENGVVSGVYVRDSHEAESAEPQLIRARRGVILACGGFEHNEQMRIKYQRAPITTEWTVGASANTGDGILAAEKLGAALDLMDDAWWGPTVPLVGKPWFALSERNSPGSIIVNMSGKRFMNESMPYVEACHHMYGGEHGQGPGPGENIPAWLVFDQRYRDRYIFAGLQPGQRIPSRWLDSGVIVQADTLAELAGKAGLPADELTATVQRFNAFARSGVDEDYHRGESAYDRYYGDPSNKPNPNLGEVGHPPYYGAKMVPGDLGTKGGIRTDVNGRALRDDGSIIDGLYAAGNVSAPVMGHTYPGPGGTIGPAMTFGYLAALHIADQAGKRE'

        res=sequtils.pairwise_blast(SeqRecord(id='Query', seq=Seq(query_str)), SeqRecord(id='Target', seq=Seq(target_str)),temp_folder=self.tempdir)
        self.assertEqual(res.differences, [564])
        self.assertEqual(res.variants, {564: 'E->P'})

    def test_pairwise_blast_2(self):
        query_str='VSAALLAADAGVPVLLAPAADAATALADASVGTVFAAR'
        target_str='mrsphrdairtarglvvkVGTTALTTPSGMFDAGRlaglaeaverrmkagsdvvivssgaiaagieplglsrrpkdlatkqaaasvgqvalvnswsaafarygrtvgqvlltahdismrvqhtnaqrtldrlrALHAVAIVNENDTVATNEIRfgdndrlsalvahlvgadalvllsdidglydcdprKTADATFIPEVSGPADLDGVVAGRsshlgtggmaskvaaallaadagvpvllapaadaataladasvgtvfaarparlsarrfwvrYAAEATGALTLDAGAVRavvrqrrsllaagitavsgrfcggdvvelrapdaamvargvvaydaselatmvgrstselpgelrRPVVHADDLVAVSAKqakqv'
        target_str = target_str.upper()

        res=sequtils.pairwise_blast(SeqRecord(id='Query', seq=Seq(query_str)), SeqRecord(id='Target', seq=Seq(target_str)),temp_folder=self.tempdir)
        print(target_str[225])
        self.assertEqual(res.differences, [226])
        self.assertEqual(res.variants, {226: 'A->S'})

    def test_pairwise_blast_3(self):
        query_str='RFRPGGGLTVPRTDNDSWAITESVGATALGVAAARAAETESDNPLINDPFARIFVDAAGDGIWSMYTNRTLLAGATDLDPDLRAPIQQMIDFMAARTAFFDEYFLATADAGVRQVVILASGLDSRAWRLPWPDGTVVYELDQPKVLEFKSATLRQHGAQPASQLVNVPIDLRQDWPKALQKAGFDPSKPCAWLAEGLVRYLPARAQDLLFERIDALSRPGSWLASNVPGAGFLDPERMRRQRADMRRMRAAAAKLVETEISDVDDLWYAEQRTAVAEWLRERGWDVSTATLPELLARYGRSIPHSGEDSIPPNLFVSAQRATS*'
        target_str='HRHHENQREHRMSAMRTHDDTWDIKTSVGATAVMVAAARAVETDRPDPLIRDPYARLLVTNAGAGAIWEAMLDPTLVAKAAAIDAETAAIVAYLRSYQAVRTNFFDTYFASAVAAGIRQVVILASGLDSRAYRLDWPAGTIVYEIDQPKVLSYKSTTLAENGVTPSAGRREVPADLRQDWPAALRDAGFDPTARTAWLAEGLLMYLPAEAQDRLFTQVGAVSVAGSRIAAETAPVHGEERRAEMRARFKKVADVLGIEQTIDVQELVYHDQDRASVADWLTDHGWRARSQRAPDEMRRVGRWVEGVPMADDPTAFAEFVTAERL*'
        target_str = target_str.upper()
        res=sequtils.pairwise_blast(SeqRecord(id='Query', seq=Seq(query_str)), SeqRecord(id='Target', seq=Seq(target_str)),temp_folder=self.tempdir)

class frameshift_peptides(unittest.TestCase):

    def setUp(self):
        self.tempdir=tempfile.mkdtemp() 
    
    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_frameshift_peptide_identification(self):
        table=11
        frame='ATCCAAGTGTGCTCCCAGTGCGGAACTGGCTGGAACGTCCGTGAGCGGCAACGCGTGTGGTGTCCACGCTGTCGGGGGATGTTGCTGGCGCCGTTGGCGGATATGCCGGCCGAGGCGCGCTGGCGTACACCGGCGCGCCCGCAGGTGCCGACGGCCTCCGACACGCGGCGCACACCGCCGCGGCTTCCCCCAGGTTTTCGGTGGATAGCGGTGCGACCCGGGGCGGCACCGCCGCCACGGCACGGCCCACGGTTACGTGGGCCTACTCCCCGCTACGCCGGGATACCGCGGTGGGGGCTAACCGACCATGTCGACCAGGCTCCCGTGCCGGCCTCGGCGAAAGCAGGGCCATCGCCGGCGGCGGTGCGCACCACGCTGTTGGTGAGCCTGCTGGTGTTCAGCATCGCGGTCGTGGTGTTTGTGGTGCGGTATGTGTTGCTGGTCATCAACCGAAACACATTGTTGAACTCGGTGGTGGCCAGCGCCTCGGTCTGGCTGGGGGTTTTGGTCAGCTTGGCGGCGATTGCGGCGGCCGGCACTACCATTGTCCTGTTGGTCCGGTGGCTGGTCGCCCGTCGGGCCGCCGCGTTTATGCATCAAGGCTTGCCGGAGCGGCGTTCCGCCCGTGAGTTATGGGCCGGCTGCCTATTGCCGATGGTCAATCTGCTGTGGGCTCCGCTGTACGTCATCGAGTTGGCGCTGGTCGAGGACCGCTACACGCGGCTGCGCAGGCCGATCGTGGTGTGGTGGATCGTGTGGATCGTCAGCAACGCGATATCGATGTTCGCGTTCGCCACCAGCTGGGTCACCGACGCTCAGGGCATCGCCAACAACACCACCATGATGGTGCTGGCGTATCTGTGTGCGGCGGCCGCGGTGGCCGCTGCTGCGCGGGTCTTCGAGGGGTTCGAGCAAAAGCCGGTCGAACGCCCAGCGCATCGCTGGGTGGTGGTGAACACAGACGGGCGTTCCGCGCCGGCATCTTCTGTTGCGGTTGAGTTGGACGGGCAGGAACCGGCAGCATAG'
        peptide_nucs=frame[90:180]  # Get a 90 aa chunk
        peptide= str(translate(Seq(peptide_nucs), table=table, cds=False))

        frame1=SeqRecord(seq=translate(Seq(frame[:135]  + 'G' + frame[135:-1]),table=table, cds=False), id='frame1', description='None')
        frame2=SeqRecord(seq=translate(Seq(frame[1:135] + 'G' + frame[135:]), table=table, cds=False), id='frame2', description='None')
        frame_original=SeqRecord(seq=translate(Seq(frame), table=table, cds=False), id='frame_orig', description='None')
        
        peptides = [peptide]
        frames = [frame1, frame2, frame_original]

        res = sequtils.frameshift_peptides(frames, peptides, self.tempdir)
        self.assertEqual(res.frameshifts, {'frame1 frame2': ['PLADMPAEARWRTPARPQVPTASDTRRTPP']})
        
    def test_frameshift_peptide_identification_1(self):
        
        contigs = SeqIO.to_dict(list(SeqIO.parse(testData + '/S5527_comb_assmbly_18_03_16.fasta','fasta')))
        contig = contigs['scaffold8|size164766']

        peptide='FLSGQSPTTIVAPPAAK'

        frame1_id='S5527_scaffold8_size164766|S5527_scaffold8_size164766_recno_5843.0|(-)34640:38587'
        frame1_seq='PATGAAHRAVGDRTSSRGTRKVPVRTVPDDHCGAPAAKTVELSVQVPVPDMANLTDNTIWPDVEARLVDLIESHNSTIVFANSRRLAERLTARLNEIHAARCGIELAPDTNQQVAGGAPAHIMGSGQTFGAPPVLARAHHGSISKEQRAVVEEDLKRGQLKAVVATSSLELGIDMGAVDLVIQVQAPPSVASGLQRIGRAGHQVGEISRGVLFPKHRTDLLGCAVSVQRMLAGEIETMRVPANPLDILAQHTVAAAALEPLDADAWFDTVRRAAPFATLPRSLFEATLDLLSGTYPSTEFAELRPRLVYDRDTGTLTARPGAQRLAVTSGGAIPDRGLFAVYLATERPSRVGELDEEMVYESRPGDVISLGATSWRITEITHDRVLVIPAPGQPARLPFWRGDDAGRPAELGAALGALTGELAALDRTAFGTRCAGLGFDDYATDNLWRLLDDQRTATAVVPTDSTLLVERFRDELGDWRVILHSPYGLRVHGPLALAVGRRLRDRYGIDEKPTASDNGIVVRLPDTVSAGEDSPPGAELFVFDADEIDPIVTTEVAGSALFASRFRESAARALLLPRRHPGRRSPLWQQRQRAARLLEVARKYPDFPIVLETVRECLQDVYDVPILVELMARIAQRRVRVAEAETAKPSPFAASLLFGYVGAFMYEGDTPLAERRAAALALDGTLLAELLGRVELRELLDPDVIAATSRQLQHLAADRVARDAEGVADLLRLLGPLTEDEIAARAGAPEVSGWLDGLRAAKRALVVSFAGRSWWVAVEDMGRLRDGVGAAVPVGLPASFTEAVADPLGELLGRYARTHTPFTTAAAAARFGLGLRVTADVLGRLASDGRLVRGEFVAAAKGSAGGEQWCDAEVLRILRRRSLAALRAQAEPVSTAAYGRFLPAWQHVSAGNSGIDGLAAVIDQLAGVRIPASAIEPLVLAPRIRDYSPAMLDELLASGDVTWSGAGSISGSDGWIALHPADSAPMTLAEPAEIDFTDAHRAILASLGTGGAYFFRQLTHDGLTEAELKAALWELIWAGRVTGDTFAPVRAVLGGAGTRKRAAPAHGGHRPPRLSRYRLTHAQARNADPTVAGRWSALPLPEPDSTLRAHYQAELLLNRHGVLTKDAVAAEGVAGGFATLYKVLSAFEDAGRCQRGYFIESLGGAQFAVASTVDRLRSYLDGVDPEQPDYHAVVLAAADPANPYGAALPWPASSADGTARPGRKAGALVVLVDGELAWFLERGGRSLLTFTDDPEANHAAAIGLADLVTAGRVASILVERADGMPVLQPGGRASAALTALLAAGFVRTPRGLRRR*'

        frame1=SeqRecord(seq=Seq(frame1_seq), id=frame1_id )
        
        frame2_id='S5527_scaffold8_size164766|S5527_scaffold8_size164766_recno_3805.0|(-)38425:39180'
        frame2_seq='VRFAQPSALSRFSALTRDWFTSTFAAPTAAQASAWAAIADGDNTLVIAPTGSGKTLAAFLWALDSLAGSEPMSERPAATRVLYVSPLKALAVDVERNLRTPLAGLTRLAERQGLPAPQIRVGVRSGDTPPALRRQLVSQPPDVLITTPESLFLMLTSAARQTLTGVQTVIIDEIHAIAATKRGAHLALSLERLDDLSSRRRAQRIGLSATVRPPEELARFLSGQSPTTIVAPRPPRPLSCPCRCRCPTWPT*'
        
        frame2=SeqRecord(seq=Seq(frame2_seq), id=frame2_id)        
        
        peptides = [peptide]
        frames = [frame1, frame2]

        res = sequtils.frameshift_peptides(frames, peptides, self.tempdir)
        for r in res.frameshift_results:
            ep = res.frameshift_results[r][0]['end_pep']
            ep_ns = res.frameshift_results[r][0]['end_rec_nuc_start']
            ep_ne = res.frameshift_results[r][0]['end_rec_nuc_end']
            ep_nstrand = res.frameshift_results[r][0]['end_rec_strand']
            nucs = contig.seq[ep_ns -1 : ep_ne]
            if ep_nstrand == '-':
                nucs = nucs.reverse_complement()
            self.assertEqual(str(translate(nucs, cds=False, table=11)), ep )

            sp = res.frameshift_results[r][0]['start_pep']
            sp_ns = res.frameshift_results[r][0]['start_rec_nuc_start']
            sp_ne = res.frameshift_results[r][0]['start_rec_nuc_end']
            sp_nstrand = res.frameshift_results[r][0]['start_rec_strand']
            nucs = contig.seq[sp_ns -1 : sp_ne]
            if sp_nstrand == '-':
                nucs = nucs.reverse_complement()
            self.assertEqual(str(translate(nucs, cds=False, table=11)), sp )
        self.assertEqual(res.frameshifts, {'{} {}'.format(frame2_id, frame1_id): [peptide]})
         

class pep2genome(unittest.TestCase):

    def setUp(self):
        self.tempdir=tempfile.mkdtemp() 
    
    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_reverse_strand_correct(self):
        nucs='ATCCAAGTGTGCTCCCAGTGCGGAACTGGCTGGAACGTCCGTGAGCGGCAACGCGTGTGGTGTCCACGCTGTCGGGGGATGTTGCTGGCGCCGTTGGCGGATATGCCGGCCGAGGCGCGCTGGCGTACACCGGCGCGCCCGCAGGTGCCGACGGCCTCCGACACGCGGCGCACACCGCCGCGGCTTCCCCCAGGTTTTCGGTGGATAGCGGTGCGACCCGGGGCGGCACCGCCGCCACGGCACGGCCCACGGTTACGTGGGCCTACTCCCCGCTACGCCGGGATACCGCGGTGGGGGCTAACCGACCATGTCGACCAGGCTCCCGTGCCGGCCTCGGCGAAAGCAGGGCCATCGCCGGCGGCGGTGCGCACCACGCTGTTGGTGAGCCTGCTGGTGTTCAGCATCGCGGTCGTGGTGTTTGTGGTGCGGTATGTGTTGCTGGTCATCAACCGAAACACATTGTTGAACTCGGTGGTGGCCAGCGCCTCGGTCTGGCTGGGGGTTTTGGTCAGCTTGGCGGCGATTGCGGCGGCCGGCACTACCATTGTCCTGTTGGTCCGGTGGCTGGTCGCCCGTCGGGCCGCCGCGTTTATGCATCAAGGCTTGCCGGAGCGGCGTTCCGCCCGTGAGTTATGGGCCGGCTGCCTATTGCCGATGGTCAATCTGCTGTGGGCTCCGCTGTACGTCATCGAGTTGGCGCTGGTCGAGGACCGCTACACGCGGCTGCGCAGGCCGATCGTGGTGTGGTGGATCGTGTGGATCGTCAGCAACGCGATATCGATGTTCGCGTTCGCCACCAGCTGGGTCACCGACGCTCAGGGCATCGCCAACAACACCACCATGATGGTGCTGGCGTATCTGTGTGCGGCGGCCGCGGTGGCCGCTGCTGCGCGGGTCTTCGAGGGGTTCGAGCAAAAGCCGGTCGAACGCCCAGCGCATCGCTGGGTGGTGGTGAACACAGACGGGCGTTCCGCGCCGGCATCTTCTGTTGCGGTTGAGTTGGACGGGCAGGAACCGGCAGCATAG'

        protein = str(translate(Seq(nucs).reverse_complement(), table = 11, cds= False))
        start = 1
        end = len(nucs)
        strand = '-'
        peptide = protein[30:60]
        contig='test'

        pg = sequtils.pep2genome(peptide, protein, strand, start, end, contig)
        
        genomic_start = pg.genomic_start
        genomic_end = pg.genomic_end

        pepnucs = Seq(nucs[genomic_start-1 : genomic_end]).reverse_complement()
        
        self.assertEqual(peptide, str(translate(pepnucs, cds=False, table=11)))
        self.assertEqual(len(str(pepnucs)) % 3, 0)

    def test_forward_strand_correct(self):
        nucs='ATCCAAGTGTGCTCCCAGTGCGGAACTGGCTGGAACGTCCGTGAGCGGCAACGCGTGTGGTGTCCACGCTGTCGGGGGATGTTGCTGGCGCCGTTGGCGGATATGCCGGCCGAGGCGCGCTGGCGTACACCGGCGCGCCCGCAGGTGCCGACGGCCTCCGACACGCGGCGCACACCGCCGCGGCTTCCCCCAGGTTTTCGGTGGATAGCGGTGCGACCCGGGGCGGCACCGCCGCCACGGCACGGCCCACGGTTACGTGGGCCTACTCCCCGCTACGCCGGGATACCGCGGTGGGGGCTAACCGACCATGTCGACCAGGCTCCCGTGCCGGCCTCGGCGAAAGCAGGGCCATCGCCGGCGGCGGTGCGCACCACGCTGTTGGTGAGCCTGCTGGTGTTCAGCATCGCGGTCGTGGTGTTTGTGGTGCGGTATGTGTTGCTGGTCATCAACCGAAACACATTGTTGAACTCGGTGGTGGCCAGCGCCTCGGTCTGGCTGGGGGTTTTGGTCAGCTTGGCGGCGATTGCGGCGGCCGGCACTACCATTGTCCTGTTGGTCCGGTGGCTGGTCGCCCGTCGGGCCGCCGCGTTTATGCATCAAGGCTTGCCGGAGCGGCGTTCCGCCCGTGAGTTATGGGCCGGCTGCCTATTGCCGATGGTCAATCTGCTGTGGGCTCCGCTGTACGTCATCGAGTTGGCGCTGGTCGAGGACCGCTACACGCGGCTGCGCAGGCCGATCGTGGTGTGGTGGATCGTGTGGATCGTCAGCAACGCGATATCGATGTTCGCGTTCGCCACCAGCTGGGTCACCGACGCTCAGGGCATCGCCAACAACACCACCATGATGGTGCTGGCGTATCTGTGTGCGGCGGCCGCGGTGGCCGCTGCTGCGCGGGTCTTCGAGGGGTTCGAGCAAAAGCCGGTCGAACGCCCAGCGCATCGCTGGGTGGTGGTGAACACAGACGGGCGTTCCGCGCCGGCATCTTCTGTTGCGGTTGAGTTGGACGGGCAGGAACCGGCAGCATAG'

        protein = str(translate(Seq(nucs), table = 11, cds= False))
        start = 1
        end = len(nucs)
        strand = '+'
        peptide = protein[30:60]
        contig='test'

        pg = sequtils.pep2genome(peptide, protein, strand, start, end, contig)
        
        genomic_start = pg.genomic_start
        genomic_end = pg.genomic_end
        pepnucs = Seq(nucs[genomic_start-1 : genomic_end])
        self.assertEqual(peptide, str(translate(pepnucs, cds=False, table=11)))
        self.assertEqual(len(str(pepnucs)) % 3, 0)

class proteogenomics(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
   
    def test_methionine_cleavage_start_site_annotation_validated(self):
        orf_id ='S5527_scaffold1|S5527_scaffold1_recno_6834.0|(+)262383:262553'
        orf_sequence='GCTCGCGTCATGACCCAACCCACCGCATGGGAGTACGCCACGGTCCCGCTGTTGACGCACGCCACCAAACAGATCCTCGACCAGTGGGGAGCCGATGGCTGGGAGCTGGTGGCGGTGCTGCCCGGGCCCACCGGTGAGCAGCACGTCGCTTACCTGAAGCGCCCGAAGTAG'
        peptides = ['TQPTAWEYATVPLLTHATK']
        orf_rec=SeqRecord(seq=Seq(orf_sequence), id=orf_id)
        ref_id='I6X824_MYCTU'
        ref_sequence='MTQPTAWEYATVPLLTHATKQILDQWGADGWELVAVLPGPTGEQHVAYLKRPK'
        ref_sequence = ref_sequence.upper()
        ref_rec = SeqRecord(id=ref_id, seq=Seq(ref_sequence))
        pg = sequtils.proteogenomics(peptides, orf_rec, ref_rec)
        self.assertEqual(pg.tss_peptides, peptides)
        self.assertEqual(pg.met_ap_peptides, peptides)
        self.assertEqual(len(pg.identified_tss_sites),1)
        validation_type='Identified TSS by Non-tryptic N-terminus peptide TQPTAWEYATVPLLTHATK following an ATG start-codon (initiator Methionine cleavage). Annotated TSS validated.'
        self.assertEqual(list(pg.variant_sequences)[0].description, validation_type)
        self.assertEqual(pg.annotation_type, ['Annotated TSS validated'])
    
    def test_methionine_cleavage_start_site_annotation_upstream(self):
        orf_id ='S5527_scaffold1|S5527_scaffold1_recno_6834.0|(+)262383:262553'
        orf_sequence='GCTCGCGTCATGACCCAACCCACCGCATGGGAGTACGCCACGGTCCCGCTGTTGACGCACGCCACCAAACAGATCCTCGACCAGTGGGGAGCCGATGGCTGGGAGCTGGTGGCGGTGCTGCCCGGGCCCACCGGTGAGCAGCACGTCGCTTACCTGAAGCGCCCGAAGTAG'
        peptides = ['TQPTAWEYATVPLLTHATK']
        orf_rec=SeqRecord(seq=Seq(orf_sequence), id=orf_id)
        ref_id='I6X824_MYCTU'
        ref_sequence='MTAWEYATVPLLTHATKQILDQWGADGWELVAVLPGPTGEQHVAYLKRPK'
        ref_sequence = ref_sequence.upper()
        ref_rec = SeqRecord(id=ref_id, seq=Seq(ref_sequence))
        pg = sequtils.proteogenomics(peptides, orf_rec, ref_rec)
        self.assertEqual(pg.tss_peptides, peptides)
        self.assertEqual(pg.met_ap_peptides, peptides)
        self.assertEqual(len(pg.identified_tss_sites),1)
        validation_type='Identified TSS by Non-tryptic N-terminus peptide TQPTAWEYATVPLLTHATK following an ATG start-codon (initiator Methionine cleavage). Upstream TSS identified.'
        self.assertEqual(list(pg.variant_sequences)[0].description, validation_type)
        self.assertEqual(pg.annotation_type, ['Upstream TSS identified'])
   
    def test_methionine_cleavage_start_site_annotation_downstream(self):
        orf_id ='S5527_scaffold1|S5527_scaffold1_recno_6834.0|(+)262383:262553'
        orf_sequence='GCTCGCGTCATGACCCAACCCACCGCATGGGAGTACGCCACGGTCCCGCTGTTGACGCACGCCACCAAACAGATCCTCGACCAGTGGGGAGCCGATGGCTGGGAGCTGGTGGCGGTGCTGCCCGGGCCCACCGGTGAGCAGCACGTCGCTTACCTGAAGCGCCCGAAGTAG'
        peptides = ['TQPTAWEYATVPLLTHATK']
        orf_rec=SeqRecord(seq=Seq(orf_sequence), id=orf_id)
        ref_id='I6X824_MYCTU'
        #ref_sequence='MTAWEYATVPLLTHATKQILDQWGADGWELVAVLPGPTGEQHVAYLKRPK'
        ref_sequence='MMTQPTAWEYATVPLLTHATKQILDQWGADGWELVAVLPGPTGEQHVAYLKRPK'
        ref_sequence = ref_sequence.upper()
        ref_rec = SeqRecord(id=ref_id, seq=Seq(ref_sequence))
        pg = sequtils.proteogenomics(peptides, orf_rec, ref_rec)
        self.assertEqual(pg.tss_peptides, peptides)
        self.assertEqual(pg.met_ap_peptides, peptides)
        self.assertEqual(len(pg.identified_tss_sites),1)
        validation_type='Identified TSS by Non-tryptic N-terminus peptide TQPTAWEYATVPLLTHATK following an ATG start-codon (initiator Methionine cleavage). Downstream TSS identified.'
        self.assertEqual(list(pg.variant_sequences)[0].description, validation_type)
        self.assertEqual(pg.annotation_type, ['Downstream TSS identified'])

    def test_methionine_cleavage_start_site_annotation_downstream_and_validation(self):
        orf_id ='S5527_scaffold1|S5527_scaffold1_recno_6834.0|(+)262383:262553'
        orf_sequence='GCTCGCGTGATGACCCAACCCACCGCATGGGAGTACGCCACGGTCCCGCTGTTGACGCACGCCACCAAACAGATCCTCGACCAGTGGGGAGCCGATGGCTGGGAGCTGGTGGCGGTGCTGCCCGGGCCCACCGGTGAGCAGCACGTCGCTTACCTGAAGCGCCCGAAGTAG'
        peptides = ['TQPTAWEYATVPLLTHATK', 'MMTQPTAWEYATVPLLTHATK']
        orf_rec=SeqRecord(seq=Seq(orf_sequence), id=orf_id)
        ref_id='I6X824_MYCTU'
        #ref_sequence='MTAWEYATVPLLTHATKQILDQWGADGWELVAVLPGPTGEQHVAYLKRPK'
        ref_sequence='MMTQPTAWEYATVPLLTHATKQILDQWGADGWELVAVLPGPTGEQHVAYLKRPK'
        ref_sequence = ref_sequence.upper()
        ref_rec = SeqRecord(id=ref_id, seq=Seq(ref_sequence))
        pg = sequtils.proteogenomics(peptides, orf_rec, ref_rec)
        self.assertEqual(pg.tss_peptides, peptides)
        self.assertEqual(pg.met_ap_peptides, [peptides[0]])
        self.assertEqual(len(pg.identified_tss_sites),2)
        #validation_type='Identified TSS by Non-tryptic N-terminus peptide TQPTAWEYATVPLLTHATK following an ATG start-codon (initiator Methionine cleavage). Downstream TSS identified.'
        self.assertEqual(len(list(pg.variant_sequences)), 2)
        self.assertEqual(set(pg.annotation_type),set(['Downstream TSS identified', 'Annotated TSS validated']))


    def test_downstream_start_site_annotation(self):

        contigs = SeqIO.to_dict(list(SeqIO.parse(testData + '/S507_comb_assmbly_18_03_16.fasta','fasta')))
        
        contig = contigs['scaffold2|size353395']
        orf_id='S507_scaffold2|S507_scaffold2_recno_6883.0|(-)329053:329808'

        orf_sequence='GCTGGTGTCGCTTCCCCGCCGTGTTGGAGGCGACCTGTACGGAAGCTTCCTGCCCGGCACTCAGCAGAACTTGGAGCGTGCGCTGGACGGCTTGCTAGAGCTGCTCCCTGCGGGCGCTTGGCTAGATCACACCTCAGATCACGCACAAGCCTCCTCCCGAGGCTGACCCCTCACATCTCCGCTACGACTTCAGAAAGGGACGCCATGGTGGACCCGCCGGGCAACGACGACGACCACGGTGATCTCGACGCCCTCGATTTCTCCGCCGCCCACACCAACGAGGCGTCGCCGCTGGACGCCTTAGACGACTATGCGCCGGTGCAGACCGATGACGCCGAAGGCGACCTGGACGCCCTCCATGCGCTCACCGAACGCGACGAGGAGCCGGAGCTGGAGTTGTTCACGGTGACCAACCCTCAAGGGTCGGTGTCGGTCTCAACCCTGATGGACGGCAGAATCCAGCACGTCGAGCTGACGGACAAGGCGACCAGCATGTCCGAAGCGCAGCTGGCCGACGAGATCTTCGTTATTGCCGATCTGGCCCGCCAAAAGGCGCGGGCGTCGCAGTACACGTTCATGGTGGAGAACATCGGTGAACTGACCGACGAAGACGCAGAAGGCAGCGCCCTGCTGCGGGAATTCGTGGGGATGACCCTGAATCTGCCGACGCCGGAAGAGGCTGCCGCAGCCGAAGCCGAAGTGTTCGCCACCCGCTACGATGTCGACTACACCTCCCGGTACAAGGCCGATGACTGA'

        peptides = ['MVENIGELTDEDAEGSALLR']
        
        orf_rec=SeqRecord(seq=Seq(orf_sequence), id=orf_id)
       
        ref_id='sp|O69732|ESPH_MYCTU'
        ref_sequence='mvdppgndddhgdldaldfsaahtneaspldalddyapvqtddaegdldalhalterdeepelelftvtnpqgsvsvstlmdgriqhveltdkatsmseaqladeifviadlarqkarasqytfMVENIGELTDEDAEGSALLREFVGMTLNLPTPEEAAAAEAEVFATRydvdytsrykadd'
        ref_sequence = ref_sequence.upper()
        ref_rec = SeqRecord(id=ref_id, seq=Seq(ref_sequence))

        pg = sequtils.proteogenomics(peptides, orf_rec, ref_rec)
       
        #print(pg.start)
        #print(pg.end)
        #print(pg.strand)
        #print(pg.novel_peptides)
        #print(pg.annotated_peptides)
        #print(pg.tss_peptides)
        #print(pg.met_ap_peptides)
        #print(pg.identified_tss_sites)
        #print(pg.peptide_starts)
        #print(pg.other_peptidase)
        
        for rec in pg.variant_sequences:
            #print(rec.format('fasta'))
            strand = rec.id.split('|')[2].split(')')[0].split('(')[1]
            start = int(rec.id.split(')')[1].split(':')[0])
            end = int(rec.id.split(')')[1].split(':')[1])
            if strand == '+':
                temp_orf = Seq(str(contig.seq)[start-1:end])
                translated_temp_orf = str(translate(temp_orf, table=11, cds=True))
            elif strand == '-':
                temp_orf = Seq(str(contig.seq)[start-1:end]).reverse_complement()
                translated_temp_orf = str(translate(temp_orf, table=11, cds=True))
            self.assertEqual(translated_temp_orf[1:], str(rec.seq)[1:])
        
        #print('\n'.join(pg.variant_sequences_trie))
        #print(dir(pg.pairwise_blast.hsps[0].expect))
        #print(pg.pairwise_blast.hsps[0].sbjct_start)
        #print(pg.pairwise_blast.hsps[0].sbjct_end)
        #print(pg.pairwise_blast.hsps[0].query_start)
        #print(pg.pairwise_blast.hsps[0].query_end)
        #print(pg.reference_sequence)
        self.assertEqual(pg.tss_peptides, ["MVENIGELTDEDAEGSALLR"])
        self.assertEqual(len(pg.identified_tss_sites),1)
        self.assertEqual(pg.annotation_type,['Downstream TSS identified'])

        #print(pg.reference.format('fasta'))

    def test_upstream_start_site_annotation(self):

        contigs = SeqIO.to_dict(list(SeqIO.parse(testData + '/S507_comb_assmbly_18_03_16.fasta','fasta')))
        contig = contigs['scaffold1|size371796']

        peptides = ['AQSAQYAFILDRMSQQVDADEHRVALLR',
                    'KTVGETWGLPSPEEAAAAEAEVFATR',
                    'AQSAQYAFILDR',
                    'TVGETWGLPSPEEAAAAEAEVFATR',
                    'MSVSTLMDGRIDHVELSARVAWMSESQLASEILVIADLAR',
                    'TVSVSTLMDGRIDHVELSARVAWMSESQLASEILVIADLAR',
                    'VSVSTLMDGRIDHVELSARVAWMSESQLASEILVIADLAR',
                    'VAWMSESQLASEILVIADLAR' ]
        orf_sequence='CTGCCCACAATGCCCTGGGCTCGTCCTTGCATACGGCCGGTGTCGATCTCGCCAAAAGTCTTCGAATTGCGGCGAAGATATATAGCGAGGCCGACGAAGCGTGGCGCAAGGCTATCGACGGGTTGTTTACCTGACCACGTTTGCTGCCCGCAGTGCAGGCCACAGCGTCTTCCCAACGACCTGTTCGGACTGACCACGCCAGCTGCCCAGGCCGACCCTTCCCGGGTGGCAATGAATTCCGAAGGGACGGTGGACTTGCCCGGAAATGACTTTGACAGCAACGATTTCGACGCCGTGGATCTCTGGGGTGCCGACGGCGCGGAGGGCTGGACTGCGGATCCGATTATTGGCGTCGGGTCGGCGGCGACCCCGGACACCGGACCCGACCTGGACAATGCCCACGGTCAGGCGGAGACGGACACCGAACAAGAGATCGCGCTTTTTACCGTGACGAATCCCCCACGCACGGTGTCGGTATCGACGCTGATGGACGGCCGGATTGACCATGTCGAGCTGTCGGCCAGGGTGGCCTGGATGAGTGAGTCGCAGCTCGCTTCTGAGATCCTGGTGATTGCCGACCTGGCGCGGCAGAAGGCGCAGTCGGCCCAGTACGCCTTCATCCTTGACAGGATGAGTCAACAGGTCGATGCAGATGAACACCGCGTCGCACTGCTACGTAAGACCGTGGGCGAAACCTGGGGGTTACCATCGCCGGAAGAAGCCGCGGCAGCAGAAGCTGAGGTGTTCGCGACGCGCTACAGCGACGATTGTCCAGCACCAGACGACGAGAGCGATCCATGGTGA'
        orf_id='S507_scaffold1_size371796|S507_scaffold1_size371796_recno_1764.0|(+)264478:265281'
        orf_rec=SeqRecord(seq=Seq(orf_sequence), id=orf_id)
       
        ref_id='sp|P9WJD5|ESPD_MYCTU'
        ref_sequence='MDLPGNDFDSNDFDAVDLWGADGAEGWTADPIIGVGSAATPDTGPDLDNAHGQAETDTEQEIALFTVTNPPRTVSVSTLMDGRIDHVELSARVAWMSESQLASEILVIADLARQKAQSAQYAFILDRMSQQVDADEHRVALLRKTVGETWGLPSPEEAAAAEAEVFATRYSDDCPAPDDESDPW'
        ref_rec = SeqRecord(id=ref_id, seq=Seq(ref_sequence))

        pg = sequtils.proteogenomics(peptides, orf_rec, ref_rec)
       
        #print(pg.start)
        #print(pg.end)
        #print(pg.strand)
        #print(pg.novel_peptides)
        #print(pg.annotated_peptides)
        #print(pg.tss_peptides)
        #print(pg.met_ap_peptides)
        #print(pg.identified_tss_sites)
        #print(pg.peptide_starts)
        #print(pg.other_peptidase)
        
        for rec in pg.variant_sequences:
            #print(rec.format('fasta'))
            strand = rec.id.split('|')[2].split(')')[0].split('(')[1]
            start = int(rec.id.split(')')[1].split(':')[0])
            end = int(rec.id.split(')')[1].split(':')[1])
            if strand == '+':
                temp_orf = Seq(str(contig.seq)[start-1:end])
                translated_temp_orf = str(translate(temp_orf, table=11, cds=True))
            elif strand == '-':
                temp_orf = Seq(str(contig.seq)[start-1:end]).reverse_complement()
                translated_temp_orf = str(translate(temp_orf, table=11, cds=True))
            self.assertEqual(translated_temp_orf[1:], str(rec.seq)[1:])
        
        #print('\n'.join(pg.variant_sequences_trie))
        #print(dir(pg.pairwise_blast.hsps[0].expect))
        #print(pg.pairwise_blast.hsps[0].sbjct_start)
        #print(pg.pairwise_blast.hsps[0].sbjct_end)
        #print(pg.pairwise_blast.hsps[0].query_start)
        #print(pg.pairwise_blast.hsps[0].query_end)
        #print(pg.reference_sequence)
        #print(pg.annotation_type) 

        #print(pg.reference.format('fasta'))

class mapping2peptides(unittest.TestCase):

    def setUp(self):
        self.tempdir=tempfile.mkdtemp() 
    
    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_most_upstream_trie_correct(self):
        genome = [contig]
        peptides_list = contig_peptides
        translation_table=11
        g2p = sequtils.peptides2genome(genome, assembly_name = 'test', translation_table=translation_table, peptides_list=peptides_list, outdir = self.tempdir, threads=8)
        mapping = sequtils.mapping2peptides(g2p.peptides, translation_table)
        self.assertEqual(len(mapping.orf_mapping), 45)
        peps = mapping.get_peptides(['test_scaffold1_size371796|test_scaffold1_size371796_recno_3846.0|(-)4379:4996'])
        self.assertEqual(len(peps), 14)

class peptides2genome(unittest.TestCase):

    def setUp(self):
        self.tempdir=tempfile.mkdtemp() 
    
    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_most_upstream_trie_correct(self):
        genome = [contig]
        peptides_list = contig_peptides
        translation_table=11
        
        g2p = sequtils.peptides2genome(genome, assembly_name = 'test', translation_table=translation_table, peptides_list=peptides_list, outdir = self.tempdir, threads=8)
        
        mapped_orfs=45
        mapped_peptides = 540
        self.assertEqual(len(g2p.peptides), mapped_peptides)
        self.assertEqual(len(g2p.mapped_orfs), 45)
        self.assertEqual(len(g2p.mapped_trans_orfs), 45)

def main():
    unittest.main()

if __name__ == '__main__':
    main()

