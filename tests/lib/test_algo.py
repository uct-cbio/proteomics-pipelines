#!/usr/bin/env python3

import shutil
import algo
import Bio; from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO
import unittest
import os
from Bio.Seq import Seq, translate
import collections
import tempfile
import time
import pandas as pd
import re

testData = os.getcwd() + '/testdata'

contig  = SeqIO.read(testData + '/genome_contig.fasta','fasta')
contigs = list(SeqIO.parse(testData + '/S507_comb_assmbly_18_03_16.fasta','fasta'))
peptides = pd.read_csv(testData + '/peptides.txt', sep=None, engine='python')
peptides = peptides[(peptides['Potential contaminant'] != '+') & (peptides['Reverse'] != '+')]
peptides = peptides['Sequence'].tolist()
translated_orfs = list(SeqIO.parse(testData + '/S507_S5527_combined_allstarts_proteins.fasta','fasta'))


class trie(unittest.TestCase):
    def setUp(self):
        self.peptides=['ADGAAIDSLLGPEMK', 'AEGLLAVTGDGAHAAR','AEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVER', 'AGEVNMVINTPYGNSGPR', 'AGFWTAPDPDGGIEEALTR', 'ATQLSPEHPVLVDR', 'DDVLYVLEANPR', 'DLGIAILR', 'ETVAELGLPVVVRPSFTMGGLGSGIAYSTDEVDR', 'EVGVDTGGCNIQFAVNPR','FAFEKFPGADPTLTTTMK','FPGADPTLTTTMK','GAFGDLLSAAGLPAPK','HFEPAQPGRPTMSAVDAIR','HSGLSDHQIASLRPELAGEAGVR','IAEEIGYPVLVR','IAEEIGYPVLVRPSYVLGGR','IMLGATIAQLR','LADAGVPIVGTPPEAIDLAEDR','LADLGFR','LGIHPVYK','LIVIEMNPR','LYDIELALR','MAGAGLAASPSANVLIEESIYGWK','NELVAAPVLNAELLR','NFVEALGK','RADGAAIDSLLGPEMK','SLVFPVK','SLVFPVKR','SQTAAYGSLPAQGTVFVSVANR','SVGEAMSLGR','VAEASGVDPWFIAQINELVNLR','VLATEGTAEMLR','VLILGSGPNR','YGTATTFSQAR','YGVELIGADFDAIQR']
        self.Trie = algo.Trie(peptides)
        self.protein='VPRRTDLHHVLVIGSGPIVIGQACEFDYSGTQACRVLRAEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVERVIAQQAERGNKIDALLATLGGQTALNTAVALYESGVLEKYGVELIGADFDAIQRGEDRQRFKDIVAKAGGESARSRVCFTMAEVRETVAELGLPVVVRPSFTMGGLGSGIAYSTDEVDRMAGAGLAASPSANVLIEESIYGWKEFELELMRDGHDNVVVVCSIENVDPMGVHTGDSVTVAPAMTLTDREYQRMRDLGIAILREVGVDTGGCNIQFAVNPRDGRLIVIEMNPRVSRSSALASKATGFPIAKIAAKLAIGYTLDEIVNDITGETPACFEPTLDYVVVKAPRFAFEKFPGADPTLTTTMKSVGEAMSLGRNFVEALGKVMRSLETTRAGFWTAPDPDGGIEEALTRLRTPAEGRLYDIELALRLGATVERVAEASGVDPWFIAQINELVNLRNELVAAPVLNAELLRRAKHSGLSDHQIASLRPELAGEAGVRSLRVRLGIHPVYKTVDTCAAEFEAQTPYHYSSYELDPAAETEVAPQTERPKVLILGSGPNRIGQGIEFDYSCVHAATTLSQAGFETVMVNCNPETVSTDYDTADRLYFEPLTFEDVLEVYHAEMESGSGGPGVAGVIVQLGGQTPLGLAHRLADAGVPIVGTPPEAIDLAEDRGAFGDLLSAAGLPAPKYGTATTFSQARRIAEEIGYPVLVRPSYVLGGRGMEIVYDEETLQGYITRATQLSPEHPVLVDRFLEDAVEIDVDALCDGAEVYIGGIMEHIEEAGIHSGDSACALPPVTLGRSDIAKVRKATEAIAHGIGVVGLLNVQYALKDDVLYVLEANPRASRTVPFVSKATAVPLAKACARIMLGATIAQLRAEGLLAVTGDGAHAARNAPIAVKEAVLPFHRFRRADGAAIDSLLGPEMKSTGEVMGIDRDFGSAFAKSQTAAYGSLPAQGTVFVSVANRDKRSLVFPVKRLADLGFRVLATEGTAEMLRRNGIPCDDVRKHFEPAQPGRPTMSAVDAIRAGEVNMVINTPYGNSGPRIDGYEIRSAAVAGNIPCITTVQGASAAVQGIEAGIRGDIGVRSLQELHRVIGGVER'
    def tearDown(self):
        del self.Trie
        del self.protein

    def test_trie_upper(self):
        TrieMatch=algo.TrieMatch(self.Trie, self.protein) 
        trie_upper = TrieMatch.trie_upper()

    def test_trie_postions(self):
        pep_positions = []
        TrieMatch=algo.TrieMatch(self.Trie, self.protein) 
        trie_positions = TrieMatch.trie_matching()

        for pep in self.peptides:
            pos=self.protein.find(pep)
            
            pos = [m.start() for m in re.finditer('(?={})'.format(pep), self.protein)] 
            pep_positions += pos
        self.assertEqual(set(pep_positions), set(trie_positions))

    def test_trie_postions_2(self):
        #trie = algo.trie_graph(peptides)
        s = 'MSFTFARFLKFIVTEIFPGGRLPSIPMVQECASANGFTVTRVQSLQPHYAKTLDLWSAALQANKGQAIALQSEEVYERYMKYLTGCAEMFRIGYIDVNQFTCQK'
        s_peps=['FIVTEIFPGGR', 'GQAIALQSEEVYER', 'IGYIDVNQFTCQK', 'TLDLWSAALQANK']

        TrieMatch=algo.TrieMatch(self.Trie, s) 
        trie_positions = TrieMatch.trie_matching()
        
        pep_positions = []

        for pep in s_peps:
            if pep in s:
                pos = [m.start() for m in re.finditer('(?={})'.format(pep), s)] 
                pep_positions += pos
        self.assertEqual(set(pep_positions), set(trie_positions))

def main():
    unittest.main()

if __name__ == '__main__':
    main()


