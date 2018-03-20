#!/usr/bin/env python

import unittest
import pandas as pd
import algo
import re

#rosalind.info/problems/fibo
class fibonacci(unittest.TestCase): 
    def setUp(self):
        pass
    def tearDown(self):
        pass
    def test_fibonacci(self):
        self.assertEqual(algo.fibonacci(6, a=0, b=1), 8)
        self.assertEqual(algo.fibonacci(20, a=0, b=1), 6765)

class Vertex(unittest.TestCase):
    def setUP(self):
        pass
    def tearDown(self):
        pass
    def test_Vertex(self):
        vt = algo.Vertex(0)
        vt.addNeighbour(1,'a')
        self.assertEqual(list(vt.connectedTo.keys()), [1])
        self.assertEqual(vt.getId(), 0)
        self.assertEqual(vt.getEdge(1), 'a')
        self.assertEqual(list(vt.getConnections()), [1])

class Graph(unittest.TestCase):
    def setUP(self):
        pass
    def tearDown(self):
        pass
    def test_Graph(self):
        g = algo.Graph()
        
        g.addVertex(1)
        self.assertEqual(list(g.getVertices()), [1])
        
        g.addEdge(2,5,'b')
        self.assertEqual(g.vertices[2].connectedTo[5], 'b')

        self.assertEqual(g.__contains__(2), True)
        self.assertEqual(list(g.getVertices()), [1,2,5])

class Trie(unittest.TestCase):
    def setUP(self):
        pass
    def tearDown(self):
        pass
    def test_Trie(self):
        seqs = ['ATAGA','ATC','GAT']
        Trie = algo.Trie(seqs)
        vertices = Trie.trie.getVertices()
        self.assertEqual(list(vertices), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    def test_prefix_trie_match(self):
        string='AATCGGGTTCAATCGGGGT'
        words=['ATCG', 'GGGT']
        Trie = algo.Trie(words)
        TrieMatch=algo.TrieMatch(Trie, string)
        matches = TrieMatch.trie_matching()
        self.assertEqual(matches, [1,4,11,15])
    def test_trie_upper(self):
        string='AATCGGGTTCAATCGGGGTBBB'
        words=['ATCG', 'GGGT', 'BB']
        Trie = algo.Trie(words)
        TrieMatch=algo.TrieMatch(Trie, string)
        upper = TrieMatch.trie_upper()
        self.assertEqual(upper, 'aATCGGGTtcaATCGGGGTBBB')

    def test_trie_example(self):
        f = open('testdata/TrieMatching.txt')
        data = f.read()
        f.close()
        data=data.split('Input')[1].split('Output')
        d1 = data[0].split('\n')
        string = d1[1]
        words = [i.strip() for i in d1[2:] if i != '']
        positions = [int(i) for i in data[1].split()]
        Trie = algo.Trie(words)
        TrieMatch=algo.TrieMatch(Trie, string)
        matches = TrieMatch.trie_matching()
        self.assertEqual(matches, positions)
    def test_trie_example(self):
        f = open('testdata/rosalind_ba9b.txt')
        data = f.readlines()
        f.close()
        string = data[0].rstrip()
        words = [i.strip() for i in data[1:] if i != '']
        Trie = algo.Trie(words)
        TrieMatch=algo.TrieMatch(Trie, string)
        matches = TrieMatch.trie_matching()
        self.assertEqual(matches,[85, 143, 204, 211, 254, 261, 297, 321, 328, 443, 459, 466, 495, 533, 613, 626, 633, 672, 704, 711,742, 780, 827, 879])
    def test_trie_upper_and_coverage(self):
        peptides=['AA','ADGAAIDSLLGPEMK', 'AEGLLAVTGDGAHAAR','AEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVER', 'AGEVNMVINTPYGNSGPR', 'AGFWTAPDPDGGIEEALTR', 'ATQLSPEHPVLVDR', 'DDVLYVLEANPR', 'DLGIAILR', 'ETVAE    LGLPVVVRPSFTMGGLGSGIAYSTDEVDR', 'EVGVDTGGCNIQFAVNPR','FAFEKFPGADPTLTTTMK','FPGADPTLTTTMK','GAFGDLLSAAGLPAPK','HFEPAQPGRPTMSAVDAIR','HSGLSDHQIASLRPELAGEAGVR','IAEEIGYPVLVR','IAEEIGYPVLVRPSYVLGGR','IMLGATIAQLR','LADAGVPIVGTPPEAIDLAEDR','LADLGFR','LGIHPVYK','LIVIEMNPR','LYDIELALR','MAGAGLAASPSANVLIEESIYGWK','NELVAAPVLNAELLR','NFVEALGK','RADGAAIDSLLGPEMK','SLVFPVK','SLVFPVKR','SQTAAYGSLPAQGTVFVSVANR',    'SVGEAMSLGR','VAEASGVDPWFIAQINELVNLR','VLATEGTAEMLR','VLILGSGPNR','YGTATTFSQAR','YGVELIGADFDAIQR']
        protein='AAAVPRRTDLHHVLVIGSGPIVIGQACEFDYSGTQACRVLRAEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVERVIAQQAERGNKIDALLATLGGQTALNTAVALYESGVLEKYGVELIGADFDAIQRGEDRQRFKDIVAKAGGESARSRVCFTMAEVRETVAELGLPVVVRPSFTMGGLGSGIAYSTDEVDRMAGAGLAASPSANVLIEESIYGWKEFELELMRDGHDNVVVVCSIENVDPMGVHTGDSVTVAPAMTLTDREYQRMRDLGIAILREVGVDTGGCNIQFAVNPRDGRLIVIEMNPRVSRSSALASKATGFPIAKIAAKLAIGYTLDEIVNDITGETPACFEPTLDYVVVKAPRFAFEKFPGADPTLTTTMKSVGEAMSLGRNFVEALGKVMRSLETTRAGFWTAPDPDGGIEEALTRLRTPAEGRLYDIELALRLGATVERVAEASGVDPWFIAQINELVNLRNELVAAPVLNAELLRRAKHSGLSDHQIASLRPELAGEAGVRSLRVRLGIHPVYKTVDTCAAEFEAQTPYHYSSYELDPAAETEVAPQTERPKVLILGSGPNRIGQGIEFDYSCVHAATTLSQAGFETVMVNCNPETVSTDYDTADRLYFEPLTFEDVLEVYHAEMESGSGGPGVAGVIVQLGGQTPLGLAHRLADAGVPIVGTPPEAIDLAEDRGAFGDLLSAAGLPAPKYGTATTFSQARRIAEEIGYPVLVRPSYVLGGRGMEIVYDEETLQGYITRATQLSPEHPVLVDRFLEDAVEIDVDALCDGAEVYIGGIMEHIEEAGIHSGDSACALPPVTLGRSDIAKVRKATEAIAHGIGVVGLLNVQYALKDDVLYVLEANPRASRTVPFVSKATAVPLAKACARIMLGATIAQLRAEGLLAVTGDGAHAARNAPIAVKEAVLPFHRFRRADGAAIDSLLGPEMKSTGEVMGIDRDFGSAFAKSQTAAYGSLPAQGTVFVSVANRDKRSLVFPVKRLADLGFRVLATEGTAEMLRRNGIPCDDVRKHFEPAQPGRPTMSAVDAIRAGEVNMVINTPYGNSGPRIDGYEIRSAAVAGNIPCITTVQGASAAVQGIEAGIRGDIGVRSLQELHRVIGGVER'
        Trie = algo.Trie(peptides)
        TrieMatch=algo.TrieMatch(Trie, protein)
        trie_upper = TrieMatch.trie_upper()
        indices =set()
        prot_indices = set()
        for p in peptides:
            positions = [m.start() for m in re.finditer('(?={})'.format(p), protein)]
            for position in positions:
                for i in range(len(p)):
                    pos = i + position
                    indices.add(pos)            
        new_sequence_list = []
        for i in range(len(protein)):
            prot_indices.add(i)
            symbol = protein[i]
            if i not in indices:
                new_sequence_list.append(symbol.lower())
            else:
                new_sequence_list.append(symbol)
        new_sequence = ''.join(new_sequence_list)
        self.assertEqual(new_sequence, trie_upper)
    
        trie_coverage = TrieMatch.trie_coverage()
    
        coverage = len(indices)/float(len(prot_indices)) * 100
        self.assertEqual(trie_coverage, coverage)

    def test_trie_export(self):
        peptides=['AA', 'ADGAAIDSLLG' ,'ADGAAIDSLLGPEMK','AIDSLLGPEMK', 'AEGLLAVTGDGAHAAR','AEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVER', 'AGEVNMVINTPYGNSGPR', 'AGFWTAPDPDGGIEEALTR', 'ATQLSPEHPVLVDR', 'DDVLYVLEANPR', 'DLGIAILR', 'ETVAE    LGLPVVVRPSFTMGGLGSGIAYSTDEVDR', 'EVGVDTGGCNIQFAVNPR','FAFEKFPGADPTLTTTMK','FPGADPTLTTTMK','GAFGDLLSAAGLPAPK','HFEPAQPGRPTMSAVDAIR','HSGLSDHQIASLRPELAGEAGVR','IAEEIGYPVLVR','IAEEIGYPVLVRPSYVLGGR','IMLGATIAQLR','LADAGVPIVGTPPEAIDLAEDR','LADLGFR','LGIHPVYK','LIVIEMNPR','LYDIELALR','MAGAGLAASPSANVLIEESIYGWK','NELVAAPVLNAELLR','NFVEALGK','RADGAAIDSLLGPEMK','SLVFPVK','SLVFPVKR','SQTAAYGSLPAQGTVFVSVANR',    'SVGEAMSLGR','VAEASGVDPWFIAQINELVNLR','VLATEGTAEMLR','VLILGSGPNR','YGTATTFSQAR','YGVELIGADFDAIQR']
        protein='AAAVPRRTDLHHVLVIGSGPIVIGQACEFDYSGTQACRVLRAEGLQVSLVNSNPATIMTDPEFADHTYVEPITPAFVERVIAQQAERGNKIDALLATLGGQTALNTAVALYESGVLEKYGVELIGADFDAIQRGEDRQRFKDIVAKAGGESARSRVCFTMAEVRETVAELGLPVVVRPSFTMGGLGSGIAYSTDEVDRMAGAGLAASPSANVLIEESIYGWKEFELELMRDGHDNVVVVCSIENVDPMGVHTGDSVTVAPAMTLTDREYQRMRDLGIAILREVGVDTGGCNIQFAVNPRDGRLIVIEMNPRVSRSSALASKATGFPIAKIAAKLAIGYTLDEIVNDITGETPACFEPTLDYVVVKAPRFAFEKFPGADPTLTTTMKSVGEAMSLGRNFVEALGKVMRSLETTRAGFWTAPDPDGGIEEALTRLRTPAEGRLYDIELALRLGATVERVAEASGVDPWFIAQINELVNLRNELVAAPVLNAELLRRAKHSGLSDHQIASLRPELAGEAGVRSLRVRLGIHPVYKTVDTCAAEFEAQTPYHYSSYELDPAAETEVAPQTERPKVLILGSGPNRIGQGIEFDYSCVHAATTLSQAGFETVMVNCNPETVSTDYDTADRLYFEPLTFEDVLEVYHAEMESGSGGPGVAGVIVQLGGQTPLGLAHRLADAGVPIVGTPPEAIDLAEDRGAFGDLLSAAGLPAPKYGTATTFSQARRIAEEIGYPVLVRPSYVLGGRGMEIVYDEETLQGYITRATQLSPEHPVLVDRFLEDAVEIDVDALCDGAEVYIGGIMEHIEEAGIHSGDSACALPPVTLGRSDIAKVRKATEAIAHGIGVVGLLNVQYALKDDVLYVLEANPRASRTVPFVSKATAVPLAKACARIMLGATIAQLRAEGLLAVTGDGAHAARNAPIAVKEAVLPFHRFRRADGAAIDSLLGPEMKSTGEVMGIDRDFGSAFAKSQTAAYGSLPAQGTVFVSVANRDKRSLVFPVKRLADLGFRVLATEGTAEMLRRNGIPCDDVRKHFEPAQPGRPTMSAVDAIRAGEVNMVINTPYGNSGPRIDGYEIRSAAVAGNIPCITTVQGASAAVQGIEAGIRGDIGVRSLQELHRVIGGVER'
        matched_peptides = []
        for i in peptides:
            if i in protein:
                matched_peptides.append(i)
        Trie = algo.Trie(peptides)
        TrieMatch=algo.TrieMatch(Trie, protein)
        words = TrieMatch.trie_export()
        self.assertEqual(words, set(matched_peptides))
        


