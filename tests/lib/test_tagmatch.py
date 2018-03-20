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
import tagmatch
from tagmatch import missed_cleavages
from tagmatch import valid_cleavage

class peptide_mass(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_peptide_mass(self):
        query='TSEYPCMQVLR'
        mass = np.round(tagmatch.peptide_mass(query, fixed_modifications=[]), 3)

        self.assertEqual(mass, 1325.611)

        query='EPIRQLPCPVDDLRTSEYPCMQVLRASLSDHDNQEQNFDNLFHNYGYLYPLWKQLHDGLIIASRWHGDIYQGHIDICSHFQCYSHDLLMRRYMCGIGDGKSLVWKLNGWRYFQGLAAHTCGMTRDMNQREHNRLNPNHNPYVDGQWQCVYVKLEGMCRKTPLKQANVVISMISIGNIQIKTAAEGHHNGAFDHIWASACVGMNNMIQKQEVGNGIPAVYFGSAVIQSAIVLMCCTAPYEAHPLACDQCIHEQTYVMHFNDEVWYISHAAIATTDTCSSDMEALEKMNENHPPCCFAERQTYGNYNDNMCEIEHMCKQIIKNQMPWIISCLNIKLWQWMFQYLRTMTKFVMIERCNQNVHQTSGDILTYYGDIGMDKGKNKKNELVIHYKCCNMWVQQFGENEMAMHEFLRLYWWSWKYWGTFDYADNMRQYQACLEWVCYYIQKWNDLMQRRLFHLINHHWRLTHSNTGRCRRGWMPNCCGTWIIRVRDIGVEQAHGIKNMNPRGFSGPHACASTKKHLAWPYKDKLSDNGVEMTEFSMMDLRKCNITMWNLNMADSVTFGWHWEGAYQPPECERVQFHSRTTYIWCPTHCKQLDGTGKELGISFTYWRLDELILKTFEFPEFFFSCIMPRFFMYLDNKYRAGSSHPFCPNVIHGIFMQSFDYDSCLPTEAKNNISESLPCYTYTSWKAEYNQIPETFSHCAEPYSCETREICCDYMMKEQIWSHHNTCHMESAMWDMHQGGHSYKRIWESEIDIDFTAVWKRGNAVAPEVCKAKGAEDQPVAMNEFKSEAFPPWYRAEKHGNYIVFRTTDEFFGNAPVWISRVSHFPNVENFTRIMMWSNGRLVQPYDNEPLSVFRASYSKCCDIKSYSRHFYILIRNTMPCQYPMEFGVKQFWAINEVFDAVNLICAFMEHAEEM'
        mass = np.round(tagmatch.peptide_mass(query, fixed_modifications=[]), 3)
        self.assertEqual(mass, 107877.398)

class TagMatch(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_positions(self):
        query='LVLFLPLELMR'
        query_mz=672.409
        query_charge=2
        query_mw=tagmatch.mz2mw(query_mz, query_charge)
        target='KLVLFLPLELMRY'
        tm = tagmatch.TagMatch(query, [query_mw], target)
        for coord in tm.query_positions:
            self.assertEqual(query, target[coord[0]:coord[1]])
        self.assertEqual(tm.query_positions, [(1,12)])
    def test_positions(self):
        query='LVLFLPLELMR'
        query_mz=672.409
        query_charge=3
        query_mw=tagmatch.mz2mw(query_mz, query_charge)
        #print(query_mw)
        target='KLVLFLPLELMRY'
        
        tm = tagmatch.TagMatch(query, [(128.09496 , query_mw, 163.06333)], target, max_missed_cleavages=3)
        for coord in tm.query_positions:
            self.assertEqual(query, target[coord[0]:coord[1]])
        self.assertEqual(tm.query_positions, [(1,12)])
        self.assertEqual(tm.validated_peptides, []) # move to own function

    def test_get_tryptic_peptide(self):
        query='FLPLE'
        query_mz=672.409
        query_charge=2
        query_mw=tagmatch.mz2mw(query_mz, query_charge)
        target='VFGKLVLFLPLELMRCAALP'
        ngap=325.236538
        cgap=400.225665

        precursor = (ngap,query_mw,cgap)
        tm = tagmatch.TagMatch(query, [precursor], target, max_missed_cleavages=0)
        self.assertEqual(list(tm.peptides.values()), ['LVLFLPLELMR']) 
        self.assertEqual(tm.validated_peptides, ['LVLFLPLELMR']) 
        self.assertEqual(tm.validated, True)
        self.assertEqual(tm.validated_precursor, [precursor])

        precursor = (ngap + 100,query_mw,cgap)
        tm = tagmatch.TagMatch(query, [precursor], target)
        #self.assertEqual(list(tm.peptides.values()), ['LVLFLPLELMR']) 
        self.assertEqual(tm.validated_peptides, []) 
        self.assertEqual(tm.validated, False)
        self.assertEqual(tm.validated_precursor, [])

        precursor = (ngap,query_mw,cgap + 100)
        tm = tagmatch.TagMatch(query, [precursor], target)
        #self.assertEqual(list(tm.peptides.values()), ['LVLFLPLELMR']) 
        self.assertEqual(tm.validated_peptides, []) 
        self.assertEqual(tm.validated, False)
        self.assertEqual(tm.validated_precursor, [])
        
        precursor = (ngap,query_mw + 100,cgap)
        tm = tagmatch.TagMatch(query, [precursor], target)
        #self.assertEqual(list(tm.peptides.values()), ['LVLFLPLELMR']) 
        self.assertEqual(tm.validated_peptides, []) 
        self.assertEqual(tm.validated, False)
        self.assertEqual(tm.validated_precursor, [])

        target='FLPLELMRCAALP'
        tm = tagmatch.TagMatch(query, [(1000,query_mw,1000)], target)
        #self.assertEqual(set(tm.peptides.values()), set(['FLPLELMR', 'FLPLELMRCAALP']))
        self.assertEqual(list(tm.validated_peptides), []) 
        
        target='DFFLPLEL'
        tm = tagmatch.TagMatch(query, [(1000,query_mw,1000)], target)
        #self.assertEqual(list(tm.peptides.values()), ['DFFLPLEL']) 
        
        target='FLPLE'
        tm = tagmatch.TagMatch(query, [(1000,query_mw,1000)], target)
        #self.assertEqual(list(tm.peptides.values()), ['FLPLE']) 

        target='FFLPLEFLPLELPLE'
        tm = tagmatch.TagMatch(query, [(1000,query_mw,1000)], target)
        #self.assertEqual(list(tm.peptides.values()), ['FFLPLEFLPLELPLE'])

    def test_tag_precursor_against_own_sequence(self):
        query='TSIA'
        query_mz=672.409
        query_charge=2
        query_mw=tagmatch.mz2mw(query_mz, query_charge)
        target='TSIA'
        tm = tagmatch.TagMatch(query, [(0, query_mw, 0)], target)
        self.assertEqual(list(tm.peptides.values()), ['TSIA']) 
        self.assertEqual(tm.validated_peptides, [])
        self.assertEqual(tm.validated, False)
        self.assertEqual(tm.validated_precursor, [])
        
        query='SSFSGTELFGK'
        query_mz=580.286
        query_charge=2
        query_mw=tagmatch.mz2mw(query_mz, query_charge)
        target='SSFSGTELFGK'
        tm = tagmatch.TagMatch(query, [(0,query_mw,0)], target)
        
        self.assertEqual(list(tm.peptides.values()), [ 'SSFSGTELFGK']) 
        self.assertEqual(tm.validated_peptides, ['SSFSGTELFGK'])
        self.assertEqual(tm.validated, True)
        self.assertEqual(tm.validated_precursor, [(0,query_mw,0)])
    
    def test_tag_against_protein(self):
        query='ISPPLPPAPPS' 
        tag_mass_list = [(57.02146, 1189.6706738937305, 61.052766)] 
        target='MPPWPIRPADPPAPPAPPRPFGLVTVPPAPPVPPLPMSCPPAPPAPPAPPPWPRISPPLPPAPPSPISQALPPAPPAPLIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSYSPRLCLQPQFAAAPAGPAAARQAAGRRRRCAGPAVAS' 
        prec_tol=0.02 
        specificity='semi-specific'

        tm = tagmatch.TagMatch(query, tag_mass_list, target=target, prec_tol=prec_tol, specificity = specificity)
        
        self.assertEqual(tm.validated_peptides, [])

    def test_tag_against_peptide_example_1(self):
        query='TLMEYLE' 
        tag_mass_list = [(258.1, 1269.5434917, 114.0)] 
        target='PPEETLMEYLENPK'
        prec_tol=0.02 
        specificity='unspecific'

        tm = tagmatch.TagMatch(query, tag_mass_list, target=target, prec_tol=prec_tol, specificity = specificity)
        
        self.assertEqual(tm.validated_peptides, ['EETLMEYLEN'])

    def test_tag_against_peptide_example_2(self):
        query='IPDT' 
        tag_mass_list=[(928.822, 1572.1766276365938, 199.132), 
                (928.822, 1572.176787132922, 199.13299999999998), 
                (1697.89, 2341.2397817277724, 199.132), 
                (470.27, 1967.016318972657, 1052.52), 
                (469.274, 1966.013737230306, 1052.52), 
                (928.822, 1572.1766574494618, 199.13299999999998), 
                (0.0, 2257.9430767851572, 1884.01), 
                (1130.54, 1773.8975445585932, 199.13299999999998)] 
        target='MGAGSGVGNVGTLGSGVLNLGSGISGFYNTSVLPFGTPAAVSGIGNLGQQLSGVSAAGTTLRSMLAGNLGLANVGNFNTGFGNVGDVNLGAANIGGHNLGLGNVGDGNLGLGNIGHGNLGFANLGLTAGAAGVGNVGFGNAGINNYGLANMGVGNIGFANTGTGNIGIGLVGDHRTGIGGLNSGIGNIGLFNSGTGNVGFFNSGTGNFGIGNSGRFNTGIGNSGTASTGLFNAGSFSTGIANTGDYNTGSFNAGDTNTGGFNPGGINTGWFNTGHANTGLANAGTFGTGAFMTGDYSNGLLWRGGYEGLVGVRVGPTISQFPVTVHAIGGVGPLHVAPVPVPAVHVEITDATVGLGPFTVPPISIPSLPIASITGSVDLAANTISPIRALDPLAGSIGLFLEPFRLSDPFITIDAFQVVAGVLFLENIIVPGLTVSGQILVTPTPIPLTLNLDTTPWTLFPNGFTIPAQTPVTVGMEVANDGFTFFPGGLTFPRASAGVTGLSVGLDAFTLLPDGFTLDTVPATFDGTILIGDIPIPIIDVPAVPGFGNTTTAPSSGFFNTGGGGGSGFANVGAGTSGWWNQGHDVLAGAGSGVANAGTLSSGVLNVGSGISGWYNTSTLGAGTPAVVSGIGNLGQQLSGFLANGTVLNRSPIVNIGWADVGAFNTGLGNVGDLNWGAANIGAQNLGLGNLGSGNVGFGNIGAGNVGFANSGPAVGLAGLGNVGLSNAGSNNWGLANLGVGNIGLANTGTGNIGIGLVGDYQTGIGGLNSGSGNIGLFNSGTGNVGFFNTGTGNFGLFNSGSFNTGIGNSGTGSTGLFNAGNFNTGIANPGSYNTGSFNVGDTNTGGFNPGDINTGWFNTGIMNTGTRNTGALMSGTDSNGMLWRGDHEGLFGLSYGITIPQFPIRITTTGGIGPIVIPDTTILPPLHLQITGDADYSFTVPDIPIPAIHIGINGVVTVGFTAPEATLLSALKNNGSFISFGPITRSNIDIPPMDFTLGLPVLGPITGQLGPIHLEPIVVAGIGVPLEIEPIPLDAISLSESIPIRIPVDIPASVIDGISMSEVVPIDASVDIPAVTITGTTISAIPLGFDIRTSAGPLNIPIIDIPAAPGFGNSTQMPSSGFFNTGAGGGSGIGNLGAGVSGLLNQAGAGSLVGTLSGLGNAGTLASGVLNSGTAISGLFNVSTLDATTPAVISGFSNLGDHMSGVSIDGLIAILTFPPAESVFDQIIDAAIAELQHLDIGNALALGNVGGVNLGLANAGEFNLGAGNVGNINVGAGNLGGSNLGLGNVGTGNLGFGNIGAGNFGFGNAGLTAGAGAWAMWGWVTPAAAAGGWPTWVWAISGWPTPAPATSGSG'

        prec_tol=0.02 
        specificity='semi-specific'

        tm = tagmatch.TagMatch(query, tag_mass_list, target, prec_tol=prec_tol, specificity = specificity)  
        self.assertEqual(tm.validated_peptides, [])
    def test_tag_against_peptide_example_3(self):
        query='TLMEYLE' 
        tag_mass_list = [(258.1, 1155.5005647, 0)] 
        target='PPEETLMEYLE'
        prec_tol=0.02 
        specificity='semi-specific'

        tm = tagmatch.TagMatch(query, tag_mass_list, target=target, prec_tol=prec_tol, specificity = specificity)
        
        self.assertEqual(tm.validated_peptides, ['EETLMEYLE'])

    def test_tag_unspecific(self):
        query='SFSGTELF'
        query_mw=943.4287177
        target='KSSFSGTELFGK'
        tm = tagmatch.TagMatch(query, [(0,query_mw,57.021464)], target, specificity='unspecific')
        
        #self.assertEqual(set(tm.peptides.values()),{'SFSGTELFGK','SSFSGTELF', 'SFSGTELF','SSFSGTELFGK', 'SSFSGTELFG', 'SFSGTELFG'}) 
        self.assertEqual(tm.validated_peptides, ['SFSGTELFG'])
        self.assertEqual(tm.validated, True)
        self.assertEqual(tm.validated_precursor, [(0,query_mw,57.021464)])
    
    def test_tag_peptide_carbamido(self):
        #NH2-FC<cmm>LRRLKARRPR-COOH
        query='FCLRRLKARRPR'
        query_mw = 1627.968151
        fixed=['Carbamidomethylation of C']
        var = []
        target='FCLRRLKARRPR'
        tm = tagmatch.TagMatch(query, [(0,query_mw,0)], target, enzymes=['Trypsin'], specificity='specific', max_missed_cleavages=10)
        self.assertEqual(tm.validated_peptides, [query])
       
        # In case mass gaps are not strictly residue mass
        tm = tagmatch.TagMatch(query, [(18.01056,query_mw,18.01056)], target, enzymes=['Trypsin'], specificity='specific', max_missed_cleavages=10, gap_tol=18.2)
        self.assertEqual(tm.validated_peptides, [query])

        # artificial tag
        tm = tagmatch.TagMatch(query[1:-1], [(147.068414,query_mw,156.101111)], target, enzymes=['Trypsin'], specificity='specific', max_missed_cleavages=10)
        self.assertEqual(tm.validated_peptides, [query])

    def test_enzymes(self):
        query='SFSGTELF'
        query_mw=943.4287177
        target='KSSFSGTELFGK'
        tag_mass_list=[(0,query_mw,57.021464)] 
        
        ######################
        # Trypsin, no P rule #
        ######################
        # Trypsin, no P rule - specific
        enzymes=['Trypsin, no P rule']
        specificity='specific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(tagmatch.valid_cleavage('K', 'P', 'K', '', tm.enzymes, tm.specificity),   True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'L', 'K', '', tm.enzymes, tm.specificity),   True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'L', 'L', '', tm.enzymes, tm.specificity),   True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'L', 'L', 'K', tm.enzymes, tm.specificity),  False )
        # Trypsin, no P rule - semi-specific
        enzymes=['Trypsin, no P rule']
        specificity='semi-specific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(tagmatch.valid_cleavage('K', 'P', 'K', '' , tm.enzymes, tm.specificity),  True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'L', 'K', '' , tm.enzymes, tm.specificity),  True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'L', 'L', '' , tm.enzymes,  tm.specificity),  True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'L', 'L', 'K', tm.enzymes, tm.specificity),  True  )
        self.assertEqual(tagmatch.valid_cleavage('K', 'P', 'L', 'L', tm.enzymes, tm.specificity),  True  )
        self.assertEqual(tagmatch.valid_cleavage('L', 'L', 'L', 'L', tm.enzymes, tm.specificity),  False )
        
        # Trypsin, no P rule - missed cleavages
        enzymes=['Trypsin, no P rule']
        specificity='specific'
        peptide = 'PRKPLP'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(missed_cleavages(peptide, tm.enzymes), 2)

        ###########
        # Trypsin #
        ###########
        # Trypsin - specific
        enzymes=['Trypsin']
        specificity='specific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(valid_cleavage('K', 'P', 'K', '', tm.enzymes, tm.specificity ),  False )
        self.assertEqual(valid_cleavage('K', 'L', 'K', '' , tm.enzymes, tm.specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', '' , tm.enzymes, tm.specificity ),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', 'K' , tm.enzymes, tm.specificity),  False )
        self.assertEqual(valid_cleavage('', 'T', 'R', 'L' , tm.enzymes, tm.specificity),  True )
        
        # Trypsin - semi-specific
        enzymes=['Trypsin']
        specificity='semi-specific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(valid_cleavage('K', 'P', 'K', '', tm.enzymes, tm.specificity ),  True  )
        self.assertEqual(valid_cleavage('K', 'P', 'L', 'L', tm.enzymes, tm.specificity), False  )
        self.assertEqual(valid_cleavage('K', 'L', 'K', '' , tm.enzymes, tm.specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', '' , tm.enzymes, tm.specificity ),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', 'K', tm.enzymes, tm.specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', 'L', tm.enzymes, tm.specificity),  True  )
        self.assertEqual(valid_cleavage('L', 'L', 'L', 'L', tm.enzymes, tm.specificity),  False )
        # Trypsin - unspecific
        enzymes=['Trypsin']
        specificity='unspecific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(valid_cleavage('K', 'P', 'K', '' , enzymes, specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'P', 'L', 'L', enzymes, specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'K', '' , enzymes, specificity ),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', '' , enzymes, specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', 'K' , enzymes, specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', 'L' , enzymes, specificity),  True  )
        self.assertEqual(valid_cleavage('L', 'L', 'L', 'L', enzymes, specificity ),  True  )
        
        # Trypsin - missed cleavages
        enzymes=['Trypsin']
        specificity='specific'
        peptide = 'PRKPLP'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(missed_cleavages(peptide, enzymes), 1)
        
        peptide = 'PPLP'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(missed_cleavages(peptide, enzymes), 0)

        #################
        # Whole protein #
        #################
        # Whole protein - specific
        enzymes=['Whole protein']
        specificity='specific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(valid_cleavage('K', 'P', 'K', '' , enzymes, specificity),  False )
        self.assertEqual(valid_cleavage('', 'L', 'L', '' , enzymes, specificity ),  True  )
        self.assertEqual(valid_cleavage('', 'L', 'L', 'K', enzymes, specificity),  False )
        # Whole protein - semi-specific
        enzymes=['Whole protein']
        specificity='semi-specific'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(valid_cleavage('K', 'P', 'K', '', enzymes, specificity ),  True  )
        self.assertEqual(valid_cleavage('', 'P', 'L', 'L', enzymes, specificity),   True  )
        self.assertEqual(valid_cleavage('', 'L', 'K', '' , enzymes, specificity),  True  )
        self.assertEqual(valid_cleavage('K', 'L', 'L', 'P' , enzymes, specificity),  False  )

        
        # Whole protein - missed cleavages
        enzymes=['Whole protein']
        specificity='specific'
        peptide = 'PRKPLP'
        tm = tagmatch.TagMatch(query,tag_mass_list,target,enzymes=enzymes,specificity=specificity)
        self.assertEqual(missed_cleavages(peptide, tm.enzymes), 0)

class consensus_strip(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_character_strip(self):
        consensus='XXXXXXXXX'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'')

        consensus='XXXTVGXXXX'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'TVG')
        
        consensus='XXXTVGXGXX'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'TVGXG')
        
        consensus='XXXTVGXG'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'TVGXG')

        consensus='TVGXGX'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'TVGXG')

        consensus='TVGXG'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'TVGXG')
        
        consensus='X'
        self.assertEqual(tagmatch.character_strip(consensus,'X'),'')

class gap_sequence(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass
    def test_single_amino_gap(self):
        gap=147.068414
        gs = tagmatch.gap_sequence(gap, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02)
        gs.sequences()
        self.assertEqual(gs.gap_sequences, ['F'])

    def test_two_amino_gap(self):
        gap=278.108904
        gs = tagmatch.gap_sequence(gap, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02)
        gs.sequences()
        self.assertEqual(gs.gap_sequences, ['DY', 'FM', 'MF', 'YD'])
    
    def test_three_amino_gap(self):
        gap=285.11471400000005
        gs = tagmatch.gap_sequence(gap, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02)
        gs.sequences()
        self.assertEqual(len(gs.gap_sequences), 50)
    
    def test_300_amino_gap(self):
        gap= 300.16
        gs = tagmatch.gap_sequence(gap, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02)
        #gs.sequences()
        #self.assertEqual(len(gs.gap_sequences), 58)

    def test_valid_amino_gap(self):
        gap=278.108904
        gs = tagmatch.gap_sequence(gap, fixed_modifications=['Carbamidomethylation of C'], variable_modifications=[], prec_tol=0.02)
        gs.validate()
        self.assertEqual(gs.valid, True)
    

class blast_tag(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_blast_tags(self):
        
        query1='FECMNDNDLR'
        query2 = query1
        recs =[SeqRecord(seq=Seq(query1), id='test.mgf;scan=1;mw=1312.517;ngap=0;cgap=0')]
        #match='+DNDLR'
        subject='DDNDLR'
        
        matches = tagmatch.blast_tags(subject, recs)
        
        # Longest substring
        self.assertEqual(matches.substrings('FECMNDNDLR', 'DNDLR'), ['DNDLR', 'NDLR', 'DLR', 'LR', 'R'])
        self.assertEqual(matches.substrings('FECMNDNDLR', 'DNDLRR'), ['DNDLR', 'NDLR', 'DLR', 'LR', 'R', 'R'])
        self.assertEqual(matches.substrings('FECMNDNDLR', 'PPPPPPPPDNDLRR'), ['DNDLR', 'NDLR', 'DLR', 'LR', 'R', 'R'])
        self.assertEqual(matches.substrings('FECMNDNDLR', 'PPPPPPPPDNDLRRTDND'), ['DNDLR', 'NDLR', 'DLR', 'LR', 'R', 'R', 'DND', 'ND', 'D'])
        self.assertEqual(matches.substrings('FECMNDNDLR', 'PPPPPPPPDND'), ['DND', 'ND', 'D'])
        self.assertEqual(matches.substrings('FECMNDNDLR', 'D'), ['D'])

        # Substring coordinates
        self.assertEqual(matches.substring_coordinates('FECMNDNDLR', 'D'), [(5, 6), (7, 8)])
        self.assertEqual(matches.substring_coordinates('FECMNDNDLR', 'DN'), [(5, 7)])
        self.assertEqual(matches.substring_coordinates('FECMNDNDLR', 'ND'), [(4, 6), (6,8)])

        # cmass gaps
        self.assertEqual(matches.cmass_gap('FECMNDNDLR', [(4,6), (6,8)]),[498.255045,269.185175]) 
        
        # nmass gaps
        self.assertEqual(matches.nmass_gap('FECMNDNDLR', [(4,6), (6,8)]),[567.182146, 796.2520159999999]) 

    def test_segment_matching(self):
        tagseq  ='LSCFAV'
        rec =SeqRecord(seq=Seq(tagseq), id='test.mgf;scan=1;mw=695.3312457;ngap=0.0;cgap=0.0')

        blastseq='EACFAV' 
        matches = tagmatch.blast_tags(blastseq, [rec])
        
        # Test segment rearrangements
        self.assertEqual(matches.segment_matching('LSCFAV', 'EACFAV', IL_equivalence=False), set([('EACFAV', 0.036387000000000003, 0, 'CFAV')]))
        
        self.assertEqual(matches.segment_matching('LSCPAV', 'EACFAV', IL_equivalence=False), set([]))

        self.assertEqual(matches.segment_matching('LSCLAVK', 'EACIAVK', IL_equivalence=False), set([('EACIAVK', 0.036387000000000003, 0, 'C-AVK')]))


        self.assertEqual(matches.segment_matching('LSCLAVKK','EACIAVKK', IL_equivalence=False), set([('EACIAVKK', 0.036387000000000003, 0, 'C-AVKK')]))

        self.assertEqual(matches.segment_matching('AII','ALI', min_identities=2, IL_equivalence=False), set([('ALI', 0, 0, 'A-I')]))

        self.assertEqual(matches.segment_matching('LLI','ILI', min_identities=2, IL_equivalence=False), set([('ILI', 0, 0, 'LI')]))
        
        self.assertEqual(matches.segment_matching('LLIP','ILIP', min_identities=2, IL_equivalence=False), set([('ILIP', 0, 0, 'LIP')]))

        self.assertEqual(matches.segment_matching('LLIPI','ILIPL', min_identities=2, IL_equivalence=False), set([('ILIPL', 0, 0, 'LIP')]))
        self.assertEqual(matches.segment_matching('LSCLAVI','EACIAVL', min_identities=3, IL_equivalence=False), set([('EACIAVL', 0.036387000000000003, 0, 'C-AV')]))

        self.assertEqual(matches.segment_matching('ARPKWTPTLVMPSR','KVPQVSTPTLVEVSR', min_identities=3, IL_equivalence=False), set([('KVPQVSTPTLVEVSR', -0.025145000000000001, 0, 'P-TPTLV-SR')]))
        
        self.assertEqual(matches.segment_matching('LWI','ILI', min_identities=2, IL_equivalence=False), set())
        
        self.assertEqual(matches.segment_matching('LSCWAV', 'EACFAV', IL_equivalence=False), set())
        self.assertEqual(matches.segment_matching('KVPIISTPIIVEIIR', 'KVPLLSTPLLVELLR', IL_equivalence=False), set([('KVPLLSTPLLVELLR', 0, 0, 'KVP-STP-VE-R')]))
        self.assertEqual(matches.segment_matching('KVPIISTPIIVEIIR', 'KVPLLSTPLLVELLR', IL_equivalence=True), set([('KVPLLSTPLLVELLR', 0, 0, 'KVPLLSTPLLVELLR')]))
        self.assertEqual(matches.segment_matching('KVPWWSTPWWVEWWR', 'KVPLLSTPLLVELLR', min_identities=4, IL_equivalence=False), set())
        self.assertEqual(matches.segment_matching('KVPWWSTPWWVEWWR', 'KVPLLSTPLLVELLR', min_identities=3, IL_equivalence=False), set([('KVPLLSTPLLVELLR', 145.99049199999999, 291.98098400000003, 'STP'), ('KVPLLSTPLLVELLR', 0, 437.9714760000006, 'KVP')])) # manually validated mass gaps
        
        self.assertEqual(matches.segment_matching('LSCLAVLPK', 'EACIAVIWK', IL_equivalence=True), set([('EACIAVIWK', 0.036387000000000003, -89.02654999999999, 'CLAVL')  ]))
        
        self.assertEqual(matches.segment_matching('LSCLAVWPK', 'EACIAVIWK', IL_equivalence=True), set([('EACIAVIWK', 0.036387000000000003, -16.031303999999977, 'CLAV')  ])) # manually validated mass gaps
    
    def test_fix_negative_gaps(self):

        tagseq  ='LSCFAV'
        rec =SeqRecord(seq=Seq(tagseq), id='test.mgf;scan=1;mw=695.3312457;ngap=0.0;cgap=0.0')

        blastseq='EACFAV' 
        matches = tagmatch.blast_tags(blastseq, [rec])
        
        self.assertEqual(matches.fix_negative_gaps('PPPEACFAV', -291.15828, 0), ('EACFAV',0.0,0))
        self.assertEqual(matches.fix_negative_gaps('PPPEACFAVWWWP', -291.15828, -655.29069), ('EACFAV',0.0,0.0))
        self.assertEqual(matches.fix_negative_gaps('CFAV', 0, 0), ('CFAV', 0, 0))
        self.assertEqual(matches.fix_negative_gaps('PCFAVP', 0.05, 0.05), ('PCFAVP', 0.05, 0.05))
        self.assertEqual(matches.fix_negative_gaps('PCFAVP', -0.05, -0.05), ('PCFAVP', -0.05, -0.05))



        


if __name__ == '__main__': unittest.main()
