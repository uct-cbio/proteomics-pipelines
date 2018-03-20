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
import mqparse
import os
import shutil
import pickle
import yaml
import pandas as pd

class mq_txt(unittest.TestCase):

    def setUp(self):
        self.config = 'testdata/example_mq/config.yaml'
    
    def tearDown(self):
        pass

    def test_a_load_tables(self):
        outpath='testdata/example_mq/txt/analysis'
        if os.path.exists(outpath): 
            shutil.rmtree(outpath)
        self.mq_txt = mqparse.mq_txt(self.config)
        self.assertEqual(len(self.mq_txt.reference_peptides), 773)
        self.assertEqual(len(self.mq_txt.non_reference_peptides), 185)
        self.assertEqual(os.path.exists(outpath + '/gsea/summary.txt'), True)
        self.assertEqual(os.path.exists(outpath + '/diff/peptide_diff/summary.txt'), True)
        
    def test_b_mq_diff(self):
        exp_design = 'testdata/example_mq/txt/analysis/diff/peptide_experimental_design.R'
        infile='testdata/example_mq/txt/analysis/unipept/pept2lca_taxon_sc.csv'
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        mqparse.mq_txt.diff(None, exp_design, infile, outpath)

    
    def test_b_human_variant_peptides(self):
        pass
        #design='testdata/example_mq/txt/analysis/diff/peptide_experimental_design.R'
        #infile='testdata/example_mq/txt/analysis/diff/peptide_normalization/msnbase/normalized.csv'
        #outpath='testdata/example_mq/txt/analysis/testing'
   
    def test_b_differential_peptides(self):

        infile='testdata/example_mq/txt/analysis/unipept/pept2lca_peptides.csv'
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath + '/peptide_experimental_design.txt'
        quant='Intensity'
        with open(self.config) as f:
            config = yaml.load(f.read())
        mqparse.mq_txt.create_R_parameters(None, config, outfile, quant, 'LEVEL_1')
        design=outfile
        mqparse.mq_txt.diff(None, design, infile, outpath)

        # Check that the files exist - piecharts
        self.assertEqual(os.path.exists(outpath +'/summed_intensity_all.jpeg'), True)
        self.assertEqual(os.path.exists(outpath +'/summed_intensity_After.jpeg'), True)
        self.assertEqual(os.path.exists(outpath +'/summed_intensity_Before.jpeg'), True)

        # Group variance comparisons 
        self.assertEqual(os.path.exists(outpath +'/group_variance'), True)
        self.assertEqual(os.path.exists(outpath +'/group_variance/kw.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/group_variance/group_variance.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/group_variance/dunn_bh.csv'), True)

        # Group limma comparisons
        self.assertEqual(os.path.exists(outpath +'/limma'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/limma_After-Before_intensity.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/limma_After-Before_intensity_qval_0.05.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/limma_After-Before_intensity_qval_0.05_list.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/volcano_After-Before_intensity.jpeg'), True)
        
        # SPLOM
        self.assertEqual(os.path.exists(outpath +'/splom'), True)
        self.assertEqual(os.path.exists(outpath +'/splom/splom.png'), True)
        
        # Hierarchical clustering of replicates
        self.assertEqual(os.path.exists(outpath +'/replicate_hclust'), True)
        self.assertEqual(os.path.exists(outpath +'/replicate_hclust/replicate_cluster_ward.D2.jpeg'), True)

        # Heatmaps
        self.assertEqual(os.path.exists(outpath +'/heatmaps'), True)
        self.assertEqual(os.path.exists(outpath +'/heatmaps/heatmap_intensity_all.jpeg'), True)
        
        # PCA
        self.assertEqual(os.path.exists(outpath +'/PCA'), True)
        self.assertEqual(os.path.exists(outpath +'/PCA/all_identified_pca.png'), True)
        self.assertEqual(os.path.exists(outpath +'/PCA/all_identified_pca_labelled.png'), True)
    
    def test_summarize_diff(self):

        diff='testdata/example_mq/txt/analysis/diff/peptide_diff'
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath +'/summary.txt'
        mqparse.mq_txt.summarize_diff(None, diff, outfile)
    
    def test_summarize_gsea(self):
        gsea='testdata/example_mq/txt/analysis/gsea/LEVEL_1/'
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath +'/summary.txt'
        with open(self.config) as f:
            config = yaml.load(f.read())
        mqparse.mq_txt.summarize_gsea(None, gsea, outfile, config, 'LEVEL_1')

    def test_b_peptide_normalization(self):
        
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath + '/peptide_experimental_design.txt'
        quant='Intensity'
        with open(self.config) as f:
            config = yaml.load(f.read())
        mqparse.mq_txt.create_R_parameters(None, config, outfile, quant, 'LEVEL_1')
        infile = 'testdata/example_mq/txt/analysis/peptides/target_peptides.txt'
        design = outfile
        mqparse.mq_txt.normalize(None, 'Intensity.', infile, outpath)
        self.assertEqual(os.path.exists(outpath + '/msnbase/normalized.csv'), True)

    def test_b_peptide_experimental_design(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath + '/peptide_experimental_design.txt'
        quant='Intensity'
        with open(self.config) as f:
            config = yaml.load(f.read())
        mqparse.mq_txt.create_R_parameters(None, config, outfile, quant, 'LEVEL_1')
        self.assertEqual(os.path.exists(outfile), True)


    def test_b_protein_experimental_design(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath + '/protein_experimental_design.txt'
        quant='iBAQ'
        with open(self.config) as f:
            config = yaml.load(f.read())
        mqparse.mq_txt.create_R_parameters(None, config, outfile, quant, 'LEVEL_1')
        self.assertEqual(os.path.exists(outfile), True)
    
    def test_b_protein_normalization(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        infile = 'testdata/example_mq/txt/analysis/proteins/target_proteins.txt'
        design = 'testdata/example_mq/txt/analysis/diff/protein_experimental_design_LEVEL_1.R'
        mqparse.mq_txt.normalize(None, 'iBAQ.', infile, outpath)
        self.assertEqual(os.path.exists(outpath + '/msnbase/normalized.csv'), True)

    def test_b_proteingroup_identifier(self):
        proteins = pd.read_csv('testdata/example_mq/txt/proteinGroups.txt',sep='\t')
        proteins = mqparse.mq_txt.create_protein_group_identifier(None, proteins)
        self.assertIn('Identifier', proteins.columns.tolist())
    
    def test_b_differential_proteins(self):
        design='testdata/example_mq/txt/analysis/diff/protein_experimental_design.R'
        infile='testdata/example_mq/txt/analysis/diff/protein_normalization/msnbase/normalized.csv'
        outpath='testdata/example_mq/txt/analysis/testing'
        
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        
        mqparse.mq_txt.diff(None, design, infile, outpath)

        # Check that the files exist - piecharts
        self.assertEqual(os.path.exists(outpath +'/summed_intensity_all.jpeg'), True)
        self.assertEqual(os.path.exists(outpath +'/summed_intensity_After.jpeg'), True)
        self.assertEqual(os.path.exists(outpath +'/summed_intensity_Before.jpeg'), True)

        # Group variance comparisons 
        self.assertEqual(os.path.exists(outpath +'/group_variance'), True)
        self.assertEqual(os.path.exists(outpath +'/group_variance/kw.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/group_variance/group_variance.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/group_variance/dunn_bh.csv'), True)

        # Group limma comparisons
        self.assertEqual(os.path.exists(outpath +'/limma'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/limma_After-Before_intensity.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/limma_After-Before_intensity_qval_0.05.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/limma_After-Before_intensity_qval_0.05_list.txt'), True)
        self.assertEqual(os.path.exists(outpath +'/limma/volcano_After-Before_intensity.jpeg'), True)
        
        # SPLOM
        self.assertEqual(os.path.exists(outpath +'/splom'), True)
        self.assertEqual(os.path.exists(outpath +'/splom/splom.png'), True)
        
        # Hierarchical clustering of replicates
        self.assertEqual(os.path.exists(outpath +'/replicate_hclust'), True)
        self.assertEqual(os.path.exists(outpath +'/replicate_hclust/replicate_cluster_ward.D2.jpeg'), True)

        # Heatmaps
        self.assertEqual(os.path.exists(outpath +'/heatmaps'), True)
        self.assertEqual(os.path.exists(outpath +'/heatmaps/heatmap_intensity_all.jpeg'), True)
        
        # PCA
        self.assertEqual(os.path.exists(outpath +'/PCA'), True)
        self.assertEqual(os.path.exists(outpath +'/PCA/all_identified_pca.png'), True)
        self.assertEqual(os.path.exists(outpath +'/PCA/all_identified_pca_labelled.png'), True)

    def test_b_protein_id_lists(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        with open(self.config) as f:
            config = yaml.load(f.read())
        proteingroups = pd.read_csv('testdata/example_mq/txt/analysis/proteins/target_proteins.txt',sep='\t')
        outfile = outpath +'/protein_ids.txt'
        mqparse.mq_txt.protein_id_lists(None, proteingroups, outfile)
        self.assertEqual(os.path.exists(outfile), True)
    
    def test_b_export_pg_fasta(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        outfile = outpath + '/pg.fasta'
        with open(self.config) as f:
            config = yaml.load(f.read())
        ids = 'testdata/example_mq/txt/analysis/proteins/protein_ids.txt'
        fasta_file = config['search_fasta']
        mqparse.mq_txt.export_pg_fasta(None, fasta_file, ids, outfile)
        self.assertEqual(os.path.exists(outfile), True)
    
    def test_b_ips_fasta(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        infile = 'testdata/example_mq/txt/analysis/fasta/leading_proteins.fasta'
        mqparse.mq_txt.ips_fasta(None, infile, outpath)
        #self.assertEqual(os.path.exists(outfile), True)
    
    def test_b_ips_genesets(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.mkdir(outpath)
        ipr='testdata/example_mq/txt/analysis/fasta/leading_proteins.fasta.tsv'
        proteingroups = pd.read_csv('testdata/example_mq/txt/analysis/proteins/target_proteins.txt',sep='\t')
        mqparse.mq_txt.ips_genesets(None, ipr, proteingroups, outpath)
        self.assertEqual(os.path.exists(outpath +'/ipr_target_proteins.txt'),True)
        self.assertEqual(os.path.exists(outpath +'/iprset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/ccset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/bpset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/mfset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/keggset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/metacycset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/reactomeset.Rdata'),True)
        self.assertEqual(os.path.exists(outpath +'/ecset.Rdata'),True)

        # gage.R --outdir $outpath --keggid $kegg_id  || rm -rf ${outpath}/gsea 
    
    def test_b_ips_gsea(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        design='testdata/example_mq/txt/analysis/diff/protein_experimental_design.R'
        table='testdata/example_mq/txt/analysis/diff/protein_normalization/msnbase/normalized.csv'
        t = pd.read_csv(table)
        
        t['Mygene.entrez'] = "5236;85469"
        kocol='Leading.Protein.Kegg.Orthology.ID'
        genecol='Mygene.entrez'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        shutil.copytree('testdata/example_mq/txt/analysis/gsea', outpath) 
        
        t.to_csv(outpath + '/table.csv')
        pval=1000
        mqparse.mq_txt.ips_gsea(None, outpath, outpath +'/LEVEL_1',  design, outpath + '/table.csv', genecol, kocol,  pval=pval)
        self.assertEqual(os.path.exists(outpath +'/LEVEL_1/After_Before/REACTOME.up.csv'), True)

    def test_b_up2ko(self):
        up='Q05025'
        ko = mqparse.mq_txt.up2ko(None, up)
        self.assertEqual(ko.ko,'K00134')
        self.assertEqual(ko.name,'glyceraldehyde 3-phosphate dehydrogenase [EC:1.2.1.12]')

        up='H2Q028'
        ko = mqparse.mq_txt.up2ko(None, up)
        self.assertEqual(ko.ko,'K21127')
        self.assertEqual(ko.name,'protein S100-A8') 
    
    def test_b_leading_protein_mygene(self):
        proteingroups = pd.read_csv('testdata/example_mq/txt/analysis/proteins/target_proteins.txt',sep='\t')
        ids = proteingroups['Identifier']
        #self.assertEqual(len(ids), len(set(ids)))
        #proteingroups['Leading Protein']= 'Q05025'
        proteingroups = mqparse.mq_txt.leading_protein_mygene(None, proteingroups)
        cols = proteingroups.columns.tolist()
        self.assertIn('Mygene.entrez', cols)
        
    def test_b_leading_protein_ko(self):
        proteingroups = pd.read_csv('testdata/example_mq/txt/analysis/proteins/target_proteins.txt',sep='\t')
        #proteingroups['Leading Protein']= 'Q05025'
        #proteingroups = mqparse.mq_txt.leading_protein_mygene(None, proteingroups)
        #proteingroups = proteingroups[proteingroups['Mygene.entrezgene'].notnull()]
        p = mqparse.mq_txt.leading_protein_ko(None, proteingroups)
        cols = p.columns.tolist()
        self.assertIn('Leading Protein Kegg Orthology ID', cols)
        self.assertIn('Leading Protein Kegg Orthology Name', cols)
        #print(p['Leading Protein Kegg Orthology Name']) 

    def test_b_aggregate_quant(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        table='testdata/example_mq/txt/analysis/diff/protein_normalization/msnbase/normalized.csv'
        design='testdata/example_mq/txt/analysis/diff/protein_experimental_design.R'
        df = pd.read_csv(table)
        agg_col='X_ec.term.union'
        quant_prefix ='iBAQ.'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        mqparse.mq_txt.aggregate_quant(None, df, agg_col, quant_prefix, outpath)
        self.assertEqual(os.path.exists(outpath +'/' + agg_col + '.csv'), True)

    def test_b_aggregate_quant_ko(self):
        outpath='testdata/example_mq/txt/analysis/testing'
        table='testdata/example_mq/txt/analysis/diff/protein_normalization/msnbase/normalized.csv'
        design='testdata/example_mq/txt/analysis/diff/protein_experimental_design.R'
        df = pd.read_csv(table)
        agg_col='Leading.Protein.Kegg.Orthology'
        quant_prefix ='iBAQ.'
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        mqparse.mq_txt.aggregate_quant(None, df, agg_col, quant_prefix, outpath)
        self.assertEqual(os.path.exists(outpath +'/' + agg_col + '.csv'), True)
        res = pd.read_csv(outpath +'/' + agg_col + '.csv')

    def test_b_parse_id(self):
        ids = 'tr|A6QNW7|A6QNW7_BOVIN;CON__ENSEMBL:ENSBTAP00000031360'
        correct = 'A6QNW7;CON__ENSEMBL:ENSBTAP00000031360'
        new_ids = mqparse.mq_txt.parse_ids(None, ids)
        self.assertEqual(correct, new_ids)

def main():
    unittest.main()

if __name__ == '__main__':
    main()

