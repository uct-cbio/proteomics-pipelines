#!/usr/bin/env python3
import json
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import os
from collections import defaultdict
import pickle
import yaml
import pgfunctions
from Bio import SeqIO
from collections import Counter
import pgfunctions
import functools
import numpy as np

config = yaml.load(open(sys.argv[1]).read(),  Loader=yaml.Loader)

output = sys.argv[2]

annotation_folder = output +'/annotations/'

output_folder = output + '/stats/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def get_fasta(path):
    recs = list(SeqIO.parse(path,'fasta'))
    rec_dict = {}
    for _ in recs:
        rec_dict[_.id] = _
    return rec_dict

strain_analysis = {}

analysis = pd.DataFrame()
annotations_summary = pd.DataFrame()

def add_protein_names(df):
    df = df.reset_index()
    if 'index' in df.columns:
        del df['index']

    df['Gene Name' ] = df['BLASTP'].map(proteome_genes)
    df['Protein Name' ] = df['BLASTP'].map(proteome_names)
    return df

def get_var_length(df):
    
    varId = df['VarId']
    if pd.notnull(varId) or varId is None:
        #print(varId)
        coord  = varId.split('|')[-1]
        coord = coord.split(')')[1]
        coords = coord.split(':')
        c1 = int(coords[0])
        c2 = int(coords[1])
        nuc_len = c2+1-c1
        amino_len = nuc_len/3
        return amino_len
strains=[]
global_group_acetyl = defaultdict(set)
global_group_identified = defaultdict(set)

for reference in config['reference']:
    combined = pd.read_csv(output+ '/{}_combined.csv'.format(reference), sep='\t')
    group_acetyl = defaultdict(set)
    group_identified = defaultdict(set)

    group_columns = [i for i in combined.columns if i.startswith('All peptides group ')]
    for col in group_columns:
        group = col.split('All peptides group ' )[1].split()[0] 
        group_peps = combined[col].dropna().tolist()
        group_peps = functools.reduce(lambda x, y: x + y.split("\n"), group_peps, [])
        for pep in group_peps:
            pep_seq = ''.join([char for char in pep if char.isupper()])
            if pep.startswith('_(ac)'):
                group_acetyl[group].add(pep_seq)
            group_identified[group].add(pep_seq)
    
    global_group_acetyl.update(group_acetyl)
    global_group_identified.update(group_identified)

    up_id = config['reference'][reference]['proteome_id']
    up_fasta = output + '/uniprot/{}/{}.fasta'.format(up_id, up_id)
    proteome = SeqIO.to_dict(SeqIO.parse(up_fasta,'fasta'))
    proteome_names = {}
    proteome_genes = {}

    for key in proteome:
        desc = proteome[key].description
        p_id = proteome[key].id
        pname = desc.split( p_id + ' ')[1].split(' OS=')[0] 
        gname = desc.split(p_id + ' ')[1].split(' GN=')[1].split(' PE=')[0] 
        p_id = p_id.split('|')[1]
        proteome_names[p_id] = pname
        proteome_genes[p_id] = gname
    strain_peptides = defaultdict(set)
    specific_strain_peptides = defaultdict(set)
    predicted_strain_peptides = defaultdict(set)
    #pep2genome_holder[strain] = {}
    peptide_dict = {}
    pep2orf = {} 

    variant_holder = {}
    trie_holder = {}
    proteome_id = config['reference'][reference]['proteome_id']
    with open(output + '/blast/orfs2proteins/{}_mapping.json'.format(proteome_id)) as f:
        ref_mapping =  json.loads(f.read())
    globally_specific_peptides = []
    reference_entries_identified = set() 
    
    translated_orfs = {}

    ################################
    # create the annotation folder #
    ################################
    export_folder = annotation_folder + '/export/{}/'.format(reference)
    if not os.path.exists(export_folder):
        os.makedirs(export_folder)
    
    
    combined_nterm_acetylation_all = []
    combined_nterm_acetylation_differences = []
    combined_nterm_acetylation_targets = []
    
    for strain in config['strains']:
        print(strain)
        strains.append(strain)
        spes_peps = combined["Specific peptides strain {}".format(strain)].dropna().tolist()
        spes_peps = functools.reduce(lambda x, y: x + y.split("\n"), spes_peps, [])
        specific_strain_peptides[strain].update(spes_peps)

        nspes_peps = combined["All peptides strain {}".format(strain)].dropna().tolist()
        nspes_peps = functools.reduce(lambda x, y: x + y.split("\n"), nspes_peps, [])
        strain_peptides[strain].update(nspes_peps)

        peptide_path = annotation_folder+'/{}_{}_peptide_annotations.csv'.format(reference, strain)
        trie_path = annotation_folder+'/{}_{}_variants_trie.fasta'.format(reference, strain)
        variant_path = annotation_folder+'/{}_{}_variants.fasta'.format(reference, strain)
        
        variant_holder[strain]= get_fasta(variant_path)
        trie_holder[strain]= get_fasta(trie_path)

        pep2genome = pd.read_csv(output + '/strains/{}/{}_mapped_peptides.csv'.format(strain, strain))
        #pep2genome_holder[strain] = pep2genome 
        orfmap = pep2genome.set_index('ORF_id')['ORF_translation'].to_dict()
        translated_orfs.update(orfmap)

        predicted_strain_peptides[strain].update(set(pep2genome['Peptide_sequence'].tolist()))
        
        if not os.path.exists(peptide_path):
            continue
        peptides = pd.read_csv(annotation_folder+'/{}_{}_peptide_annotations.csv'.format(reference, strain))
        peptides['BLASTP'] = peptides['ORF_id'].apply(lambda x : x.split('|')[1]).map(ref_mapping)
        
        
        def clean(x):
            try:
                #return ';'.join(x[0])
                return x[0]
            except:
                return x
        peptides['BLASTP'] = peptides['BLASTP'].apply(lambda x : clean(x))
        peptides = add_protein_names(peptides)
        
        


        #print(peptides['BLASTP'])
        del peptides['Unnamed: 0']
        peptides['VarLength'] = peptides.apply(get_var_length, axis=1)
        peptides = peptides.drop_duplicates(['ORF_id','Ref_id','PeptideSequence'])
        
        group_peptides = peptides[['PeptideSequence']]
        group_peptides = group_peptides.drop_duplicates()
        for col in group_acetyl:
            group_peptides['Identified {}'.format(col)] = peptides['PeptideSequence'].apply(lambda x : x in group_identified[col])
            group_peptides['Acetylated {}'.format(col)] = peptides['PeptideSequence'].apply(lambda x : x in group_acetyl[col])

        peplen = len(peptides)
        peptides  = pd.merge(peptides, group_peptides)
        peptides.to_csv(annotation_folder+'/{}_{}_peptide_annotations_groups.csv'.format(reference, strain))
        #print(peplen)
        #print(len(peptides)) 
        #print(len(group_peptides)) 
        assert len(peptides) == peplen


        peptide_dict[strain] = peptides
        peps = peptides['PeptideSequence'].tolist()

        def get_peptide_map(df, dct):
            sequence = df['PeptideSequence']
            orf_id = df['ORF_id']
            dct[sequence].add(orf_id)

        #m= defaultdict(set)
        #peptides.apply(get_peptide_map,m, axis=1)
        #print(m)
        res = pep2genome.groupby('Peptide_sequence')['ORF_id'].apply(list).to_dict()
        pep2orf[strain] = res

    strain_exclusive_tss_ = []
    strain_exclusive_tss_fasta_ = []
    
    annotation_peptides = []
    annotation_variants = {}
    def process_peptide_df(df, strain_annot):
        df['Cluster'] = df['BLASTP'] + '_' + df['PeptidePosition'].apply(str)
        df['Strain'] =  strain_annot
        return df
    for strain in peptide_dict.keys():
        peptides = peptide_dict[strain] 
        peptides = process_peptide_df(peptides, strain)
        orf_peptides = peptides[peptides['Ref_id'].isnull()]
        orf_peptides = orf_peptides.drop_duplicates(['ORF_id','PeptideSequence'], keep='first')
        tss_peptides = orf_peptides[orf_peptides['TSS'] == True] # only count TSSs 
        
        strain_orfs = peptides['ORF_id'].tolist()
        #for orf in strain_orfs:
        #    print(orf)
        #    orf_peps = peptides[peptides['ORF_id'] == orf]
        #    peps = list(set(orf_peps['PeptideSequence'].tolist()))
        
        #####################
        # Strain acetylated #
        #####################
        strain_acetylated = peptides[peptides['Strain_Nterm_Acetylated'] == True ]  
        combined_nterm_acetylation_all.append(strain_acetylated)
        strain_acetylated.to_csv(export_folder + '/{}_nterm_acetylation.csv'.format(strain), index=False)


        def process_group(df, pep2orf, strain, other_strain, other_strain_peptides):
            '''DataFrame of peptide annotations df, mapping of peptide to orfs
            of all strains as a dict from genome_mapping file as pep2orf, strain
            name as strain, other strain name as other_strain, DataFrame of 
            other strain annotations as other_strain_peptides'''
            other_orfs  = set()
            if len(df[df['Ref_id'].notnull()]) > 0:
                df = df[df['Ref_id'].notnull()] 
            #df = df.drop_duplicates()
            peps = list(set(df['PeptideSequence'].tolist()))
            

            ref_ids = set(df['Ref_id'].dropna().tolist())
            orf_novel = False

            novel_peps = []
            
            var_lengths = df.set_index('VarId')['VarLength'].to_dict()
            
            results = []
            orf_id = df['ORF_id'].tolist()[0]
            
            try:
                blasted = ';'.join(ref_mapping[orf_id.split('|')[1]])
            except:
                blasted = None
            for pep in peps:
                if pep in pep2orf[other_strain]:
                    other_orfs.update(set(pep2orf[other_strain][pep]))
                    #assert orf_novel == True  
                if not pep in strain_peptides[other_strain]:
                    #other_orfs.update(set(other_strain_peptides[other_strain_peptides['PeptideSequence']==pep]['ORF_id'].tolist()))
                    novel_peps.append(pep)
        
                #other_stain_orfs = predicted_strain_peptides[other_strain]['PeptideSequ:
            
            other_strain_peps = set() # get peptides of all mapped orfs in other strain
            for orf in other_orfs:
                other_strain_peps.update(other_strain_peptides[other_strain_peptides['ORF_id']==orf]['PeptideSequence'].tolist())
                #print(other_strain_peps, novel_peps)
           
            pep_intersect = other_strain_peps & strain_peptides[other_strain]
            #pep_intersect = set(peps) & strain_peptides[other_strain]
            if len(list(pep_intersect)) == 0:
                print("NOVEL ORF", reference, strain, ref_ids, novel_peps)
                orf_novel=True
            for pep in novel_peps:
                pepdf = df[df['PeptideSequence'] == pep]        
                pep_orf_ids = set()
                pep_orf_ids.update(pep2orf[strain][pep]) 
                if pep in pep2orf[other_strain]:
                    pep_orf_ids.update(pep2orf[other_strain][pep])
                pep_orf_ids = list(pep_orf_ids)
                pep_orf_ids.sort()
                pep_orf_ids = ';'.join(pep_orf_ids)

                # alternative tss
                #assert len(tss_peps) <= 1
                other_strain_identified = False
                if pep in strain_peptides[other_strain]:
                    other_strain_identified = True
                other_strain_predicted  = False
                if pep in predicted_strain_peptides[other_strain]:
                    other_strain_predicted = True
                tss = False
                tss_peps = pepdf[(pepdf['TSS'] == True) ] #& (pepdf['Ref_id'].isnull())]
                if len(tss_peps) > 0:

                    tss = True
                    varIds=tss_peps["VarId"].tolist()
                    assert len(varIds) > 0
                    for ref_var in varIds: 
                        var_length = int(var_lengths[ref_var] )
                        found = False
                        for orf in other_orfs:
                            orf_df = other_strain_peptides[other_strain_peptides['ORF_id'] == orf]
                            orf_df = orf_df.drop_duplicates('VarId')
                            other_var_lengths = orf_df.set_index('VarId')['VarLength'].to_dict()
                            otherVarIds= orf_df['VarId'].dropna().tolist()
                            for other_var in otherVarIds:
                                if var_lengths[ref_var] == other_var_lengths[other_var]:
                                    found = True
                        result = orf_id, ref_var, var_length, pep, other_strain_identified, other_strain_predicted, tss, orf_novel,  ','.join(ref_ids), pep_orf_ids, blasted 
                        if (found == False) or (other_strain_predicted == False):
                            results.append(result)
                elif (other_strain_identified == False) :
                    result = orf_id, None, None, pep, other_strain_identified, other_strain_predicted, tss, orf_novel, ','.join(ref_ids), pep_orf_ids, blasted 
                    results.append(result)
                    #print(result) 
                #if orf_novel == True:
                #    print(result)
                #    #print(pepdf.stack())
            
                # alternative pos-processing

                # alternative peptide sequence
            return results
        

        for other_strain in peptide_dict.keys():
            if other_strain == strain:
                continue
            print('Other strain: ', other_strain)
            #pepmap = pep2orf[other_strain]
            res = peptides.groupby('ORF_id').apply((lambda x : process_group(x, \
                    pep2orf,\
                    strain, \
                    other_strain, peptide_dict[other_strain])))
            results = [item for row in res for item in row]
            results = pd.DataFrame(results)
            results.columns = ['ORF_id', 'VarId','VarLength','PeptideSequence', 'OtherStrainIdentified', 'OtherStrainPredicted', 'TSS', 'SpecificPeptidesExclusive',  'Ref_Ids','ORF_ids', 'BLASTP']
            missing_cols = []
            for col in peptides.columns:
                if col.startswith('Identified ') or col.startswith('Acetylated '):
                    missing_cols.append(col)
            missing_cols.append('PeptideSequence')
            result_len = len(results)
            results = pd.merge(results, peptides[missing_cols].drop_duplicates(), how='left')
            assert len(results) == result_len
            results.to_csv(annotation_folder +'/{}_peptides_{}_vs_{}.csv'.format(reference, strain, other_strain), index=False)            
            strain_exclusive_tss=results[results['VarId'].notnull()]
            strain_exclusive_tss.to_csv(annotation_folder + '/{}_exclusive_tss_{}_vs_{}.csv'.format(reference, strain, other_strain))
            strain_exclusive_tss.loc[:, 'Comparison'] = '{} vs {}'.format(strain, other_strain)
            strain_exclusive_tss_.append(strain_exclusive_tss)
            


            for varId in strain_exclusive_tss['VarId'].tolist():
                strain_exclusive_tss_fasta_.append(variant_holder[strain][varId])
            
            
            ###############################
            # strain acetylation analysis #
            ###############################
            

            other_strain_peptides = peptide_dict[other_strain]
            other_strain_peptides = process_peptide_df(other_strain_peptides, other_strain) 
            other_strain_clusters = set(other_strain_peptides['Cluster'].tolist())
            other_strain_acetylated = other_strain_peptides[other_strain_peptides['Strain_Nterm_Acetylated'] == True ]  
            

            #other_strain_acetylated['SequenceExclusive'] = other_strain_acetylated['PeptideSequence'].apply(lambda x : x not in strain_peptides[strain])
            
            strain_acetylated_peptides = set(strain_acetylated['PeptideSequence'].tolist())
            other_strain_acetylated_peptides = set(other_strain_acetylated['PeptideSequence'].tolist())
            # get acetylated positions where the same position was identified by any peptide in the other strain
            strain_clusters_filt = set(strain_acetylated['Cluster'].tolist()) & other_strain_clusters
            
            # remove the position that are acetylated in the other strain
            strain_clusters_filt_exclusive = strain_clusters_filt - set(other_strain_acetylated['Cluster'].tolist())
            
            # get a dataframe of clusters where any acetylation was not seen in other strain
            strain_only_clusters = set(strain_acetylated[strain_acetylated['Cluster'].apply(lambda x : x not in set(other_strain_acetylated['Cluster'].tolist()))]['Cluster'].tolist()) 
            strain_common_peptides = set(strain_acetylated[strain_acetylated['PeptideSequence'].apply(lambda x : x  in set(other_strain_acetylated['PeptideSequence'].tolist()))]['Cluster'].tolist()) # remove any sequences that are also acetylated in other strain by peptide sequence
            nterm_acetylation_diff_filt = strain_only_clusters - strain_common_peptides
            nterm_acetylated_diff = strain_acetylated[strain_acetylated['Cluster'].apply(lambda x : x in nterm_acetylation_diff_filt)]
            nterm_acetylated_diff['PeptideIdentifiedInOtherStrain'] = nterm_acetylated_diff['PeptideSequence'].apply(lambda x : x in strain_peptides[other_strain])
            nterm_acetylated_diff['ClusterIdentifiedInOtherStrain'] = nterm_acetylated_diff['Cluster'].apply(lambda x : x in other_strain_clusters)
    
            # get a dataframe of clusters where the same cluster (ref and peptide position in ORF) was not acetylated in the other strain
            #exclusive_clusters = strain_acetylated[strain_acetylated['Cluster'].apply(lambda x : x in strain_clusters_filt_exclusive)] 
            #strain_acetylated_peptides_common = strain_acetylated_peptides & strain_peptides[other_strain]
            #filt_set = strain_acetylated_peptides_common - other_strain_acetylated_peptides
            
            # get a dataframe of peptides where the same peptide was not acetylated in other strain
            #exclusive_peptides = strain_acetylated[strain_acetylated['PeptideSequence'].apply(lambda x : x in filt_set)]
            #acetylated_targets =pd.concat([exclusive_clusters, exclusive_peptides])
            acetylated_targets =nterm_acetylated_diff[(nterm_acetylated_diff['PeptideIdentifiedInOtherStrain'] ==True) | (nterm_acetylated_diff['ClusterIdentifiedInOtherStrain']==True)] 

            #combined_acetylated = pd.concat([strain_acetylated,other_strain_acetylated])
            #combined_acetylated = combined_acetylated.drop_duplicates(['Strain', 'PeptideSequence'])
            #combined_acetylated=combined_acetylated.sort_values(['BLASTP', 'PeptidePosition'])
            #combined_acetylated=combined_acetylated[['PeptideSequence','PeptidePosition', 'BLASTP', 'ORF_id','Strain']]
            #combined_acetylated['Cluster'] = combined_acetylated['BLASTP'] + '_' + combined_acetylated['PeptidePosition'].apply(str)
            #filt = strain_acetylated.drop_duplicates(['Cluster'])
            #filt['AcetylationId'] = np.arange(1, len(filt) + 1)
            #acetyl_group = filt.set_index("Cluster")['AcetylationId'].to_dict()
            #strain_acetylated['AcetylationId'] = strain_acetylated['Cluster'].map(acetyl_group)
            
            #acetylated_targets['AcetylationId'] = acetylated_targets['Cluster'].map(acetyl_group)
             
            acetylated_targets['Comparison'] = '{}_vs_{}'.format(strain, other_strain)
            acetylated_targets.to_csv(export_folder + '/{}_vs_{}_nterm_acetylation_targets.csv'.format(strain, other_strain), index=False)
            
            nterm_acetylated_diff['Comparison'] = '{}_vs_{}'.format(strain, other_strain)
            nterm_acetylated_diff.to_csv(export_folder + '/{}_vs_{}_nterm_acetylation_differences.csv'.format(strain, other_strain), index=False)
            
            combined_nterm_acetylation_targets.append(acetylated_targets)
            combined_nterm_acetylation_differences.append(nterm_acetylated_diff)


        start_count = pgfunctions.start_site_count_map(tss_peptides) 
        peptides['ORF_start_count'] = peptides['ORF_id'].map(start_count)
        peptides['ORF_start_count'] = peptides['ORF_start_count'].fillna(0) 
        peptides = peptides.sort_values("ORF_start_count", ascending=False)
        
        tss_counts = peptides[['ORF_id','ORF_start_count']].drop_duplicates()
        tss_counts['ORF_start_count'] = tss_counts['ORF_start_count'].astype(int)
        tss_counts.to_csv(output_folder + 'tss_counts_{}.csv'.format(strain))
        #print(tss_counts)
        #print(peptides.head())

        start_count_df = peptides[['ORF_id','ORF_start_count']]
        start_count_df = start_count_df.drop_duplicates(["ORF_id", "ORF_start_count"])
        
        gene_model_count = Counter(start_count_df['ORF_start_count'])
        for key in gene_model_count.keys():
            col = 'ORFs with ' + str(int(key)) + ' TSSs'
            analysis.loc[strain, col] = gene_model_count[key]
        
        #codon_table = orf_peptides.drop_duplicates(['ORF_id','StartPosition'])

        specific_peptides = orf_peptides.drop_duplicates('PeptideSequence')['PeptideSequence'].tolist()

        #specific_strain_peptides[strain] = specific_peptides
        globally_specific_peptides += specific_peptides

        analysis.loc[strain, 'GloballySpecificPeptides']   = len(specific_peptides)
        

        #tss_peptides = orf_peptides[orf_peptides['TSS'] == True]
        
        #print(Counter(orf_peptides['StartType']))

        #peptides = peptides[peptides['Ref_id'].notnull()]
        #peptides  = peptides.drop_duplicates(['StartPosition', 'ORF_id','Ref_id'])
        #peptides = peptides[peptides['StartPosition'].notnull()]
        #peptides = peptides.sort_values('ORF_start_count')
        #print(peptides.tail())
        #print(orf_peptides.head(1).stack())
       
        
        def process_annot_df(peptides, strain, reference):
            global analysis
            global annotations_summary 
            global translated_orfs

            df = peptides.copy()
            assert 'Enzymatic N-terminus non-TSS peptide' in df['StartType'].tolist()
            # strain refs
            strain_reference_entries = df['BLASTP'].dropna().unique().tolist()
            reference_entries_identified.update(strain_reference_entries)
            analysis.loc[strain, ' {} Reference entries identified'.format(reference)] = len(strain_reference_entries)
            refs = df[df['BLASTP'].apply(lambda x : x in reference_entries_identified)]
            refs = refs.sort_values(['BLASTP','ORF_id'])
            filt_refs = refs.drop_duplicates(['BLASTP'])
            filt_refs['AnnotationId'] = np.arange(1, len(filt_refs) + 1)
            filt_refs_map = filt_refs.set_index('BLASTP')['AnnotationId'].to_dict()
            refs['AnnotationId'] = refs['BLASTP'].map(filt_refs_map)
            refs = refs[['ORF_id','BLASTP', 'AnnotationId']].drop_duplicates()
            refs.to_csv(export_folder + '/{}_Annotated_ORFs.csv'.format(strain), index=False)

            novel = df[~df['BLASTP'].apply(lambda x : x in reference_entries_identified)]
            novel = novel.sort_values(['ORF_id'])
            filt_novel = novel.drop_duplicates(['ORF_id'])
            filt_novel['AnnotationId'] = np.arange(1, len(filt_novel) + 1)
            filt_novel_map = filt_novel.set_index('ORF_id')['AnnotationId'].to_dict()
            novel['AnnotationId'] = novel['ORF_id'].map(filt_novel_map)
            novel = novel[['ORF_id', 'AnnotationId']].drop_duplicates()
            novel.to_csv(export_folder + '/{}_Novel_ORFs.csv'.format(strain), index=False)
            
            novel_export = []
            for orf_id in novel['ORF_id'].unique().tolist():
                orf_seq = translated_orfs[orf_id]
                rec =Bio.SeqRecord.SeqRecord(id=orf_id, seq=Bio.Seq.Seq(orf_seq),description='Novel ORF')
                novel_export.append(rec)
            SeqIO.write( novel_export,  export_folder + '/{}_Novel_ORFs.fasta'.format(strain), 'fasta')
            

            
            df = df.sort_values(['Ref_id','VarLength'],ascending=False)
            nr_df = df.drop_duplicates(['VarId'], keep='first')
             


            filt_df =  nr_df.drop_duplicates(['Ref_id','VarLength'], keep='first')
            # where the peptide starts in the same position in different strain ORFs, but the stop codon may be different leading to different variant sequence lengths.
            filt_df = filt_df.drop_duplicates(['Ref_id','StartPosition','AnnotationType'], keep='first')
            filt_df =  filt_df[filt_df['Ref_id'].notnull()]
            
            filt_vars = set(filt_df['VarId'].tolist())

            annotations = Counter(filt_df['AnnotationType'].dropna())
            annot_summaries = []
            variants = trie_holder[strain]
            for annot in annotations:
                analysis.loc[strain,annot + ' - ' + reference] = annotations[annot]
                if strain == 'Combined':
                    annotations_summary.loc[reference, annot ] = annotations[annot]
                annot_proteins = Counter(nr_df[nr_df['AnnotationType']==annot]['Ref_id'])
                
               
                variant_df=df[df['AnnotationType']==annot]
                filt = variant_df[variant_df['VarId'].apply(lambda x : x in filt_vars)]
                filt = filt.drop_duplicates(['VarId'])
                filt['AnnotationId'] = np.arange(1, len(filt) + 1)
                annot_group = filt.set_index("VarId")['AnnotationId'].to_dict()

                variant_ids=variant_df['VarId'].unique().tolist()
                variant_df['AnnotationId'] = variant_df['VarId'].map(annot_group)
                variant_df['AnnotationId'] = variant_df['AnnotationId'].ffill()
                variant_df['AnnotationId'] = variant_df['AnnotationId'].astype(int)

                export = [variants[i] for i in variant_ids]
                name = '_'.join(annot.split())
                SeqIO.write( export,  export_folder + '/{}_{}.fasta'.format(strain, name), 'fasta')
                variant_df.to_csv(export_folder + '/{}_{}.csv'.format(strain, name))

                annot_proteins_df = pd.DataFrame(annot_proteins.values(),index=annot_proteins.keys(), columns=[annot], dtype=int)
                annot_summaries.append( annot_proteins_df)
                num_proteins = len(annot_proteins.keys())
                analysis.loc[strain,annot + ' (proteins) - ' + reference] = num_proteins

            annot_summaries = pd.concat(annot_summaries, axis=1)
            annot_summaries = annot_summaries.fillna(0)
            annot_summaries = annot_summaries.astype(int)
            
            annot_summaries = annot_summaries.reset_index().rename(columns={'index':'Protein'})
            annot_summaries.to_csv(export_folder+ '/{}_Annotation_Counts.csv'.format(strain),index=False)
            return peptides, variants
        annotation_df, variants = process_annot_df(peptides, strain, reference)
        annotation_peptides.append(annotation_df)
        annotation_variants.update(variants)
        #print(Counter(peptides['Ref_id']))
        #tries = get_fasta(trie_path)
   
    
    ###################################################################################
    # make a combined dataframe off all nterm acetylatyion protein IDs across strains #
    ###################################################################################
    def process_acetylation(df):
        strains = list(set(df['Strain'].tolist()))
        strains.sort()
        return ';'.join(strains)
    

    combined_acetylated = pd.concat(combined_nterm_acetylation_all)
    acetylated_proteins = combined_acetylated[combined_acetylated["BLASTP"].notnull()].groupby('BLASTP').apply((lambda x : process_acetylation(x)))
    if not isinstance(acetylated_proteins, pd.DataFrame):
        acetylated_proteins = pd.DataFrame(acetylated_proteins)
    acetylated_proteins.rename({0: "Strains"}, axis=1, inplace=True)
    
    acetylated_proteins = add_protein_names(acetylated_proteins)
    acetylated_proteins.to_csv(export_folder + '/Combined_Proteins_Nterm_Acetylation.csv')
    combined_acetylated.to_csv(export_folder + '/Combined_Nterm_Acetylation.csv', index=False)
    
    combined_targets = pd.concat(combined_nterm_acetylation_targets)
    target_proteins = combined_targets[combined_targets["BLASTP"].notnull()].groupby('BLASTP').apply((lambda x : process_acetylation(x)))
   # Convert to DataFrame if it's not already
    if not isinstance(target_proteins, pd.DataFrame):
        target_proteins = pd.DataFrame(target_proteins)
    target_proteins.rename({0: "Strains"}, axis=1, inplace=True)
    target_proteins = add_protein_names(target_proteins)
    target_proteins.to_csv(export_folder + '/Combined_Proteins_Nterm_Acetylation_Targets.csv')
    combined_targets.to_csv(export_folder + '/Combined_Nterm_Acetylation_Targets.csv', index=False)
                    
    combined_differences = pd.concat(combined_nterm_acetylation_differences)
    diff_proteins = combined_differences[combined_differences["BLASTP"].notnull()].groupby('BLASTP').apply((lambda x : process_acetylation(x)))
   # Convert to DataFrame if it's not already
    if not isinstance(diff_proteins, pd.DataFrame):
        diff_proteins = pd.DataFrame(diff_proteins)
    diff_proteins.rename({0: "Strains"}, axis=1, inplace=True)
    assert 'Strains' in diff_proteins.columns
    diff_proteins = add_protein_names(diff_proteins)
    diff_proteins.to_csv(export_folder + '/Combined_Proteins_Nterm_Acetylation_Differences.csv')
    combined_differences.to_csv(export_folder + '/Combined_Nterm_Acetylation_Differences.csv', index=False)





    analysis.loc['Combined', 'GloballySpecificPeptides'] = len(set(globally_specific_peptides))
    #analysis.loc['Combined', '{} Reference entries identified'.format(reference)] = len(reference_entries_identified)
    annotation_df = pd.concat(annotation_peptides)
    trie_holder['Combined'] = annotation_variants
    annotations_summary.loc[reference, 'Reference entries identified'] = len(reference_entries_identified)
    annotation_df = process_annot_df(annotation_df, 'Combined', reference)
    #annotation_df = annotation_df.drop_duplicates(['Ref_id','VarLength'], keep='first')
    # where the peptide starts in the same position in different strain ORFs, but the stop codon may be different leading to different variant sequence lengths.
    #annotation_df = annotation_df.drop_duplicates(['Ref_id','StartPosition','AnnotationType'], keep='first')
    #annotation_df = annotation_df[annotation_df['Ref_id'].notnull()]
    #temp = annotation_df[annotation_df['AnnotationType']=='Annotated TSS validated']
    #filt = temp[temp['Ref_id'].duplicated(keep=False)]
    #filt.to_csv(export_folder+ '/debug.csv')
    #annotations = Counter(annotation_df['AnnotationType'].dropna())
    #annot_summaries = []
    #for annot in annotations:
    #    analysis.loc['Combined', annot + ' - ' + reference] = annotations[annot]
    #    annotations_summary.loc[reference, annot ] = annotations[annot]
    #    annot_proteins = Counter(annotation_df[annotation_df['AnnotationType']==annot]['Ref_id'])
    #    annot_proteins_df = pd.DataFrame(annot_proteins.values(),index=annot_proteins.keys(), columns=[annot], dtype=int)
    #    annot_summaries.append( annot_proteins_df)
    #    num_proteins = len(annot_proteins.keys())
    #    analysis.loc['Combined', annot + ' (proteins) - ' + reference] = num_proteins
    #    #annotations_summary.loc[reference, annot + ' (proteins)'] = num_proteins

    #annot_summaries = pd.concat(annot_summaries, axis=1)
    #annot_summaries = annot_summaries.fillna(0)
    #annot_summaries = annot_summaries.astype(int)
    
    #annot_summaries = annot_summaries.reset_index().rename(columns={'index':'Protein'})
    #annot_summaries.to_csv(export_folder+ '/annotation_counts.csv',index=False)


    combined_exclusive_tss = pd.concat(strain_exclusive_tss_)
    
    combined_exclusive_tss.to_csv(annotation_folder + '/{}_exclusive_tss.csv'.format(reference),index=False)
    SeqIO.write( strain_exclusive_tss_fasta_,  annotation_folder + '/{}_exclusive_tss.fasta'.format(reference), 'fasta')
        # export table, var_fasta, trie fasta
    
    def process_diff_tss(df):
        comparisons = list(set(df['Comparison'].tolist()))
        comparisons.sort()
        return ';'.join(comparisons)

    diff_tss = combined_exclusive_tss[combined_exclusive_tss["BLASTP"].notnull()].groupby('BLASTP').apply((lambda x : process_diff_tss(x)))
   # Convert to DataFrame if it's not already
    if not isinstance(diff_tss, pd.DataFrame):
        diff_tss = pd.DataFrame(diff_tss)
    diff_tss.rename({0: "Comparisons"}, axis=1, inplace=True)
    assert 'Comparisons' in diff_tss.columns
    diff_tss = add_protein_names(diff_tss)
    diff_tss.to_csv(export_folder + '/Combined_Proteins_TSS_Differences.csv')

# export summary table
tss_cols = []
other_cols = []

for col in analysis.columns:
    if col.startswith('ORFs with '):
        tss_cols.append(col)
        analysis[col]  = analysis[col].fillna(0)
        analysis[col] = analysis[col].astype(int)
    else:
        other_cols.append(col)
        try:
            analysis[col] = analysis[col].astype(int)
        except:
            pass
strains = list(set(strains))
analysis.loc['Combined',tss_cols] = analysis.loc[strains, tss_cols].sum()
tss_cols.sort(reverse=False, key=lambda x : int(x.split("ORFs with ")[1].split()[0]))
combined_cols = other_cols + tss_cols
analysis['GloballySpecificPeptides'] = analysis['GloballySpecificPeptides'].apply(lambda x : int(x))
analysis = analysis[combined_cols]
analysis = analysis.T
analysis = analysis.reset_index()
analysis  = analysis.rename(columns ={'index':''})
analysis.to_csv(annotation_folder + '/summary.csv', index=False)

annotations_summary = annotations_summary.T
annotations_summary = annotations_summary.astype('int')
annotations_summary = annotations_summary.reset_index()
annotations_summary = annotations_summary.rename(columns={'index':''})
annotations_summary.to_csv(annotation_folder + '/annotations.csv', index=False)



#######################
# Dump JSON to a file #
#######################

for group in global_group_identified.keys():
    global_group_identified[group] = list(global_group_identified[group])
    global_group_acetyl[group] = list(global_group_acetyl[group])

with open(annotation_folder+'/global_group_acetyl.json','w') as file:
    json.dump(dict(global_group_acetyl), file)
with open(annotation_folder+'/global_group_identified.json','w') as file:
    json.dump(dict(global_group_identified), file)
