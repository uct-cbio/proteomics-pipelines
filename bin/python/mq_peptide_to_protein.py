#!/usr/bin/env python3

import multiprocessing
import numpy as np
import pandas as pd
import sys
import importlib.machinery
import Bio; from Bio import SeqIO
import sequtils
import shutil
import algo
import os
import subprocess
from collections import defaultdict
import json
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
from collections import Counter
from io import StringIO
import uniprot
import pickle
import yaml

# to do : blast orfs between strains to identify ICDS's without protein evidence
# check if there are orfs on either side of a fs , that dont have peptides - inactivating fs relative to other strain

config = yaml.load(open(sys.argv[1]), Loader=yaml.Loader)
output = sys.argv[2]

pg = pd.read_csv(config['mq_txt'] +'/proteinGroups.txt',sep='\t')
pg =pg[(pg['Potential contaminant'].isnull()) & (pg['Reverse'].isnull())]
peptides=pd.read_csv(config['mq_txt'] +'/peptides.txt',sep='\t')
peptides=peptides[(peptides['Potential contaminant'].isnull()) & (peptides['Reverse'].isnull())]

evidence=pd.read_csv(config['mq_txt'] +'/evidence.txt',sep='\t')

evidence=evidence[(evidence['Potential contaminant'].isnull()) & (evidence['Reverse'].isnull())]



for reference in config['reference']:
    var_fasta = defaultdict(list)
    trie_fasta = defaultdict(list)
    peptide_annotations = defaultdict(list)
    
    proteome_id = config['reference'][reference]['proteome_id']
    with open(output +'/mapping/{}_peptides.json'.format(proteome_id)) as f:
        reference_peptides = json.loads(f.read())


    reference_proteome = SeqIO.to_dict(list(SeqIO.parse(output +'/uniprot/{}/{}.fasta'.format(proteome_id, proteome_id),'fasta')))

    ref_entry = reference_peptides[list(reference_peptides.keys())[0]][0].split('|')[1]

    #gi_operons = json.loads(open( output + '/mapping/operons.json').read() )

    idmapping = pd.read_csv(output + '/uniprot/{}/{}.tab'.format(proteome_id, proteome_id),sep='\t')
    blast_mapping = json.loads(open(output +'/blast/orfs2proteins/{}_mapping.json'.format(proteome_id)).read())
    #taxon_peptides=pickle.load(open(output +'/mapping/pep2entry.p','rb'))
    #taxon_data=pickle.load(open(output +'/mapping/entrydata.p','rb'))

    reference_features = sequtils.gff3(output +'/uniprot/{}/{}.gff3'.format(proteome_id,proteome_id))

    #orf_features = sequtils.gff3(output +'/fasta/proteins.fasta.gff3')
    #with open(output + '/fasta/id_mapping.json') as f:
    #    orf_feature_mapping = json.loads(f.read())
    #orf_features.expand_table(orf_feature_mapping)

    #orf_features = pd.read_csv(output +'/fasta/proteins.fasta.tsv',sep='\t')
    orf_features = sequtils.gff3(output +'/fasta/proteins.fasta.tsv')
    
    design = pd.read_csv(config['design'])


    design = design[design['exclude'] != '+']

    if 'rename' in design.columns:
        samples =list(set(design['rename'].tolist()))
        rename = design.set_index('rename')['sample'].to_dict()
    else:
        samples =list(set(design['sample'].tolist()))
        design['rename'] = design['sample']
        rename = design.set_index('rename')['sample'].to_dict()
    # Get the mapping of samples to strain
    strain_dct = design.set_index('rename')['Strain'].to_dict()
    strainGroup_dct = design.set_index('rename')['StrainCondition'].to_dict()

    try:
        shutil.rmtree(output+'/clustalw')
        os.mkdir(output+'/clustalw')
        os.mkdir(output+'/clustalw/orf')
        os.mkdir(output+'/clustalw/prot')
    except: 
        os.mkdir(output+'/clustalw')
        os.mkdir(output+'/clustalw/orf')
        os.mkdir(output+'/clustalw/prot')

    global_non_specific_peptides=[]
    global_specific_peptides=[]
    global_specific_paralog_peptides = []

    strain_map={}
    strain_sets={}

    for strain in config['strains']:
        strain_sets[strain]=set()
        if config['strains'][strain]['assembly'] != None:
            st =  pd.read_csv(output +'/strains/{}/{}_mapped_peptides.csv'.format(strain, strain))
            peptide_map = sequtils.mapping2peptides(st, config['translation_table'])
            global_non_specific_peptides += peptide_map.non_specific()
            global_specific_peptides += peptide_map.specific()
            global_specific_paralog_peptides += peptide_map.paralogous_peptides
            strain_map[strain]  = peptide_map
            strain_sets[strain] = set(peptide_map.mapping['Peptide_sequence'].tolist())

    global_specific_peptides = set(global_specific_peptides) - set(global_non_specific_peptides)
    global_non_tryptic_nterm = global_specific_peptides.copy()
    global_non_tryptic_cterm = global_specific_peptides.copy()
    global_non_atg_starts    = global_specific_peptides.copy()
    global_specific_paralog_peptides = set(global_specific_paralog_peptides) - set(global_non_specific_peptides)
    global_non_specific_peptides = set(global_non_specific_peptides)


    strain_sets_nterm_acetylated = {}

    for strain in config['strains']:
        if config['strains'][strain]['assembly'] != None:
            peptide_map = strain_map[strain]
            st = peptide_map.mapping
            nterm_acetylated = set(st[st['Strain_Nterm_Acetylated'] == '+']['Peptide_sequence'].tolist())
            strain_sets_nterm_acetylated[strain] = nterm_acetylated

    for strain in config['strains']:
        if config['strains'][strain]['assembly'] != None:
            peptide_map = strain_map[strain]
            st = peptide_map.mapping
            strain_non_tryptic_nterm = set(st[st['Peptide_tryptic_nterm'].apply(str) == 'False' ]['Peptide_sequence'])
            strain_non_tryptic_cterm = set(st[st['Peptide_tryptic_cterm'].apply(str) == 'False' ]['Peptide_sequence'])
            strain_non_atg_starts    = peptide_map.non_atg_m()
            global_non_tryptic_nterm = global_non_tryptic_nterm & strain_non_tryptic_nterm
            global_non_tryptic_cterm = global_non_tryptic_cterm & strain_non_tryptic_cterm
            global_non_atg_starts    = global_non_atg_starts    & strain_non_atg_starts
                
    # Annotation status
    annotated_peptides = set(reference_peptides.keys())
    novel_peptides = global_specific_peptides - annotated_peptides

    # Get strain_exclusive
    strain_sets_exclusive={}
    for strain in strain_sets:
        tmp = strain_sets[strain] 
        for comparison_strain in strain_sets:
            if comparison_strain != strain:
                comp = strain_sets[comparison_strain]
                tmp = tmp - comp  # remove all peptides in the another strain
        strain_sets_exclusive[strain] = tmp - annotated_peptides # make sure all are novel

    def peptide_list_blast(peptides, targets, features):
        count=1
        blasted = []
        
        feature_dict = {}
        blast_dict = {}

        for p in peptides:
            feature_dict[p] = {}
            blast_dict[p] = {}

            for target in targets:
                pairwise_blast = []
                pairwise_features = []

                id = 'peptide_{}_{}'.format(strain, str(count))
                seq = Seq(p)
                rec = SeqRecord(id = id, seq = seq)
                
                out = sequtils.pairwise_blast(rec, target, output)
                if len(out.differences) > 0:
                    header = '***Results for peptide {} in record {}***'.format(p, target.id) 
                    pairwise_blast.append(header)
                    pairwise_blast.append(out.results)
                    pairwise_blast.append('\nIdentified polymorphism positions: {}'.format('; '.join(['{} ({})'.format(str(i), out.variants[i]) for i in out.differences])))
                    #try:
                    overlap = out.feature_overlap(features)
                    #except:
                    #    # help this code please 
                    #    overlap = None


                    if overlap != None: 
                        pairwise_features.append(header)
                        pairwise_features.append(overlap)
                    
                    count += 1

                    pairwise_blast = '\n\n'.join(pairwise_blast)
                    pairwise_features = '\n\n'.join(pairwise_features)
                    
                    feature_dict[p][target.id] = pairwise_features
                    blast_dict[p][target.id] = pairwise_blast

        return feature_dict, blast_dict


    #pg = pg[pg['id']==349] #mutant
    #pg = pg[pg['id']==934] #mutant
    #pg = pg[(pg['id']==1528) | (pg['id']==934)] #mutants
    #pg = pg[pg['id'] ==  410]  two frames in 5527

    #pg = pg[pg['id'] == 934]
    #pg = pg[pg['id'] == 2362]
    #pg = pg[pg['id'] ==246]

    #pg = pg[pg['id'] == 35]

    strain_columns = defaultdict(list)

    strain_peptides = defaultdict(set)
    strain_evidences = defaultdict(set)

    orf_start_codons = {}
    for strain in config['strains']:
        
        orf_start_codons[strain] = defaultdict(list)

        samplecols = []
        
        groupcols = []

        for sample in samples:
            if strain_dct[sample] == strain:
                sample_peps = set(evidence[evidence['Experiment'] == rename[sample]]['Sequence'].tolist()) 
                sample_evs = set(evidence[evidence['Experiment'] == rename[sample]]['Modified sequence'].tolist()) 
                strain_peptides[strain].update(sample_peps)
                strain_evidences[strain].update(sample_evs)

                samplecols.append("All peptides sample {}".format(sample))       
                groupcols.append("All peptides group {}".format(strainGroup_dct[sample]))

        scols = sorted(list(set(samplecols))) + list(set(groupcols))
        scols.append("All peptides strain {}".format(strain))
        scols.append("Specific peptides strain {}".format(strain))
        scols.append("Exclusive peptides strain {}".format(strain))       
        scols.append("Exclusive peptides reference BLAST strain {}".format(strain))
        scols.append("Exclusive peptides polymorphism reference feature overlap strain {}".format(strain))
        scols.append("Novel specific peptides strain {}".format(strain))
        scols.append("Annotated specific peptides strain {}".format(strain))
        scols.append("Non-genomic peptides strain {}".format(strain))
        scols.append("Translated orfs strain {}".format(strain))
        scols.append('Variant orfs strain {}'.format(strain))
        scols.append('Annotation type strain {}'.format(strain))
        
        scols.append('Strain ICDS intrastrain evidence {}'.format(strain))
        scols.append('Strain ICDS interstrain exclusive evidence {}'.format(strain))
        scols.append('Strain ICDS intrastrain peptide validation {}'.format(strain))
        scols.append('Strain ICDS interstrain exclusive peptide validation {}'.format(strain))

        scols.append("Best orf-reference blast evalue strain {}".format(strain))
        scols.append("Best orf-reference blast match {}".format(strain))
        scols.append("Reference BLAST strain {}".format(strain))
        scols.append('Identifier')
        scols.append('Gene_OrderedLocusName')
        #scols.append('GI')
        scols.append('GeneID')
        scols.append('UniProtKB-ID')
        scols.append("Reference BLAST - all strains")
        strain_columns[strain] = scols

    #combined= pg[:100]

    #pg = pg.head(10)

    combined = pg
    
    
    for row in combined.iterrows():
        combined_blasted = set()
        row_paralogous = '-'

        peptide_ids = row[1]['Peptide IDs'].split(';')  # Get the ids of peptides in the row
        evidence_ids = row[1]['Evidence IDs'].split(';') # Get the ids of envidence.txt entries
        row_peps = peptides[peptides['id'].apply(str).isin(peptide_ids)]
        row_evs = evidence[evidence['id'].apply(str).isin(evidence_ids)]

        all_row_peptides = set(row_peps['Sequence'].tolist())
        row_paralog_peptides  = all_row_peptides & global_specific_paralog_peptides
        row_specific  = all_row_peptides - global_non_specific_peptides  #- row_paralog_peptides
        row_novel     = row_specific & novel_peptides
        row_annotated = row_specific & annotated_peptides
        row_non_tryptic_cterm = row_specific & global_non_tryptic_cterm  
        row_non_tryptic_nterm = row_specific & global_non_tryptic_nterm  
        row_non_atg_starts = row_specific & global_non_atg_starts 
        
        #print('specific peptides : ', len(list(row_specific)))
        if len(row_paralog_peptides) > 0:
            row_paralogous = '+'

        #row_ref_strain_peps = set()

        strain_peps=defaultdict(list)
        
        group_evs=defaultdict(list)
        
        for sample in samples:
            strain=strain_dct[sample]
            group=strainGroup_dct[sample]
            #sample_peps = list(set(row_peps[row_peps['Experiment {}'.format(sample)]>=1]['Sequence'].tolist()))
            sample_evs = list(set(row_evs[row_evs['Experiment'] == rename[sample]]['Modified sequence'].tolist()))
            sample_peps = list(set(row_evs[row_evs['Experiment'] == rename[sample]]['Sequence'].tolist()))
            sample_column = "All peptides sample {}".format(sample)
            combined.loc[row[0], sample_column] = '\n'.join(sorted(list(sample_evs)))
            
            strain_peps[strain] += sample_peps
            
            group_evs[group] += sample_evs
            
            #if strain in config['reference_strains']:
            #    row_ref_strain_peps.update(sample_peps)
        
        #row_ref_strain_peps = row_ref_strain_peps & row_specific
        
        for group in group_evs:
            combined.loc[row[0],"All peptides group {}".format(group)]='\n'.join(sorted(list(set(group_evs[group]))))
        
        combined.loc[row[0], 'Row paralogous specific peptides'] = '\n'.join(sorted(list(row_paralog_peptides)))
        combined.loc[row[0], 'Row non-paralogous specific peptides'] = '\n'.join(sorted(list(row_specific)))

        all_orfs_fasta=[]
        all_prot_fasta=[]
        fasta_holder = {}
        peptide_holder = {}
        id_holder={}
        frameshifts =[]
        
        strains_exclusive=[]
             
        mapped_reference = set()
        
        for peptide in row_annotated:
            mapped_reference.update(reference_peptides[peptide]) 
        
        refps = []
        
        for mr in mapped_reference:
            refps.append(reference_proteome[mr])
        
        assert len(refps) == len(mapped_reference)

        reftrie  =  sequtils.list_trie_upper(refps, all_row_peptides)
        
        variant_dict = {}
        
        min_eval_ref_blast = config['eval_cutoff']
        best_eval_strain = None
        
        row_variant_features_ref, row_variant_blast_ref = peptide_list_blast(row_specific, refps, reference_features)
        
        row_variant_features_ref_list = []
        
        row_variant_blast_ref_list = []
        
        for pep in row_variant_blast_ref:
            row_variant_blast_ref_list += list(row_variant_blast_ref[pep].values())
            if pep in row_variant_features_ref:
                row_variant_features_ref_list += list(row_variant_features_ref[pep].values()) 
        
        mapped_orfs = []
        
        row_specific_strain_peps = {}
        
        all_orf_ids = []

        for strain in config['strains']:
            
            strain_blasted = set()
            
            x = []
            y = []
            
            if strain in strain_map:
                peptide_mapping = strain_map[strain]
                x,y = peptide_mapping.return_fasta( row_specific )
            orf_ids = [] 

            for rf in x:
                orf_ids.append(rf.id)
                rf_id = rf.id.split('|')[1]
                if rf_id in blast_mapping:
                    rfmap=set(blast_mapping[rf_id])
                    strain_blasted.update(rfmap)
                    combined_blasted.update(rfmap)
            
            all_orf_ids += orf_ids

            orfs_peps = set(strain_map[strain].get_peptides(orf_ids)) & strain_peptides[strain]  # Combined all identified peptides for the orf and strain
            
            #st_peps = orfs_peps
            st_peps = set(strain_peps[strain]) | orfs_peps # add other peptides for the identified ORFs to the same row. Some of the identified specific peptides may not belong to mapped ORFs (if predicted by ORFs of other strains but identified in the sample - these are termed genome exclusive
            specific_strain_peps = st_peps & global_specific_peptides
            
            # Dict to hold the specific peptides per strain
            
            row_specific_strain_peps[strain] = specific_strain_peps
            genome_unmapped = specific_strain_peps - strain_sets[strain]
            other_unmapped = row_specific - specific_strain_peps   # remove peptides identified in the same strain
            other_unmapped = other_unmapped - strain_sets[strain]  # remove peptides identified in other strains that can be mapped to the genome

            strain_exclusive_novel = specific_strain_peps & strain_sets_exclusive[strain]
            
            strain_n_term_acetylated= specific_strain_peps & strain_sets_nterm_acetylated[strain]
            
            all_orfs_fasta += x
            all_prot_fasta += y
            
            trie = sequtils.list_trie_upper(y, st_peps)
            
            annotation_trie = []
            annotation_type = []
            

            for _ in x:
                # compare the ORF to itself, get TSS positions if no other annotation is available
                pgomics = sequtils.proteogenomics(specific_strain_peps, _, None, n_term_acetylated=strain_n_term_acetylated) # all peptides for the strain.
                #orf_start_codons[strain][_.id][ _.id] = pgomics.start_codons_found
                var_fasta[strain] += pgomics.variant_sequences
                trie_fasta[strain] += pgomics.variant_sequences_trie
                peptide_annotations[strain] += pgomics.get_peptide_info()
                
                for refp in refps:
                    #if len(refps) > 0:
                    pgomics = sequtils.proteogenomics(specific_strain_peps, _, refp, n_term_acetylated=strain_n_term_acetylated) # all peptides for the strain.
                    #orf_start_codons[strain][_.id][ refp.id] = pgomics.start_codons_found
                    annotation_trie += pgomics.variant_sequences_trie
                    annotation_type += pgomics.annotation_type
                    
                    peptide_annotations[strain] += pgomics.get_peptide_info()
                    var_fasta[strain] += pgomics.variant_sequences
                    trie_fasta[strain] += pgomics.variant_sequences_trie

            fasta_holder['{}_prots'.format(strain)] = y
            fasta_holder['{}_nucs'.format(strain)] = x
            peptide_holder['{}_genome_unmapped_same'.format(strain)]=genome_unmapped
            peptide_holder['{}_genome_unmapped_other'.format(strain)]=other_unmapped
            id_holder[strain] =[rec.id for rec in y]
            
            mapped_orfs += y
            
            for comparison_strain in strain_peps:
                if comparison_strain != strain:
                    comp_peps=set(strain_peps[comparison_strain])
                    strain_exclusive_novel = strain_exclusive_novel - comp_peps
            
            if len(strain_exclusive_novel) > 0:
                strains_exclusive.append(strain)
            
            combined.loc[row[0], "All peptides strain {}".format(strain)]='\n'.join(sorted(st_peps))
            combined.loc[row[0], "Specific peptides strain {}".format(strain)]='\n'.join(sorted(specific_strain_peps))
            combined.loc[row[0], "Exclusive peptides strain {}".format(strain)]='\n'.join(sorted(strain_exclusive_novel) )       
            combined.loc[row[0], "Novel specific peptides strain {}".format(strain)]='\n'.join(sorted(specific_strain_peps - annotated_peptides))
            combined.loc[row[0], "Annotated specific peptides strain {}".format(strain)]='\n'.join(sorted(specific_strain_peps & annotated_peptides))
            combined.loc[row[0], "Reference BLAST strain {}".format(strain)]='\n'.join(strain_blasted)
            variant_dict[strain] = strain_exclusive_novel
            blasted = []
            features = []
            count=1
            
            for p in strain_exclusive_novel:
                if p in row_variant_blast_ref:
                    blasted += list(row_variant_blast_ref[p].values())
                    if p in row_variant_features_ref:
                        features += list(row_variant_features_ref[p].values())

            #combined.loc[row[0], "Exclusive peptides reference BLAST strain {}".format(strain)]='\n'.join(blasted) 
            #combined.loc[row[0], "Exclusive peptides polymorphism reference feature overlap strain {}".format(strain)]='\n'.join(features) 
            combined.loc[row[0], "Non-genomic peptides strain {}".format(strain)]='\n'.join(genome_unmapped) 
            combined.loc[row[0], "Translated orfs strain {}".format(strain)] = '\n'.join( [i.format('fasta') for i in trie])
            combined.loc[row[0], 'Variant orfs strain {}'.format(strain)] = '\n'.join(list(set([i.format('fasta') for i in annotation_trie])))
            combined.loc[row[0], 'Annotation type strain {}'.format(strain)] = '\n'.join(list(set(annotation_type)))
            
            
            ref_blast = sequtils.reference_mapping_blast(y, refps, output)
            
            if ref_blast[0] < min_eval_ref_blast:
                min_eval_ref_blast = ref_blast[0]
                best_eval_strain = strain

            combined.loc[row[0], "Best orf-reference blast evalue strain {}".format(strain)] = ref_blast[0]
            combined.loc[row[0], "Best orf-reference blast match {}".format(strain)] = ref_blast[1]

        for strain in id_holder:
            y = fasta_holder['{}_prots'.format(strain)]
            genome_unmapped = peptide_holder['{}_genome_unmapped_same'.format(strain)]
            other_unmapped =  peptide_holder['{}_genome_unmapped_other'.format(strain)]

            ids = id_holder[strain]
            if len(ids) > 1:
                reference_orfs = [i for i in mapped_orfs if not i.id in ids]
                fs_st = sequtils.frameshift_peptides(y, reference_orfs, genome_unmapped, output)
                fs_other = sequtils.frameshift_peptides(y, reference_orfs, other_unmapped, output)
            
                if len(fs_st.icds_pairs) or len(fs_other.icds_pairs) > 0:
                    combined.loc[row[0], 'Strain ICDS intrastrain evidence {}'.format(strain)] = fs_st.frameshift_report 
                    combined.loc[row[0], 'Strain ICDS interstrain exclusive evidence {}'.format(strain)] = fs_other.frameshift_report 
                    
                    intrastrain_evidence = '-'
                    interstrain_evidence = '-'
                    if len(fs_st.frameshift_peptides) > 0:
                        intrastrain_evidence = '+'
                    if len(fs_other.frameshift_peptides) > 0:
                        interstrain_evidence = '+'

                    combined.loc[row[0], 'Strain ICDS intrastrain peptide validation {}'.format(strain)] = intrastrain_evidence
                    combined.loc[row[0], 'Strain ICDS interstrain exclusive peptide validation {}'.format(strain)] = interstrain_evidence

                    frameshifts.append(strain) 
                    
                    

        #combined.loc[row[0], "All ORF IDs"] = ';'.join(sorted(all_orf_ids))

        #####################
        # STRAIN VARIANTS   #
        #####################
        #var_dict = defaultdict(list)
        #orf_dict = {}
        
        for strain in row_specific_strain_peps:
            strain_variants = set()
            specific_peps = row_specific_strain_peps[strain]
            strain_variant_features_orfs_list = []
            strain_variant_blast_orfs_list = []
            for comp_strain in row_specific_strain_peps:
                comp_peps = strain_sets[comp_strain]
                var_peps = specific_peps - comp_peps
                comp_orfs = fasta_holder['{}_prots'.format(comp_strain)]
                strain_variant_features_orfs, strain_variant_blast_orfs = peptide_list_blast(var_peps, comp_orfs, orf_features)
                for pep in strain_variant_blast_orfs:
                    strain_variants.add(comp_strain)
                    strain_variant_blast_orfs_list += list(strain_variant_blast_orfs[pep].values())
                    if pep in strain_variant_features_orfs:
                        strain_variant_features_orfs_list += list(strain_variant_features_orfs[pep].values()) 
            
            combined.loc[row[0], "Variant peptide BLAST strain {}".format(strain)]='\n'.join(strain_variant_blast_orfs_list) 
            combined.loc[row[0], "Variant peptide feature overlap strain {}".format(strain)]='\n'.join(strain_variant_features_orfs_list) 
            combined.loc[row[0], "Variant peptide strains strain {}".format(strain)]='\n'.join(strain_variants) 
            

        ####################### 
        # ORF feature overlap #
        #######################
        
        #row_variant_features_orfs, row_variant_blast_orfs = peptide_list_blast(row_specific, mapped_orfs, orf_features)
        #row_variant_features_orfs_list = []
        #row_variant_blast_orfs_list = []
        #for pep in row_variant_blast_orfs:
        #    row_variant_blast_orfs_list += list(row_variant_blast_orfs[pep].values())
        #    if pep in row_variant_features_orfs:
        #        row_variant_features_orfs_list += list(row_variant_features_orfs[pep].values()) 
        #blasted = []
        #features = []
        #for p in row_specific:
        #    if p in row_variant_blast_orfs:
        #        blasted += list(row_variant_blast_orfs[p].values())
        #        if p in row_variant_features_orfs:
        #            features += list(row_variant_features_orfs[p].values())
        
        #combined.loc[row[0], "Specific peptide variants"]='\n'.join(blasted) 
        #combined.loc[row[0], "Identified polymorphisms feature overlap"]='\n'.join(features) 
        combined.loc[row[0], "ICDS"] = '\n'.join(list(frameshifts))
        combined.loc[row[0], "Exclusive peptide strains"] = '\n'.join(list(strains_exclusive))
        combined.loc[row[0], "Specific peptides - all strains"] ='\n'.join(list(row_specific))
        combined.loc[row[0], "Specific annotated peptides - all strains"] ='\n'.join(row_annotated)
        combined.loc[row[0], "Specific novel peptides - all strains"] ='\n'.join(row_novel)
        combined.loc[row[0], "ORF ids - all strains"] = '\n'.join(rec.id for rec in mapped_orfs)
        combined.loc[row[0], "Reference proteins mapped - all strains"]= '\n'.join([i.format('fasta') for i in reftrie])
        combined.loc[row[0], "Reference proteins mapped count - all strains"]=len(mapped_reference)
        combined.loc[row[0], "Reference entries mapped - all strains"] = ';'.join([i.split('|')[1] for i in mapped_reference])
        combined.loc[row[0], "Variant peptides reference blast - all strains"]='\n'.join(row_variant_blast_ref_list)
        combined.loc[row[0], "Variant peptides reference features - all strains"]='\n'.join(row_variant_features_ref_list)
        combined.loc[row[0], "Non-tryptic nterm peptides - all strains"] = '\n'.join(row_non_tryptic_nterm)
        combined.loc[row[0], "Non-tryptic cterm peptides - all strains"] = '\n'.join(row_non_tryptic_cterm)
        combined.loc[row[0], "Non-ATG Methionine start peptides - all strains"] = '\n'.join(row_non_atg_starts)
        combined.loc[row[0], "Reference BLAST - all strains"] = '\n'.join(combined_blasted)
        
        if len(mapped_reference) > 0:
            fs_ref = sequtils.frameshift_peptides(refps, mapped_orfs,  row_novel, output)
            if len(fs_ref.icds_pairs) > 0: 
                combined.loc[row[0], 'Reference ICDS evidence'] = fs_ref.frameshift_report
                if len(set(fs_ref.frameshift_peptides)) > 0:
                    combined.loc[row[0], 'Reference ICDS peptide validation'] = '+' 

        ##########################
        # Mapping of protein ids #
        ##########################
        #protein_ids =[  i.split('|')[1] for i in mapped_reference ]
        protein_ids = combined_blasted
        rowGI = []
        rowGene_OrderedLocusName= []
        rowGeneID = []
        rowUPKBID = []
        
        for i in protein_ids:
            #print(idmapping[i])
            try:
                gis = idmapping[i]['GI']
                rowGI += gis
            except:
                pass

            try:
                goln = idmapping[i]['Gene_OrderedLocusName']
                rowGene_OrderedLocusName += goln
            except:
                pass

            try:
                gid = idmapping[i]['GeneID']
                rowGeneID += gid
            except:
                pass

            try:
                upkbid = idmapping[i]['UniProtKB-ID']
                rowUPKBID += upkbid
            except:
                pass

        identifier = rowGene_OrderedLocusName.copy()
        identifier.append('(protein group {})'.format(str(row[0])))
        combined.loc[row[0], 'Identifier'] = ' '.join(identifier)
        
        combined.loc[row[0], 'Gene_OrderedLocusName'] = ';'.join(rowGene_OrderedLocusName)
        combined.loc[row[0], 'GI'] = ';'.join(rowGI)
        combined.loc[row[0], 'GeneID'] = ';'.join(rowGeneID)
        combined.loc[row[0], 'UniProtKB-ID'] = ';'.join(rowUPKBID)

        rowOperons = []

        #for gi in rowGI:
        #    gi = str(gi)
        #    if gi in gi_operons:
        #        rowOperons.append(str(int(gi_operons[gi])))

        #combined.loc[row[0], 'DOOR2_Operons'] = ';'.join(rowOperons)

    #for strain in config['strains']:
    #    scols = strain_columns[strain]
    #    scols = [i for i in scols if i in combined.columns]
    #    sdf = combined[scols]
    #    if not os.path.exists(output +'/strains/{}'.format(strain)):
    #        os.mkdir(output +'/strains/{}'.format(strain))
    #
    #    print(sdf[:1].stack())
    #    
    #    specific_peptides = sdf['Specific peptides strain {}'.format(strain)].apply(lambda x : x.split('\n'))
    #    specific_peptides  = list(set([val for sublist in specific_peptides for val in sublist if val != '']))
    #
    #    novel_specific_peptides = sdf['Novel specific peptides strain {}'.format(strain)].apply(lambda x : x.split('\n'))
    #    novel_specific_peptides  = list(set([val for sublist in novel_specific_peptides for val in sublist if val != '']))
    #
    #    annotated_specific_peptides = sdf['Annotated specific peptides strain {}'.format(strain)].apply(lambda x : x.split('\n'))
    #    annotated_specific_peptides  = list(set([val for sublist in annotated_specific_peptides for val in sublist if val != '']))
    #        
    #    w = open(output +'/strains/{}/{}_novel_specific_peptides.txt'.format(strain, strain),'w')
    #    w.write('\n'.join(novel_specific_peptides))
    #    w.close()
    #    w = open(output +'/strains/{}/{}_annotated_specific_peptides.txt'.format(strain, strain),'w')
    #    w.write('\n'.join(annotated_specific_peptides))
    #    w.close()
    #    w = open(output +'/strains/{}/{}_specific_peptides.txt'.format(strain, strain),'w')
    #    w.write('\n'.join(specific_peptides))
    #    w.close()
    #    
    #    pgdf = sdf[sdf['Specific peptides strain {}'.format(strain)].notnull()]
    #    pgs = pgdf['Identifier'].tolist()
    #    w = open(output +'/strains/{}/{}_proteinGroups_list.txt'.format(strain, strain),'w')
    #    w.write('\n'.join(pgs))
    #    w.close()
    #
    #    sdf.to_csv(output +'/strains/{}/{}.csv'.format(strain, strain))
    #    #break

    #combined = combined.sort_values('iBAQ', ascending=False).drop_duplicates('All ORF IDs', keep='first')
    
    combined.to_csv(output+'/{}_combined.csv'.format(reference), sep='\t')
    
    #break
    # Save the start codons found
    #orf_start_codons = json.dumps(orf_start_codons)
    annotation_path = output + '/annotations/'
    if not os.path.exists(annotation_path):
        #shutil.rmtree(annotation_path)
        os.mkdir(annotation_path)

    for strain in peptide_annotations:
        annotations = pd.DataFrame(peptide_annotations[strain])
        annotations.to_csv(annotation_path +'/{}_{}_peptide_annotations.csv'.format(reference,strain))
        strain_var_fasta = var_fasta[strain]
        SeqIO.write(strain_var_fasta, annotation_path+'/{}_{}_variants.fasta'.format(reference,strain),'fasta')
        strain_trie_fasta = trie_fasta[strain]
        SeqIO.write(strain_trie_fasta, annotation_path+'/{}_{}_variants_trie.fasta'.format(reference,strain),'fasta')

    #with open(output +'/{}_start_codons_found.json'.format(reference),'w') as w:
    #    w.write(orf_start_codons)
    

