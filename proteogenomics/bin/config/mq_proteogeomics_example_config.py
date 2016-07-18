#!/usr/bin/env python

# Threads for concurrency

threads=50

# Proteogenomics
translation_table=11

# Specify the full path to the maxquant output 'txt' file
mq_txt='/researchdata/fhgfs/ptgmat003/BLACKBURNLAB/HH_PROTEOGENOMICS/17_SIX_FRAME_3_missed/combined/txt/'

# Specify the reference genome to be used to identify possible new annotations (corresponding to the Reference proteome)
reference_genome='/researchdata/fhgfs/ptgmat003/HYPOHYPER/NCBI_2_6_2016/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna'

# UniProt reference proteome to use for sequence mapping (Reference proteome or pan proteome - should not contain redundant sequences - this is for functional analysis of sequences across strains)
reference_proteome='/researchdata/fhgfs/ptgmat003/HYPOHYPER/proteomes/UP000001584_14_03_2016.fasta'

# UniProt mapping table (required fields: 'Organism ID','Sequence Length'), for functional and othr annotation of BLAST results (select appropiate taxa level eg. taxonomy:"Mycobacterium [1763]")
# The accompanying fasta is used to create the BLAST db (the two files need to correspond, and are selected by changing the uniprot download format. Gene set enrichment analysis info using here.
mapping_databaase='/home/ptgmat003/hypohyper/S507_S5527_hexdata/tax_1763_export_16_07_2016.tab'
mapping_fasta='/home/ptgmat003/hypohyper/S507_S5527_hexdata/tax_1763_export_16_07_2016_can_iso.fasta'

# The combined fasta database from SF_Pipeline (protein database, translated six frame db with alternative start sites included 
six_frame_translated ='/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/six_frame_databases/S507_S5527_combined_allstarts_proteins.fasta'

# A dictionary of all the strains to include in the search, please include the Reference genome again if using a Reference strain.
# sf_genome: The genome used to generate the six frame database, such as de novo assembly, can be the same as the variant_genome or else the Reference genome
# variant_genome: A modified Reference genome based on VCF data, specify 'None' here if none is available.
# orfs_raw: the stop-to-stop six frame nucleotide sequences
# orfs_trans: the translated six frame sequences with alternative posible TSS sites included
strains={ 'S5527' : {'sf_genome':'/researchdata/fhgfs/ptgmat003/HYPOHYPER/genomes/S5527_comb_assmbly_18_03_16.fasta',
                      'variant_genome':'/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/vcf_variants/S5527_variant.fasta',
                      'orfs_raw'  :'/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/six_frame_databases/S5527_sixframes.fasta',
                      'orfs_trans':'/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/six_frame_databases/S5527_allstarts_proteins.fasta'},

            'S507' : {'sf_genome':'/researchdata/fhgfs/ptgmat003/HYPOHYPER/genomes/S507_comb_assmbly_18_03_16.fasta',
                      'variant_genome':'/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/vcf_variants/S507_variant.fasta',
                      'orfs_raw'  :'/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/six_frame_databases/S507_sixframes.fasta',
                      'orfs_trans':'/researchdata/fhgfs/ptgmat003/OUTPUT/S507_S5527_proteogenomics/six_frame_databases/S507_allstarts_proteins.fasta'},
        }
# A dictionary of samples included in the experiment, with strain_name values corresponding to the keys in the strains dictionary above
# an error will be generated if there is a mismatch. Sample will correspond to the experiment name in MQ given when entering the experimental design template
# for readability make the names as concise as possible
# check that you have your strains right!!
samples={ '507_ML_1':     {'STRAIN':'S507', 'GROUP':'S507ML'},
          '507_ML_2':     {'STRAIN':'S507', 'GROUP':'S507ML'},
          '507_ML_3':     {'STRAIN':'S507', 'GROUP':'S507ML'},
          '507_MLexp_1':  {'STRAIN':'S507', 'GROUP':'S507MLexp'},
          '507_MLexp_2':  {'STRAIN':'S507', 'GROUP':'S507MLexp'},
          '507_MLexp_3':  {'STRAIN':'S507', 'GROUP':'S507MLexp'},
          '507_ST_1':     {'STRAIN':'S507', 'GROUP':'S507ST'},
          '507_ST_2':     {'STRAIN':'S507', 'GROUP':'S507ST'},
          '5527_ML_1':    {'STRAIN':'S5527', 'GROUP':'S5527ML'},
          '5527_ML_2':    {'STRAIN':'S5527', 'GROUP':'S5527ML'},
          '5527_ML_3':    {'STRAIN':'S5527', 'GROUP':'S5527ML'},
          '5527_MLexp_1': {'STRAIN':'S5527', 'GROUP':'S5527MLexp'},
          '5527_MLexp_3': {'STRAIN':'S5527', 'GROUP':'S5527MLexp'},
          '5527_ST_1':    {'STRAIN':'S5527', 'GROUP':'S5527ST'},
          '5527_ST_2':    {'STRAIN':'S5527', 'GROUP':'S5527ST'},
          }

# Decide on the comparisons based on LFQ and give the order eg. ['B','A']. FC will be calculated log2(B/A). A list of tuples is required.
# LIMMA will be used for LFQ, edgeR will be used for spectral counts
comparisons= [('S507ML','S507MLexp'),
              ('S507ML','S5527ML'),
              ('S5527ML','S5527MLexp'),
              ('S507MLexp','S5527MLexp'),
              ('S507ML','S507ST'),
              ('S5527ML','S5527ST'),
              ('S507ST','S5527ST')]


