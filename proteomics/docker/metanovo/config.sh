#####################
# Full paths please #
#####################

MGF_FOLDER="/home/thys/bio/tagmatch/sample"
FASTA_FILE="/home/thys/bio/tagmatch/UP000001584_83332.fasta"
#FASTA_FILE="/home/thys/bio/uniprot/uniprot_sprot.fasta"
OUTPUT_FOLDER='/home/thys/bio/tagmatch/sample_output'
CHUNKSIZE=200000 # size to split fasta
THREAD_LIMIT=4
JVM_Xmx=3500M
JVM_Xms=1024M

#######################
# MetaNovo parameters #
#######################

mn_specificity='specific'
mn_enzymes='Trypsin'
mn_max_missed_cleavages=3
mn_filter_database=1
mn_search_database=1
mn_search_fdr_value=1  
mn_search_fdr_level="accession"

######################
# DeNovoGUI Settings #
######################

dg_pepnovo=0
dg_pnovo=0
dg_novor=0
dg_directag=1

#################################
# Identification Parameters CLI #
#################################

# Spectrum matching parameters
prec_tol=0.02
prec_ppm=0
frag_tol=0.02
frag_ppm=0
digestion=0
enzyme='Trypsin'
specificity=0
mc=2
fixed_mods="Carbamidomethylation of C"
variable_mods="Oxidation of M, Acetylation of protein N-term"
min_charge=2
max_charge=4
fi='b'
ri='y'
min_isotope='0'
max_isotope='1'

# Spectrum Annotation
annotation_level=0.75
annotation_high_resolution=1

# Sequence Matching
sequence_index_type=0
sequence_matching_type=2
sequence_matching_x=0.25

# Import Filters
import_peptide_length_min=8
import_peptide_length_max=30
import_precursor_mz_ppm=0
exclude_unknown_ptms=1

# PTM Localization
ptm_score=1
score_neutral_losses=0
ptm_sequence_matching_type=1
ptm_alignment=1

# Gene Annotation
useGeneMapping=1
updateGeneMapping=1

# Protein Inference
simplify_groups=1
simplify_score=1
simplify_enzymaticity=1
simplify_evidence=1
simplify_uncharacterized=1

# Validation Levels
psm_fdr=1
peptide_fdr=1
protein_fdr=1
group_psms=1
group_peptides=1
merge_subgroups=1

# Fraction Analysis
protein_fraction_mw_confidence='95.0'

# XTandem advanced parameters
# MyriMatch advanced parameters
# MS Amanda advanced parameters
# MS-GF advanced parameters
# OMSSA advanced parameters
# Comet advanced parameters
# Tide advanced parameters

#PepNovo advanced parameters
pepnovo_hitlist_length=1
pepnovo_estimate_charge=1
pepnovo_correct_prec_mass=1
pepnovo_discard_spectra=1
pepnovo_fragmentation_model='CID_IT_TRYP'
pepnovo_generate_blast=0

#DirecTag
directag_tic_cutoff=85
directag_max_peak_count=400
directag_intensity_classes=3
directag_adjust_precursor=0
directag_min_adjustment='-2.5'
directag_max_adjustment='2.5'
directag_adjustment_step='0.1'
directag_charge_states='3'
directag_ms_charge_state='0'
directag_duplicate_spectra='1'
directag_deisotoping='0'
directag_isotope_tolerance='0.25'
directag_complement_tolerance='0.5'
directag_tag_length='3'
directag_max_var_mods='2'
directag_max_tag_count='30'
directag_intensity_weight='1.0'
directag_fidelity_weight='1.0'
directag_complement_weight='1.0'

#Novor
novor_fragmentation=HCD
novor_mass_analyzer=Trap


############
# MSGFPlus #
############

msgfplus_t="0.02Da"
msgfplus_ti="0,1"
msgfplus_tda=1
msgfplus_m=3 
msgfplus_inst=3 
msgfplus_e=1
msgfplus_protocol=0 
msgfplus_ntt=2 
msgfplus_minLength=6 
msgfplus_maxLength=40 
msgfplus_minCharge=2 
msgfplus_maxCharge=3 
msgfplus_n=1
msgfplus_addFeatures=0
msgfplus_ccm=1.00727649










