######################
# Full paths to data #
######################

MGF_FOLDER=${HOME}/my_metanovo_project/mgf_files
FASTA_FILE=${HOME}/my_metanovo_project/uniprot_sprot.fasta
OUTPUT_FOLDER=${HOME}/my_metanovo_project
######################
# Processing Control #
######################

CHUNKSIZE=100000 # size to split fasta for paralellel processing
THREAD_LIMIT=2   # How many threads to use per node
JVM_Xmx=10000M   # Maximum memory allocated to each Java thread
JVM_Xms=1024M    # Minimum memory allocated to each Java thread

#######################
# MetaNovo parameters #
#######################

mn_specificity='specific'      # specific|semi-specific|unspecific   (Enzyme specificity)
mn_enzymes='Trypsin'           # 'Trypsin, no P rule'|'Trypsin'|'Whole protein' (Enzyme rule)
mn_max_missed_cleavages=3      # Number of enzymatic missed cleavages
mn_filter_database=0           # Wether  to filter the database using MetaNovo algorithm (1=yes, 0=no)
mn_search_database=1           # Wether to run an !X Tandem search. if 'mn_filter_database=0', then the search will run against the original database without MetaNovo (1=yes, 0=no)
mn_prot_fdr_value=1            # Protein level FDR for !X Tandem post-processing
mn_pep_fdr_value=1             # Peptide level FDR for !X Tandme post-processing

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
directag_tag_length='4'
directag_max_var_mods='2'
directag_max_tag_count='5'
directag_intensity_weight='1.0'
directag_fidelity_weight='1.0'
directag_complement_weight='1.0'

#Novor
novor_fragmentation=HCD
novor_mass_analyzer=Trap







