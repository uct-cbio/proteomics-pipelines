# MetaNovo 
A probabilistic database export tool for peptide identification of complex mass spectrometry data.
## MetaNovo with PBS on a cluster
### 1. Clone the repository
`git clone https://github.com/uct-cbio/proteomics-pipelines.git`
### 2. Install metanovo dependencies
#### 2.1 The following dependencies will be installed (skip this step if they are already installed):
##### http://genesis.ugent.be/maven2/com/compomics/utilities/4.11.19/utilities-4.11.19.zip
##### http://www.proteoannotator.org/datasets/releases/ProteoAnnotator-1.7.86.zip
##### http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz
##### http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz
##### https://cache.ruby-lang.org/pub/ruby/2.2/ruby-2.2.5.tar.gz
~~~~
cd proteomics-pipelines/bin/bash/
export METANOVO_DEPENDENCIES=${HOME}/software/ # Change the path if needed
./install_metanovo.sh
~~~~
#### 2.2 Install python dependencies (if not available on your cluster)
Using a dedicated python virtual environment is highly recommended.
~~~~
cd proteomics-pipelines/lib 
pip3 install -r requirements.txt
~~~~
#### 2.3 Install R dependencies (Optional to enable post-processing with !X Tandem)
Required dependencies can be installed with a script:
~~~~
cd proteomics-pipelines/bin/R
./install_R_modules.R
~~~~
#### 2.4 Install UniPept (Optional, to allow for taxonomic characterization of !X Tandem results)
https://unipept.ugent.be/clidocs (MetaNovo is tested with UniPept version 1.1.1)
### 3. Create the project folder 
Create a project folder to run MetaNovo. Please change to a directory for data storage on your cluster.
~~~~
cd .. # Change to the folder for data storage here 
mkdir my_metanovo_project && cd my_metanovo_project
~~~~
Pull an example FASTA file from UniProt. Create your own by combining multiple species proteomes or use the whole of UniProt. For our example we will use only curated sequences avaiable in SwissProt.
~~~~
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
~~~~
Obtain some example mgf files from PRIDE (or use your own!). We will use files obtained from patients with colorectal cancer (https://www.ebi.ac.uk/pride/archive/projects/PXD00046).
~~~~
mkdir mgf_files && cd mgf_files
wget https://www.ebi.ac.uk/pride/data/archive/2014/08/PXD000467/YJC_100327SNOCRC_B11b372_N01.mgf
wget https://www.ebi.ac.uk/pride/data/archive/2014/08/PXD000467/YJC_100327SNOCRC_B11b372_N02.mgf
~~~~
### 4. Create MetaNovo config file for your analysis
Copy the metanovo config file in bin/config/metanovo_config.sh to the project folder and edit.
~~~~
cd ..
cp proteomics-pipelines/bin/config/metanovo_config.sh my_metanovo_project/config.sh
~~~~
Configure the full system paths to the data to be analyzed. A folder containing the MGF files, a FASTA file to search, and and output folder needs to be specified. Please modify the paths for your specific system and data storage requirements. Data files do not all need to be in the same location, but are done so here for simplicity.
~~~~
MGF_FOLDER=${HOME}/my_metanovo_project/mgf_files
FASTA_FILE=${HOME}/my_metanovo_project/uniprot_sprot.fasta
OUTPUT_FOLDER=${HOME}/my_metanovo_project
~~~~
Configure resource use parameters. The THREAD_LIMIT multiplied by JVM_Xmx should be less than the total RAM available per node, and THREAD_LIMIT should be one less than the total number of available cores. "Out of Memory" issues can be corrected by reducing CHUNKSIZE, assuming the THREAD_LIMIT and JVM_Xmx are correct.
~~~~
CHUNKSIZE=100000 # size to split fasta for paralellel processing
THREAD_LIMIT=2   # How many threads to use per node
JVM_Xmx=10000M   # Maximum memory allocated to each Java thread
JVM_Xms=1024M    # Minimum memory allocated to each Java thread
~~~~
Configure general parameters for the MetaNovo pipeline.
~~~~
mn_specificity='specific'      # specific|semi-specific|unspecific   (Enzyme specificity)
mn_enzymes='Trypsin'           # 'Trypsin, no P rule'|'Trypsin'|'Whole protein' (Enzyme rule)
mn_max_missed_cleavages=3      # Number of enzymatic missed cleavages
mn_filter_database=1           # Wether  to filter the database using MetaNovo algorithm (1=yes, 0=no)
mn_search_database=1           # Wether to run an !X Tandem search. if 'mn_filter_database=0', then the search will run against the original database without MetaNovo (1=yes, 0=no)
mn_prot_fdr_value=1            # Protein level FDR for !X Tandem post-processing
mn_pep_fdr_value=1             # Peptide level FDR for !X Tandme post-processing
~~~~
Configure wich sequencing engines to use with DeNovoGUI. Currently only DirecTag is supported by MetaNovo. http://compomics.github.io/projects/denovogui/wiki/denovocli.html
~~~~
dg_pepnovo=0
dg_pnovo=0
dg_novor=0
dg_directag=1
~~~~~
Configure the identification parameters for DenovoGUI. https://github.com/compomics/compomics-utilities/wiki/IdentificationParametersCLI
~~~~~
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
~~~~

### 3. Define the PBS job parameters
3) Copy the pbs script in bin/pbs/metanovo.pbs to project folder and edit
4) Run the metanovo.pbs script "qsub metanovo.pbs", the job will end quickly, edit the X!tandem config files (default is ok)
5) Restart the script "qsub metanovo.pbs"
## MetaNovo with Docker
Ensure docker is installed, and allocate at least 2 cores and 4 GB of RAM to the docker engine for this example.

# Run the tests!
~~~~
git clone https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git
~~~~
