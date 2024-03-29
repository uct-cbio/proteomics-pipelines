# MetaNovo 
MetaNovo is an open-source pipeline that integrates existing tools with a custom algorithm to produce targeted databases for mass spectrometry analysis. As inputs it requires raw mass spectrometry data in MGF format, and a UniProt FASTA file to search. 

MetaNovo uses open-source tools to match raw spectra to database entries in a parallalized and scalable pipeline that can be installed in a cluster or standalone on a linux machine. 

Database matches are stored in an SQLite database, and used to estimate the abundance of proteins using sequence tag matches, taking into account protein length using a normalized spectral abundance factor estimation for each protein in each sample. The MetaNovo algorithm uses the protein abundance level estimation to rank the data, producing a parsimonious list of protein identifiers that can explain all the database matches (such that each spectral match maps to at least one protein in the non-redundant list). Taxonomic representation in this list is calculated using the UniProt FASTA header "OS" entry, and a score for each organism is obtained. Database proteins are re-ranked based on the combined scores for spectral and organism abundance, and a database is exported.

#### Recomended use - MetaNovo with Singularity. For this you only need to follow steps 3,4 and 6 below.

## MetaNovo installation

### 1. Clone the repository

`git clone https://github.com/uct-cbio/proteomics-pipelines.git`

### 2. Install metanovo dependencies
#### 2.1 The following dependencies will be installed (skip this step if they are already installed):
##### http://genesis.ugent.be/maven2/com/compomics/utilities/4.12.0/utilities-4.12.0.zip
##### http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz
##### http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz

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

### 3. Create the project folder 
Create a project folder to run MetaNovo. Please change to a directory for data storage on your cluster.
~~~~
cd .. # Change to the folder for data storage here 
mkdir my_metanovo_project && cd my_metanovo_project
~~~~
Pull an example FASTA file from UniProt. Create your own by combining multiple species proteomes or use the whole of UniProt. For our example we will use only curated sequences available in SwissProt. For demonstration purposes, use only uniprot_sprot.fasta. For accurate results, please use more MGF (>10) files, and use the whole of UniProt. The suggested method is to concatenate UniProt trembl, sprot and sprot_varsplic databases into a single FASTA file.
~~~~
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz
gunzip *.fasta.gz
cat *.fasta > combined.fasta
~~~~
Obtain some example mgf files from PRIDE (or use your own!). We will use files obtained from https://www.ebi.ac.uk/pride/archive/projects/PXD029490. The more MGF files you use, the more accurate the database should be - but the longer the analysis may take.
~~~~
mkdir mgf_files && cd mgf_files
wget http://ftp.pride.ebi.ac.uk/pride/data/archive/2021/11/PXD029490/20180326_CAINO_01_band02.mgf
wget http://ftp.pride.ebi.ac.uk/pride/data/archive/2021/11/PXD029490/20180326_CAINO_01_band03.mgf
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
FASTA_FILE=${HOME}/my_metanovo_project/combined.fasta # or just use swissprot for demonstration only purpose. We recommend all UniProt
OUTPUT_FOLDER=${HOME}/my_metanovo_project/output
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
mn_enzymes='Trypsin, no P rule'           # 'Trypsin, no P rule'|'Trypsin'|'Whole protein' (Enzyme rule)
mn_max_missed_cleavages=2      # Number of enzymatic missed cleavages

~~~~
Configure wich sequencing engines to use with DeNovoGUI. Currently only DirecTag is supported by MetaNovo. http://compomics.github.io/projects/denovogui/wiki/denovocli.html
~~~~
dg_pepnovo=0
dg_pnovo=0
dg_novor=0
dg_directag=1
~~~~~
Configure the identification parameters for DenovoGUI and PeptideMapper. The DirecTag 'directag_max_tag_count' can be reduced to 1 for very large datasets. Support for other search engines is in progress. The default settings below should be sufficient for most applications.  https://github.com/compomics/compomics-utilities/wiki/IdentificationParametersCLI
~~~~~
# Spectrum matching parameters
prec_tol=0.02
prec_ppm=0
frag_tol=0.02
frag_ppm=0
digestion=0
enzyme='Trypsin (no P rule)'
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

# PepNovo advanced parameters
pepnovo_hitlist_length=1
pepnovo_estimate_charge=1
pepnovo_correct_prec_mass=1
pepnovo_discard_spectra=1
pepnovo_fragmentation_model='CID_IT_TRYP'
pepnovo_generate_blast=0

# DirecTag
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

# Novor
novor_fragmentation=HCD
novor_mass_analyzer=Trap
~~~~~
### 5. Define the PBS job parameters
#### 5.1 Copy the pbs script in bin/pbs/metanovo.pbs to project folder
~~~~
cd ..
cp proteomics-pipelines/bin/pbs/metanovo.pbs my_metanovo_project/metanovo.pbs
nano my_metanovo_project/metanovo.pbs
~~~~
#### 5.2 Define the PBS job settings 
The below are only examples, please adjust accordingly.
~~~~
#PBS -P CBBI0825 
#PBS -M matthys.potgieter@gmail.com 
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -l walltime=48:00:00
#PBS -N my_metanovo_example
#PBS -q smp
#PBS -W group_list=largeq
#PBS -m be
~~~~

#### 5.3 Load any required dependencies
eg.
~~~~
module add chpc/gnu/parallel-20160422
~~~~

#### 5.4 Define paths to software dependencies installed above
~~~~
export SG_PATH=${HOME}/software/SearchGUI/SearchGUI-3.2.20
export DG_PATH=${HOME}/software/DeNovoGUI/DeNovoGUI-1.15.11
export CU_PATH=${HOME}/software/utilities/utilities-4.12.0
~~~~

#### 5.5 Start the pipeline
If the pipeline fails at any step, simply restart and the pipeline will continue where it left off.
~~~~
qsub metanovo.pbs
~~~~

#### 6. MetaNovo with Singularity (Recommended)
Ensure singularity is installed, and allocate at least 2 cores and 4 GB of RAM to the docker engine for this example.

~~~
wget https://github.com/uct-cbio/proteomics-pipelines/archive/refs/tags/v1.9.4.zip                                                        
cd proteomics-pipelines-1.9.4/singularity/metanovo 
./create_image.sh
~~~
Proceed with Step 3 and 4 above.

~~~
singularity shell metanovo_v1.9.4.img
Singularity metanovo_v1.9.4.img:~/my_metanovo_project> metanovo.sh config.sh 
~~~
This process can be performed on local workstatiuon, and Singularity compatible High Performance clusters. 

The database is saved as <output_folder>/metanovo/metanovo.fasta

#### 7. MetaNovo with Docker
Proceed with Step 3 and 4 above. Then:
~~~
cd proteomics-pipelines/docker/metanovo 
./run.sh ${HOME}/my_metanovo_project/config.sh 
~~~

The database is saved as <output_folder>/metanovo/metanovo.fasta

# Run the tests!
~~~~
git clone https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git
~~~~
