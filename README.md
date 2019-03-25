
# Test suite
Access the proteomics test suite and data from public repository available here:
https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git
# MetaNovo 
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
#### 2.3 Install R dependencies (Optional to enable post-processing with !XTandem)
Required dependencies can be installed with a script:
~~~~
cd proteomics-pipelines/bin/R
./install_R_modules.R
~~~~
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
### 2. Create MetaNovo config file for your analysis
Copy the metanovo config file in bin/config/metanovo_config.sh to the project folder and edit
### 3. Define the PBS job parameters
3) Copy the pbs script in bin/pbs/metanovo.pbs to project folder and edit
4) Run the metanovo.pbs script "qsub metanovo.pbs", the job will end quickly, edit the X!tandem config files (default is ok)
5) Restart the script "qsub metanovo.pbs"
## MetaNovo with Docker
Ensure docker is installed, and allocate at least 2 cores and 4 GB of RAM to the docker engine for this example.
