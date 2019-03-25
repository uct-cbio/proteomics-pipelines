
# Test suite

Access the proteomics test suite and data from public repository available here:
https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git

# MetaNovo 

## MetaNovo with PBS on a cluster

### 1. Clone the repository

`git clone https://github.com/uct-cbio/proteomics-pipelines.git`

### 2. Install metanovo dependencies


The following dependencies will be installed (skip this step if they are already installed):


http://genesis.ugent.be/maven2/com/compomics/utilities/4.11.19/utilities-4.11.19.zip

http://www.proteoannotator.org/datasets/releases/ProteoAnnotator-1.7.86.zip

http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.2.20/SearchGUI-3.2.20-mac_and_linux.tar.gz

http://genesis.ugent.be/maven2/com/compomics/denovogui/DeNovoGUI/1.15.11/DeNovoGUI-1.15.11-mac_and_linux.tar.gz

https://cache.ruby-lang.org/pub/ruby/2.2/ruby-2.2.5.tar.gz


`cd proteomics-pipelines/bin/bash/`

`export METANOVO_DEPENDENCIES=${HOME}/software/ # Change the path if needed`

`./install_metanovo.sh`

### 3. Create the project folder 

Create a project folder to run MetaNovo. Please change to a directory for data storage on your cluster.

`cd .. # Change to the folder for data storage here `

`mkdir my_metanovo_project && cd my_metanovo_project` 

Pull an example fasta file from UniPro. Create your own by combining multiple species proteomes or use the whole of UniProt:

ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/


### 2. Create MetaNovo config file for your analysis
Copy the metanovo config file in bin/config/metanovo_config.sh to the project folder and edit

### 3. Define the PBS job parameters
3) Copy the pbs script in bin/pbs/metanovo.pbs to project folder and edit
4) Run the metanovo.pbs script "qsub metanovo.pbs", the job will end quickly, edit the X!tandem config files (default is ok)
5) Restart the script "qsub metanovo.pbs"
