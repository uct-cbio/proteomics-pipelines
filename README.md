
# Test suite

Access the proteogenomics test suite and data from public repository available here:
https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git

# MetaNovo 

## MetaNovo with PBS on a cluster

### 1 . Install metanovo dependencies
`cd proteomics-pipelines/bin/bash/`
`./install_metanovo.sh`

### Create MetaNovo config file for your analysis
Copy the metanovo config file in bin/config/metanovo_config.sh to the project folder and edit

### Create the PB
3) Copy the pbs script in bin/pbs/metanovo.pbs to project folder and edit
4) Run the metanovo.pbs script "qsub metanovo.pbs", the job will end quickly, edit the X!tandem config files (default is ok)
5) Restart the script "qsub metanovo.pbs"
