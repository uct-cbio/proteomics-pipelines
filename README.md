
Test suite

Access the proteogenomics test suite and data from public repository available here::
https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git


This will be moved to the script for peptideshaker pipelines:
# To convert Q exactive Thermo RAW files to mgf
#$ msconvert.exe * --64 --zlib --filter "peakPicking true 1-" --mgf

## MetaNovo

1) Install metanovo dependencies, run the script: bin/bash/install_metanovo.sh 
2) Copy the metanovo config file in bin/config/metanovo_config.sh to the project folder and edit
3) Copy the pbs script in bin/pbs/metanovo.pbs to project folder and edit
4) Run the metanovo.pbs script "qsub metanovo.pbs", the job will end quickly, edit the X!tandem config files (default is ok)
5) Restart the script "qsub metanovo.pbs"
