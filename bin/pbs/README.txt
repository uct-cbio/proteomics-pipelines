
############
# METANOVO #
############

Use interactive mode or qsub on CHPC
1) Install software dependencies:
    ./cbio-pipelines/proteomics/bin/bash/install_all.sh
2) Install R dependencies:
    ./install_R_modules.R
3) Edit example pbs script (path to config file and software dependencies):
    nano cbio-pipelines/proteomics/bin/pbs/metanovo.pbs 
4) Edit example config file:
    nano cbio-pipelines/proteomics/bin/config/metanovo_test_config.sh
5) Submit the job:
    qsub metanovo.pbs



