1. First create the singularity image:

~~~
cd proteomics-pipelines/singularity/mqmetaproteomics
./create_image.sh
~~~

This pulls the latest version from DockerHub, and saves it to the current directory.

2. Create an interactive shell with all the required dependencies, excluding InterProScan, wich needs to be installed separately (https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload)

~~~
singularity shell mqmetaproteomics_v1.1.img
~~~
