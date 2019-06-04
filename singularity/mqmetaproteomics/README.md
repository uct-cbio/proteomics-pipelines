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

3. Export the path to the folder where interproscan.sh is located.

~~~
PATH=$PATH::/users/ptgmat003/ceph/cbio/users/ptgmat003/interproscan/5.35-74.0/interproscan-5.35-74.0
~~~

4. Clone the cbio-proteogenomics-tests repo
~~~
git clone https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git
~~~

5. Go to the MaxQuant example output ('txt') file included int the tests repo

~~~
cd cbio-proteogenomics-tests/lib/testdata/example_mq
~~~

6. Look at the example config file
~~~ 
nano singularity_config.yaml
~~~

7. Run the mqmetaproteomics pipeline
~~~
mqmetaproteomics.sh singularity_config.yaml
~~~

