# MQMetaproteomics Pipeline

A pipeline for metaproteomics analysis of MaxQuant output.

1. First create the singularity image:

~~~
cd proteomics-pipelines/singularity/mqmetaproteomics
./create_image.sh
~~~

This pulls the latest version from DockerHub, and saves it to the current directory.

2. Create an interactive shell with all the required dependencies
~~~
singularity shell mqmetaproteomics_v2.6.img
~~~

3. Clone the cbio-proteogenomics-tests repo
~~~
git clone https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git
~~~

4. Go to the MaxQuant example output ('txt') file included int the tests repo

~~~
cd cbio-proteogenomics-tests/lib/testdata/example_mq
~~~

5. Look at the example config file
~~~ 
nano singularity_config.yaml
~~~

6. Run the mqmetaproteomics pipeline
~~~
mqmetaproteomics.sh singularity_config.yaml
~~~

The output will appear as per the 'outdir' parameter specified in the config file.

