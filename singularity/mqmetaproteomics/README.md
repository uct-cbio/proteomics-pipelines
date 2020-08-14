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
singularity shell mqmetaproteomics_v2.7.2.img
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

7. Contribute!

Download both repos eg:
~~~
${HOME}/repos/proteomics-pipelines
${HOME}/repose/cbio-proteogenomics-tests
~~~
Load the singularity environment (step 2 above)

To adjust the paths to the git repo, and overwrite the existing binaries in the singularity image, run:
~~~
source ${HOME}/repos/proteomics-pipelines/bin/bash/dev.sh
~~~
Run the python unittests:
~~~
cd ${HOME}/repose/cbio-proteogenomics-tests/lib/ 
python3 test_mqparse.py # verify that all tests pass
~~~
Add more unit-tests to test_mqparse.py and modify the source code in proteomics-pipelines. 

Submit a pull request




