
# MQProteogenomics

MQProteogenomics is an open-source proteogenomics pipeline for multistrain genome annotation, variant and frameshift detection between strains and relative to a reference strain.

### 1. Clone the repository
`git clone https://github.com/uct-cbio/proteomics-pipelines.git`

#### 2. MQProteogenomics with Singularity (Recommended)
Ensure singularity is installed, and make sure at least 2 cores and 4 GB of RAM are available.

~~~
cd proteomics-pipelines/singularity/mqproteogenomics
./create_image.sh
~~~

#### 3. Run the example 
Clone the example data:
~~~
git clone https://thys_potgieter@bitbucket.org/thys_potgieter/cbio-proteogenomics-tests.git
~~~

Create an interactive singularity session:
~~~
cd proteomics-pipelines/singularity/mqproteogenomics
singularity shell mqproteogenomics_v1.5.img 
~~~

Run the example:
~~~
cd cbio-proteogenomics-tests/lib/testdata/proteogenomics_sample
mqproteogenomics.sh mq_proteogeomics_test.yml
~~~
