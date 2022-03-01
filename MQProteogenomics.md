
# MQProteogenomics

MQProteogenomics is an open-source proteogenomics pipeline for multistrain genome annotation, variant and frameshift detection between strains and relative to a reference strain.

### 1. Clone the repository
`git clone https://github.com/uct-cbio/proteomics-pipelines.git`

#### 2. MQProteogenomics with Singularity (Recommended)
Ensure singularity is installed, and allocate at least 2 cores and 4 GB of RAM to the docker engine for this example.

~~~
cd proteomics-pipelines/singularity/mqproteogenomics
./create_image.sh
~~~
