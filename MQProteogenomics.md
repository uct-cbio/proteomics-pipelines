
# MQProteogenomics

MQProteogenomics is an open-source proteogenomics pipeline for multistrain genome annotation, variant and frameshift detection between strains and relative to a reference strain.

### 1. Clone the repository
`git clone https://github.com/uct-cbio/proteomics-pipelines.git`

#### 2. MQProteogenomics with Singularity (Recommended)
Ensure singularity is installed, and make sure at least 2 cores and 4 GB of RAM are available.

~~~
cd proteomics-pipelines/singularity/mqproteogenomics
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
singularity shell mqproteogenomics_v2.0.4.img 
~~~

Run the example:
~~~
cd cbio-proteogenomics-tests/lib/testdata/proteogenomics_sample
mqproteogenomics.sh mq_proteogeomics_test.yml
~~~

Now adjust the config file for our own data.

#### 4. Notes
MaxQuant protein search database entries are mapped to STOP to STOP genome six frame translated open reading frame sequences for each strain. Alligned ORFs with a evalue of 0 are considered mapped. ORFs are mapped to UniProt proteomes using BLAST, with the top hit per ORF for each proteome selected with an evalue cutoff of 0.0001 used for annotation. The mapped proteins from the Reference proteome are used for functional analysis using limma t-test and gene set enrichment analysis.
