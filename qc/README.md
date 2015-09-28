### How to run
1. First create your config file
2. Then run the fastqc batch script:
```/home/gerrit/projects/cbio-pipelines/qc/fastqc.batch.sh /home/gerrit/projects/cbio-pipelines/qc/config.txt```
3. Once 2 has completed reports can then be combined: ```/home/gerrit/projects/cbio-pipelines/qc/fastqc_combine.batch.sh /home/gerrit/projects/cbio-pipelines/qc/config.txt```

The scripts can be called from any directory. The base path of the pipeline scripts needs to be defined in the config.txt ($scripts_dir).

```sh
$ pwd
/home/gerrit/projects/uct/medmicbio/project08/qc
$ ls
config.txt  fastq.list
$ /home/gerrit/projects/cbio-pipelines/qc/fastqc.batch.sh config.txt
```

The test data that is used in the fastq.list contains
* 2 BEI control samples
* 1 sample with a mixture of Staph species
* 20 stool samples

The directory /scratch/DB/bio/test_data/16S/raw/mix with the raw data is readable to all on the hex cluster. 

The example output generated from the scripts in this repository are accesible in /scratch/DB/bio/test_data/16S/process. See output paths in config.txt.

If you would want to run the scripts on the test data, use the config.txt as an example but change the output paths to where you have write permissions to.


