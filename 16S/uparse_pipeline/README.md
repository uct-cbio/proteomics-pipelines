### How to run
1. First create your config file (an example config and list file are included)
2. Then run the fastqc batch script:
```/home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/uparse_merge.batch.sh /home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/config.txt```
3. Once 2 has completed the merge reads can be filtered: ```/home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/uparse_filter.batch.sh /home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/config.txt```

The sid.fastq.pairs.list file needs to have per line entries which contains the sample id and the full paths of the forward and reverse pairs. The pair entries needs to be seperated by a TAB. For example the file can be created like this:

```paste -d '\t' <(for i in `ls -1 /scratch/DB/bio/test_data/16S/raw/mix/*R1_001.fastq | sort`; do basename=$(basename $i}); echo ${basename%_R*}; done) <(ls -1 /scratch/DB/bio/test_data/16S/raw/mix/*R1_001.fastq | sort) <(ls -1 /scratch/DB/bio/test_data/16S/raw/mix/*R2_001.fastq | sort) > sid.fastq_pair.list```

