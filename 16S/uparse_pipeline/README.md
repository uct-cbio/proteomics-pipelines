### How to run
1. First create your config file (an example config and list file are included)
2. Then run the fastqc batch script:
```./uparse_merge.batch.sh config.txt```
3. Once 2. is done the merge reads can be filtered: ```./uparse_filter.batch.sh config.txt```


The sample_id.fastq.pairs.list file needs to have per line entries which contains the sample id and the full paths of the forward and reverse pairs. The pair entries needs to be seperated by a TAB. For example the file can be created like this:

paste -d '\t' <(for i in `ls -1 /scratch/DB/bio/test_data/16S/*R1.fastq | sort`; do basename=$(basename $i}); echo ${basename%_R*}; done) <(ls -1 /scratch/DB/bio/test_data/16S/*R1.fastq | sort) <(ls -1 /scratch/DB/bio/test_data/16S/*R2.fastq | sort) > sample_id.fastq_pair.list

