### How to run
1. First create your config file and mapping file with primer sequences (example config, mapping and list files are included)
2. Then run the fastqc batch script:
```/home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/uparse_merge.batch.sh /home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/config.txt```
3. Once 2 has completed the merge reads can be filtered: ```/home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/uparse_filter.batch.sh /home/gerrit/projects/cbio-pipelines/16S/uparse_pipeline/config.txt```
4. For the final UPARSE downstream steps run ```uparse_downstream.qsub```. This script is not written in the same fashion as the other scipts in this folder. It is up to you to change the parameter and database setting in the script it selve (it is not driven by a config file). There is also a qiime activation script that needs to be downloaded from here: http://hex.uct.ac.za/~gerrit/software/qiime/hex/activate_qiime.sh , replace the section ```source /home/gerrit/activate_qiime.sh``` with where you will be saving the file. The file just loads paths necessary for some of the QIIME steps (it actually needs to be revised if QIIME activation is needed in the UPARSE step).

The sid.fastq.pairs.list file needs to have per line entries which contains the sample id and the full paths of the forward and reverse pairs. The pair entries needs to be seperated by a TAB. For example the file can be created like this:

```paste -d '\t' <(for i in `ls -1 /scratch/DB/bio/test_data/16S/raw/mix/*R1_001.fastq | sort`; do basename=$(basename $i}); echo ${basename%_R*}; done) <(ls -1 /scratch/DB/bio/test_data/16S/raw/mix/*R1_001.fastq | sort) <(ls -1 /scratch/DB/bio/test_data/16S/raw/mix/*R2_001.fastq | sort) > sid.fastq_pair.list```

Something extra

- Currenlty the filter step creates .fasta files so that it can me stripped from primers and trimmed```. If it is however .fastq files we can then run FastQC on them and see what the effect of the UPARSE QC was. For now I've created a script that just create the .fastq files (after filtering). We later need to add the strip primer and trim steps, once that is replaced we can replace ```uparse_filter.batch.sh``` and ```uparse_filter.single.sh```. To run ```uparse_filter.fastq.batch.sh config.txt```

