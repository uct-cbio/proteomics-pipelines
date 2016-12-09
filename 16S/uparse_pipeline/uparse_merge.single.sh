#!/bin/bash
. $config

echo -n "" > $cmds_log

cmd="$fastx_base/fastx_renamer -Q33 -n COUNT -i $fastq_r1 -o $fastq_r1_tmp"
echo $cmd >> $cmds_log
eval $cmd

cmd="sed \"s/^\(@\|+\)\([0-9]*\)$/\1$sid\2;barcodelabel=$sid\/1/\" $fastq_r1_tmp > $fastq_r1_renamed"
echo $cmd >> $cmds_log
eval $cmd

cmd="$fastx_base/fastx_renamer -Q33 -n COUNT -i $fastq_r2 -o $fastq_r2_tmp"
echo $cmd >> $cmds_log
eval $cmd

cmd="sed \"s/^\(@\|+\)\([0-9]*\)$/\1$sid\2;barcodelabel=$sid\/2/\" $fastq_r2_tmp > $fastq_r2_renamed"
echo $cmd >> $cmds_log
eval $cmd

cmd="rm -f $fastq_r1_tmp $fastq_r2_tmp"

# All stdout and stderr will go to individual per job *.o and *.e, so will not be piping to a log file.
cmd="$uparse_base/usearch9 -fastq_mergepairs $fastq_r1_renamed -reverse $fastq_r2_renamed -fastq_maxdiffs $fastq_maxdiffs -fastqout $merged_fastq"

echo $cmd >> $cmds_log
eval $cmd
