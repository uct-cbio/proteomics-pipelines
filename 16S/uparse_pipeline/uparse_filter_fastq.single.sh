#!/bin/bash
. $config

echo -n "" > $cmds_log

# All stdout and stderr will go to individual per job *.o and *.e, so will not be piping to a log file.
#Quality filter based on E scores
if [ ! -f $filtered_1_fastq ]
then
  cmd="$uparse_base/usearch9 -fastq_filter $merged_fastq  -fastq_maxee $fastq_maxee -fastqout $filtered_1_fastq"
  echo $cmd >> $cmds_log
  eval $cmd
fi

#Now strip primers

#Now do absolute length truncation
