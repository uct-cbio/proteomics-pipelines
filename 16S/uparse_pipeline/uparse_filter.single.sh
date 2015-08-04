#!/bin/bash
. $config

echo "" > $cmds_log

# All stdout and stderr will go to individual per job *.o and *.e, so will not be piping to a log file.
cmd="$uparse_base/usearch7 -fastq_filter $merged_fastq -fastq_trunclen $fastq_trunclen -fastq_maxee $fastq_maxee -fastaout $filtered_fastq"

echo $cmd >> $cmds_log
`eval $cmd`
