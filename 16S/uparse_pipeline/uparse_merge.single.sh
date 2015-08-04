#!/bin/bash
. $config

echo "" > $cmds_log

# All stdout and stderr will go to individual per job *.o and *.e, so will not be piping to a log file.
cmd="$uparse_base/usearch7 -fastq_mergepairs $fastq_r1 -reverse $fastq_r2 -fastq_maxdiffs $fastq_maxdiffs -fastqout $merged_fastq"

echo $cmd >> $cmds_log
`eval $cmd`
