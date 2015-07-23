#!/bin/bash
. $config

echo "" > $cmds_log

cmd="$fastqc_base/fastqc $fastq -f fastq -o $out_dir -t $fastqc_threads"

echo $cmd >> $cmds_log
`eval $cmd`

