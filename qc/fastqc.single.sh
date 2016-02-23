#!/bin/bash
. $config

echo -n "" > $cmds_log

cmd="$fastqc_base/fastqc --extract $fastq -f fastq -o $out_dir -t $fastqc_threads"

echo $cmd >> $cmds_log
eval $cmd
