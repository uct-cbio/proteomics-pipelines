#!/bin/bash
. $config
$fastqc_base/fastqc $fastq -f fastq -o $out_dir -t $fastqc_threads
