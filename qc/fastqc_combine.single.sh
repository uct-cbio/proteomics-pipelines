#!/bin/bash
. $config

echo "" > $cmds_log

$fastqc_combine_base/fastqc_combine.pl -v --out $out_dir --skip --files "$out_dir/*_fastqc"

echo $cmd >> $cmds_log
`eval $cmd
