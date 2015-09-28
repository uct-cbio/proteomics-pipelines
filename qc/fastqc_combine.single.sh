#!/bin/bash
. $config

echo -n "" > $cmds_log

cmd="$fastqc_combine_base/fastqc_combine.pl -v --out $out_dir --skip --files \"$out_dir/*_fastqc\""

echo $cmd >> $cmds_log
eval $cmd
