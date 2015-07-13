#!/bin/bash
. $config
$fastqc_combine_base/fastqc_combine.pl -v --out $out_dir --skip --files "$out_dir/*_fastqc"
