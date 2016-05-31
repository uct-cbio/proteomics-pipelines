#!/bin/bash
. $config

echo -n "" > $cmds_log

# All stdout and stderr will go to individual per job *.o and *.e, so will not be piping to a log file.
#Quality filter based on E scores
if [ ! -f $filtered_1_fasta ]
then
  cmd="$uparse_base/usearch7 -fastq_filter $merged_fastq  -fastq_maxee $fastq_maxee -fastaout $filtered_1_fasta"
  echo $cmd >> $cmds_log
  eval $cmd
fi

#Now strip primers
if [ ! -f $filtered_2_fasta ]
then
  module add python/anaconda-python-2.7
  cmd="python /opt/exp_soft/qiime/packages/other/primer_stripping/strip_primers_hex.py $mapping_file $filtered_1_fasta $filtered_2_fasta $primer_strip_log" 
  echo $cmd >> $cmds_log
  eval $cmd
fi

#Now do absolute length truncation
if [ ! -f $filtered_3_fasta ]
then
  cmd="/opt/exp_soft/qiime/python-2.7.3/bin/python  /opt/exp_soft/qiime/packages/other/primer_stripping/truncate_seq_lens.py $filtered_2_fasta 250 260 250 $filtered_3_fasta"
  echo $cmd >> $cmds_log
  eval $cmd
fi
