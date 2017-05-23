#!/bin/bash

set -x

#if the first column is empty it prints the message "Please provide a config"
if [ "$1" == "" ]
then
  echo "Please provide a config"
  echo
  exit
else
 config=`readlink -m $1`
  . $config
fi

in_dir=$uparse_merge_out_dir
out_dir=$uparse_filter_fastq_out_dir
log_dir=$uparse_filter_fastq_log_dir

if [ ! -d $out_dir ];
 then
   mkdir -p $out_dir
fi

log_dir=$log_dir"/uparse_filter_fastq."`date +"%y%m%d%H%M%S"`

if [ ! -d $log_dir ];
 then
   mkdir -p $log_dir
fi

count=1

while read sid_fastq_pair;
do 
  sid=`echo $sid_fastq_pair | awk -F ' ' '{print $1}'` # It is actually not reading it as a TAB. Maybe the read operation replaces the TABs with spaces.
#KL: the above line gets the sample ID from the sid_fastq_pair file
  merged_fastq=$in_dir/$sid.merged.fastq #here the sid created in the above line is concatenated with .merged.fastq at the end (a '.' seperates variables from other stuff)
  filtered_1_fastq=$out_dir/$sid.merged.filtered_1.fastq #this is for the first filtering step (max_e filtering) before primer stripping
  cmds_log=$log_dir/uparse_filter_fastq.$sid.$count.cmds

    #KL changed to:
  qsub="qsub -N uparse_filter_fastq.$sid.$count -o $log_dir/uparse_filter_fastq.$sid.$count.o -e $log_dir/uparse_filter_fastq.$sid.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$uparse_filter_fastq_fastq_threads -l walltime=$uparse_filter_fastq__walltime -v config=$config,merged_fastq=$merged_fastq,filtered_1_fastq=$filtered_1_fastq,fastq_maxee=$uparse_filter_fastq_maxee,cmds_log=$cmds_log $scripts_dir/uparse_filter_fastq.single.sh"
  echo $qsub > $log_dir/uparse_filter_fastq.$sid.$count.qsub
  cat $scripts_dir/uparse_filter_fastq.single.sh > $log_dir/uparse_filter_fastq.$sid.$count.sh

  if [ $uparse_filter_fastq_DEBUG -eq 1 ]
   then
     echo $qsub
   else
     jobid=`$qsub`
     echo $jobid
     echo $jobid > $log_dir/uparse_filter_fastq.$sid.$count.jobid
   fi

   (( count+=1 ))


done < $sid_fastq_pair_list

