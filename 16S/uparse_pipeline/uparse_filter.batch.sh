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
out_dir=$uparse_filter_out_dir
log_dir=$uparse_filter_log_dir

if [ ! -d $out_dir ];
 then
   mkdir -p $out_dir
fi

log_dir=$log_dir"/uparse_filter."`date +"%y%m%d%H%M%S"`

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
  #filtered_fasta=$out_dir/$sid.merged.filtered.fasta
  #KL changed to:
  filtered_1_fasta=$out_dir/$sid.merged.filtered_1.fasta #this is for the first filtering step (max_e filtering) before primer stripping
  filtered_2_fasta=$out_dir/$sid.merged.filtered_2.fasta #this is for the second filtering step: primer stripping
  filtered_3_fasta=$out_dir/$sid.merged.filtered_3.fasta #this is for the third filtering step: length trunctation
  
  #if [ ! -f $filtered_fasta ]
  #KL changed to:
  if [ ! -f $filtered_3_fasta ]
  then

    cmds_log=$log_dir/uparse_filter.$sid.$count.cmds
   # qsub="qsub -N uparse_filter.$sid.$count -M $pbs_status_mail_to -m $pbs_status_mail_events  -o $log_dir/uparse_filter.$sid.$count.o -e $log_dir/uparse_filter.$sid.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$uparse_filter_threads -l walltime=$uparse_filter_walltime -v config=$config,merged_fastq=$merged_fastq,filtered_fasta=$filtered_fasta,fastq_trunclen=$uparse_filter_fastq_trunclen,fastq_maxee=$uparse_filter_fastq_maxee,cmds_log=$cmds_log $scripts_dir/uparse_filter.single.sh"
    #KL changed to:
    qsub="qsub -N uparse_filter.$sid.$count -o $log_dir/uparse_filter.$sid.$count.o -e $log_dir/uparse_filter.$sid.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$uparse_filter_threads -l walltime=$uparse_filter_walltime -v config=$config,merged_fastq=$merged_fastq,filtered_1_fasta=$filtered_1_fasta,fastq_maxee=$uparse_filter_fastq_maxee,filtered_2_fasta=$filtered_2_fasta,mapping_file=$mapping_file,primer_strip_log=$log_dir/SP_log.$sid.$count.txt,filtered_3_fasta=$filtered_3_fasta,fastq_trunclen=$uparse_filter_fastq_trunclen,cmds_log=$cmds_log $scripts_dir/uparse_filter.single_KL_primer_strip.sh"
    echo $qsub > $log_dir/uparse_filter.$sid.$count.qsub
    cat $scripts_dir/uparse_filter.single_KL_primer_strip.sh > $log_dir/uparse_filter.$sid.$count.sh

   if [ $uparse_filter_DEBUG -eq 1 ]
    then
      echo $qsub
    else
      jobid=`$qsub`
      echo $jobid
      echo $jobid > $log_dir/uparse_filter.$sid.$count.jobid
    fi

    (( count+=1 ))

   else
    echo "UPARSE filter ouput for "$sid" exists already. Skipping this sample." 
   fi

done < $sid_fastq_pair_list

