#!/bin/bash

if [ "$1" == "" ]
then
  echo "Please provide a config"
  echo
  exit
else
 config=`readlink -m $1`
  . $config
fi

echo $out_dir
echo $log_dir

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

  merged_fastq=$out_dir/$sid.merged.fastq
  filtered_fastq=$out_dir/$sid.merged.filtered.fasta

  echo $merged_fastq

  if [ ! -f $filtered_fastq ]
  then

    cmds_log=$log_dir/uparse_filter.$sid.$count.cmds
    qsub="qsub -N uparse_filter.$sid.$count -M $pbs_status_mail_to -m $pbs_status_mail_events  -o $log_dir/uparse_filter.$sid.$count.o -e $log_dir/uparse_filter.$sid.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$uparse_filter_threads -l walltime=$uparse_filter_walltime -v config=$config,merged_fastq=$merged_fastq,filtered_fastq=$filtered_fastq,fastq_trunclen=$uparse_filter_fastq_trunclen,fastq_maxee=$uparse_filter_fastq_maxee,cmds_log=$cmds_log ./uparse_filter.single.sh"

    echo $qsub > $log_dir/uparse_filter.$sid.$count.qsub
    cat uparse_filter.single.sh > $log_dir/uparse_filter.$sid.$count.sh

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

