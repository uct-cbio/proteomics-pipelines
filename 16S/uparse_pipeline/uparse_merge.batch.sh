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

log_dir=$log_dir"/uparse_merge."`date +"%y%m%d%H%M%S"`

if [ ! -d $log_dir ];
 then
   mkdir -p $log_dir
fi

count=1

while read sid_fastq_pair;
do 
  sid=`echo $sid_fastq_pair | awk -F ' ' '{print $1}'` # It is actually not reading it as a TAB. Maybe the read operation replaces the TABs with spaces.
  fastq_r1=`echo $sid_fastq_pair | awk -F ' ' '{print $2}'` # It is actually not reading it as a TAB. Maybe the read operation replaces the TABs with spaces.
  fastq_r2=`echo $sid_fastq_pair | awk -F ' ' '{print $3}'`
  echo $sid
  echo $fastq_r1 
  echo $fastq_r2 

  merged_fastq=$out_dir/$sid.merged.fastq

  echo $merged_fastq

  if [ ! -f $merged_fastq ]
  then

    cmds_log=$log_dir/uparse_merge.$sid.$count.cmds
    qsub="qsub -N uparse_merge.$sid.$count -M $pbs_status_mail_to -m $pbs_status_mail_events  -o $log_dir/uparse_merge.$sid.$count.o -e $log_dir/uparse_merge.$sid.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$uparse_merge_threads -l walltime=$uparse_merge_walltime -v config=$config,merged_fastq=$merged_fastq,fastq_r1=$fastq_r1,fastq_r2=$fastq_r2,fastq_maxdiffs=$uparse_merge_fastq_maxdiffs,cmds_log=$cmds_log ./uparse_merge.single.sh"

    echo $qsub > $log_dir/uparse_merge.$sid.$count.qsub
    cat uparse_merge.single.sh > $log_dir/uparse_merge.$sid.$count.sh

   if [ $uparse_merge_DEBUG -eq 1 ]
    then
      echo $qsub
    else
      jobid=`$qsub`
      echo $jobid
      echo $jobid > $log_dir/uparse_merge.$sid.$count.jobid
    fi

    (( count+=1 ))

   else
    echo "UPARSE merged ouput for "$sid" exists already. Skipping this sample." 
   fi

done < $sid_fastq_pair_list

