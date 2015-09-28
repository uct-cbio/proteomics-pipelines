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

out_dir=$uparse_merge_out_dir
log_dir=$uparse_merge_log_dir

tmp_dir=$tmp_dir
rename_dir=$uparse_merge_rename_dir

if [ ! -d $out_dir ];
 then
   mkdir -p $out_dir
fi

if [ ! -d $tmp_dir ];
 then
   mkdir -p $tmp_dir
fi

if [ ! -d $rename_dir ];
 then
   mkdir -p $rename_dir
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
  
  fastq_r1_tmp=$tmp_dir"/"$(basename $fastq_r1)
  fastq_r2_tmp=$tmp_dir"/"$(basename $fastq_r2)

  fastq_r1_renamed=$uparse_merge_rename_dir"/"$(basename $fastq_r1)
  fastq_r2_renamed=$uparse_merge_rename_dir"/"$(basename $fastq_r2)

  merged_fastq=$out_dir/$sid.merged.fastq

  if [ ! -f $merged_fastq ]
  then

    cmds_log=$log_dir/uparse_merge.$sid.$count.cmds
    qsub="qsub -N uparse_merge.$sid.$count -M $pbs_status_mail_to -m $pbs_status_mail_events  -o $log_dir/uparse_merge.$sid.$count.o -e $log_dir/uparse_merge.$sid.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$uparse_merge_threads -l walltime=$uparse_merge_walltime -v config=$config,merged_fastq=$merged_fastq,fastq_r1=$fastq_r1,fastq_r2=$fastq_r2,fastq_r1_tmp=$fastq_r1_tmp,fastq_r2_tmp=$fastq_r2_tmp,fastq_r1_renamed=$fastq_r1_renamed,fastq_r2_renamed=$fastq_r2_renamed,sid=$sid,fastq_maxdiffs=$uparse_merge_fastq_maxdiffs,cmds_log=$cmds_log $scripts_dir/uparse_merge.single.sh"

    echo $qsub > $log_dir/uparse_merge.$sid.$count.qsub
    cat $scripts_dir/uparse_merge.single.sh > $log_dir/uparse_merge.$sid.$count.sh

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

