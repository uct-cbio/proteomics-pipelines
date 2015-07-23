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

echo $in_dir
echo $out_dir
echo $log_dir

log_dir=$log_dir"/fastqc_combine."`date +"%y%m%d%H%M%S"`

if [ ! -d $log_dir ];
 then
   mkdir -p $log_dir
fi

count=1

cmds_log=$log_dir/fastqc_combine.$count.cmds
qsub="qsub -N fastqc_combine.$count -M $pbs_status_mail_to -m $pbs_status_mail_events -o $log_dir/fastqc_combine.$count.o -e $log_dir/fastqc_combine.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$fastqc_combine_threads -l walltime=$fastqc_combine_walltime  -v config=$config,cmds_log=$cmds_log ./fastqc_combine.single.sh"

echo $qsub > $log_dir/fastqc_combine.$count.qsub
cat fastqc_combine.single.sh > $log_dir/fastqc_combine.$count.sh

if [ $fastqc_combine_DEBUG -eq 1 ]
then
  echo $qsub
  else
    jobid=`$qsub`
    echo $jobid
    echo $jobid > $log_dir/fastqc_combine.$count.jobid
fi
