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

qsub="qsub -N fastqc_combine.$count -M gerrit.botha@uct.ac.za -m abe  -o $log_dir/fastqc_combine.$count.o -e $log_dir/fastqc_combine.$count.e -d $out_dir -q UCTlong -S /bin/bash -l nodes=1:series600:ppn=$fastqc_combine_threads -v config=$config ./fastqc_combine.single.sh"

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
