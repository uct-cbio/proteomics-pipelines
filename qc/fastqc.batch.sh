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
echo $fastqc_base

if [ ! -d $out_dir ];
 then
   mkdir -p $out_dir
fi

log_dir=$log_dir"/fastqc."`date +"%y%m%d%H%M%S"`

if [ ! -d $log_dir ];
 then
   mkdir -p $log_dir
fi

count=1

while read fastq;
do 
  basename=$(basename $fastq)
  dirname=$(dirname $fastq)
  ext=${basename##*.}
  filename=${basename%.*}
  fastqc_output_zip=$out_dir/$filename"_fastqc.zip"

  echo $fastqc_output_zip

  if [ ! -f $fastqc_output_zip ]
  then
    cmds_log=$log_dir/fastqc.$filename.$count.cmds
    qsub="qsub -N fastqc.$filename.$count -M $pbs_status_mail_to -m $pbs_status_mail_events  -o $log_dir/fastqc.$filename.$count.o -e $log_dir/fastqc.$filename.$count.e -d $out_dir -q $pbs_queue -S /bin/bash -l nodes=1:$pbs_series:ppn=$fastqc_threads -l walltime=$fastqc_walltime -v config=$config,fastq=$fastq,cmds_log=$cmds_log ./fastqc.single.sh"

    echo $qsub > $log_dir/fastqc.$filename.$count.qsub
    cat fastqc.single.sh > $log_dir/fastqc.$filename.$count.sh

    if [ $fastqc_DEBUG -eq 1 ]
    then
      echo $qsub
    else
      jobid=`$qsub`
      echo $jobid
      echo $jobid > $log_dir/fastqc.$filename.$count.jobid
    fi
    
    (( count+=1 ))
   
   else
    echo "FastQC ouput for "$filename" exists already. Skipping this sample." 
   fi

done < $fastq_list


