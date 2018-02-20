#!/usr/bin/env bash


function submit {

     local cmd=$1;
     local previous=$2
     
     local deps=$( echo $( IFS=$':'; echo "${previous[*]}" ) )
     

     while (( $( qstat | grep $pbs_user | grep $q | wc -l ) >= $pbs_queue_limit ))
         do
             sleep 5
         done

     if [ -z "$previous" ]; then
        job=$(echo "$cmd" | qsub -P $P -q $q -l $l )
        echo $job
     else
        job=$(echo "$cmd" | qsub -P $P -q $q -l $l -W depend=afterok:${deps})
        echo $job;
     fi
 }
