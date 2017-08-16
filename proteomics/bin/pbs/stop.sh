#!/usr/bin/env bash

function stop {

     local log_file=$1

     if [ -f $log_file ]; then
         while IFS= read -r line
         do
             qdel $line
             qsig -s SIGINT $line
         done < $log_file
         rm $log_file
     fi
 }
