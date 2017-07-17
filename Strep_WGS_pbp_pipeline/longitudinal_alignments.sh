#PURPOSE OF THIS SCRIPT: 
#1. FOR EACH PATIENT COLLATE ALL LONGITUDINAL SAMPLES (FOR A GIVEN PBP) INTO ONE FILE

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

out_dir=$out_dir

if [ ! -d $out_dir ];
 then
   mkdir -p $out_dir
fi

tmp_dir=$tmp_dir

if [ ! -d $tmp_dir ];
 then
   mkdir -p $tmp_dir
fi

while read accession; #read in pbp gene accession
do
   ID=$accession #extract accession number from query.list 
   pid_record=`awk 'FNR==1 {print $2}' $pid_list` #initialize to first patient ID, record patient id (pid) use later as check in loop
   cmds_log=$log_dir/$ID.cmds
   echo -n "" > $cmds_log #in case file was made previously first clear it before continuing..   
   #****************************************
   #Initialisation for 1st patient in loop which has to be done here due to loop which runs over same patient multiple times for different samples (see below)
   #For each PBP gene lets make a summary file describing whether each of the samples for a given patient had the pbp gene seq we were looking for
   #just for the first patient to start loop (this is just to make output summary file more readable
   echo Sample details for PID $pid_record: > $out_dir/$ID.summary 
   #For first patient check that the longitudinal file has not been created in a previous run (the rest will be overwritten in the loop below if they have)
   echo -n "" > $out_dir/$pid_record.$ID.longitudinal.fasta #create file and empty in case it was previously created
   cmd="awk '/>/{p=0} /$ID/{p=1} p' /home/kviljoen/strep_pbp_project/KL_Strep_pneumoniae_NCBI_random_reference_pbp_seqs.fasta > $out_dir/$pid_record.$ID.longitudinal.fasta" #add the relevant PBP reference sequence to be aligned against
   echo $cmd >> $cmds_log
   echo $'\n' >> $cmds_log #add newline
   eval $cmd
   #****************************************
   while read file;
   do
   pid=`echo $file | awk -F ' ' '{print $2}'` #extract patient ID from pid_list (column 2)
   sid=`echo $file | awk -F ' ' '{print $1}'` #extract sample ID from pid_list (column 1) 
   #-----------------------------------
   #First we append sample IDs to fasta headers so we can keep track of where each entry comes from in the collated file later
   cmd="awk '/^>$ID/ {\$0=\$0 \"_sid_$sid\"}1' $out_dir/$sid.$ID.matches.fa > $out_dir/$sid.$ID.matches_2.fa"  
   echo $cmd >> $cmds_log
   echo $'\n' >> $cmds_log #add newline
   eval $cmd
   #Next we remove all the gap '-' characters from the fasta files that were introduced during BLAST alignments
   cmd="awk '{gsub(\"-\",\"\")}1' $out_dir/$sid.$ID.matches_2.fa > $out_dir/$sid.$ID.matches_3.fa" 
   echo $cmd
   echo $cmd >> $cmds_log
   echo $'\n' >> $cmds_log #add newline
   eval $cmd
   #-----------------------------------
   #make a a file to capture relevant seqs
   if [ $pid == $pid_record ]; #if the pid is still the same, keep collating to the same file
    then 
     cmd="cat $out_dir/$sid.$ID.matches_3.fa >> $out_dir/$pid.$ID.longitudinal.fasta"
     echo $cmd
     echo $cmd >> $cmds_log
     echo $'\n' >> $cmds_log #add newline
     eval $cmd
     if [ -s $out_dir/$sid.$ID.matches_3.fa ]; #this checks whether the files exists AND has content
     then
       echo Suitable contig found for sample $sid PBP $ID, will be added to longitudinal file >> $out_dir/$ID.summary
     else
       echo No suitable contigs found for sample $sid PBP $ID >> $out_dir/$ID.summary
     fi
   else #if we've moved on to the next pid
   echo Sample details for PID $pid: >> $out_dir/$ID.summary 
     #First we check if the fasta file in question has content, if not we leave a message in the output file
     if [ -s $out_dir/$sid.$ID.matches_3.fa ]; #this checks whether the files exists AND has content
     then
       echo Suitable contig found for sample $sid PBP $ID, will be added to longitudinal file >> $out_dir/$ID.summary
     else
       echo No suitable contigs found for sample $sid PBP $ID >> $out_dir/$ID.summary
     fi
     cmd="awk '/>/{p=0} /$ID/{p=1} p' $pbp_refs > $out_dir/$pid.$ID.longitudinal.fasta" #add the relevant PBP reference sequence to be aligned against test seqs. #note here the file will be overwritten if it already exists (> vs >>)
     echo $cmd >> $cmds_log
     echo $'\n' >> $cmds_log #add newline
     eval $cmd
     cmd="cat $out_dir/$sid.$ID.matches_3.fa >> $out_dir/$pid.$ID.longitudinal.fasta"
     echo $cmd >> $cmds_log
     echo $'\n' >> $cmds_log #add newline
     eval $cmd
   fi
   pid_record=$pid
   #((count+=1))
   done < $pid_list
done < $query_list




