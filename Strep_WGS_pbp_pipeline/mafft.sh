#PBS -N mafft_alignment_test
#PBS -S /bin/bash
#PBS -q UCTlong
#PBS -l nodes=srvslshpc617:series600:ppn=4
#PBS -V
#PBS -d /home/kviljoen/strep_pbp_project/

#PURPOSE OF THIS SCRIPT: 
#1.Take the files for each individual with longitudinal pbp seqs and do multiples sequence alignments using MAFFT
#2. Do phylogenetic tree of 1.

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


#get list of longitudinal alignment files:
cmd="ls $out_dir/*longitudinal.fasta > $tmp_dir/longitudinal_fastas_filenames.list"
cmds_log=$log_dir/mafft_alignments.cmds
echo $cmd > $cmds_log
eval $cmd
while read file; #read in pbp gene accession
do
	b=$(basename $file) #get just the filename without directory
	echo $b
	pid_id=$(echo $b | awk 'sub(/longitudinal.fasta/,//)') #now get just the PID and pbp ID section from existing files to recycle in filename for mafft alignments
	echo $pid_id
	cmd="mafft --adjustdirection $file > $out_dir/$pid_id.mafft_alignment.fasta"
	echo $cmd >> $cmds_log
	eval $cmd
done < $tmp_dir/longitudinal_fastas_filenames.list
