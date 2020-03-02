#PBS -P CBBI0825 
#PBS -M matthys.potgieter@gmail.com 
#PBS -l select=15:ncpus=24:nodetype=haswell_reg
#PBS -l walltime=48:00:00
#PBS -N adultstool_blast
#PBS -q large
#PBS -W group_list=largeq
#PBS -m be

###########W##@@@#
#  MetaTaxa v1   #
############W#####

##############
# Parameters #
##############

output_folder='/home/mpotgieter1/lustre/andrew_stool_out/target_helminths/output'
input_fasta='/home/mpotgieter1/lustre/andrew_stool_out/target_helminths/all_identified_helminth_targets.fasta'

bp_python_chunksize=1
bp_blast_evalue=200000
bp_blast_matrix='PAM30'
bp_blast_gapopen=9
bp_blast_word_size=2
bp_blast_gapextend=1
bp_blast_num_threads=1
bp_blast_outfmt=5
bp_blast_max_target_seqs=10 # ignored during BLAST, applied during export
bp_blast_max_hsps=1
bp_blast_fasta='/home/mpotgieter1/lustre/uniprot/uniprot_sept/uniprot_sept.fasta'
bp_blast_comp_based_stats=0
bp_blast_window_size=15
bp_blast_threshold=16
bp_blast_word_size=2
bp_blast_seg='yes' # mask regions of low complexity (yes|no|...)

bp_gnu_parallel_j=24

############
# Pipeline #
############

set -e
module add chpc/gnu/parallel-20160422

# Create output folder if needed
if [ ! -d $output_folder ] ; then
    mkdir $output_folder
fi
script="$(readlink -f $0)"
echo $script

cd ${PBS_O_WORKDIR}

# Check that config does not exist or is unchanged
if [ ! -f $output_folder/pipeline.pbs ] ; then
    cp $script $output_folder/pipeline.pbs
else
    cmp --silent $script $output_folder/pipeline.pbs && echo "$script unchanged."|| { echo "'$script' has changed, please delete '$output_folder' or replace '$script' with the contents of pipeline.sh in "${output_folder}; exit 1; }
fi

#########
# BLAST #
#########

if [ ! -d $output_folder/blast ]; then
    mkdir $output_folder/blast
fi

if [ ! -d $output_folder/db ] ; then
    mkdir $output_folder/db && cp $bp_blast_fasta $output_folder/db && cd $output_folder/db && makeblastdb -in "$(basename $bp_blast_fasta)" -dbtype 'prot' -out "$(basename $bp_blast_fasta)"
fi

if [ ! -d $output_folder/blast/fasta ]; then
    mkdir $output_folder/blast/fasta && chunkfasta.py $input_fasta $bp_python_chunksize ${output_folder}/blast/fasta/ || { rm -rf ${output_folder}/blast/fasta ; echo "Error in fasta splitting in blast section"; exit 1; }    
fi

wait
NODES=$(cat ${PBS_NODEFILE} | sort | uniq)
BLASTDB="$( basename $bp_blast_fasta )"

# copy blast databases to ram disk
#for node in ${NODES}
#  do
#      ssh ${node} "mkdir -p /dev/shm/${USER}/BLAST && cp -r ${output_folder}/db/${BLASTDB}* /dev/shm/${USER}/BLAST && echo 'successfully added DBs on ${node}' || exit 1" &
#  done

wait  # wait for parallel copies to finish
  
fasta_count=$( find "${output_folder}/blast/fasta" -name "*.fasta" | wc -l  ) 

cmd="blastp -query {} -outfmt '$bp_blast_outfmt' -seg $bp_blast_seg -word_size '$bp_blast_word_size' -out {}.xml -comp_based_stats '$bp_blast_comp_based_stats' -window_size '$bp_blast_window_size' -threshold '$bp_blast_threshold' -db ${output_folder}/db/${BLASTDB} -max_hsps $bp_blast_max_hsps -num_threads '$bp_blast_num_threads' -evalue '$bp_blast_evalue' -matrix '$bp_blast_matrix' -gapopen '$bp_blast_gapopen' -gapextend '$bp_blast_gapextend' && blast_XML_to_csv.py {}.xml ${input_fasta} {}.csv $bp_blast_max_target_seqs  &> {}.log && gzip --best {}"

if [ "${fasta_count}" -ne "0" ]; then
    ls ${output_folder}/blast/fasta/*.fasta | parallel -j $bp_gnu_parallel_j -u --sshloginfile ${PBS_NODEFILE} "cd ${PBS_O_WORKDIR}; ${cmd}"
fi
wait 

# clean up ram disk
for node in ${NODES}
  do
    ssh ${node} "rm -rf /dev/shm/${USER}/BLAST && echo 'successfully deleted DBs on ${node}' || exit 1" &
done
wait

fasta_count=$( find "${output_folder}/blast/fasta" -name "*.fasta" | wc -l  ) 
if [ "${fasta_count}" -ne "0" ]; then
    echo 'There are unprocessed fasta files'
    exit 1
fi

if [ ! -d $output_folder/blast/combined.csv ] ; then
   cat $output_folder/blast/fasta/*.csv > $output_folder/blast/combined.csv || rm -rf $output_folder/blast/combined.csv

fi

if [ ! -d $output_folder/blast/ranked_matches.csv ] ; then
   bp_blast_taxa.py $output_folder/blast/combined.csv $output_folder/ranked_matches.csv || rm -rf $output_folder/blast/ranked_matches.csv

fi



