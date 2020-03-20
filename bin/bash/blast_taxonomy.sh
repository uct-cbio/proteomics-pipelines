###########W##@@@#
#  MetaTaxa v1   #
############W#####

##############
# Parameters #
##############

config=$1
source $config

output_folder=${output_folder}
input_fasta=${input_fasta}
bp_python_chunksize=${bp_python_chunksize}
bp_blast_evalue=${bp_blast_evalue}
bp_blast_matrix=${bp_blast_matrix}
bp_blast_gapopen=${bp_blast_gapopen}
bp_blast_word_size=${bp_blast_word_size}
bp_blast_gapextend=${bp_blast_gapextend}
bp_blast_num_threads=${bp_blast_num_threads}
bp_blast_outfmt=${bp_blast_outfmt}
bp_blast_max_target_seqs=${bp_blast_max_target_seqs} # ignored during BLAST, applied during export
bp_blast_max_hsps=${bp_blast_max_hsps}
bp_blast_fasta=${bp_blast_fasta}
bp_blast_comp_based_stats=${bp_blast_comp_based_stats}
bp_blast_window_size=${bp_blast_window_size}
bp_blast_threshold=${bp_blast_threshold}
bp_blast_word_size=${bp_blast_word_size}
bp_blast_seg=${bp_blast_seg} # mask regions of low complexity (yes|no|...)
bp_gnu_parallel_j=${bp_gnu_parallel_j}

echo $output_folder
############
# Pipeline #
############

set -e

# Create output folder if needed
if [ ! -d $output_folder ] ; then
    mkdir $output_folder
fi
script="$(readlink -f $0)"
echo $script

#cd ${PBS_O_WORKDIR}
####################################################
# Check that config does not exist or is unchanged #
####################################################

if [ ! -f $output_folder/config.env ] ; then
    cp $config $output_folder/config.env
else
    cmp --silent $config $output_folder/config.env && echo "$config unchanged."|| { echo "'$config' has changed, please delete '$output_folder' or replace '$config' with the contents of config.env in "${output_folder}; exit 1; }
fi

# Check that pipeline does not exist or is unchanged
if [ ! -f $output_folder/pipeline.sh ] ; then
    cp $script $output_folder/pipeline.sh
else
    cmp --silent $script $output_folder/pipeline.sh && echo "$script unchanged."|| { echo "'$script' has changed, please delete '$output_folder' or replace '$script' with the contents of pipeline.sh in "${output_folder}; exit 1; }
fi


echo Proceeding...
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

echo Done splitting query fasta...

wait
#NODES=$(cat ${PBS_NODEFILE} | sort | uniq)
BLASTDB="$( basename $bp_blast_fasta )"

echo Preparing for BLAST...
# copy blast databases to ram disk
#for node in ${NODES}
#  do
#      ssh ${node} "mkdir -p /dev/shm/${USER}/BLAST && cp -r ${output_folder}/db/${BLASTDB}* /dev/shm/${USER}/BLAST && echo 'successfully added DBs on ${node}' || exit 1" &
#  done

wait  # wait for parallel copies to finish
  
fasta_count=$( find "${output_folder}/blast/fasta" -name "*.fasta" | wc -l  ) 

echo Need to process $fasta_count fasta files...

cmd="blastp -query {} -outfmt '$bp_blast_outfmt' -seg $bp_blast_seg -word_size '$bp_blast_word_size' -out {}.xml -comp_based_stats '$bp_blast_comp_based_stats' -window_size '$bp_blast_window_size' -threshold '$bp_blast_threshold' -db ${output_folder}/db/${BLASTDB} -max_hsps $bp_blast_max_hsps -num_threads '$bp_blast_num_threads' -evalue '$bp_blast_evalue' -matrix '$bp_blast_matrix' -gapopen '$bp_blast_gapopen' -gapextend '$bp_blast_gapextend' && blast_XML_to_csv.py {}.xml ${input_fasta} {}.csv $bp_blast_max_target_seqs  &> {}.log && gzip --best {}"

if [ "${fasta_count}" -ne "0" ]; then
    #ls ${output_folder}/blast/fasta/*.fasta | parallel -j $bp_gnu_parallel_j -u --sshloginfile ${PBS_NODEFILE} "cd ${PBS_O_WORKDIR}; ${cmd}"
    ls ${output_folder}/blast/fasta/*.fasta | parallel -j $bp_gnu_parallel_j -u "${cmd}"
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



