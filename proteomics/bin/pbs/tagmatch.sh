#PBS -P CBBI0825 
#PBS -M matthys.potgieter@gmail.com 
#PBS -l select=30:ncpus=24:nodetype=haswell_reg
#PBS -l walltime=48:00:00
#PBS -N adultstool_specific
#PBS -q large
#PBS -W group_list=largeq
#PBS -m be

##############
#  TAGS v8   #
##############

# Parameters
output_folder='/mnt/lustre/users/mpotgieter1/andrew_stool_out/tagmatch_v8_stool_specific_top'
tag_matches='/home/mpotgieter1/lustre/andrew_stool_out/Adult_Stool_denOVO_Exports/tag_matches.txt'
bp_python_chunknumber=5000 # chunknumber = cores available, set bp_blast_num_threads=1
bp_minimum_tag_export_length=4

bp_blast=1 # 1-on, 0-off
bp_blast_denovo_options_handler=consensus # concatenated|consensus|multiple
bp_blast_evalue=200000
bp_blast_matrix='PAM30'
bp_blast_gapopen=9
bp_blast_word_size=2
bp_blast_gapextend=1
bp_blast_num_threads=1
bp_blast_outfmt=5
bp_blast_max_target_seqs=1
bp_blast_max_hsps=2
bp_blast_fasta='/home/mpotgieter1/lustre/uniprot/uniprot_current.fasta'
bp_blast_comp_based_stats=0
bp_blast_window_size=15
bp_blast_threshold=16
bp_blast_word_size=2
bp_blast_seg='no' # mask regions of low complexity (yes|no|...)

blast_tag_min_identities=6

bp_gnu_parallel_j=24
headnode_user_ip=mpotgieter1@lengau.chpc.ac.za  #headnode user account - NB for derivative qsub jobs

bp_sum_pept2lca=0
bp_unipept_prec_tol=0.02 # Validate unmodified sequences against prec ion mass.

unipept_q='serial'
unipept_l='select=1:ncpus=24:mpiprocs=24'
unipept_P='CBBI0825'

#Analysis based on subject sequences from BLAST output
#bp_analysis_hsp_score_cutoff=0

# Trie Graph analysis
bp_tagmatch_prec_tol=0.02
bp_tagmatch_specificity="unspecific"
bp_tagmatch_enzymes="Whole protein" 
bp_tagmatch_max_missed_cleavages=0 
bp_tagmatch_gap_tol=0.5 
bp_tagmatch_fixed_modifications="Carbamidomethylation of C"
bp_tagmatch_variable_modifications=""

bp_analysis_q='smp'
bp_analysis_l='select=1:ncpus=24:mpiprocs=24'
bp_analysis_P='CBBI0825'

bp_analysis_gnu_parallel_j=6
############
# Pipeline #
############

set -e
module add chpc/gnu/parallel-20160422

# Create output folder if needed
if [ ! -d $output_folder ] ; then
    mkdir $output_folder
fi

cd ${PBS_O_WORKDIR}

# Check that config does not exist or is unchanged
if [ ! -f $output_folder/pipeline.pbs ] ; then
    cp "$(readlink -f $0)" $output_folder/pipeline.pbs
else
    cmp --silent "$(readlink -f $0)" $output_folder/pipeline.pbs && echo "'$(readlink -f $0)' unchanged."|| { echo "'$(readlink -f $0)' has changed, please delete '$output_folder' or replace '$(readlink -f $0)' with the contents of pipeline.sh in "${output_folder}; exit 1; }
fi

tagdb=$output_folder/sqlite3.db

if [ ! -d $output_folder/tags ] ; then
    mkdir $output_folder/tags
fi

if [ ! -f $output_folder/tags/tags.fasta ] ; then
    tagfasta.py $tag_matches $bp_minimum_tag_export_length $output_folder/tags/tags.fasta && bp_tagfasta2sqlite3.py $output_folder/tags/tags.fasta $tagdb 'denovo' || { echo "tag export error" ; rm -rf $output_folder/tags/tags.fasta ; exit 1 ; }
fi


#if [ ! -f $output_folder/tags/combined_sequences.fasta ] ; then 
#    bp_fasta.py ${output_folder}/tags/tags.fasta $output_folder/tags/combined_sequences.fasta
#fi
#tag_search_db=${output_folder}/tags/combined_sequences.fasta


if [ ! -d $output_folder/tags/tags ] ; then
    mkdir $output_folder/tags/tags && bp_consensus_groups.py $output_folder/tags/tags.fasta 100 $output_folder/tags/tags || { rm -rf $output_folder/tags/tags ; exit 1; } 
fi

if [ ! -d $output_folder/tags/$bp_blast_denovo_options_handler ] ; then
    mkdir $output_folder/tags/$bp_blast_denovo_options_handler && ls ${output_folder}/tags/tags/*.fasta | parallel -j $bp_gnu_parallel_j -u --sshloginfile ${PBS_NODEFILE} "cd ${PBS_O_WORKDIR}; bp_denovogui_combine_seqs.py {} ${output_folder}/tags/$bp_blast_denovo_options_handler $bp_blast_denovo_options_handler" || { rm -rf $output_folder/tags/$bp_blast_denovo_options_handler ; exit 1 ; }
fi 

if [ ! -f $output_folder/tags/${bp_blast_denovo_options_handler}_tags.fasta ] ; then
    cat $output_folder/tags/${bp_blast_denovo_options_handler}/*.fasta > $output_folder/tags/${bp_blast_denovo_options_handler}_tags.fasta && bp_combinedtags2sqlite3.py $output_folder/tags/${bp_blast_denovo_options_handler}_tags.fasta $tagdb $bp_blast_denovo_options_handler || { rm -rf $output_folder/tags/${bp_blast_denovo_options_handler}_tags.fasta; exit 1 ; }
fi

if [ "$bp_blast" -eq 1 ] ; then
    
    if [ ! -d $output_folder/blast ]; then
        mkdir $output_folder/blast
    fi
    
    if [ ! -d $output_folder/db ] ; then
        mkdir $output_folder/db && cp $bp_blast_fasta $output_folder/db && cd $output_folder/db && makeblastdb -in "$(basename $bp_blast_fasta)" -dbtype 'prot' -out "$(basename $bp_blast_fasta)"
    fi
    
    if [ ! -d $output_folder/blast/fasta ]; then
       
        bp_fasta.py ${output_folder}/tags/${bp_blast_denovo_options_handler}_tags.fasta $output_folder/blast/combined_sequences.fasta
        mkdir $output_folder/blast/fasta && chunkfasta.py $output_folder/blast/combined_sequences.fasta "$(expr $(grep -c '>' $output_folder/blast/combined_sequences.fasta ) / $bp_python_chunknumber)" ${output_folder}/blast/fasta/ || { rm -rf ${output_folder}/blast/fasta ; echo "Error in fasta splitting in blast section"; exit 1; }    
    fi
    
    wait
    NODES=$(cat ${PBS_NODEFILE} | sort | uniq)
    BLASTDB="$( basename $bp_blast_fasta )"
    # copy blast databases to ram disk
    for node in ${NODES}
      do
          ssh ${node} "mkdir -p /dev/shm/${USER}/BLAST && cp -r ${output_folder}/db/${BLASTDB}* /dev/shm/${USER}/BLAST && echo 'successfully added DBs on ${node}' || exit 1" &
      done

    wait  # wait for parallel copies to finish
      
    fasta_count=$( find "${output_folder}/blast/fasta" -name "*.fasta" | wc -l  ) 

    cmd="blastp -query {} -outfmt '$bp_blast_outfmt' -seg $bp_blast_seg -word_size '$bp_blast_word_size' -out {}.xml -comp_based_stats '$bp_blast_comp_based_stats' -window_size '$bp_blast_window_size' -threshold '$bp_blast_threshold' -db /dev/shm/${USER}/BLAST/${BLASTDB} -max_target_seqs $bp_blast_max_target_seqs -max_hsps $bp_blast_max_hsps -num_threads '$bp_blast_num_threads' -evalue '$bp_blast_evalue' -matrix '$bp_blast_matrix' -gapopen '$bp_blast_gapopen' -gapextend '$bp_blast_gapextend' && blast_XML_to_csv.py {}.xml ${output_folder}/blast/combined_sequences.fasta {}.csv $bp_blast_max_target_seqs  &> {}.log && bp_blast_tag_export.py {}.csv ${output_folder}/sqlite3.db $blast_tag_min_identities {}.tag.export.fa && gzip --best {}"

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
    
    if [ ! -f $output_folder/blast/blast_tags.fasta ]; then 
        cat ${output_folder}/blast/fasta/*.tag.export.fa > ${output_folder}/blast/blast_tags.fasta && bp_tagfasta2sqlite3.py ${output_folder}/blast/blast_tags.fasta ${output_folder}/sqlite3.db 'blast' || { rm -rf $output_folder/blast/blast_tags.fasta ; exit 1; }
    fi

    #if [ ! -f ${output_folder}/blast/modified_combined_sequences.fasta ] ; then
    #    cat ${output_folder}/blast/blast_tags.fasta \
    #        ${output_folder}/tags/tags.fasta > ${output_folder}/blast/denovo_blast_tags.fasta  && bp_fasta.py $output_folder/blast/denovo_blast_tags.fasta $output_folder/blast/modified_combined_sequences.fasta
    #fi
    #tag_search_db=${output_folder}/blast/modified_combined_sequences.fasta
fi

if [ ! -d ${output_folder}/analysis ] ; then
    mkdir ${output_folder}/analysis
fi

if [ ! -f ${output_folder}/analysis/trie.p ] ; then
    create_trie.py $output_folder/sqlite3.db 'tags' 'seqstr' ${output_folder}/analysis/trie.p 
fi

if [ ! -d ${output_folder}/analysis/db ] ; then
    mkdir ${output_folder}/analysis/db && chunkfasta.py $bp_blast_fasta "$(expr $(grep -c '>' $bp_blast_fasta) / $bp_python_chunknumber) " ${output_folder}/analysis/db/ || { rm -rf ${output_folder}/analysis/db ; echo "Error in database splitting in analysis section"; exit 1; }
fi

if [ ! -d ${output_folder}/analysis/tables ] ; then
    mkdir ${output_folder}/analysis/tables
fi


db_count=$( find "${output_folder}/analysis/db" -name "*.fasta" | wc -l  ) 
if [ "${db_count}" -ne "0" ]; then 
    cmd="bp_peptide2db.py {} '$bp_tagmatch_prec_tol' '${bp_tagmatch_specificity}' '$bp_tagmatch_enzymes' '$bp_tagmatch_max_missed_cleavages' '$bp_tagmatch_gap_tol' '$bp_tagmatch_fixed_modifications' '$bp_tagmatch_variable_modifications' ${output_folder}/analysis/trie.p ${output_folder}/analysis/tables $output_folder/sqlite3.db && gzip --best {}"
    ls ${output_folder}/analysis/db/*.fasta | parallel -j $bp_analysis_gnu_parallel_j -u --sshloginfile ${PBS_NODEFILE} "cd ${PBS_O_WORKDIR}; ${cmd}"
fi

db_count=$( find "${output_folder}/analysis/db" -name "*.fasta" | wc -l  ) 
if [ "${db_count}" -ne "0" ]; then
    echo 'There are unprocessed database fasta files'
    exit 1
fi

if [ ! -f "${output_folder}/analysis/analysis.csv" ] ; then
    csvcat.py ${output_folder}/analysis/tables > ${output_folder}/analysis/analysis.csv &&  cd ${output_folder} && bp_NSAF.py ${output_folder}/analysis/analysis.csv ${output_folder}/analysis/nsaf_ranked.csv
    #cmd="csvcat.py ${output_folder}/analysis/tables > ${output_folder}/analysis/analysis.csv &&  cd ${output_folder} && bp_NSAF.py ${output_folder}/analysis/analysis.csv ${output_folder}/analysis/nsaf_ranked.csv"
    #ssh $headnode_user_ip 'cd '$output_folder' && echo "'$cmd'" | qsub -N bp_analysis -P '$bp_analysis_P' -q '$bp_analysis_q' -l '$bp_analysis_l' -l walltime=04:00:00'
fi

# Unipept LCA analysis
if [ "$bp_sum_pept2lca" -eq 1 ] ; then
    if [ ! -d $output_folder/unipept ] ; then
       cmd="cd ${output_folder} && mkdir unipept && bp_exclude_tags.py $bp_input_fasta $bp_unipept_prec_tol ${output_folder}/unipept/tags_excluded.fasta && fasta2unipept.py ${output_folder}/unipept/tags_excluded.fasta ${output_folder}/unipept/"
       echo ${cmd}
       ssh $headnode_user_ip 'cd '$output_folder' && echo "'$cmd'" | qsub -N unipept -P '$unipept_P' -q '$unipept_q' -l '$unipept_l' -l walltime=48:00:00'
    fi
fi





