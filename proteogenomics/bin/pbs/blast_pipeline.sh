#!/usr/bin/env bash

set -e

source stop.sh 
source submit.sh

config=$1
eval "$(egrep '^#|^[^ ]*=[^;&]*'  "$config")"

# Create output folder if needed
if [ ! -d $output_folder ] ; then
    mkdir $output_folder
fi

# Check that config does not exist or is unchanged
if [ ! -f $output_folder/config.cfg ] ; then
    cp $config $output_folder/config.cfg
else
    cmp --silent $config $output_folder/config.cfg && echo ${config}' is unchanged.'|| { echo ${config}' has changed, please delete '$output_folder' or replace '$config' with the contents of config.cfg in '$output_folder; exit 1; }
fi

# BLAST path

export BLASTDB=$bp_blast_BLASTDB

# PBS parameters (global)

P=$bp_pbs_P

# Stop running jobs from same config

stop.sh $bp_log_file

jobs=()

# blast_pipeline.pbs part 1: Split fasta with python

q=$bp_python_pbs_q
l=$bp_python_pbs_l

if [ ! -d $output_folder/blast ]; then
    mkdir $output_folder/blast
fi

if [ ! -d $output_folder/blast/fasta ]; then
    mkdir $output_folder/blast/fasta
    cmd="chunkfasta.py $bp_input_fasta $bp_python_chunknumber $output_folder/blast/fasta"
    job=$(submit "$cmd" $jobs); jobs+=($job)
    echo $job >> $bp_log_file
else
    echo "$output_folder/blast/fasta already exists - stopping blast_pipeline.pbs part 1"  
fi


# blast_pipeline.pbs part 2: RUN BLAST

q=$bp_blast_q
l=$bp_blast_l

newlist=()

if [ ! -d $output_folder/blast/results ]; then
    mkdir $output_folder/blast/results
fi

jobs=()
echo "Waiting for all query fasta files"

for i in $(seq 1 $bp_python_chunknumber) ; do
    file=$(basename $bp_input_fasta)
    infile=$output_folder/blast/fasta/${file}.part.${i}.fasta
    while [ ! -f $infile ]
    do 
        echo "Waiting for $infile"
        sleep 5
    done

    outfile=$output_folder/blast/results/${file}.part.${i}.fasta
    if [ ! -f ${outfile}.out ]; then
       cmd="blastp -query '$infile' -outfmt '$bp_blast_outfmt' -out '${outfile}.temp' -db UniProtCurrent -num_threads '$bp_blast_num_threads' -evalue '$bp_blast_evalue' -matrix '$bp_blast_matrix' -gapopen '$bp_blast_gap_open' -gapextend '$bp_blast_gap_extend' -word_size '$bp_blast_word_size' && mv '${outfile}.temp' '${outfile}.out'";
       job=$(submit "$cmd"  $jobs) ; newlist+=($job);
       echo $job >> $bp_log_file
       echo $job
    fi
done

echo 'Waiting for BLAST to finish processing the chunks...'


# Make sure everythin has processed correctly
completed='FALSE'
while [ "$completed" = 'FALSE' ]
do 
  sleep 5
  count=$( ls  "${output_folder}/blast/results" | wc -l ) 
  
  if [ "$count" -eq "$bp_python_chunknumber" ] ; then
      temp=$( ls "${output_folder}/blast/results/" | grep '.fasta.temp' | wc -l )
      if [ "$temp" -eq 0  ] ; then
          completed='TRUE' 
      fi
  fi
done


echo 'Starting BLAST export'

# Export a table of all blast results
q=$bp_xml_pbs_q
l=$bp_xml_pbs_l
jobs=()

if [ ! -d $output_folder/blast/tables ] ; then
    mkdir $output_folder/blast/tables
fi

jobs=()
for i in $(seq 1 $bp_python_chunknumber) ; do
    file=$(basename $bp_input_fasta)
    infile=$output_folder/blast/results/${file}.part.${i}.fasta.out
    outfile=$output_folder/blast/tables/${file}.part.${i}.fasta.csv

    if [ ! -f ${outfile} ]; then
        cmd="blast_XML_to_csv.py $infile $bp_input_fasta $outfile"
        job=$(submit "$cmd"  $jobs) ; 
        echo $job >> $bp_log_file
    fi
done



# Make sure everything has processed correctly

echo 'Waiting for BLAST export to complete'
completed='FALSE'

while [ "$completed" = 'FALSE' ]
do 
  sleep 5
  count=$( ls  "${output_folder}/blast/tables" | wc -l )   
  if [ "$count" -eq "$bp_python_chunknumber" ] ; then
      completed='TRUE' 
  fi
done

# Create summary table from blast output 'tables' folder
q=$bp_sum_pbs_q
l=$bp_sum_pbs_l
jobs=()
echo "Starting filtering alignments"

if [ ! -d "${output_folder}/blast/filtered" ] ; then
    mkdir "${output_folder}/blast/filtered"
fi

for i in $(seq 1 $bp_python_chunknumber) ; do
    infile=$output_folder/blast/tables/${file}.part.${i}.fasta.csv
    outfile=$output_folder/blast/filtered/${file}.part.${i}.fasta.csv
    if [ ! -f ${outfile} ]; then
        cmd="bp_blast_filter.py $infile $bp_sum_aln_cutoff $outfile"
        job=$(submit "$cmd"  $jobs) ; 
        echo $job >> $bp_log_file
    fi
done


if [ ! -f "${output_folder}/blast/combined_blast.tsv" ] ; then
    cmd="csvcat.py '${output_folder}/blast/filtered' > '${output_folder}/blast/combined_blast.tsv'"
    job=$(submit "$cmd"  $jobs) ; 
    echo $job >> $bp_log_file
fi


# Unipept LCA analysis
if [ "$bp_sum_pept2lca" -eq 1 ] ; then
    if [ ! -d $output_folder/unipept ] ; then
        mkdir $output_folder/unipept
    fi

    if [ ! -f $output_folder/unipept/merged_blast_unipept.csv ] ; then
        
    cmd="csv2unipeptlca.py '${output_folder}/blast/combined_blast.tsv' 'hsp.sbjct' '${output_folder}/unipept/'"
    job=$(submit "$cmd"  $jobs) ; 
    echo $job >> $bp_log_file
    fi
fi



echo 'DONE'
