
#PURPOSE OF THIS SCRIPT: 
#1. CREATE LOCAL BLAST DB FROM FASTA FILE CONTAINING TARGET PBP GENE SEQS
#2. FOR A GIVEN SAMPLE, BLAST THE CONTIGS FILE CONTAINING CONTIGS THAT MAP TO STREPTOCOCCUS AGAINST THE DB CREATED IN 1. --> OUTPUT = TABULAR FILE SUMMARY OF HITS
#3. Get the contig IDs from #2 to filter contigs file
#4. Now filter the original contigs file using matching IDs from #3.
#5. Perform BLAST on filtered file (saves compuational time) to extract the aligned part of these contigs to file so that we can use these downstream
#6. For reference lets also do the plain BLAST alignments file: so same as steps 2 and 3, just default output
#7. Now the file from 6 needs to be split into subfiles one for each accession number in the query database
#8. Now use the IDs from step 7. to filter the fasta file to get all seqs matching specific accession number

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


count=1


while read file;
do 
sid=`echo $file | awk -F ' ' '{print $1}'` #extract sample ID from file.list
cmds_log=$log_dir/Strep_pbp.$sid.$count.cmds
echo -n "" > $cmds_log #in case file was made previously first clear it before continuing..

#2. BLAST Strep contigs against local Strep pbp DB --> output = tabular summary of hits

cmd="$blast_base/blastn -db $scripts_dir/KL_Strep_pneumoniae_NCBI_random_reference_pbp_seqs.fsa -query $in_dir/$sid.03_B_strep_ctgs.fasta -num_threads $num_threads -evalue $evalue -num_alignments $num_output -out $out_dir/$sid.Strep_contigs_against_pbp_seqs_tabular.blastn.report -outfmt 6"

echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd

#3. Get the contig IDs from #2 to filter contigs file
cmd="cut -f1 $out_dir/$sid.Strep_contigs_against_pbp_seqs_tabular.blastn.report > $tmp_dir/$sid.Strep_contigs_against_pbp_seqs_IDs.txt"
echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd
#4. Now filter the original contigs file using matching IDs from #3.
cmd="source $scripts_dir/activate_qiime.sh"
echo $cmd >> $cmds_log
eval $cmd
cmd="filter_fasta.py -f $in_dir/$sid.03_B_strep_ctgs.fasta -o $tmp_dir/$sid.strep_ctgs_pbs_filtered.fa -s $tmp_dir/$sid.Strep_contigs_against_pbp_seqs_IDs.txt"
echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd
#5. Perform BLAST on filtered file (saves compuational time) to extract the aligned part of these contigs to file so that we can use these downstream
cmd="$blast_base/blastn -db $scripts_dir/KL_Strep_pneumoniae_NCBI_random_reference_pbp_seqs.fsa -query $tmp_dir/$sid.strep_ctgs_pbs_filtered.fa -num_threads $num_threads -evalue $evalue -num_alignments $num_output -out $tmp_dir/$sid.Strep_contigs_against_pbp_seqs_output_query_seqs.blastn.report -outfmt '6 sseqid qseqid qseq'"
echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd
#6. For reference lets also do the plain BLAST alignments file: so same as steps 2 and 3, just default output
cmd="$blast_base/blastn -db $scripts_dir/KL_Strep_pneumoniae_NCBI_random_reference_pbp_seqs.fsa -query $tmp_dir/$sid.strep_ctgs_pbs_filtered.fa -num_threads $num_threads -evalue $evalue -num_alignments $num_output -out $out_dir/$sid.Strep_contigs_against_pbp_seqs_alignments.blastn.report"
echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd
#5. Now take the aligned part of contigs (output from step 3.) and convert this tabular file to fasta format so we can filter it based on seq. length
cmd="$python_base $scripts_dir/tab2fasta.py $tmp_dir/$sid.Strep_contigs_against_pbp_seqs_output_query_seqs.blastn.report 3 1 2 > $tmp_dir/$sid.hits.fa"
echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd
#6. Now filter output from step 5 by sequence length (See http://itrylinux.com/use-awk-to-filter-fasta-file-by-minimum-sequence-length/)
cmd="awk '!/^>/ { next } { getline seq } length(seq) >= $min_len { print \$0 \"\n\" seq }' $tmp_dir/$sid.hits.fa > $tmp_dir/$sid.hits.length_filtered.fa"
echo $cmd >> $cmds_log
echo $'\n' >> $cmds_log #add newline
eval $cmd
#7. Now the file from 6 needs to be split into subfiles one for each accession number in the query database
   while read accession;
   do
   ID=$accession #extract accession number from query.list
   cmd="grep $ID $tmp_dir/$sid.hits.length_filtered.fa > $tmp_dir/$sid.$ID.matches.txt" #fasta header in which the given accession occurs
   echo $cmd >> $cmds_log
   echo $'\n' >> $cmds_log #add newline
   eval $cmd
   cmd="awk 'sub(/^>/, \"\")' $tmp_dir/$sid.$ID.matches.txt > $tmp_dir/$sid.$ID.matches2.txt" #remove the > from the fasta header
   echo $cmd >> $cmds_log
   echo $'\n' >> $cmds_log #add newline
   eval $cmd
#8. Now use the IDs from step 7. to filter the fasta file to get all seqs matching specific accession numbeir 
   cmd="filter_fasta.py -f $tmp_dir/$sid.hits.length_filtered.fa -o $out_dir/$sid.$ID.matches.fa -s $tmp_dir/$sid.$ID.matches2.txt" 
   echo $cmd >> $cmds_log
   echo $'\n' >> $cmds_log #add newline
   eval $cmd
   done < $query_list
 
(( count+=1 ))

done < $file_list

#NOTE: to make file.list from the list of files do something like:
#paste -d '\t' <(for i in `ls -1 /researchdata/fhgfs/clinton.moodley/JCVI_data_correct/*.03_B* | sort`; do basename=$(basename $i}); echo ${basename%.03_*}; done) <(ls -1 /researchdata/fhgfs/clinton.moodley/JCVI_data_correct/*.03_B* | sort) > /home/kviljoen/strep_pbp_project/file.list

