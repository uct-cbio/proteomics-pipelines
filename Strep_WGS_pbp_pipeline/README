This folder contains all the necessary scripts and config file to: 
1. Extract PBP sequences from the Strep enriched contig files (based on specified reference sequences for different PBP genes, obtained from NCBI - see KL_Strep_pneumoniae_NCBI_random_reference_pbp_seqs.fasta, modify as required) and save sequences for each PBP in seperate files for each sample (e.g. in this case each sample will have 3 fasta files (if they all had matches), one for each of the 3 PBP genes in the reference fasta file. Extrction paramters, e.g. expect scores and min matched length can be modified in the config file. Coded in script: automation_script_v2.sh which need to be executed on UCT hex by running the .qsub file like: qsub submit_script.qsub (modify file paths in submit script as necessary).
Output files:
-Commands are recorded in .cmds files e.g. Strep_pbp.107767.65.cmds for sample '107767' in the /log_dir
-BLAST results (contigs vs. NCBI PBP ref seqs) are summarized in tabular format (e.g. 107767.Strep_contigs_against_pbp_seqs_tabular.blastn.report) and as alignments (e.g. 107767.Strep_contigs_against_pbp_seqs_alignments.blastn.report)
-Fasta files of the part of the contig that aligned to the pbp ref seq e.g. 107767.DQ056780.1.matches.fa  

2. For each individual (PID) create a fasta file (for each PBP gene) containing all longitudinal samples available for that individual, for that specific gene and start the file by listing the sequence for the relevant reference gene e.g. pbp2b. Create a summary file (one for each PBP gene) listing for each PID which samples had suitable contig data for the gene in question and which did not (output file name will be e.g. DQ056780.1.summary) Coded in script: longitudinal_alignemtns.sh which needs to be executed on UCT hex by running the .qsub file like qsub submit_script_2.qsub (modify file paths as required)
Output files:
-Commands are recorded for each PBP gene in .cmds files e.g. DQ056780.1.cmds in the log_dir/
-Summary file (one for each PBP gene) listing for each PID which samples had suitable contig data for the gene in question and which did not e.g DQ056780.1.summary
-For each individual (PID) create a fasta file (for each PBP gene) containing all longitudinal samples available for that individual, for that specific gene e.g. 1248.1.1.NPIF.DQ056780.1.longitudinal.fasta

3. The final step is to perform multiple sequence alignements for each of the longitudinal fasta files created in 2. Coded in script: mafft.sh which needs to be executed on UCT hex by running the .qsub file like qsub submit_script_3.qsub (modify file paths as required). This script requires mafft software which is installed on hex.
Output files:
-Commands are recordedin the file mafft_alignments.cmds
-For each longitudinal file from 2. we do a mafft alignment e.g 1248.1.1.NPIF.DQ056780.1.1.mafft_alignment.fasta - this file is used for downstream analysis e.g view alignments and translate to amino acids in AliViewer and export as nexus file from AliViewer for import in FigTree for phylogenetic trees.


#NOTE: to make file.list from the list of contig files to be used you'll need to do something like:
#paste -d '\t' <(for i in `ls -1 /researchdata/fhgfs/clinton.moodley/JCVI_data_correct/*.03_B* | sort`; do basename=$(basename $i}); echo ${basename%.03_*}; done) <(ls -1 /researchdata/fhgfs/clinton.moodley/JCVI_data_correct/*.03_B* | sort) > /home/kviljoen/strep_pbp_project/file.list 
