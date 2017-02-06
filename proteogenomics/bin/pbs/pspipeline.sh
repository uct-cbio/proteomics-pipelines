#PBS -P CBBI0825
#PBS -M matthys.potgieter@gmail.com
#PBS -l select=60:ncpus=24:nodetype=haswell_reg
#PBS -l walltime=48:00:00
#PBS -N for_nyari_EpicPilot
#PBS -q large
#PBS -W group_list=largeq
#PBS -m be

###################################
# PeptideShaker/MSnID pipeline v5 #
###################################

##############################################################################################################################
##############################################################################################################################
# NB the number of select= above, must be set to the number of unprocessed mgf files (or lower) but not more!                #
# Otherwise system resources will be wasted.                                                                                 #
# Adjust the parameter if restarting the job to number of remaining uncompressed mgf files in the mgf folder...              #
# Large queue is needed if more than 10 files need to processed at a time.                                                   #
# See PeptideShaker website for details on parameters: http://compomics.github.io/searchgui/wiki/searchcli.html              #
##############################################################################################################################
##############################################################################################################################

##########
# Config #
##########

experiment='EpicPilotSemiTryptic'
output_folder="/home/mpotgieter1/lustre/nyari_out/EpicPilotTest_db1"
spectrum_files="/home/mpotgieter1/lustre/data/for_nyari/mgf"
target_fasta="/home/mpotgieter1/lustre/nyari_out/nsaf_ranked_compactedLatest_Uniprot_human_proteome_Feb2016.fasta"
contaminant_fasta="/mnt/lustre/users/mpotgieter1/blackburn/hypohyper/S507_S5527_hexdata/HYPOHYPER/proteomes/gpm_crap_2016_07_03.fasta"
ps_folder='/home/mpotgieter1/software/PeptideShaker/PeptideShaker-1.14.6'
sg_folder='/home/mpotgieter1/software/SearchGUI/SearchGUI-3.2.5'      

########################
# SearchGUI parameters #
########################

########################################## 
# Optional output compression parameters #
##########################################

output_data='1'     # default 0

##############################
# Optional common parameters #
##############################

xtandem='1' 
myrimatch='1' 
ms_amanda='0'
msgf='1'
omssa='1'
comet='1'
tide='1'
andromeda='0'        # not on linux!

################################
# Optional advanced parameters #
################################

threads=24
mgf_splitting='1000'

################################
# Spectrum matching parameters #
################################

prec_tol=10
prec_ppm=1
frag_tol=0.02
frag_ppm=0
enzyme="Trypsin"
specificity=1    # 0: Specific, 1: Semi-Specific, 2: N-term Specific or 3: C-term Specific
digestion=0      # The type of digestion to consider: 0: Enzyme, 1: Unspecific or 2: Whole Protein. Default is 0.
fixed_mods="Carbamidomethylation of C"
variable_mods="Oxidation of M, Acetylation of protein N-term"
min_charge=2
max_charge=4
mc=2 #Number of allowed missed cleavages, default is '2'.
fi='b'
ri='y'

#######################
# Spectrum annotation #
#######################

annotation_level=0.75

##################
# Import filters #
##################

import_peptide_length_min=8
import_peptide_length_max=30
psm_fdr=1
peptide_fdr=1
protein_fdr=1

#################################
# MyriMatch advanced parameters #
#################################

myrimatch_min_pep_length=8
myrimatch_max_pep_length=30

#############################
# MS-GF advanced parameters #
#############################

msgf_instrument=3
msgf_min_pep_length=8
msgf_max_pep_length=30

############################
# Tide advanced parameters #
############################

tide_min_pep_length=8            # default 8
tide_max_pep_length=30           # default 30

####################
# PTM localization #
####################

ptm_score=1                      # default 1
score_neutral_losses=0           # default 0
ptm_sequence_matching_type=1     # default 1
ptm_alignment=1                  # default 1

###################
# Gene Annotation #
###################

useGeneMapping=0                 # default 1
updateGeneMapping=0              # default 1

######################
# MzidCLI parameters #
######################

contact_first_name='Matthys'
contact_last_name='Potgieter'
contact_email='matthys.potgieter@gmail.com'
contact_address='Same as organization adress'
organization_name='University of Cape Town'
organization_email='organization@email.com'
organization_address='Anzio Road, Observatory'
contact_url='http://www.cbio.uct.ac.za/'
organization_url='http://www.cbio.uct.ac.za'

#################
# Follow-up cli #
#################

psm_type=0
recalibrate=0                # 0=None, 1=True (but twice as likely to fail)

###############
# FDR Control #
###############

MSnID_FDR_value=1            #FDR to control global identifications (%)
MSnID_FDR_level="accession"  # options are 'PSM','peptide','accession'

#################
# gnu paralllel #
#################

ps_gnu_parallel_j=1          # Number of runs per node, default is 1

###################
# Derivative jobs #
###################

headnode_user_ip=mpotgieter1@scp.chpc.ac.za  #headnode user account - NB for derivative qsub jobs
d_q='smp'
d_l='select=1:ncpus=24:mpiprocs=24'
d_P='CBBI0825'

################
# Output files #
################



############
# Pipeline #
############

set -e

JVM_ARGS="-d64 -Xms1024M -Xmx153600M -server"

source compomics.sh

if [ ! -d ${output_folder} ] ; then
    mkdir ${output_folder}
fi

mgf_file_count=$( find "${spectrum_files}" -name "*.mgf" | wc -l  )
# Check that config does not exist or is unchanged
if [ ! -f $output_folder/pipeline.pbs ] ; then
    cp "$(readlink -f $0)" $output_folder/pipeline.pbs
else
    cmp --silent "$(readlink -f $0)" $output_folder/pipeline.pbs && echo "'$(readlink -f $0)' unc    hanged."|| { echo "'$(readlink -f $0)' has changed, please delete '$output_folder' or replace '$(    readlink -f $0)' with the contents of pipeline.sh in "${output_folder}; exit 1; }
fi

if [ ! -d $output_folder/mgf ] ; then
    ps_prepare
fi

cd ${PBS_O_WORKDIR}
module add chpc/gnu/parallel-20160422

mgf_count=$( find "${output_folder}/mgf" -name "*.mgf" | wc -l  )
source `which env_parallel.bash`

if [ "${mgf_count}" -ne "0" ] ; then
    ls ${output_folder}/mgf/*.mgf | env_parallel --timeout 200% -j $ps_gnu_parallel_j -u --sshloginfile ${PBS_NODEFILE} "cd ${PBS_O_WORKDIR}; search $output_folder {}"
fi
wait

mzid_count=$( find "${output_folder}/mzIdentMLs" -name "*.mzid" | wc -l  )
if [ "${mzid_count}" -ne "${mgf_file_count}" ]; then
       echo 'There are unprocessed mgf files'
       exit 1
fi

if [ ! -d $output_folder/mzIdentMLS/analysis ]; then
   cmd="cd ${output_folder} && MSnIDshake.R -i mzIdentMLs/ -v ${MSnID_FDR_value} -l ${MSnID_FDR_level} && cd mzIdentMLs/analysis && unipept pept2lca -i peptides_cleaned.txt -e -o pept2lca.csv && unipept.R"
   ssh $headnode_user_ip 'cd '$output_folder' && echo "'$cmd'" | qsub -N MSnidFDRControlUnipept -P '$d_P' -q '$d_q' -l '$d_l' -l walltime=48:00:00'
    
fi

