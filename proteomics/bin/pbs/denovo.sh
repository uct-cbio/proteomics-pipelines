#PBS -P CBBI0825
#PBS -M matthys.potgieter@gmail.com
#PBS -l select=10:ncpus=24:nodetype=haswell_reg
#PBS -l walltime=48:00:00
#PBS -q normal
#PBS -m be

###############################
# Mandatory script parameters #
###############################
output_folder="/mnt/lustre/users/mpotgieter1/suereta_stool_out/baby_stool_universal"
spectrum_files="/mnt/lustre/users/mpotgieter1/blackburn/suereta_baby_stool/MGF"  #this must be a folder containing mgf files

dg_folder=/home/mpotgieter1//software/DeNovoGUI/DeNovoGUI-1.12.3

# DenovoGUI parameters
pepnovo=1
directag=1
pnovo=0
novor=1

# Spectrum matching parameters
fixed_mods="Carbamidomethylation of C"
variable_mods="Oxidation of M, Acetylation of protein N-term"
frag_tol=0.02 # in Da, NB for High Res qexactive


############
# Pipeline #
############

set -e

JVM_ARGS="-d64 -Xms1024M -Xmx15360M -server"

if [ ! -d ${output_folder}]; then
  mkdir ${output_folder}
fi

temp_folder=${output_folder}"/temp"
if [ ! -d ${temp_folder} ]; then
    mkdir ${temp_folder}
fi

log_folder=${output_folder}"/log"
if [ ! -d ${log_folder} ]; then
    mkdir ${log_folder}
fi

cp -R ${dg_folder} ${output_folder}/dg
wait

# Create search parameters
java $JVM_ARGS -cp ${output_folder}/dg/DeNovoGUI-*.jar com.compomics.denovogui.cmd.IdentificationParametersCLI -out ${output_folder}/identification.par -fixed_mods "${fixed_mods}" -variable_mods     "${variable_mods}" -frag_tol ${frag_tol}
wait
java $JVM_ARGS -cp ${output_folder}/dg/DeNovoGUI-*.jar com.compomics.denovogui.cmd.PathSettingsCLI -temp_folder ${temp_folder} -log ${log_folder}






