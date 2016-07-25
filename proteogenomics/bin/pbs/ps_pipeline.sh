#!/usr/bin/env bash

N="MyJob"
q="UCTlong"
l="nodes=1:ppn=10:series600"
M="matthys.potgieter@gmail.com"
m="ae"

config='/home/ptgmat003/cbio-pipelines/proteogenomics/bin/config/PS_Pipeline_example.cfg'

######## PIPELINE RUNS FROM HERE #########

#### DELETE PENDING JOBS FOR THIS SCRIPT
if [ -f $config.joblist ]; then
    while IFS= read -r line 
    do
        echo $line
        qdel $line 
        qsig -s SIGINT $line
    done < $config.joblist
    rm $config.joblist 
fi
rm -rf MyJob*
eval "$(egrep '^#|^[^ ]*=[^;&]*'  "$config")"


function dependencies {
    arr=$1
    arr=$(printf ":%s" "${arr[@]}")
    arr=${arr:1}
    echo $arr
}

############################################
# CREATE PARAMETERS
############################################

newruns=()
oldruns=()
if [ ! -d $output_folder ]; then
    ps_prepare=$(echo "PS_prepare.sh $config" | qsub -N $N.PS_prepare.sh -q $q -l $l -M $M -m $m)
    echo $ps_prepare >> $config.joblist
    newruns+=($ps_prepare)

else
    cmp --silent $config $output_folder/config.cfg && echo ${config}' is unchanged.'|| { echo ${config}' has changed, please delete '$output_folder' or replace '$config' with the contents of config.cfg in '$output_folder; exit 1; }   
fi

############################################
# SEARCH SAMPLES WITH PEPTIDESHAKER
############################################

oldruns+=$newruns
newruns=()
for spectrum_file in ${spectrum_files}/*.mgf; do
    sample=$(basename $spectrum_file)
    if [ $recalibrate -eq 1 ]; then
        sample=${sample%.mgf}_recalibrated.mgf
    fi
    if [ ! -f $output_folder/mzIdentMLs/$sample.mzid ]; then	
	if [ -z $oldruns ]; then
 	    ps_sample=$(echo "PS_sample.sh '$output_folder' '$spectrum_file'" | qsub -N $N -q $q -l $l -M $M -m $m)
    	else
	    ps_sample=$(echo "PS_sample.sh '$output_folder' '$spectrum_file'" | qsub -N $N -q $q -l $l -M $M -m $m -W depend=afterok:$(dependencies $oldruns))
        fi
	echo $ps_sample
        echo $ps_sample >> $config.joblist
        newruns+=($ps_sample)
    fi
done

############################################
############################################

oldruns+=$newruns
newruns=()
if [ ! -d $output_folder/mzIdentMLS/analysis ]; then
    cmd="cd ${output_folder} && MSnIDshake.R -i mzIdentMLs/ -v ${MSnID_FDR_value} -l ${MSnID_FDR_level}"
    if [ -z $oldruns ]; then 
        msnid_shake=$(echo "${cmd}" | qsub -N $N -q $q -l $l -M $M -m $m )
    else
	msnid_shake=$(echo "${cmd}" | qsub -N $N -q $q -l $l -M $M -m $m -W depend=afterok:$(dependencies $oldruns))
    fi
    echo $msnid_shake
    newruns+=($msnid_shake)
fi

############################################
############################################
oldruns+=$newruns
newruns=()

