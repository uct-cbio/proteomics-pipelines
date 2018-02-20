
output_folder=$(readlink -f $1)
genome_inpath=$(readlink -f $2)
provider=$3
source_url=$4
genome_name=$5
genome_version=$6

BSCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
RSCRIPTPATH=$(dirname $BSCRIPTPATH)/R_scripts


line=$(head -n 1 $genome_inpath)
line="${line%%,*}"
line=$(sed -e 's#.*>\(\)#\1#' <<< $line)


genus=$(echo $line | cut -d " " -f2)
species=$(echo $line | cut -d " " -f3)
strain=$(echo $line | cut -d " " -f4-)

#name=$(echo $line | cut -d " " -f2-)

name=$(echo "$line" | tr ' ' '_')
name=$(echo "$name" | tr '.' '_')

echo $name
#genome_outname=$( echo $(basename $genome_inpath) | sed 's/\.[^.]*$//' ).fa

genome_outpath=$output_folder/chr1.fa

rm -rf ${output_folder}
mkdir ${output_folder} && cd ${output_folder}

cp $genome_inpath $genome_outpath

Rscript --vanilla $RSCRIPTPATH/BSgenome.R \
-g $genome_outpath \
-p $provider \
-s $source_url \
-n $genome_name \
-v $genome_version

R CMD build 'BSgenome.'$genome_name'.NCBI.'$genome_version
R CMD check 'BSgenome.'$genome_name'.NCBI.'$genome_version
R CMD INSTALL --library=${HOME}/R/x86_64-unknown-linux-gnu-library/3.2/ *.tar.gz
