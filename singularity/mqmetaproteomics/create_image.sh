version=v2.7.1
export SINGULARITY_CACHEDIR=/cbio/users/ptgmat003/singularity # comment this out if unsure 
export SINGULARITY_TMPDIR=/cbio/users/ptgmat003/singularity/tmp

singularity build mqmetaproteomics_${version}.img docker://thyscbio/mqmetaproteomics:${version}
