version=v4.0.1
export SINGULARITY_CACHEDIR=/cbio/users/ptgmat003/singularity # comment this out if unsure 
export SINGULARITY_TMPDIR=/cbio/users/ptgmat003/singularity/tmp
singularity build mqproteogenomics_${version}.img docker://thyscbio/mqproteogenomics:${version}
