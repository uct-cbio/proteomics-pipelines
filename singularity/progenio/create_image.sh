version=v1.0.0
export SINGULARITY_CACHEDIR=/cbio/users/ptgmat003/singularity # comment this out if unsure 
export SINGULARITY_TMPDIR=/cbio/users/ptgmat003/singularity/tmp
singularity build progenio_${version}.img docker://thyscbio/progenio:${version}
