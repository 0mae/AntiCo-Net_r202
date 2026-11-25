#PBS -S /bin/bash
#PBS -l ncpus=40
#PBS -V
cd $HOME/crest/cpnet
singularity exec --bind $(pwd):/temp cpnet_latest.sif bash /temp/scripts/bash/count_og_displacement_2.sh cog archaea 40
