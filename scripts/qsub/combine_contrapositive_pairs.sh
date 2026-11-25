source /etc/profile.d/modules.sh 
cd $PBS_O_WORKDIR 

bash scripts/bash/combine_contrapositive_pairs.sh
