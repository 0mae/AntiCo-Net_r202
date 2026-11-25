source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

bash scripts/bash/join_proteins_raw_genomes_analyzed.sh
