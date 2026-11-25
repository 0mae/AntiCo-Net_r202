source /etc/profile.d/modules.sh
module load CheckM/1.1.3
cd $PBS_O_WORKDIR

checkm analyze -x fna -t 12 raw/gtdb/cpr_43_markers.hmm data/gtdb/cpr_genomes_reps data/gtdb/cpr_checkm
