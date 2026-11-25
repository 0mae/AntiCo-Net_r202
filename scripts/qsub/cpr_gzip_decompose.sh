source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

find data/gtdb/cpr_genomes_reps/ -name "*.fna.gz" | xargs -I {} -P 8 -n 1 sh -c "gzip -d {}"
