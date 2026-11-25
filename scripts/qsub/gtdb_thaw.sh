source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
tar -zxvf db/gtdb_release202/genomic_files_reps/gtdb_proteins_aa_reps_r202.tar.gz -C raw/gtdb/
tar -zxvf db/gtdb_release202/genomic_files_reps/gtdb_proteins_nt_reps_r202.tar.gz -C raw/gtdb/
tar -zxvf db/gtdb_release202/genomic_files_reps/gtdb_genomes_reps_r202.tar.gz -C raw/gtdb/
