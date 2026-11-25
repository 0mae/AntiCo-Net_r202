source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
cp db/gtdb_release202/genomic_files_reps/ar122_ssu_reps_r202.tar.gz raw/gtdb/
cp db/gtdb_release202/genomic_files_reps/bac120_ssu_reps_r202.tar.gz raw/gtdb/
tar -zxvf raw/gtdb/ar122_ssu_reps_r202.tar.gz -C raw/gtdb/
tar -zxvf raw/gtdb/bac120_ssu_reps_r202.tar.gz -C raw/gtdb/
