source /etc/profile.d/modules.sh 
cd $PBS_O_WORKDIR 

perl scripts/perl/split_fasta.pl 100 raw/gtdb/ar122_ssu_reps_r202.fna data/gtdb_ssu_seqs 3
perl scripts/perl/split_fasta.pl 100 raw/gtdb/bac120_ssu_reps_r202.fna data/gtdb_ssu_seqs 3
