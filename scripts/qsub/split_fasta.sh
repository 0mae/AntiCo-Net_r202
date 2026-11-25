source /etc/profile.d/modules.sh 
cd $PBS_O_WORKDIR 

perl scripts/perl/split_fasta.pl 421293 data/seqs/proteins_raw.faa data/seqs_split 3
