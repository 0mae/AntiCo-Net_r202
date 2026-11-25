source /etc/profile.d/modules.sh
module load seqkit/0.16.1
cd $PBS_O_WORKDIR

seqkit stats data/seqs/proteins_raw.faa
cat data/seqs_split/* | seqkit stats
