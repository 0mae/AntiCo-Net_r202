source /etc/profile.d/modules.sh
module load seqkit/0.16.1
cd $PBS_O_WORKDIR

cat data/tsv/proteins_raw.tsv | tail -n +2 | cut -f 2,13 | seqkit tab2fx > data/seqs/proteins_raw.faa
