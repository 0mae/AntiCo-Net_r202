source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

cat data/prokatlas/bac120_ssu_reps_r202.tsv.* data/prokatlas/ar122_ssu_reps_r202.tsv.* | sed -E 's/seq_.*barcodelabel=//g' | awk '$10/$6 >= 0.9 && $9/$10 >= 0.97 {print}' | sed '1s/^/qseqid\tsseqid\tqlen\tqstart\tqend\tslen\tsstart\tsend\tnident\tlength\n/' > data/tsv/gtdb_ssu_prokatlas.tsv
