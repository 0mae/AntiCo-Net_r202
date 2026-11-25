source /etc/profile.d/modules.sh
module load blast+/2.12.0
cd $PBS_O_WORKDIR

seq -w 0 15 | sed 's/^/0/' | xargs -P 40 -n 1 -I {} sh -c "blastn -db db/prokatlas/ProkAtlas.fa -query data/gtdb_ssu_seqs/ar122_ssu_reps_r202.fna.{} -out data/prokatlas/ar122_ssu_reps_r202.tsv.{} -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send nident length' -evalue 1e-10 -max_target_seqs 361474 -num_threads 1"

seq -w 0 313 | xargs -P 40 -n 1 -I {} sh -c "blastn -db db/prokatlas/ProkAtlas.fa -query data/gtdb_ssu_seqs/bac120_ssu_reps_r202.fna.{} -out data/prokatlas/bac120_ssu_reps_r202.tsv.{} -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send nident length' -evalue 1e-10 -max_target_seqs 361474 -num_threads 1"
