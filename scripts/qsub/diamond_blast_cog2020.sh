source /etc/profile.d/modules.sh
module load diamond/2.0.11
cd $PBS_O_WORKDIR

seq -w 0 249 | xargs -I {} -P 36 -n 1 sh -c "diamond blastp -q data/seqs_split/proteins_raw.faa.{} -d db/diamond_cog2020/cog-20.dmnd -o data/diamond_blast_cog2020/proteins_raw.tsv.{} -p 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp slen scovhsp 2>qsub_out/diamond_blast_cog2020/{}.err"
