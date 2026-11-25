source /etc/profile.d/modules.sh 
cd $PBS_O_WORKDIR 
source scripts/bash/functions.sh 

cat data/tsv/genomes_analyzed_bac.tsv data/tsv/genomes_analyzed_ar.tsv | cut -f 1 | sed '/accession/d' | gtdbacc2path | xargs -I {} -P 8 -n 1 sh -c "cp raw/gtdb/gtdb_genomes_reps_r202/{}*_genomic.fna.gz data/genomes"
