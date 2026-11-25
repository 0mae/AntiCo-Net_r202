source /etc/profile.d/modules.sh
module load seqkit/0.16.1
cd $PBS_O_WORKDIR
source scripts/bash/functions.sh

(cat data/tsv/genomes_analyzed_bac.tsv | cut -f 1 | tail -n +2 | proteins2tsv_bac; cat data/tsv/genomes_analyzed_ar.tsv | cut -f 1 | tail -n +2 | proteins2tsv_ar) | protein_names2tsv | sed '1s/^/accession\tprodigal_acc\tprodigal_start\tprodigal_end\tprodigal_strand\t\tprodigal_ID\tprodigal_partial\tprodigal_start_type\tprodigal_rbs_motif\tprodigal_rbs_spacer\tprodigal_gc_cont\tprotein_seq\n/' > data/tsv/proteins_raw.tsv
