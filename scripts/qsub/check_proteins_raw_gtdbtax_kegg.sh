source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

# File size
du -hs data/tsv/proteins_raw_gtdbtax_ko.tsv
# All lines
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | wc -l
# Number of gene_id
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 1 | sort | uniq | wc -l
# Number of genomes
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 2 | sort | uniq | wc -l
# Number of ko_id
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 27 | sort | uniq | wc -l
# Number of ko_ids
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 28 | sort | uniq | wc -l
# Number of genes assigned to ko_cog
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 31 | grep -v -w "NA" | wc -l
# Number of assigned ko_cog
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 31 | grep -v -w "NA" | sort | uniq | wc -l
# Number of genes assigned to ko_ec
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 32 | grep -v -w "NA" | wc -l
# Number of assigned ko_ec
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 32 | grep -v -w "NA" | sort | uniq | wc -l
# Number of genes assigned to ko_path
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 33 | grep -v -w "NA" | wc -l
# Number of assigned ko_path
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 33 | grep -v -w "NA" | sort | uniq | wc -l
# Number of genes assigned to ko_rn
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 34 | grep -v -w "NA" | wc -l
# Number of assigned ko_rn
tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 34 | grep -v -w "NA" | sort | uniq | wc -l
