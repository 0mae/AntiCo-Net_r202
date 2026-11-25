source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR

# File size
du -hs data/tsv/proteins_raw_gtdbtax_cog2020.tsv
# All lines
tail -n +2 data/tsv/proteins_raw_gtdbtax_cog2020.tsv | wc -l
# Number of gene_id
tail -n +2 data/tsv/proteins_raw_gtdbtax_cog2020.tsv | cut -f 1 | sort | uniq | wc -l
# Number of genomes
tail -n +2 data/tsv/proteins_raw_gtdbtax_cog2020.tsv | cut -f 2 | sort | uniq | wc -l
# Number of cog_id
tail -n +2 data/tsv/proteins_raw_gtdbtax_cog2020.tsv | cut -f 27 | sort | uniq | wc -l
# Number of cog_ids
tail -n +2 data/tsv/proteins_raw_gtdbtax_cog2020.tsv | cut -f 28 | sort | uniq | wc -l
