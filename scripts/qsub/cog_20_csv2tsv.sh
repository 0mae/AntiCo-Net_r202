source /etc/profile.d/modules.sh 
cd $PBS_O_WORKDIR 

#perl scripts/perl/cog_20_csv2tsv.pl db/cog2020_fixed/cog-20.wrong_names.seq_csv.tsv db/cog2020/cog-20.cog.csv > db/cog2020_fixed/cog-20.cog.tsv
perl scripts/perl/cog_20_csv2tsv.pl db/cog2020_fixed/cog-20.wrong_names.seq_csv.tsv db/cog2020/cog-20.cog.csv db/cog2020_fixed/cog-20.wrong_columns.geneid_ncol.tsv > db/cog2020_fixed/cog-20.cog.tsv
