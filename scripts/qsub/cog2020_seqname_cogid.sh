source /etc/profile.d/modules.sh
module load R/4.1.2
cd $PBS_O_WORKDIR

Rscript --vanilla --slave scripts/r/cog2020_seqname_cogid.R
