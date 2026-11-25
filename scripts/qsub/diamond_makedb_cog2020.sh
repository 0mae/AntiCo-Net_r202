source /etc/profile.d/modules.sh
module load diamond/2.0.11
cd $PBS_O_WORKDIR

diamond makedb --in db/diamond_cog2020/cog-20.fa -d db/diamond_cog2020/cog-20 -p 12
