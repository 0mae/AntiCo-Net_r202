source /etc/profile.d/modules.sh
cd $PBS_O_WORKDIR
source scripts/bash/functions.sh

ls data/pastml_cog2020_ar/*.txt | perl -ne 's/^.*\///;s/.txt$//;print' | xargs -I {} sh -c "bash scripts/bash/check_pastml_results.sh pastml_cog2020_ar {}" > data/txt/check_pastml_cog2020_ar.txt
ls data/pastml_cog2020_bac/*.txt | perl -ne 's/^.*\///;s/.txt$//;print' | xargs -I {} sh -c "bash scripts/bash/check_pastml_results.sh pastml_cog2020_bac {}" > data/txt/check_pastml_cog2020_bac.txt
ls data/pastml_ko_ar/*.txt | perl -ne 's/^.*\///;s/.txt$//;print' | xargs -I {} sh -c "bash scripts/bash/check_pastml_results.sh pastml_ko_ar {}" > data/txt/check_pastml_ko_ar.txt
ls data/pastml_ko_bac/*.txt | perl -ne 's/^.*\///;s/.txt$//;print' | xargs -I {} sh -c "bash scripts/bash/check_pastml_results.sh pastml_ko_bac {}" > data/txt/check_pastml_ko_bac.txt
