echo -ne "$2\t"
tail -n +2 data/$1/$2.txt | cut -f 2 | sort | uniq | perl -ne 's/\n/\t/; print' | perl -ne 's/\t$/\n/;print'
