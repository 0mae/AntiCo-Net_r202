cut -f 2,13-19,27 $1 | \
  sed '1s/^/00_/; 1s/\t/\t00_/g' | \
  sort -k 1,1 -k 9,9 | \
  uniq -c | \
  sed 's/^ *//;1s/1 /count\t/; 1s/00_//g; s/ /\t/' \
  > $2
