
b='f'
for OUTPUT in $(ls *.ff)
do
  d=${OUTPUT%?}
  echo $d
  perl -pe '$_= uc($_)' $OUTPUT > $d
done
