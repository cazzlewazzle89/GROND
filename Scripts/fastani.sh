
ls Genomes_Complete/*.fna.gz > genomelist.txt
split -l 2000 genomelist.txt

for i in $(ls xa*)
do
    fastANI --ql ${i} --rl genomelist.txt -t 20 -o fastani_complete_${i}.tsv
done

cat fastani_complete_xa*.tsv > fastani_complete.tsv
rm -f fastani_complete_xa*.tsv
rm -f xa*

sed -i 's/Genomes_Complete\///g ; s/\.fna\.gz//g' fastani_complete.tsv

