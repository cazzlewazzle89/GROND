# testing full-length primer sequences to see which pair gives the highest % of amplicons

for i in "RefSeq" "GTDB"
do
    in_silico_PCR.pl \
        -s Database_${i}_207/Full/full.fna.gz \
        -p primerpairs.txt \
        -l 10000 \
        > FilesForPlotting_207/${i}_primertest_description_full.txt \
        2> Fangorn_${i}_207/${i}_primertest_amplicons_full.fna

    in_silico_PCR.pl \
        -s Database_${i}_207/NR/nrRep.fna.gz \
        -p primerpairs.txt \
        -l 10000 \
        > FilesForPlotting_207/${i}_primertest_description_nrRep.txt \
        2> Fangorn_${i}_207/${i}_primertest_amplicons_nrRep.fna
done

mv FilesForPlotting_207/GTDB_primertest_description_full.txt FilesForPlotting_207/gtdb_primertest_amplicons_full.txt
mv FilesForPlotting_207/GTDB_primertest_description_nrRep.txt FilesForPlotting_207/gtdb_primertest_amplicons_nrRep.txt
mv FilesForPlotting_207/RefSeq_primertest_description_full.txt FilesForPlotting_207/refseq_primertest_amplicons_full.txt
mv FilesForPlotting_207/RefSeq_primertest_description_nrRep.txt FilesForPlotting_207/refseq_primertest_amplicons_nrRep.txt

mv Fangorn_GTDB_207/GTDB_primertest_amplicons_full.fna Fangorn_GTDB_207/gtdb_primertest_amplicons_full.fna
mv Fangorn_GTDB_207/GTDB_primertest_amplicons_nrRep.fna Fangorn_GTDB_207/gtdb_primertest_amplicons_nrRep.fna
mv Fangorn_RefSeq_207/RefSeq_primertest_amplicons_full.fna Fangorn_RefSeq_207/refseq_primertest_amplicons_full.fna
mv Fangorn_RefSeq_207/RefSeq_primertest_amplicons_nrRep.fna Fangorn_RefSeq_207/refseq_primertest_amplicons_nrRep.fna
