
qiime feature-classifier extract-reads \
      --i-sequences $1 \
      --p-f-primer CAGCMGCCGCGGTAA \
      --p-r-primer CCRAMCTGTCTCACGACG \
      --p-min-length 1000 \
      --o-reads reads_pairA.qza

qiime feature-classifier extract-reads \
      --i-sequences $1 \
      --p-f-primer AGRGTTTGATYHTGGCTCAG \
      --p-r-primer ACCRCCCCAGTHAAACT \
      --p-min-length 1000 \
      --o-reads reads_pairB.qza

qiime feature-classifier extract-reads \
      --i-sequences $1 \
      --p-f-primer AGRGTTTGATYHTGGCTCAG \
      --p-r-primer CCRAMCTGTCTCACGACG \
      --p-min-length 1000 \
      --o-reads reads_pairC.qza

