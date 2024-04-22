
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads reads_pairA.qza \
  --i-reference-taxonomy $1 \
  --o-classifier classifier_pairA.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads reads_pairB.qza \
  --i-reference-taxonomy $1 \
  --o-classifier classifier_pairB.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads reads_pairC.qza \
  --i-reference-taxonomy $1 \
  --o-classifier classifier_pairB.qza

