rm -rf TaxonKit/

wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

mkdir -p TaxonKit/
mv names.dmp nodes.dmp delnodes.dmp merged.dmp TaxonKit/
rm -f citations.dmp division.dmp gc.prt gencode.dmp readme.txt taxdump.tar.gz
