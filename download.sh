#!/bin/sh

DB_SIZE=$1
TAX_SYSTEM=$2
TAX_METHOD=$3

curl https://zenodo.org/records/10889037/files/gtdb207full_full.fna.gz -o gtdb207full_full.fna.gz
curl https://zenodo.org/records/10889037/files/gtdb207full_taxFull.tsv -o gtdb207full_taxFull.tsv