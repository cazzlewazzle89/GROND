suppressMessages(suppressWarnings(library('tidyverse')))

genbank_accessions <- read_tsv('genbank_gtdb_accessions.txt', col_names = F)
genbank_assemblyinfo <- read_tsv('assembly_summary_genbank.txt', skip = 1, quote = "")

assembly_summary_genbank_gtdb <- genbank_accessions %>%
    left_join(genbank_assemblyinfo, by = c('X1' = '#assembly_accession')) %>%
    group_by(X1) %>% 
    slice_head(n = 1) %>%
    filter(is.na(ftp_path) == FALSE) %>%
    rename('Assembly_Accession' = 'X1')

assembly_summary_genbank_gtdb %>% 
    filter(assembly_level == 'Complete Genome') %>%
    write_tsv('assembly_summary_genbank_gtdb_complete.txt', quote = 'none', col_names = F)
assembly_summary_genbank_gtdb %>% 
    filter(assembly_level != 'Complete Genome') %>%
    write_tsv('assembly_summary_genbank_gtdb_incomplete.txt', quote = 'none', col_names = F)

refseq_accessions <- read_tsv('refseq_gtdb_accessions.txt', col_names = F)
refseq_assemblyinfo <- read_tsv('assembly_summary_refseq.txt', skip = 1, quote = "")

assembly_summary_refseq_gtdb <- refseq_accessions %>%
    left_join(refseq_assemblyinfo, by = c('X1' = '#assembly_accession')) %>%
    group_by(X1) %>% 
    slice_head(n = 1) %>%
    filter(is.na(ftp_path) == FALSE) %>%
    rename('Assembly_Accession' = 'X1')

assembly_summary_refseq_gtdb %>% 
    filter(assembly_level == 'Complete Genome') %>%
    write_tsv('assembly_summary_refseq_gtdb_complete.txt', quote = 'none', col_names = F)
assembly_summary_refseq_gtdb %>% 
    filter(assembly_level != 'Complete Genome') %>%
    write_tsv('assembly_summary_refseq_gtdb_incomplete.txt', quote = 'none', col_names = F)