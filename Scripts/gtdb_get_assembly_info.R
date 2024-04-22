suppressMessages(suppressWarnings(library('tidyverse')))

genbank_accessions <- read.delim('genbank_gtdb_accessions.txt', header = F)
genbank_assemblyinfo <- read.delim('assembly_summary_genbank.txt', header = T, sep = '\t', skip = 1, quote = "")

assembly_summary_genbank_gtdb <- genbank_accessions %>%
    left_join(genbank_assemblyinfo, by = c('V1' = 'X..assembly_accession')) %>%
    group_by(V1) %>% 
    slice_head(n = 1) %>%
    filter(is.na(ftp_path) == FALSE) %>%
    rename('Assembly_Accession' = 'V1')

assembly_summary_genbank_gtdb %>% 
    filter(assembly_level == 'Complete Genome') %>%
    write.table('assembly_summary_genbank_gtdb_complete.txt', quote = F, row.names = F, col.names = F, sep = '\t')
assembly_summary_genbank_gtdb %>% 
    filter(assembly_level != 'Complete Genome') %>%
    write.table('assembly_summary_genbank_gtdb_incomplete.txt', quote = F, row.names = F, col.names = F, sep = '\t')

refseq_accessions <- read.delim('refseq_gtdb_accessions.txt', header = F)
refseq_assemblyinfo <- read.delim('assembly_summary_refseq.txt', header = T, sep = '\t', skip = 1, quote = "")

assembly_summary_refseq_gtdb <- refseq_accessions %>%
    left_join(refseq_assemblyinfo, by = c('V1' = 'X..assembly_accession')) %>%
    group_by(V1) %>% 
    slice_head(n = 1) %>%
    filter(is.na(ftp_path) == FALSE) %>%
    rename('Assembly_Accession' = 'V1')

assembly_summary_refseq_gtdb %>%
    filter(assembly_level == 'Complete Genome') %>%    
    write.table('assembly_summary_refseq_gtdb_complete.txt', quote = F, row.names = F, col.names = F, sep = '\t')
assembly_summary_refseq_gtdb %>%
    filter(assembly_level != 'Complete Genome') %>%    
    write.table('assembly_summary_refseq_gtdb_incomplete.txt', quote = F, row.names = F, col.names = F, sep = '\t')