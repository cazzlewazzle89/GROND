suppressMessages(suppressWarnings(library('tidyverse')))

# comparing downloaded sequence and annotation filenames for complete genomes
# identifying assemblies which have sequence data but no annotations
# these will be annotated with barrnap to ID rRNA genes
genomes_complete <- list.files(path = 'Genomes_Complete/') %>% 
    str_replace(., '\\.fna\\.gz', '')
annotation_complete <- list.files(path = 'Annotation_Complete/') %>% 
    str_replace(., '\\.gff\\.gz', '')

genomes_complete[!genomes_complete %in% annotation_complete] %>% 
    as.data.frame() %>%
    write_delim('missingannotation_complete.txt', quote = F, col_names = F)

# doing the same as above for incomplete assemblies
genomes_incomplete <- list.files(path = 'Genomes_Incomplete/') %>% 
    str_replace(., '\\.fna\\.gz', '')
annotation_incomplete <- list.files(path = 'Annotation_Incomplete/') %>% 
    str_replace(., '\\.gff\\.gz', '')

genomes_incomplete[!genomes_incomplete %in% annotation_incomplete] %>% 
    as.data.frame() %>%
    write_delim('missingannotation_incomplete.txt', quote = F, col_names = F)
