suppressMessages(suppressWarnings(library('tidyverse')))

taxonomy <- read_delim('assemblyaccession_taxid_lineage_mpa.txt',
                       col_names = c('Assembly', 'Taxonomy'), delim = '\t')

rrn_complete <- read_delim('Outputs_Complete/master_rrna.gff', delim = '\t') %>%
  select(seqid, OperonID) %>%
  separate(seqid, into = c('Assembly', 'SeqID'), sep = '__')

rrn_incomplete <- read_delim('Outputs_Incomplete/master_rrna.gff', delim = '\t') %>%
  select(seqid, OperonID) %>%
  separate(seqid, into = c('Assembly', 'SeqID'), sep = '__')

rrn_combined <- bind_rows(rrn_complete,
                 rrn_incomplete) %>%
  group_by(Assembly) %>%
  summarise(CopyNumber = n()) %>%
  left_join(taxonomy)

stats_complete <- rrn_combined %>%
  filter(Assembly %in% rrn_complete$Assembly) %>%
  mutate(Taxonomy = str_remove(Taxonomy, '\\|t__.*')) %>% 
  mutate(Lineage = Taxonomy) %>%
  separate_rows(Taxonomy, sep = '\\|') %>%
  mutate(Rank = str_remove(Taxonomy, '__.*') %>% toupper) %>%
  mutate(Rank = ordered(Rank, levels = c('K', 'P', 'C', 'O', 'F', 'G', 'S'))) %>%
  mutate(Lineage = if_else(Rank == 'K', str_remove(Lineage, '\\|p__.*'),
                           if_else(Rank == 'P', str_remove(Lineage, '\\|c__.*'),
                                   if_else(Rank == 'C', str_remove(Lineage, '\\|o__.*'),
                                           if_else(Rank == 'O', str_remove(Lineage, '\\|f__.*'),
                                                   if_else(Rank == 'F', str_remove(Lineage, '\\|g__.*'),
                                                           if_else(Rank == 'G', str_remove(Lineage, '\\|s__.*'), Lineage))))))) %>%
  group_by(Taxonomy, Rank, Lineage) %>% 
  summarise(GenomeCount_Complete = n(),
            CopyNumber_Mean = mean(CopyNumber),
            CopyNumber_Median = median(CopyNumber)) %>%
  arrange(Rank)

stats_combined <- rrn_combined %>%
  mutate(Taxonomy = str_remove(Taxonomy, '\\|t__.*')) %>% 
  mutate(Lineage = Taxonomy) %>%
  separate_rows(Taxonomy, sep = '\\|') %>%
  mutate(Rank = str_remove(Taxonomy, '__.*') %>% toupper) %>%
  mutate(Rank = ordered(Rank, levels = c('K', 'P', 'C', 'O', 'F', 'G', 'S'))) %>%
  mutate(Lineage = if_else(Rank == 'K', str_remove(Lineage, '\\|p__.*'),
                           if_else(Rank == 'P', str_remove(Lineage, '\\|c__.*'),
                                   if_else(Rank == 'C', str_remove(Lineage, '\\|o__.*'),
                                           if_else(Rank == 'O', str_remove(Lineage, '\\|f__.*'),
                                                   if_else(Rank == 'F', str_remove(Lineage, '\\|g__.*'),
                                                           if_else(Rank == 'G', str_remove(Lineage, '\\|s__.*'), Lineage))))))) %>%
  group_by(Taxonomy, Rank, Lineage) %>% 
  summarise(GenomeCount_All = n(),
            CopyNumber_Mean_All = mean(CopyNumber),
            CopyNumber_Median_All = median(CopyNumber)) %>%
  arrange(Rank)

stats_combined %>%
  left_join(stats_complete) %>%
  write_delim('stats_copynumber.tsv', quote = 'none', delim = '\t')
