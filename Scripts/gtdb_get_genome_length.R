suppressMessages(suppressWarnings(library('tidyverse')))

taxonomy <- read_delim('gtdb_taxonomy.tsv',
                       col_names = c('Assembly', 'Taxonomy'), delim = '\t') %>%
  mutate(Assembly = str_remove(Assembly, '^[RG][SB]_')) %>%
  mutate(Assembly = str_replace(Assembly, '\\.[0-9]$', '')) 

length_complete <- read_delim('Outputs_Complete/seq_length.tsv',
                          col_names = c('Assembly__SeqID', 'Length'), delim = '\t') %>%
  separate(Assembly__SeqID, into = c('Assembly', 'SeqID'), sep = '__')

length_incomplete <- read_delim('Outputs_Incomplete/seq_length.tsv',
                            col_names = c('Assembly__SeqID', 'Length'), delim = '\t') %>%
  separate(Assembly__SeqID, into = c('Assembly', 'SeqID'), sep = '__')

length_combined <- bind_rows(length_complete,
                             length_incomplete) %>%
  group_by(Assembly) %>%
  summarise(Length = sum(Length)) %>%
  left_join(taxonomy)

stats_complete <- length_combined %>%
  filter(Assembly %in% length_complete$Assembly) %>%
  mutate(Taxonomy = str_remove(Taxonomy, '\\|t__.*')) %>% 
  mutate(Lineage = Taxonomy) %>%
  separate_rows(Taxonomy, sep = ';') %>%
  mutate(Rank = str_remove(Taxonomy, '__.*') %>% toupper) %>%
  mutate(Rank = ordered(Rank, levels = c('D', 'P', 'C', 'O', 'F', 'G', 'S'))) %>%
  mutate(Lineage = if_else(Rank == 'D', str_remove(Lineage, ';p__.*'),
                           if_else(Rank == 'P', str_remove(Lineage, ';c__.*'),
                                   if_else(Rank == 'C', str_remove(Lineage, ';o__.*'),
                                           if_else(Rank == 'O', str_remove(Lineage, ';f__.*'),
                                                   if_else(Rank == 'F', str_remove(Lineage, ';g__.*'),
                                                           if_else(Rank == 'G', str_remove(Lineage, ';s__.*'), Lineage))))))) %>%
  group_by(Taxonomy, Rank, Lineage) %>% 
  summarise(GenomeCount_Complete = n(),
            GenomeLength_Mean_Complete = mean(Length),
            GenomeLength_Median_Complete = median(Length)) %>%
  arrange(Rank)

stats_combined <- length_combined %>%
  mutate(Taxonomy = str_remove(Taxonomy, '\\|t__.*')) %>% 
  mutate(Lineage = Taxonomy) %>%
  separate_rows(Taxonomy, sep = ';') %>%
  mutate(Rank = str_remove(Taxonomy, '__.*') %>% toupper) %>%
  mutate(Rank = ordered(Rank, levels = c('D', 'P', 'C', 'O', 'F', 'G', 'S'))) %>%
  mutate(Lineage = if_else(Rank == 'D', str_remove(Lineage, ';p__.*'),
                           if_else(Rank == 'P', str_remove(Lineage, ';c__.*'),
                                   if_else(Rank == 'C', str_remove(Lineage, ';o__.*'),
                                           if_else(Rank == 'O', str_remove(Lineage, ';f__.*'),
                                                   if_else(Rank == 'F', str_remove(Lineage, ';g__.*'),
                                                           if_else(Rank == 'G', str_remove(Lineage, ';s__.*'), Lineage))))))) %>%
  group_by(Taxonomy, Rank, Lineage) %>% 
  summarise(GenomeCount_All = n(),
            GenomeLength_Mean_All = mean(Length),
            GenomeLength_Median_All = median(Length)) %>%
  arrange(Rank)

stats_combined %>%
  left_join(stats_complete) %>%
  write_delim('stats_genomelength.tsv', quote = 'none', delim = '\t')
