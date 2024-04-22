suppressMessages(suppressWarnings(library('tidyverse')))

taxonomicRanks <- c('Kingdom', 'Phylum', 'Class', 'Order',
                    'Family', 'Genus', 'Species')

# IMPORT GTDB TAXONOMY INFORMATION
gtdb_taxonomy <- read.delim('gtdb_taxonomy.tsv', header = F, sep = '\t') %>%
  mutate(V1 = str_replace(V1, '^[RG][SB]_', '')) %>%
  mutate(V1 = str_replace(V1, '\\.[0-9]$', '')) %>%
  mutate(V2 = str_replace_all(V2, ';' ,'|'))

# IMPORT OPERON INFORMATION (OPERONID AND ASSEMBLY ACCESSION ARE THE IMPORTANT COLUMNS)
complete_operons <- read.delim('Outputs_Complete/master_rrna.gff', header = T, sep = '\t') %>%
  select(seqid, OperonID) %>%
  separate(seqid, into = c('Assembly', 'SeqID'), sep = '__') %>%
  select(-SeqID)
incomplete_operons <- read.delim('Outputs_Incomplete/master_rrna.gff', header = T, sep = '\t') %>%
  select(seqid, OperonID) %>%
  separate(seqid, into = c('Assembly', 'SeqID'), sep = '__') %>%
  select(-SeqID)
combined_operons <- rbind(complete_operons,
                          incomplete_operons)

# WRITE taxFull FILES
complete_operons %>%
  left_join(gtdb_taxonomy, by = c('Assembly' = 'V1')) %>%
  select(OperonID, V2) %>%
  mutate(V2 = str_replace_all(V2, ';' ,'|')) %>%
  write.table('Outputs_Complete/taxFull.tsv', row.names = F, col.names = F, quote = F, sep = '\t')

combined_operons %>%
  left_join(gtdb_taxonomy, by = c('Assembly' = 'V1')) %>%
  select(OperonID, V2) %>%
  mutate(V2 = str_replace_all(V2, ';' ,'|')) %>%
  write.table('Outputs_Combined/taxFull.tsv', row.names = F, col.names = F, quote = F, sep = '\t')

##### ASSIGNING TAXONOMY TO NR99.9% CLUSTERS FROM COMBINED OUTPUT BASED ON THREE DIFFERENT SCHEMES
##### CAN SUBSET OUTPUT COMPLETE ONLY LATER
# import vsearch outputs and taxonomy info
vsearch_centroids <- read.delim('Outputs_Combined/vsearch_centroids.tsv', header = F, sep = '\t') %>%
  left_join(combined_operons, by = c('V9' = 'OperonID'))
vsearch_hits <- read.delim('Outputs_Combined/vsearch_hits.tsv', header = F, sep = '\t') %>%
  left_join(combined_operons, by = c('V9' = 'OperonID'))
vsearch <- rbind(vsearch_centroids, vsearch_hits) %>%
  select(V2, V9, Assembly) %>%
  rename(ClusterID = V2, OperonID = V9) %>%
  arrange(ClusterID)
tax <- gtdb_taxonomy %>%
  separate(V2, sep = '\\|', into = taxonomicRanks, remove = F) %>%
  rename(Assembly = V1, Taxonomy = V2)

# 1 - taxRep: RefSeq taxonomy of the cluster representative sequence
taxRep <- vsearch_centroids %>%
  select(V2, V9, Assembly) %>%
  rename(ClusterID = V2, OperonID = V9) %>%
  left_join(tax %>%
              select(Assembly, Taxonomy))
taxRep %>%
  select(OperonID, Taxonomy) %>%
  arrange(OperonID) %>%
  write.table('Outputs_Combined/taxRep.tsv', quote = F, row.names = F, col.names = F, sep = '\t')

# 2 - taxLCA: lowest common ancestor of all sequences in the cluster
# first identify all clusters with a 100% species consensus
# then, for the remaining clusters, move up through taxonomic ranks until their LCA has been identified
vsearch_tax <- vsearch %>%
  arrange(ClusterID, OperonID) %>%
  left_join(tax %>%
              select(-Taxonomy)) %>%
  select(-Assembly)

cons_species <- vsearch_tax %>%
  select(ClusterID, Species) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(SpeciesConsensus = (Count == 1)) %>%
  filter(SpeciesConsensus == TRUE) %>%
  select(-Count, -SpeciesConsensus) %>%
  mutate(LCA_Rank = 'Species')
cons_species <- cons_species %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  distinct()

cons_genus <- vsearch_tax %>%
  filter(!ClusterID %in% cons_species$ClusterID) %>%
  select(ClusterID, Genus) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(GenusConsensus = (Count == 1)) %>%
  filter(GenusConsensus == TRUE) %>%
  select(-Count, -GenusConsensus) %>%
  mutate(LCA_Rank = 'Genus')
cons_genus <- cons_genus %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Species = 'NA') %>%
  distinct()

cons_family <- vsearch_tax %>%
  filter(!ClusterID %in% c(cons_species$ClusterID, cons_genus$ClusterID)) %>%
  select(ClusterID, Family) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(FamilyConsensus = (Count == 1)) %>%
  filter(FamilyConsensus == TRUE) %>%
  select(-Count, -FamilyConsensus) %>%
  mutate(LCA_Rank = 'Family')
cons_family <- cons_family %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Genus = 'NA',
         Species = 'NA') %>%
  distinct()

cons_order <- vsearch_tax %>%
  filter(!ClusterID %in% c(cons_species$ClusterID, cons_genus$ClusterID, cons_family$ClusterID)) %>%
  select(ClusterID, Order) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(OrderConsensus = (Count == 1)) %>%
  filter(OrderConsensus == TRUE) %>%
  select(-Count, -OrderConsensus) %>%
  mutate(LCA_Rank = 'Order')
cons_order <- cons_order %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Family = 'NA',
         Genus = 'NA',
         Species = 'NA') %>%
  distinct()

cons_class <- vsearch_tax %>%
  filter(!ClusterID %in% c(cons_species$ClusterID, cons_genus$ClusterID, cons_family$ClusterID, cons_order$ClusterID)) %>%
  select(ClusterID, Class) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(ClassConsensus = (Count == 1)) %>%
  filter(ClassConsensus == TRUE) %>%
  select(-Count, -ClassConsensus) %>%
  mutate(LCA_Rank = 'Class')
cons_class <- cons_class %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Order = 'NA',
         Family = 'NA',
         Genus = 'NA',
         Species = 'NA') %>%
  distinct()

cons_phylum <- vsearch_tax %>%
  filter(!ClusterID %in% c(cons_species$ClusterID, cons_genus$ClusterID, cons_family$ClusterID, cons_order$ClusterID, cons_class$ClusterID)) %>%
  select(ClusterID, Phylum) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(PhylumConsensus = (Count == 1)) %>%
  filter(PhylumConsensus == TRUE) %>%
  select(-Count, -PhylumConsensus) %>%
  mutate(LCA_Rank = 'Phylum')
cons_phylum <- cons_phylum %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Class = 'NA',
         Order = 'NA',
         Family = 'NA',
         Genus = 'NA',
         Species = 'NA') %>%
  distinct()

cons_kingdom <- vsearch_tax %>%
  filter(!ClusterID %in% c(cons_species$ClusterID, cons_genus$ClusterID, cons_family$ClusterID, cons_order$ClusterID, cons_class$ClusterID, cons_phylum$ClusterID)) %>%
  select(ClusterID, Kingdom) %>%
  distinct() %>%
  group_by(ClusterID) %>%
  summarise(Count = n()) %>%
  mutate(KingdomConsensus = (Count == 1)) %>%
  filter(KingdomConsensus == TRUE) %>%
  select(-Count, -KingdomConsensus) %>%
  mutate(LCA_Rank = 'Kingdom')
cons_kingdom <- cons_kingdom %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Phylum = 'NA',
         Class = 'NA',
         Order = 'NA',
         Family = 'NA',
         Genus = 'NA',
         Species = 'NA') %>%
  distinct()

cons_unknown <- vsearch_tax %>%
  filter(!ClusterID %in% c(cons_species$ClusterID, cons_genus$ClusterID, cons_family$ClusterID, cons_order$ClusterID, cons_class$ClusterID, cons_phylum$ClusterID, cons_kingdom$ClusterID)) %>%
  select(ClusterID) %>%
  distinct() %>%
  mutate(LCA_Rank = 'Unknown')
cons_unknown <- cons_unknown %>%
  left_join(vsearch_tax) %>%
  select(-OperonID) %>%
  mutate(Kingdom = 'NA',
         Phylum = 'NA',
         Class = 'NA',
         Order = 'NA',
         Family = 'NA',
         Genus = 'NA',
         Species = 'NA') %>%
  distinct()

taxLCA <- rbind(cons_species,
                cons_genus,
                cons_family,
                cons_order,
                cons_class,
                cons_phylum,
                cons_kingdom,
                cons_unknown) %>%
  left_join(vsearch_centroids %>%
              select(V2, V9),
            by = c('ClusterID' = 'V2')) %>%
  rename(OperonID = V9) %>%
  relocate(OperonID, .after = ClusterID)

taxLCA %>%
  select(-c('ClusterID', 'LCA_Rank')) %>%
  unite(Taxonomy, Kingdom:Species, sep = '\\|') %>%
  arrange(OperonID) %>%
  write.table('Outputs_Combined/taxLCA.tsv', quote = F, row.names = F, col.names = F, sep = '\t')

taxLCA %>%
  group_by(LCA_Rank) %>%
  summarise(Count = n()) %>%
  arrange(-Count) %>%
  mutate(Percentage = 100*(Count/sum(Count))) %>%
  write.table('Outputs_Combined/taxLCA_ranksummary.tsv', quote = F, row.names = F, sep = '\t')

vsearch_tax %>%
  filter(ClusterID %in% (taxLCA %>%
                           filter(LCA_Rank != 'Species') %>%
                           pull(ClusterID))) %>%
  group_by(ClusterID, Species) %>%
  summarise(OperonCount = n()) %>%
  write.table('Outputs_Combined/taxLCA_conflictingspecies.tsv', quote = F, row.names = F, sep = '\t')

# 3 - taxMaj: lowest taxonomic rank at which there is a single majority agreement of all sequences in the cluster
# first identify all clusters with a species majority
# then, for the remaining clusters, move up through taxonomic ranks until their LCA has been identified
maj_species <- vsearch_tax %>%
  select(ClusterID, Species) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Species')

maj_genus <- vsearch_tax %>%
  filter(!ClusterID %in% maj_species$ClusterID) %>%
  select(ClusterID, Genus) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Genus')

maj_family <- vsearch_tax %>%
  filter(!ClusterID %in% c(maj_species$ClusterID, maj_genus$ClusterID)) %>%
  select(ClusterID, Family) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Family')

maj_order <- vsearch_tax %>%
  filter(!ClusterID %in% c(maj_species$ClusterID, maj_genus$ClusterID, maj_family$ClusterID)) %>%
  select(ClusterID, Order) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Order')

maj_class <- vsearch_tax %>%
  filter(!ClusterID %in% c(maj_species$ClusterID, maj_genus$ClusterID, maj_family$ClusterID, maj_order$ClusterID)) %>%
  select(ClusterID, Class) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Class')

maj_phylum <- vsearch_tax %>%
  filter(!ClusterID %in% c(maj_species$ClusterID, maj_genus$ClusterID, maj_family$ClusterID, maj_order$ClusterID, maj_class$ClusterID)) %>%
  select(ClusterID, Phylum) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Phylum')

maj_kingdom <- vsearch_tax %>%
  filter(!ClusterID %in% c(maj_species$ClusterID, maj_genus$ClusterID, maj_family$ClusterID, maj_order$ClusterID, maj_class$ClusterID, maj_phylum$ClusterID)) %>%
  select(ClusterID, Kingdom) %>%
  group_by_all() %>%
  summarise(Count = n()) %>%
  slice_max(Count, n = 1) %>%
  group_by(ClusterID) %>%
  filter(n() == 1) %>%
  select(-Count) %>%
  mutate(Maj_Rank = 'Kingdom')

taxMaj <- rbind(maj_species,
                maj_genus,
                maj_family,
                maj_order,
                maj_class,
                maj_phylum,
                maj_kingdom) 

taxMaj <- vsearch_tax %>%
  unite(Taxonomy, Kingdom:Species, sep = '|') %>%
  select(OperonID, ClusterID, Taxonomy) %>%
  mutate(tempTaxonomy = if_else(ClusterID %in% maj_species$ClusterID, Taxonomy,
                                if_else(ClusterID %in% maj_genus$ClusterID, str_replace(Taxonomy, '\\|s__.*', ''),
                                        if_else(ClusterID %in% maj_family$ClusterID, str_replace(Taxonomy, '\\|g__.*', ''),
                                                if_else(ClusterID %in% maj_order$ClusterID, str_replace(Taxonomy, '\\|f__.*', ''),
                                                        if_else(ClusterID %in% maj_class$ClusterID, str_replace(Taxonomy, '\\|o__.*', ''),
                                                                if_else(ClusterID %in% maj_phylum$ClusterID, str_replace(Taxonomy, '\\|c__.*', ''),
                                                                        str_replace(Taxonomy, '\\|p__.*', '')))))))) %>%
  select(OperonID, tempTaxonomy) %>%
  filter(OperonID %in% vsearch_centroids$V9)

taxMaj %>%
  arrange(OperonID) %>%
  write.table('Outputs_Combined/taxMaj.tsv', quote = F, row.names = F, col.names = F, sep = '\t')

data.frame(Rank = taxonomicRanks[1:7],
           Count = c(dim(maj_kingdom)[1],
                     dim(maj_phylum)[1],
                     dim(maj_class)[1],
                     dim(maj_order)[1],
                     dim(maj_family)[1],
                     dim(maj_genus)[1],
                     dim(maj_species)[1])) %>%
  mutate(Percentage = 100*(Count/sum(Count))) %>%
  write.table('Outputs_Combined/taxMaj_ranksummary.tsv', quote = F, row.names = F, sep = '\t')

# combining all taxonomy systems for NR database into a single dataframe
taxCombined <- taxRep %>%
  rename(taxRep = Taxonomy) %>%
  left_join(taxLCA %>%
              unite(taxLCA, Kingdom:Species, sep = '|')) %>%
  left_join(taxMaj %>%
              rename(taxMaj = tempTaxonomy)) %>%
  arrange(OperonID)

taxCombined %>%
  write.table('Outputs_Combined/taxCombined.tsv', quote = F, row.names = F, sep = '\t')

### NOW EXTRACTING TAXONOMY INFO FOR OPERONS IN THE COMPLETE DATABASE

complete_taxCombined <- taxCombined %>%
  filter(OperonID %in% complete_operons$OperonID) %>%
  arrange(OperonID)

complete_taxCombined %>%
  write.table('Outputs_Complete/taxCombined.tsv', quote = F, row.names = F, sep = '\t')

complete_taxCombined %>%
  select(OperonID, taxRep) %>%
  write.table('Outputs_Complete/taxRep.tsv', quote = F, row.names = F, col.names = F, sep = '\t')
complete_taxCombined %>%
  select(OperonID, taxLCA) %>%
  write.table('Outputs_Complete/taxLCA.tsv', quote = F, row.names = F, col.names = F, sep = '\t')
complete_taxCombined %>%
  select(OperonID, taxMaj) %>%
  write.table('Outputs_Complete/taxMaj.tsv', quote = F, row.names = F, col.names = F, sep = '\t')


