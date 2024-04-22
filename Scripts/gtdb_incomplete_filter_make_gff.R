suppressMessages(suppressWarnings(library('tidyverse')))

# read input gff, select informative columns, and edit rRNA gene names
rna_gff <- ape::read.gff('Outputs_Incomplete/combined_rrna.gff', GFF3 = T) %>%
  select(seqid, start, end, strand, attributes) %>%
  filter(!grepl('partial=true', attributes)) %>%
  mutate(attributes = str_replace(attributes, '.*product=', '')) %>%
  filter(attributes %in% c('16S ribosomal RNA', '23S ribosomal RNA',
                           '16S Ribosomal RNA', '23S Ribosomal RNA',
                           'ribosomal RNA-16S', 'ribosomal RNA-23S')) %>%
  mutate(rrna_gene = if_else(grepl('16S', attributes), '16S',
                             if_else(grepl('23S', attributes), '23S', 'Other'))) %>%
  filter(rrna_gene != 'Other')

# loop through rna_gff and, for every 16S gene, identify the closest 23S gene with the same seqid and strand values
# if strand is + then 23S should be after, if strand is - then 23S should be before
rowstokeep <- list()
for (i in 1:nrow(rna_gff)){
  if(rna_gff[i, 'rrna_gene'] != '16S'){next} # doesn't apply loop to rows encoding 23S genes
  if(rna_gff[i, 'strand'] == '+'){
    if(i == nrow(rna_gff)){break} # stops the loop trying to look for a 23S gene beyond the bottom of the table, would cause an error
    if(rna_gff[i+1, 'rrna_gene'] != '23S'){next} # ensures the next row is a 23S gene
    if(rna_gff[i+1, 'seqid'] != rna_gff[i, 'seqid']){next} # ensures the two genes are on the same contig/scaffold
    if(rna_gff[i+1, 'strand'] != '+'){next} # ensures the two genes are on the same strand
    if(rna_gff[i+1, 'rrna_gene'] == '23S'){rowstokeep <- append(rowstokeep, c(i, i+1))} # adds the row number to the list of rows to keep for the table of rRNA operons
    }
  if(rna_gff[i, 'strand'] == '-'){ # this loop does the same as the previous one but for 16S genes on the reverse strand
    if(i == 1){next} # stops the loop trying to look for a 23S gene beyond the top of the table, would cause an error
    if(rna_gff[i-1, 'rrna_gene'] != '23S'){next}
    if(rna_gff[i-1, 'seqid'] != rna_gff[i, 'seqid']){next}
    if(rna_gff[i-1, 'strand'] != '-'){next}
    if(rna_gff[i-1, 'rrna_gene'] == '23S'){rowstokeep <- append(rowstokeep, c(i-1, i))}
  }
}
# extract these rows into a data frame
operons <- rna_gff[unlist(rowstokeep),]

# assign unique operon identifiers
# identify linkage status (unlinked operons are those with an ITS longer than 1500bp)
# identify operons that cross the start/end boundary of the chromosome (start coordinate is 1 or end coordinate is equal to or larger than the sequence length)
operons <- operons %>%
  select(-attributes) %>%
  mutate(OperonID = paste0('U_', str_pad(rep(1:(length(seqid)/2), each = 2), width = 8, pad = 0))) %>%
  pivot_longer(cols = -c('OperonID', 'seqid', 'strand', 'rrna_gene'), names_to = 'End', values_to = 'Coordinate') %>%
  pivot_wider(id_cols = c('OperonID', 'seqid', 'strand'), names_from = c('End', 'rrna_gene'), values_from = Coordinate) %>%
  mutate(ITS_Length = if_else(strand == '+', start_23S - end_16S, start_16S - end_23S)) %>%
  mutate(Linkage = if_else(ITS_Length > 1500, 'Unlinked', 'Linked')) %>%
  left_join(read.delim('Outputs_Incomplete/seq_length.tsv', header = F, sep = '\t'),
            by = c('seqid' = 'V1')) %>%
  rename(Sequence_Length = V2) %>%
  rowwise() %>%
  mutate(Operon_Start = (min(start_23S, end_23S, start_16S, end_16S))) %>%
  mutate(Operon_End = (max(start_23S, end_23S, start_16S, end_16S))) %>%
  mutate(CrossBoundary = (Operon_Start == 1 | Operon_End >= Sequence_Length))

# write unlinked and boundary-crossing operons to file
operons %>%
  filter(CrossBoundary == TRUE | Linkage == 'Unlinked') %>%
  write.table('Outputs_Incomplete/operons_unlinked_or_crossingboundary.tsv', quote = F, row.names = F, sep = '\t')

# create a 'master gff' will all the info needed for the next two steps
operons <- operons %>%
  filter(CrossBoundary == FALSE & Linkage == 'Linked') %>%
  select(-c('ITS_Length', 'Linkage', 'Sequence_Length', 'Operon_Start', 'Operon_End', 'CrossBoundary')) %>%
  mutate(source = 'Dummy') %>%
  mutate(type = 'rRNA') %>%
  mutate(score = NA) %>%
  mutate(phase = NA) %>%
  mutate(start_ITS = if_else(strand == '+', end_16S, end_23S)) %>%
  mutate(end_ITS = if_else(strand == '+', start_23S, start_16S)) %>%
  relocate(seqid, OperonID, source, type, start_ITS, end_ITS, start_16S, end_16S, start_23S, end_23S, score, strand, phase)

# write the master gff to file for plotting later
operons %>%
  write.table('Outputs_Incomplete/master_rrna.gff', quote = F, row.names = F, sep = '\t')

# write the coordinates of each 'good quality' operon's ITS region to a gff-like file
operons %>%
  mutate(seqid_operon = paste0(seqid, ' ', OperonID)) %>%
  select(seqid_operon, source, type, start_ITS, end_ITS, score, strand, phase) %>%
  rename(start = start_ITS, end = end_ITS) %>%
  write.table('Outputs_Incomplete/full_ITS.gff', quote = F, row.names = F, col.names = F, sep = '\t')

# define the start and end point of each 'good quality' operon and format similar to GFF, write to file to allow extraction of FASTA sequence using bedtools
operons %>%
  mutate(seqid_operon = paste0(seqid, ' ', OperonID)) %>%
  mutate(start = min(start_23S, end_23S, start_16S, end_16S)) %>%
  mutate(end = max(start_23S, end_23S, start_16S, end_16S)) %>%
  select(seqid_operon, source, type, start, end, score, strand, phase) %>% 
  mutate(seqid_operon = str_replace(seqid_operon, ' .*', '')) %>% 
  write.table('Outputs_Incomplete/full_operons.gff', quote = F, row.names = F, col.names = F, sep = '\t')

# create a file with the unique operon identifiers that can be substituted into the fasta file downstream
operons %>%
  mutate(seqid_operon = paste0(seqid, ' ', OperonID)) %>%
  select(seqid_operon) %>%
  write.table('Outputs_Incomplete/operon_identifiers.txt', row.names = F, col.names = F, quote = F, sep = '\t')

