

clean_orthofinder <- function(input_path,output_path){
  ### Orthologs: 
  orthogroups <- read.delim2(paste0(input_path,"/Orthogroups/Orthogroups.tsv"))
  
  # get the one to ones 
  one_to_ones <- orthogroups %>%
    filter_all(all_vars(!grepl(",", .))) %>% # remove the rows with commas (separate multiple genes in one species)
    mutate_all(~ifelse(. == "", NA, .)) %>% # change empty cells to NA 
    na.omit()
  
  # write the one to ones to file
  write.table(one_to_ones,file=paste0(output_path,'./One_to_Ones.tsv'))
  
  
  # import the number of genes each species has in each orthogroup 
  orthogroup_gene_count <- read.delim(paste0(input_path,"./Orthogroups/Orthogroups.GeneCount.tsv"))
  
  # get the two to one to ones to zeros 
  two_to_ones <- orthogroup_gene_count %>%
    filter(if_all(2:(ncol(orthogroup_gene_count) - 1), ~. <= 2)) %>%              # keep only pairs 
    filter(rowSums(select(., 2:(ncol(orthogroup_gene_count) - 1)) == 2) == 1) %>% # make sure there is only one species with 2 copies 
    filter(Total > 2)                                                             # remove duplicates without any orthologs 
  
  # add a column with the name of the species that has the duplication
  two_to_ones$duplicate_pair_species <- 
    apply(two_to_ones[, 2:(ncol(two_to_ones) - 1)], 1, function(x) {
      col_index <- which(x == 2)
      return(colnames(two_to_ones[, 2:(ncol(two_to_ones) - 1)])[col_index])})
  
  two_to_ones <- two_to_ones[c('Orthogroup','duplicate_pair_species')]
  two_to_ones$duplicate_pair_species <- gsub('_prot','',two_to_ones$duplicate_pair_species) # DELETE THIS 
  
  # merge back with the gene names 
  two_to_ones <- merge(orthogroups,two_to_ones,by='Orthogroup')
  
  write.table(two_to_ones, file=paste0(output_path,'/Dup_Pair_Orthologs.tsv'))
  
  
  
  ### Paralogs: 
  dups <- two_to_ones %>%
    rowwise() %>%
    
    # extract the duplicate pair genes
    mutate(duplicate_pair = toString(c_across(2:(ncol(two_to_ones) - 1))[grep(",", c_across(2:(ncol(two_to_ones) - 1)))])) %>%
    select(Orthogroup, duplicate_pair) %>%
    separate(duplicate_pair, into = c("dup_1", "dup_2"), sep = ", ") %>%
    
    # ordering so that dup_1 is always > dup_2 for easy error catching  
    mutate(temp_column = dup_1) %>%
    mutate(dup_1 = case_when(dup_1 < dup_2 ~ dup_1, dup_1 > dup_2 ~ dup_2)) %>%
    mutate(dup_2 = case_when(dup_1 < dup_2 ~ dup_2, dup_1 > dup_2 ~ temp_column)) %>%
    select(-temp_column)
  
  
  
  
  #######################
  # DELETE THISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  #dups <- dups[1:10,] 
  #two_to_ones <- two_to_ones[two_to_ones$Orthogroup %in% dups$Orthogroup,]
  #write.table(two_to_ones, file=paste0(output_path,'/Dup_Pair_Orthologs.tsv'))
  #########################################################################################################################
  
  
  
  
  write.table(dups, file=paste0(output_path,'/Dup_Pairs.tsv'))
  
  return(nrow(dups))
  
}

