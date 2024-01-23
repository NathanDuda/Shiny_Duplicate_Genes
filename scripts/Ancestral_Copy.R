
# this script also gets all ortholog pairs of each duplicate pair 





get_ancestral_copy <- function(input_path, output_path, expression_path){
  
  dups <- read.csv(paste0(output_path,'/Dup_Pairs.tsv'),sep='')
  orthologs <- read.csv(paste0(output_path,'/Dup_Pair_Orthologs.tsv'),sep='')
  expression <- read.csv(paste0(expression_path,'/Expression_Data.tsv'),sep='')
  colnames(expression)[1] <- 'gn'
  
  newick_tree <- ape::read.tree(paste0(input_path,'/Species_Tree/SpeciesTree_rooted.txt'))
  
  # clean the column names of the orthologs dataframe and replace empty cells with NAs
  orthologs <- orthologs %>% 
    rename_all(~gsub("_prot", "", .)) %>%  # DELETE THIS 
    mutate_all(~na_if(.,""))               # make empty cells NA 
  
  
  # create function to find the closest expressed ortholog to each duplicate pair
  find_closest_ortholog <- function(row,species,newick_tree) {
    
    # calculate phylogenetic distances between the given species and each tip  
    species_node <- which(newick_tree$tip.label == species)
    distances <- cophenetic(newick_tree)[species_node, ]
    distances <- distances[distances > 0] # remove itself from distance calculation 
    
    # set non-expressed genes to NA 
    row[!row %in% expression$gn] <- NA
    
    # remove the duplicate pair species and missing species from the possible top choices 
    exclude_species <- colnames(row)[apply(row, 2, function(x) all(is.na(x)))]
    exclude_species <- c(exclude_species, species)
    distances <- distances[setdiff(names(distances), exclude_species)]
    
    # pick the closest available tip to the species
    closest_species <- names(which.min(distances)) # find the closest species by minimum distance
    
    
    if (is.null(closest_species)) {return(NA)}
    
    # get the ortholog from that species 
    closest_gene <- row[[closest_species]]
    
    if (exists("closest_gene")) {return(closest_gene)}
    return(NA)
    
  }
  
  # apply the function to each row of the ortholog table 
  orthologs$ancestral_copy <- NA
  for (row_num in 1:nrow(orthologs)) {
    row <- orthologs[row_num,] 
    orthologs[row_num,'ancestral_copy'] <- find_closest_ortholog(row, species = row$duplicate_pair_species, newick_tree)
    if (exists("closest_gene")) {rm(closest_gene)}
  }
  
  # remove orthogroups without ancestral copies (the orthologs they had weren't expressed)
  ancestral_copy <- orthologs[c('Orthogroup','ancestral_copy')]
  ancestral_copy <- na.omit(ancestral_copy)
  
  # merge the ancestral copy with the duplicate pairs
  dups <- merge(dups,ancestral_copy,by='Orthogroup')
  
  
  # write the duplicate pairs with their ancestral copy to file
  write.table(dups,file=paste0(output_path,'Dup_Pairs_Ancestral.tsv'))

  
  
  
  
  ### all combinations of ortholog pairs
  all_ortholog_pairs <- orthologs[c(1:(ncol(orthologs) - 2))] %>%
    mutate_all(~ifelse(grepl(",", .), NA, .)) %>%
    pivot_longer(cols = c(2:(ncol(orthologs) - 2)))
  
  all_ortholog_pairs <- 
    merge(all_ortholog_pairs,all_ortholog_pairs,by='Orthogroup') %>%
    na.omit() %>%
    filter(name.x!=name.y) %>% 
  
    # keep only one of the reciprocal combination pairs  
    mutate(x_y = paste0(value.x,'_',value.y)) %>%
    mutate(y_x = paste0(value.y,'_',value.x)) %>%
    filter(!(x_y > y_x)) %>%
    select(-x_y,-y_x)
  
  
  colnames(all_ortholog_pairs) <- c('Orthogroup','species.x','gn.x','species.y','gn.y')
  
  
  # keep orthologs with expression data and are expressed in at least one tissue
  expressed_ortholog_pairs <- all_ortholog_pairs %>%
    filter((gn.x %in% expression$gn) & (gn.y %in% expression$gn))
  
  # write ortholog pairs to file
  write.table(expressed_ortholog_pairs,file = paste0(output_path,'./Ortholog_Pairs.tsv'))
}

