


functionalization <- function(output_path, expression_path) {
  
  
  dups <- read.csv(paste0(output_path,"/Dup_Pairs_Ancestral.tsv"), sep="")
  ortho_pairs <- read.csv(paste0(output_path,"/Ortholog_Pairs.tsv"), sep="")
  raw_exp <- read.csv(paste0(expression_path,"/Expression_Data.tsv"), sep="")
  colnames(raw_exp)[1] <- 'gn'
  
  
  # calculate relative expression values 
  exp <- raw_exp
  rownames(exp) <- exp$gn
  exp <- exp %>% select(-gn)
  exp <- exp / rowSums(exp)
  exp$gn <- rownames(exp)
  exp <- exp[, c("gn", setdiff(names(exp), "gn"))]
  
  # get exp for dups and ancestral
  dup_1_exp <- exp %>% rename_at(-1, ~paste('dup_1_', ., sep = ''))
  dup_anc_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'gn', by.y ='dup_1') %>% rename(dup_1 = 'gn')
  
  dup_2_exp <- exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
  dup_anc_exp <- dup_anc_exp %>% merge(dup_2_exp, ., by.x = 'gn', by.y ='dup_2') %>% rename(dup_2 = 'gn')
  
  anc_exp <- exp %>% rename_at(-1, ~paste('anc_', ., sep = ''))
  dup_anc_exp <- dup_anc_exp %>% merge(anc_exp, ., by.x = 'gn', by.y ='ancestral_copy') %>% rename(anc = 'gn')
  
  # get exp for ortholog pairs 
  ortho_x_exp <- exp %>% rename_at(-1, ~paste('ortho_x_', ., sep = ''))
  ortho_pair_exp <- ortho_pairs %>% merge(ortho_x_exp, ., by.x = 'gn', by.y ='gn.x') %>% rename(ortho_x = 'gn')
  
  ortho_y_exp <- exp %>% rename_at(-1, ~paste('ortho_y_', ., sep = ''))
  ortho_pair_exp <- ortho_pair_exp %>% merge(ortho_y_exp, ., by.x = 'gn', by.y ='gn.y') %>% rename(ortho_y = 'gn') %>%
    select(-species.x, -species.y)
  
  # get combined expression values and after adding, calculate relative expression
  dup_1_exp <- raw_exp %>% rename_at(-1, ~paste('dup_1_', ., sep = ''))
  dups_combined_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'gn', by.y ='dup_1') %>% rename(dup_1 = 'gn')
  
  dup_2_exp <- raw_exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
  dups_combined_exp <- dups_combined_exp %>% merge(dup_2_exp, ., by.x = 'gn', by.y ='dup_2') %>% rename(dup_2 = 'gn')
  
  tissue_names <- colnames(exp[2:ncol(exp)]) 
  
  for (tissue in tissue_names) {
    dups_combined_exp <- dups_combined_exp %>%
      rowwise() %>%
      mutate("{tissue}" :=  sum(c_across(ends_with(tissue)), na.rm = TRUE))
  }
  
  dups_combined_exp <- dups_combined_exp %>% 
    select(Orthogroup,all_of(tissue_names))
  
  # calculate relative expression values of combined 
  dups_combined_exp <- as.data.frame(dups_combined_exp) 
  rownames(dups_combined_exp) <- dups_combined_exp$Orthogroup
  dups_combined_exp <- dups_combined_exp %>% select(-Orthogroup)
  dups_combined_exp <- dups_combined_exp / rowSums(dups_combined_exp)
  dups_combined_exp$Orthogroup <- rownames(dups_combined_exp)
  
  # add combined expression values to the dup_anc_exp dataframe  
  dups_combined_exp <- dups_combined_exp %>% rename_at(-ncol(dups_combined_exp), ~paste('d1d2_', ., sep = ''))
  dup_anc_exp <- merge(dups_combined_exp,dup_anc_exp,by='Orthogroup')
  
  # calculate euclidean distance values 
  ed_values <- dup_anc_exp %>%
    select(-dup_1, -dup_2, -anc) %>%
    mutate(dup1_a = sqrt(rowSums((select(., starts_with('dup_1')) - select(., starts_with('anc_'))) ^ 2))) %>%
    mutate(dup2_a = sqrt(rowSums((select(., starts_with('dup_2')) - select(., starts_with('anc_'))) ^ 2))) %>%
    mutate(d1d2_a = sqrt(rowSums((select(., starts_with('d1d2_')) - select(., starts_with('anc_'))) ^ 2))) %>%
    select(Orthogroup, dup1_a, dup2_a, d1d2_a)
  
  sc_ed_values <- ortho_pair_exp %>%
    select(-ortho_x, -ortho_y) %>%
    mutate(sc_ed = sqrt(rowSums((select(., starts_with('ortho_x')) - select(., starts_with('ortho_y'))) ^ 2))) %>%
    select(Orthogroup, sc_ed)
  
  
  # calculate the cutoff ed value 
  iqr <- IQR(sc_ed_values$sc_ed) / 2
  cutoff <- median(sc_ed_values$sc_ed) + iqr
  #cutoff <- cutoff /2
  
  
  # use euclidean distance values to classify into functional groups 
  func <- ed_values %>%
    rowwise() %>%
    mutate(func = case_when((dup1_a <= cutoff) & (dup2_a <= cutoff) ~ 'Conserved',
                            (dup1_a > cutoff) & (dup2_a <= cutoff) ~ 'Neofunctionalized', # neo dup 1
                            (dup1_a <= cutoff & dup2_a > cutoff) ~ 'Neofunctionalized', # neo dup 2 
                            (dup1_a > cutoff & dup2_a > cutoff & d1d2_a <= cutoff) ~ 'Subfunctionalized',
                            (dup1_a > cutoff & dup2_a > cutoff & d1d2_a > cutoff) ~ 'Specialized'))
  
  # add dup and anc gene ids to func 
  func <- merge(dups,func,by='Orthogroup')
  
  # write functionalization results to file
  write.table(func,file=paste0(output_path,'/Dup_Functionalizations.tsv'))
  
  return(table(func$func))
}


visualize_func <- function(func_counts) {
  func_counts <- as.data.frame(func_counts)
  plot_ly(func_counts, labels = ~Var1, values = ~Freq, type = 'pie', hole = 0.6,
          marker = list(colors = c("#3B7BBD", "#E23F51", "#6CBC4D", '#F18244'))) %>%
    layout(title = "Functionalization Categories", showlegend = TRUE)
}


