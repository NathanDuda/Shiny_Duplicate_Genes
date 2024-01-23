

calculate_tau <- function(expression_path, output_path){
  
  expression <- read.csv(paste0(expression_path,'/Expression_Data.tsv'),sep='')
  colnames(expression)[1] <- 'gn'
  
  rownames(expression) <- expression$gn
  expression <- expression %>% select(-gn)
  
  tau <- calcTau(expression)
  
  tau$gn <- rownames(tau)
  tau <- tau[c('gn','tau')]
  
  # write tau to file
  mkdir(paste0(output_path,'Tau/'))
  write.table(tau, file= paste0(output_path,'Tau/All_Tau.tsv'))
  
  
  # get tau for dups 
  #dups <- read.csv(paste0(output_path,"/Dup_Pairs_Ancestral.tsv"), sep="")
  
  # get tau for orthologs 
  # 
  
}



v <- function(){
  ggplot(all_tau, aes(x=factor(category, levels = func_levels),y=tau))+ 
    # facet_grid(cols=vars(variable)) +
    ylab('Tau') +
    geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
    theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
          #axis.text.x=element_blank(),
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 20),
          strip.background =element_rect(color='black',fill="white")) +
    ylim(0,1) +
    #geom_jitter(size=0.001) + 
    stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                       ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
    geom_hline(yintercept=median((all_tau[all_tau$category=='All_Ortho',])['tau']$tau),
               linetype="dashed", color = '#778899',linewidth=0.5) 
}













