process_annotation_DB = function(background_genes, DB_name, min_pval = 0.01,
                                 min_fisher_pval = 0.05, min_hits = 10){
  
  #REQUIRES
  #a set of differentially expressed genes in an X+Y experiment
  #a file name containing an annotation database
  
  #RETURNS
  #a list of enriched annotation terms with a distribution of 
  #synergistic, antagonistic, and additive genes that deviates from 
  #the background distributions corresponding to all differentially expressed genes
  
  background_table = table(background_genes$type)
  
  #reading DB
  pathways <- readLines(DB_name)
  pathway_list = strsplit(pathways, "\t")
  
  #get annotation terms with at least 10 genes
  pathway_list = pathway_list[sapply(pathway_list, length) >= 10]
  
  #initialize a list whose entries contain info on the annotation terms
  processed_pathways = list()
  
  #parameters for hypergeometric test
  pop_size = 20000
  sample_size = nrow(background_genes)
  
  for (k in 1:length(pathway_list)){
    
    pathway_name = pathway_list[[k]][1]
    pathway_genes = pathway_list[[k]][-(1:2)]
    genes_in_pathway = background_genes[background_genes$genes %in% pathway_genes,]
    
    all_hits = length(pathway_genes)
    hits = nrow(genes_in_pathway)
    pval = phyper(hits-1, all_hits, pop_size - all_hits, 
                  sample_size, lower.tail = FALSE)
    
    
    pathway_counts = table(genes_in_pathway$type)
    pathway_table = rep(0, 3)
    names(pathway_table) = c('A', 'N', 'P')
    pathway_table[names(pathway_counts)] = pathway_counts
    
    contingency_table = rbind(background_table, pathway_table)
    norm_contingency_table = sweep(contingency_table,1,
                                   rowSums(contingency_table),`/`)
    
    chisq_pval = chisq.test(contingency_table)$p.value
    fisher_pval = fisher.test(contingency_table)$p.value
    
    if (pval < min_pval & fisher_pval < min_fisher_pval & hits >= min_hits){
      
      processed_pathway = list(pathway_name, 
                               all_hits, hits, 
                               pval, chisq_pval, genes_in_pathway, 
                               contingency_table, norm_contingency_table, 
                               fisher_pval)
      
      processed_pathways = list.append(processed_pathways, processed_pathway)
      
    }
    
  } 
  
  return(processed_pathways)
  
}