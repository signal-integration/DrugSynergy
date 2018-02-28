cluster_annotation_terms = function(deg, inter_type, dbs, n = 20){
  
  library(enrichR)
  library(ggplot2)
  
  joined = data.frame()
  gene_set = deg[deg$type == inter_type, 'genes']
  enriched_list = lapply(enrichr(gene_set, dbs), function(x) x[order(x$P.value)[1:n],])
  enriched <- do.call("rbind", enriched_list)
  enriched = enriched[,c("Term", "Overlap", "P.value", "Genes")]
  #enriched = enriched[order(enriched$P.value)[1:n], 
  #                    c("Term", "Overlap", "P.value")]
  
  enriched$size = length(gene_set)
  enriched$gene_set = 'all'
  joined = rbind(joined, enriched)
  
  for (h in c(1, 4, 5)){
    
    gene_set = deg[deg$type == inter_type & deg$case == h, 'genes']
    
    if (h == 100){
      
      gene_set = deg[deg$type == inter_type & deg$case != c(2,3), 'genes']
    }
    
    enriched_list = lapply(enrichr(gene_set, dbs), function(x) x[order(x$P.value)[1:n],])
    enriched <- do.call("rbind", enriched_list)
    enriched = enriched[,c("Term", "Overlap", "P.value", "Genes")]
    
    enriched$size = length(gene_set)
    enriched$gene_set = as.character(h)
    joined = rbind(joined, enriched)
    joined$gene_set[joined$gene_set == 100] = 'high-val'
  }
  
  joined$Term = strtrim(joined$Term, 40)
  
  #manipulate dataframe
  joined$score = -log10(joined$P.value)
  joined$hits = sapply(strsplit(joined$Overlap, '/'), function(x) as.numeric(x[[1]]))
  joined$ratio = joined$hits/joined$size
  
  joined = na.omit(joined)
  joined = joined[joined$hits >=2, ]
  joined = joined[joined$P.value <=0.05, ]
  
  joined$gene_set = relevel(factor(joined$gene_set), "all")

  return(joined)  
  
  }