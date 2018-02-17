pre_process_RNA_seq = function(combinatorial_data, design, min_count = 4){
  
  #filter by rowcounts
  isexpr = rowSums(cpm(combinatorial_data[, -(1:2)])>1) >= min_count
  combinatorial_data = combinatorial_data[isexpr, ]
  
  # define design matrix for limma
  design_matrix = model.matrix(~ 0 + design)
  colnames(design_matrix) = gsub("group","",colnames(design))
  
  # here transform counts to apply limma
  nf = calcNormFactors(combinatorial_data[-(1:2)])
  
  y = voom(combinatorial_data[-(1:2)], design_matrix, lib.size = colSums(combinatorial_data[-(1:2)])*nf)
  
  pre_processed_counts = y$E
  
  return(data.frame(combinatorial_data[,1:2], pre_processed_counts))
  
}