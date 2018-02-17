library(edgeR)
library(limma)
source("pre_process_RNA_seq.R")

combinatorial_data = read.csv("RNA_seq_SRP064561_8h.txt", sep = '\t')

samples = ncol(combinatorial_data[,-(1:2)])/4
design = factor(c(rep("CTRL", samples), 
                  rep("X", samples), 
                  rep("Y", samples), 
                  rep("YX", samples)))

#returns counts matrix which is filtered by counts and voom-transformed
combinatorial_data = pre_process_RNA_seq(combinatorial_data, design)
#
boxplot(combinatorial_data[,-c(1:2)])


#now we can apply the RF classifier
source("extract_limma_features_v1.R")
source("find_optimal_match.R")
source("compute_minimum_delta.R")
source("resolve_integration.R")
source("filter_results.R")
load("rf_model")

PROFCODES = read.table("profile_codes_v2.txt",header = TRUE,sep = "\t") 

results = resolve_integration(combinatorial_data, design, PROFCODES)

#filter results
filtered_results = filter_results(results)


#check emergent synergies
em = filtered_results[[1]][filtered_results[[1]]$prof_index == 3,]

#sort by fc
em = em[order(-em$min_fc), ]
k = 10
plot(as.numeric(em[k,3:14]) ~ design, main = em[k,2])
