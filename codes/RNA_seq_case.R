library(recount)
library(edgeR)
library(limma)

#here extract RNA-seq dataset (counts)

abstract_search("synergistic")

#GRO-seq
SRP_id = 'SRP064561'
if(!file.exists(file.path(SRP_id, "rse_gene.Rdata"))) {
  download_study(SRP_id, type = "rse-gene")
}

load(file.path(SRP_id, "rse_gene.Rdata"))



#extract samples to have standard X+Y design
samples_names = colData(rse_gene)$title
samples_names = gsub("-rep[1-3]", "", samples_names)
table(samples_names)


sample_idx = c(which(samples_names == 'MOCK_DMSO_8'),
               which(samples_names == 'UVC_DMSO_8'),
               which(samples_names == 'MOCK_TPA_8'),
               which(samples_names == 'UVC_TPA_8'))


samples = 3
design = factor(c(rep("CTRL", samples), rep("X", samples), rep("Y", samples), rep("YX", samples)))



#subset summarized experiment
rse_gene = rse_gene[, sample_idx]

#extract count_matrix
count_matrix = assay(rse_gene)
isexpr = rowSums(cpm(count_matrix)>1) >= 4
rse_gene = rse_gene[isexpr, ]

#row_ids = as(rowData(rse_gene)$symbol, "data.frame")[, c("group_name", "value")]

#filter genes with low expression

count_matrix_f = assay(rse_gene)


# define design matrix for limma
design_matrix <- model.matrix(~ 0 + design)
colnames(design_matrix)<- gsub("group","",colnames(design))


# here transform counts to apply limma
nf <- calcNormFactors(count_matrix_f)

y <- voom(count_matrix_f, design_matrix, lib.size = colSums(count_matrix_f)*nf)
counts_voom <- y$E


a = as(rowData(rse_gene)[3], "data.frame")
counts_voom = data.frame(probe = a,
                         genes = a,
                         counts_voom)

row.names(counts_voom) = NULL

#now we can apply the RF classifier
source("extract_limma_features_v1.R")
source("find_optimal_match.R")
source("compute_minimum_delta.R")
source("resolve_integration.R")

load("rf_model")

PROFCODES = read.table("profile_codes_v2.txt",header = TRUE,sep = "\t") 

results_1h_list = resolve_integration(counts_voom, design, PROFCODES)

filter_deg = results_1h_list[[1]]$adjusted_pvals < 0.05
results_1h = results_1h_list[[1]][filter_deg, ]
features_1h = results_1h_list[[2]][filter_deg, ]

min_delta = vector()
for (k in 1:sum(filter_deg)){
  
  min_delta[k] = compute_minimum_delta(as.numeric(features_1h[k, 2:5]),
                                       PROFCODES,
                                       results_1h$prof_index[k])
  
}

results_1h$min_delta = min_delta
res_f = results_1h[results_1h$min_delta > 0.3, ]
plot(as.numeric(res_f[which(res_f$prof_index == 1)[9], 3:14]) ~ design)
