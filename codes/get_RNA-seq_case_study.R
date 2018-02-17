library(recount)
library(edgeR)
library(limma)


#here extract RNA-seq dataset (counts)
abstract_search("synergistic")

SRP_id = 'SRP064561'

download_study(SRP_id, type = "rse-gene")

load(file.path(SRP_id, "rse_gene.Rdata"))


#extract samples to have standard X+Y design
samples_names = colData(rse_gene)$title
samples_names = gsub("-rep[1-3]", "", samples_names)

#check samples for each condition (choosing 8 hours), 
#there are 3 replicates per condition
table(samples_names)

sample_idx = c(which(samples_names == 'MOCK_DMSO_8'),
               which(samples_names == 'UVC_DMSO_8'),
               which(samples_names == 'MOCK_TPA_8'),
               which(samples_names == 'UVC_TPA_8'))

samples = 3
design = factor(c(rep("CTRL", samples), 
                  rep("X", samples), 
                  rep("Y", samples), 
                  rep("YX", samples)))

#subset summarized experiment
rse_gene = rse_gene[, sample_idx]

#extract count_matrix
counts_matrix = assay(rse_gene)
colnames(counts_matrix) = make.names(samples_names[sample_idx], unique = T)
symbols_ids = as(rowData(rse_gene)[3], "data.frame")

counts_df = data.frame(ensembl = row.names(counts_matrix),
                       genes = symbols_ids,
                       counts_matrix)

write.table(counts_df, file = "RNA_seq_SRP064561_8h.txt", row.names = F, sep = '\t')