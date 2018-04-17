library(plotly)
library(enrichR)
library(DT)
library(memisc)
library(data.table)
library(shinyWidgets)
library(reshape2)
library(vembedr)
library(edgeR)
library(limma)

source("resolve_integration.R")
source("filter_data_on_deltas.R")
source("find_optimal_match.R")
source("visualize_all_profiles.R")
source("fancy.frequency.plots3.R")
source("compute_bliss.R")
source("compute_minimum_delta.R")
source("filter_results.R")
source("chooser.R")
source("filter_data.R")
source("do_pca.R")
source("pre_process_RNA_seq.R")
library("shinyjs")

#load data
#load("IFN_TNF_1h")
load("rf_model")
immune_DB = read.csv('immune_selection/DB_immune_selection.txt', sep = '\t')
immune_DB = immune_DB[,1:12]
immune_DB = immune_DB[,-9]

#read profile definitions
PROFCODES = read.table("profile_codes_v2.txt", header = TRUE,sep = "\t") 

img_link = uiOutput('case_link')

#DB for enrichment analysis
dbs <- c("GO_Biological_Process_2017b", "NCI-Nature_2016",
         "WikiPathways_2016", "Reactome_2016")

options(shiny.maxRequestSize = 30 * 1024 ^ 2)
options(error = expression(NULL))
options(warn=-1) 
