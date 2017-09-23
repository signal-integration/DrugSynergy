load("constraints_vector")
source("compute_profile_means.R")
source("compute_minimum_delta.R")
source("compute_minimum_delta_analysis.R")
source("setPowerPointStyle.R")
source("simulate_from_means.R")
source("setPowerPointStyle.R")
source("match11.R")
source("fancy.frequency.plots3.R")
library(caret)
library(randomForest) 
setPowerPointStyle()

#load profiles
PROFCODES = read.table("profile_codes_v2.txt",header = TRUE,sep = "\t")
#load model
load("rf_model")

#load data
data_file = "TNF_IFN_1.csv"
my_data = read.csv(data_file,sep = '\t')

#remember to log-transform if not on log-scale

#experimental design
samples = (dim(my_data)[2]-2)/4
design = c(rep("0",samples), rep("X",samples),
           rep("Y",samples), rep("Y+X",samples))

my_data_means = t(apply(my_data[, -(1:2)], 1,
                        function(x) tapply(x, design, mean)))

my_data_deltas = apply(my_data_means, 1 ,function(x) 
  compute_minimum_delta_analysis(x, 0.5))

my_data_proc = my_data[which(my_data_deltas),]

#generate a little noise to avoid problem of zero variances
noise = replicate(12, rnorm(dim(my_data_proc[,-(1:2)])[1],sd=0.001))

#compute statistical features for each gene
features = t(apply(my_data_proc[,-(1:2)] + noise, 1, match11))

#use random forest model to compute the probability of each class for each gene
predicted_classes = predict(rf_model,
                          newdata = data.frame(features),
                          type = 'prob')


#get largest probability and corresponding class
max_match = t(apply(predicted_classes,1,function(x)  
  c(max(x),which.max(x))) )

#generate results
my_results = my_data_proc
my_results[c("prob_match", "class")] = max_match

head(my_results)


#visualization 
cases = cbind(max_match,PROFCODES[max_match[,2],11:12])

my_case=1
outcomes=rep(0,table(PROFCODES$case)[my_case])
outcomes[as.integer(names(table(cases[cases[,3] == my_case, 4])))]=table(cases[cases[, 3] == my_case, 4])
fancy.freq.plots3(my_case,outcomes)


#to generate and save all cases uncomment code below

# for (h in 1:17){
#   my_case=h
#   png(file = paste(paste("my_case",my_case,sep="_"),'.png',sep=""),
#     width = 400, height = 300)
#   outcomes=rep(0,table(PROFCODES$case)[my_case])
# 
# outcomes[as.integer(names(table(cases[cases[,3]==my_case,4])))]=table(cases[cases[,3]==my_case,4])
# 
#   fancy.freq.plots3(my_case,outcomes)
#   setPowerPointStyle()
#   dev.off()
# 
# }
