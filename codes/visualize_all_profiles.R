visualize_all_profiles = function(my_data_filtered_matched){
  
  source("fancy_frequency_plots3.R")
  
  PROFCODES = read.table("profile_codes_v2.txt",header = TRUE,sep = "\t")
  
  max_match = my_data_filtered_matched[, c('score', 'prof_index')]
  
  cases = cbind(max_match, PROFCODES[max_match[, 2], 11:12])
  
  for (h in 1:17){
    
    my_case=h
    
    outcomes=rep(0,table(PROFCODES$case)[my_case])
    
    outcomes[as.integer(names(table(cases[cases[,3]==my_case,4])))]=table(cases[cases[,3]==my_case,4])

    fancy.freq.plots3(my_case,outcomes)
    
  }
  
  
  
  
  
  
}