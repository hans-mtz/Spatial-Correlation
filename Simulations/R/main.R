# Setting up packages and loading functions
library(RStata)
# Tim's iMac
options("RStata.StataPath"="/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
options("RStata.StataVersion" = 16)
# Hans' MacBook
# options("RStata.StataPath"="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE")
# options("RStata.StataVersion" = 18)

# load("glob_vars.RData")
load("Simulations/Products/functions.RData")

# List of models

model_list <- list(model_1,model_2, model_3, model_4)
# model_list <- list(model_3, model_4)

# Simulations

rej_r <- lapply(model_list, rejection_ratio)
rej_r

# Collecting simulations results

rej_freq_table<-data.frame(Model=1:4, do.call(rbind,rej_r))
rej_freq_table

# rej_freq_table_3_4<-data.frame(Model=1:4, do.call(rbind,rej_r))
# rej_freq_table_3_4

#Saving results
save(rej_freq_table,file="Simulations/Products/table.RData")
# save(rej_freq_table_3_4,file="Simulations/Products/table_3_4.RData")
