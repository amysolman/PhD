#Friday 6th November
#Script to call simulation functions

rm(list=ls())
graphics.off()

source("NeutralModelDormancy.R")

#Read in the job number from the cluster. To do this your code should include a new variable iter
#and should start with the line:

#iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))  #this will not work on a local machine
#so to run locally we set the value of iter for ourselves

iter <- 1 

i <- iter

set.seed(i)

output_file_name = paste0("DormancySim_", i, ".rda")

RUN_MY_SIM(m=0.03, beta=0.04, DEATH=0.03, DISTURBANCE=0.04, A_Jmeta=10000, D_Jmeta=10000, A_J=100, D_J=100, v=0.001, wall_time = 10, output_file_name = output_file_name)
