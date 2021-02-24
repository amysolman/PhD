#Script to plot simulation output of active and dormant community species richness

rm(list=ls())
graphics.off()

library("ggplot2")

i = 1

#READ IN SIMULATION
load(paste0(file="DormancySim_", i, ".rda"))

g <- ggplot(species_timeseries, aes(x=1:nrow(species_timeseries)*5000)) + 
  geom_point(aes(y = Active), shape = "triangle", color="red") +
  geom_point(aes(y = Dormant), shape = "square", color="blue") +
  ylab("Species Richness")+
  xlab("Timestep")+
  ggtitle("Neutral Model with Dormancy")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

file_name = paste0("DormancySimPlot_", i, ".pdf")
pdf(file_name)
print(g)
dev.off()
