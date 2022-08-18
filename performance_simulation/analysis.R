#################################
#####  Performance analysis #####
#################################

library(tidyverse)
library(latex2exp)
library(xtable)

### Load simulated data

mean <- read_rds("output/sim_study_performance.rds") %>% 
  unnest(cols = c("perf")) %>%
  filter(type == "gamma_time" |type == "design_time" | type=="dirder_max" | type == "dual_gap") %>% 
  group_by(d,k,type) %>% 
  summarise(mean_value = mean(value))

### Write results in table

results = matrix(nrow=4,ncol = 8)
colnames(results)=c(("d=8,k=3"),("d=8,k=4"),("d=8,k=5"),("d=8,k=6"),
                    ("d=10,k=3"),("d=10,k=4"),("d=10,k=5"),("d=10,k=6"))
rownames(results)=c("Directional derivatives", "Duality gap","Gamma time","Design time")

for (i in 0:7) {
  results[4,i+1]= mean$mean_value[(4*i+1)]
  results[1,i+1]= mean$mean_value[(4*i+2)]
  results[2,i+1]= mean$mean_value[(4*i+3)]
  results[3,i+1]= mean$mean_value[(4*i+4)]
}

results
xtable(results)
xtable(results,display = c("s",rep("e",8)))

