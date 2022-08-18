############################
### Simulated parameters ###
############################

library(tidyverse)
library(latex2exp)
library(xtable)

### Load simulated data

mean <- read_rds("output/sim_study_simulated_parameters.rds") %>% 
  unnest(cols = c("perf")) %>%
  #  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  #  filter(category == "time") %>% 
  filter(type == "design_time" | type=="dirder_max" | type == "dual_gap"| type == "efficiency"| type == "support") %>% 
  group_by(sig,type) %>% 
  summarise(mean_value = mean(value))

### Write results in table

results = matrix(nrow=5,ncol = 4)
colnames(results)=c("0.5","1","1.5","2")
rownames(results)=c("Directional derivatives", "Duality gap","Efficiency","Support","Design time")

for (i in 0:3) {
  results[5,i+1]= mean$mean_value[(5*i+1)]
  results[1,i+1]= mean$mean_value[(5*i+2)]
  results[2,i+1]= mean$mean_value[(5*i+3)]
  results[3,i+1]= mean$mean_value[(5*i+4)]
  results[4,i+1]= mean$mean_value[(5*i+5)]
}

results
xtable(results,digits=4)
xtable(results,display = c("s",rep("e",4)))
