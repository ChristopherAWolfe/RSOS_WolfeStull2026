#################################################################################
# This is the 3rd of 6 scripts to complete the analyses in Wolfe and Stull 2026.#
# The following script imports posterior samples from a fit Stan model and      #
# prepares correlation plots. A second file 'response_vars.csv' contains        #
# a full list of response variables to aid in plotting.                         #
#################################################################################

# Necessary Packages
library(cmdstanr)
set_cmdstan_path("C:/cmdstan/cmdstan-2.33.1")
library(tidyverse)
library(posterior)
library(magrittr)
library(factoextra)
library(reshape2)

# Helper Functions for Correlation Matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

col <- c(rep("steelblue", 18), rep("tomato", 16), rep("goldenrod", 20))

# List of Response Variables for plotting
resp_vars <- data.frame(x = names(read.csv("response_vars.csv")))$x

# Import Posterior Draws of the Correlation Matrix
corr_draws <- read_cmdstan_csv(
  files = c("out_cop-1.csv", "out_cop-2.csv", "out_cop-3.csv", "out_cop-4.csv"),
  variables = c("corr_mat"), sampler_diagnostics = "")

# Summarise Draws
cor_sum <- summarise_draws(corr_draws$post_warmup_draws)

# Posterior Mean Correlation Matrix
cmat_mean <- matrix(cor_sum$mean, 54, 54, dimnames = list(resp_vars,resp_vars))

# Full Correlation Plot (Figure 3)
ggplot(melt(get_upper_tri(cmat_mean), na.rm = TRUE), 
       aes(Var1, Var2, fill=value)) + geom_tile(height=0.8, width=0.8, 
       color = "black") + scale_fill_gradient2(low="blue", mid="white", 
       high="red", midpoint = 0, limits = c(-1,1)) + theme_minimal() + 
       coord_equal() + labs(x="",y="",fill="Corr") + theme(axis.text.x = 
       element_text(size=5, angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0),
       colour = col, face = "bold"), axis.text.y=element_text(size=5, 
       margin = margin(0,-3,0,0), colour = col, face = "bold"), 
       panel.grid.major = element_blank(), panel.grid = 
       element_line(color="black"))


## Skeletal Growth
### Figure 5A

### Posterior Mean
cmat_diaphyseal <- cmat_mean[1:18,1:18]

ggplot(melt(get_upper_tri(cmat_diaphyseal), na.rm = TRUE), 
       aes(Var1, Var2, fill=value)) + geom_tile(height=0.8, width=0.8, 
       color = "black") + scale_fill_gradient2(low="blue", mid="white", 
       high="red", midpoint = 0, limits = c(-1,1)) + theme_minimal() + 
       coord_equal() + labs(x="",y="",fill="Corr") + theme(axis.text.x = 
       element_text(size=5, angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), 
       colour = "black", face = "bold"), axis.text.y=element_text(size=5, 
       margin=margin(0,-3,0,0), colour = "black", face = "bold"), 
       panel.grid.major=element_blank(), panel.grid = 
       element_line(color="black")) + geom_text(aes(label=round(value,2)), 
       size = 2) + scale_size(guide="none")


## Skeletal Development
### Figure 6A

### Posterior Mean
cmat_skel_dev <- cmat_mean[35:54,35:54]

ggplot(melt(get_upper_tri(cmat_skel_dev), na.rm = TRUE), 
       aes(Var1, Var2, fill=value)) +
       geom_tile(height=0.8, width=0.8, color = "black") + 
       scale_fill_gradient2(low="blue", mid="white", high="red", midpoint = 0, 
       limits = c(-1,1)) + theme_minimal() + coord_equal() + labs(x="", y="", 
       fill="Corr") + theme(axis.text.x=element_text(size=5, angle=45, vjust=1, 
       hjust=1, margin=margin(-3,0,0,0), colour = "black", face = "bold"), 
       axis.text.y=element_text(size=5, margin=margin(0,-3,0,0), 
       colour = "black", face = "bold"), panel.grid.major=element_blank(), 
       panel.grid = element_line(color="black")) + geom_text(aes(label = 
       round(value,2)), size = 2) + scale_size(guide="none")


## Dental Development
### Figure 7A

### Posterior Mean
cmat_dental_dev <- cmat_mean[19:34,19:34]

ggplot(melt(get_upper_tri(cmat_dental_dev), na.rm = TRUE), 
       aes(Var1, Var2, fill=value)) + geom_tile(height=0.8, width=0.8, 
       color = "black") + scale_fill_gradient2(low="blue", mid="white", 
       high="red", midpoint = 0, limits = c(-1,1)) + theme_minimal() +
       coord_equal() + labs(x="",y="",fill="Corr") + theme(axis.text.x = 
       element_text(size=5, angle=45, vjust=1, hjust=1, margin=margin(-3,0,0,0), 
       colour = "black", face = "bold"), axis.text.y=element_text(size=5, 
       margin = margin(0,-3,0,0), colour = "black", face = "bold"),
       panel.grid.major=element_blank(), panel.grid = 
       element_line(color="black")) + geom_text(aes(label=round(value,2)), 
       size = 2) + scale_size(guide="none")

####################################END#########################################
