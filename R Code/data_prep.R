#################################################################################
# This is the 1st of 6 scripts to complete the analyses in Wolfe and Stull 2026.#
# The following script prepares the data for modeling. The data file can be     #
# accessed from a Zenodo repository associated with the SVAD project at         #
# https://zenodo.org/records/5193208. All questions about the SVAD data can     #
# be directed to K.E.S.                                                         #
#################################################################################

# Necessary Packages
library(tidyverse)

# Import Data from SVAD Repository
dat <- read.csv("SVAD_US.csv")

# Wrangle

# In the below steps, we perform some general data wrangling in preparation for #
# the Stan model fitting and downstream analyses. We select all left-sided      #  
# variables, recode all NA values, and organize the columns.                    #

dat_sub <- dat %>% select(agey, location,CC_Oss,TC_Oss ,ends_with("_L")) %>% 
  na_if(-1)

dat_sub <- dat_sub %>% select(1:53,56:65) %>% na_if(-1)

dat_sub <- dat_sub %>% mutate(across(c(3:4,23:63), factor))

dat_sub <- dat_sub %>% mutate_at(vars(c(FH_EF_L:FBDE_EF_L, HPE_EF_L,
                                        HDE_EF_L:CT_EF_L)), 
                                     list(~recode_factor(.,'0' = "1", '1' = "2",
                                                   '12' = "3", '2' = "4",
                                                   '23' = '5', '3' = "6",
                                                   '4' = "7", .ordered = T)))

dat_sub <- dat_sub %>% mutate_at(vars(CC_Oss), 
                                 list(~recode_factor(.,'0' = "1", '1' = "2",
                                                     '2' = "3", '3' = "4",
                                                     '4' = '5', '5' = "6",
                                                     '6' = "7", '7' = "8",
                                                     '8' = "9",.ordered = T)))

dat_sub <- dat_sub %>% mutate_at(vars(TC_Oss), 
                                 list(~recode_factor(.,'0' = "1", '1' = "2",
                                                     '2' = "3", '3' = "4",
                                                     '4' = '5', '5' = "6",
                                                     '6' = "7", '7' = "8",
                                                     .ordered = T)))

dat_sub <- dat_sub %>% mutate_at(vars(ISPR_EF_L:ILIS_EF_L), 
                                 list(~recode_factor(.,'0' = "1", '1' = "2",
                                                     '2' = "3",.ordered = T)))

dat_final <- dat_sub %>% select(1:2, HH_Oss_L:HLT_Oss_L, HC_Oss_L:HLE_Oss_L,
                                PC_Oss_L, max_M1_L:man_I2_L, FH_EF_L:FBDE_EF_L, 
                                HPE_EF_L,HDE_EF_L:CT_EF_L, CC_Oss, TC_Oss, 
                                ISPR_EF_L:ILIS_EF_L, 5:22)

dat_final <- dat_final %>% mutate_at(vars(max_M1_L:man_I2_L), 
                                     list(~factor(., ordered = T)))

# write.csv(dat_final, "data.analysis.csv", row.names = F)\

####################################END#########################################
