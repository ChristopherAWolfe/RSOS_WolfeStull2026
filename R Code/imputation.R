#################################################################################
# This is the 5th of 6 scripts to complete the analyses in Wolfe and Stull 2026.#
# The following script imports posterior samples from a fit Stan model and      #
# prepares the imputation analyses. To do so, we use CmdStanR's GQ Function     #
# to complete a series of 2 imputation tasks. The Stan models can be found in   #
# the 'Stan_models' folder.                                                     #
#################################################################################

# Necessary Packages
library(cmdstanr)
set_cmdstan_path("C:/cmdstan/cmdstan-2.33.1")
library(tidyverse)
library(posterior)
library(magrittr)
library(bayesplot)

# Import Data File created from 'data_prep.R' file. Note all study data derives
# from a repository independent of the current project.
dat <- read.csv("data.analysis.csv")

# Load Necessary Stan Models
mod1 <- cmdstan_model(stan_file = "impute1.stan")
mod2 <- cmdstan_model(stan_file = "impute2.stan")

# Import Posterior Samples
draws1 <- read_cmdstan_csv(
  files = c("fitted_models/multivariate_models/combined_SD_final/out_cop-1.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-2.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-3.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-4.csv"),
            variables = c("corr_mat","FMSB_constant", "FMSB_exponent", 
                          "FMSB_offset", "FMSB_noise_intercept", 
                          "FMSB_noise_slope", "FDB_constant", "FDB_exponent", 
                          "FDB_offset", "FDB_noise_intercept", "FDB_noise_slope",
                          "FDL_constant", "FDL_exponent", "FDL_offset", 
                          "FDL_noise_intercept", "FDL_noise_slope"),
                          sampler_diagnostics = "")

draws2 <- read_cmdstan_csv(
  files = c("fitted_models/multivariate_models/combined_SD_final/out_cop-1.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-2.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-3.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-4.csv"),
            variables = c("corr_mat", "HDL_constant", "HDL_exponent", 
                          "HDL_offset", "HDL_noise_intercept", "HDL_noise_slope",
                          "man_pm1_constant", "man_pm1_t_pars", "ct_ef_constant", 
                          "ct_ef_t_pars"), sampler_diagnostics = "")

# Figure 8A, 8B, and 8C

## Stan Model Data Prep
dat1 <- dat %>% select(agey,FMSB_L, FDB_L,FDL_L) %>% na.omit()

standat1 <- list(N = nrow(dat1),
                M = 54, 
                K = 3, 
                y1 = dat1$FMSB_L,
                y2 = dat1$FDB_L,
                x = dat1$agey)

fit_gq1 <- mod1$generate_quantities(fitted_params = draws1$post_warmup_draws,
                                  data = standat1, parallel_chains = 4)

## Wrangle Samples
fdl <- bind_cols(dat1$FDL_L, summarise_draws(fit_gq1$draws("ypred"))$mean, 
                 summarise_draws(fit_gq1$draws("expectation"))$mean)
colnames(fdl) <- c("FDL_true", "FDL_impute", "FDL_impute_E")
fdl$age <- dat1$agey

## Figure 8A
fdl |> ggplot(aes(x = FDL_true, y = FDL_impute_E)) + 
  geom_point(size = 0.7) + geom_abline(col = "tomato",slope = 1, intercept = 0) +
  xlab("True FDL [mm]") + ylab("FDL Imputed [mm]") +
  theme_classic() + theme(legend.position = "none")

## Figure 8B
fdl |> ggplot(aes(x = age)) + 
  geom_point(aes(y = FDL_true, fill = "True"),size = 0.7) + geom_point(aes(y = FDL_impute, fill = "Imputed"),col = "tomato", size = 0.7) +
  xlab("Age [years]") + ylab("FDL [mm]") +
  theme_classic() + theme(legend.position = "bottom") + labs(fill = "")

# Figure 8C

## Stan Model Data Prep
dat2 <- dat %>% select(agey,HDL_L, man_PM1_L, CT_EF_L) %>% na.omit()

standat2 <- list(N = nrow(dat2),
                M = 54, 
                K = 3,
                y1 = dat2$HDL_L,
                y2 = dat2$man_PM1_L,
                x = dat2$agey,
                y2_c = 12,
                y3_c = 7)

fit_gq2 <- mod2$generate_quantities(fitted_params = draws2$post_warmup_draws,
                                  data = standat2, parallel_chains = 4)

## Wrangle Samples
ct <- bind_cols(dat2$CT_EF_L, summarise_draws(fit_gq2$draws("ypred"))$median, 
                ummarise_draws(fit_gq2$draws("expectation"))$median)
colnames(ct) <- c("CT_True", "CT_Impute", "CT_Impute_E")
ct$age <- dat2$agey

ct |> ggplot(aes(x = age)) + 
  geom_point(aes(y = CT_True, fill = "True"),size = 0.7) + 
  geom_point(aes(y = CT_Impute_E, fill = "Imputed"), col = "tomato", 
  size = 0.7) + xlab("Age [years]") + ylab("CT Score") + 
  scale_x_continuous(breaks = seq(0,16,1)) +
  scale_y_continuous(breaks = seq(1,7,1)) + theme_classic() + 
  theme(legend.position = "bottom") + labs(fill = "")

####################################END#########################################
