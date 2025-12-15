#################################################################################
# This is the 6th of 6 scripts to complete the analyses in Wolfe and Stull 2026.#
# The following script imports posterior samples from a fit Stan model and      #
# prepares the age estimation procedures. To do so, we use Slice Sampling of the#
# posterior to complete age estimation. We compare these results to the MCP     #
# results. In order to get MCP results we complete the steps at the following   #
# link: https://rpubs.com/elainechu/mcp_vignette. The data from the 'data_prep.R#
# script is used in both analyses.                                              #
#################################################################################

# Necessary Packages
library(cmdstanr)
set_cmdstan_path("C:/cmdstan/cmdstan-2.33.1")
library(posterior)
library(tidyverse)
library(magrittr)
library(bridgestan)
set_bridgestan_path("C:/bridgestan")
library(foreach)
library(doParallel)

# Necessary Functions
to_stan_json <- function(data, always_decimal = FALSE) {
  if (!is.list(data)) {
    stop("'data' must be a list.", call. = FALSE)
  }
  
  data_names <- names(data)
  if (length(data) > 0 &&
      (length(data_names) == 0 ||
       length(data_names) != sum(nzchar(data_names)))) {
    stop("All elements in 'data' list must have names.", call. = FALSE)
    
  }
  if (anyDuplicated(data_names) != 0) {
    stop("Duplicate names not allowed in 'data'.", call. = FALSE)
  }
  
  for (var_name in data_names) {
    var <- data[[var_name]]
    if (!(is.numeric(var) || is.factor(var) || is.logical(var) ||
          is.data.frame(var) || is.list(var))) {
      stop("Variable '", var_name, "' is of invalid type.", call. = FALSE)
    }
    if (anyNA(var)) {
      stop("Variable '", var_name, "' has NA values.", call. = FALSE)
    }
    
    if (is.table(var)) {
      var <- unclass(var)
    } else if (is.logical(var)) {
      mode(var) <- "integer"
    } else if (is.data.frame(var)) {
      var <- data.matrix(var)
    } else if (is.list(var)) {
      var <- list_to_array(var, var_name)
    }
    data[[var_name]] <- var
  }
  
  # unboxing variables (N = 10 is stored as N : 10, not N: [10])
  jsonlite::toJSON(
    data,
    auto_unbox = TRUE,
    factor = "integer",
    always_decimal = always_decimal,
    digits = NA,
    pretty = TRUE
  )
}

uni.slice <- function (x0, g, ..., w=1, m=0, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # From HBglm slice_sampler.R
  # Find the log density at the initial point, if not already known.
  if (is.null(gx0)) {gx0 <- g(x0)}
  # Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)
  # Find the initial interval to sample from.
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  if (is.infinite(m))  # no limit on number of steps
  { 
    repeat
    { if (L<=lower) break
      if (g(L)<=logy) break
      L <- L - w
    }
    repeat
    { if (R>=upper) break
      if (g(R)<=logy) break
      R <- R + w
    }
  }
  else if (m>1)  # limit on steps, bigger than one
  { 
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    while (J>0)
    { if (L<=lower) break
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    while (K>0)
    { if (R>=upper) break
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  # Shrink interval to lower and upper bounds.
  if (L<lower) { L <- lower}
  if (R>upper) { R <- upper}
  # Sample from the interval, shrinking it on each rejection.
  repeat
  { 
    x1 <- runif(1,L,R)
    gx1 <- g(x1)
    if (gx1>=logy) break
    if (x1>x0) { R <- x1}
    else { L <- x1}
  }
  return (x1)
}

# Compile Bridgestan Model
mod_bs <- compile_model("stan_models/cop_age_3var.stan",
                        make_args = c('BRIDGESTAN_AD_HESSIAN=true'))

# Import Inital Dataset
dat <- read.csv("data.analysis.csv")

# Prep Data for Slice Sampling
dat2 <- dat %>% select(agey,FDL_L, HPE_EF_L, max_M1_L) %>% 
  na.omit()

# Import Posterior Samples for Slice Sampling and Thin for Efficiency
draws <- read_cmdstan_csv(
  files = c("fitted_models/multivariate_models/combined_SD_final/out_cop-1.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-2.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-3.csv",
            "fitted_models/multivariate_models/combined_SD_final/out_cop-4.csv"),
            variables = c("corr_mat", "FDL_constant", "FDL_exponent", 
                          "FDL_offset", "FDL_noise_intercept", 
                          "FDL_noise_slope", "max_m1_constant", "max_m1_t_pars",
                          "hpe_ef_constant", "hpe_ef_t_pars"),
                          sampler_diagnostics = "", format = "draws_matrix")

draws_thin <- thin_draws(draws$post_warmup_draws,16)

# Wrangle in Prep for Slice Sampling

### FDL draws
fdl_draws <- draws_thin[,c("FDL_constant", "FDL_exponent", "FDL_offset",
                           "FDL_noise_intercept", "FDL_noise_slope")]

### M1 draws
max_m1 <- draws_thin[,2922:2933]

### HPE draws
hpe_ef <- draws_thin[,2934:2940]

corr_list <- list()
m1_t <- list()
hpe_t <- list()

resp_vars <- data.frame(x = names(read.csv("response_vars.csv")))$x

### Prep correlation matrices and ordinal thresholds for slice sampling
for(i in 1:1000){
  
  mats <- matrix(data = draws_thin[i,c(1:2916)], nrow = 54, ncol = 54, dimnames = list(resp_vars,resp_vars))
  mats <- mats[c("FDL","HPE_EF", "Max_M1"),c("FDL","HPE_EF","Max_M1")]
  
  corr_list[[i]] <- mats
  m1_t[[i]] <- c(max_m1[i, 2:12])
  hpe_t[[i]] <- c(hpe_ef[i, 2:7])
  
}

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define function for parallel computation
comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
}

# Complete Slice Sampling
results <- foreach(i = 1:nrow(dat2), .combine = "comb",.multicombine = T, 
                   .packages = c('bridgestan')) %dopar%{
  
  th <- 2
  ths <- matrix(nrow = 1000, ncol = 1)
  ths_unc <- matrix(nrow = 1000, ncol = 1)
  
  for(s in 1:1000){
    
    standati <- list(
      N = 1,
      M = 3,
      LB_C = 7,
      DD_C = 12,
      corr_mat = corr_list[[s]],
      FDL = as.array(dat2$FDL_L[i]),
      HPE = as.array(dat2$HPE_EF_L[i]),
      max_M1 = as.array(dat2$max_M1_L[i]),
      FDL_constant = as.numeric(fdl_draws[s,"FDL_constant"]),
      FDL_exponent = as.numeric(fdl_draws[s,"FDL_exponent"]),
      FDL_offset = as.numeric(fdl_draws[s,"FDL_offset"]),
      FDL_noise_intercept = as.numeric(fdl_draws[s,"FDL_noise_intercept"]),
      FDL_noise_slope = as.numeric(fdl_draws[s,"FDL_noise_slope"]),
      HPE_ef_constant = as.numeric(hpe_ef[s,"hpe_ef_constant"]),
      HPE_ef_t_pars = hpe_t[[s]],
      max_M1_constant = as.numeric(max_m1[s,"max_m1_constant"]),
      max_M1_t_pars = m1_t[[s]]
    )     
    suppressWarnings(modbs <- StanModel$new(lib = mod_bs,
                                            data = to_stan_json(standati),
                                            seed = 1))
    if (s==1) {
      # warmup
      for (r in 1:10) {
        th <- uni.slice(th, modbs$log_density)
      }
    }
    
    th <- try(uni.slice(th, modbs$log_density), silent = T)
    if(inherits(th, "try-error")){
      th <- NA
    }
    th <- th
    ths_unc[s,1] <- th
    ths[s,1] <- modbs$param_constrain(th)
    
    
  }
  return(list(ths_unc, ths))
}

# Comparison to MCP
## The MCP was fit per the same dataset above and using the steps
## outlined at https://rpubs.com/elainechu/mcp_vignette/. 

## To complete the continuously ranked probability score (CRPS) each model must
## be run twice. 
# CRPS_cop <- crps(cop1, cop2, dat2$agey)
# CRPS_mcp <- crps(mcp1, mcp2, dat2$agey)

## Prep Figure 11

### MCP Test Predictions

mcp <- read.csv("mcp_results.csv")
cop <- results[[2]]

df2 <- cbind(dat2$agey,mcp$agey,colMeans(cop, na.rm=T),mcp$xmean)

### Figure 11A
ggplot() + geom_point(aes(dat2$agey, (dat2$agey - colMeans(cop, na.rm=T)), 
  shape = "Copula"), size=0.5) + geom_smooth(aes(dat2$agey, 
  (dat2$agey - colMeans(cop, na.rm=T)), fill = "Copula"), se=F, 
  col = "tomato") + geom_point(aes(mcp$agey, (mcp$agey - mcp$xmean), 
  shape = "MCP"), size=0.5) + geom_smooth(aes(mcp$agey, (mcp$agey - mcp$xmean), 
  fill = "MCP"), se=F, col = "steelblue") + scale_x_continuous(name = 
  "Chronological Age [years]",breaks = seq(0,16,1)) + 
  scale_y_continuous(name = "Residuals", breaks = seq(-4,5,1)) +
  scale_fill_manual(name = "", values = c("tomato", "steelblue")) + 
  scale_shape_manual(name = "", values = c(1,2)) + theme_classic() + 
  guides(shape = guide_legend(override.aes = list(size =2))) +
  theme(legend.position = "bottom")

### Figure 11B
ggplot() + geom_point(aes(dat2$agey, abs((dat2$agey - 
  colMeans(cop, na.rm=T))), shape = "Copula"), size=0.5) +
  geom_smooth(aes(dat2$agey, abs((dat2$agey - colMeans(cop, na.rm=T))), 
  fill = "Copula"), se=F, col = "tomato") + geom_point(aes(mcp$agey, 
  abs((mcp$agey - mcp$xmean)), shape = "MCP"), size=0.5) + 
  geom_smooth(aes(mcp$agey, abs((mcp$agey - mcp$xmean)), fill = "MCP"), se=F, 
  col = "steelblue") + scale_x_continuous(name = "Chronological Age [years]", 
  breaks = seq(0,16,1)) + scale_y_continuous(name = "Absolute Residuals", 
  breaks = seq(0,5,1)) + scale_fill_manual(name = "", values = c("tomato", 
  "steelblue")) + scale_shape_manual(name = "", values = c(1,2)) + 
  theme_classic() + guides(shape = guide_legend(override.aes = list(size =2))) +
  theme(legend.position = "bottom")

### Figure 11C
ggplot(mcp) + 
  geom_errorbar(aes(x=agey, ymin=lower95, ymax=upper95), col = "steelblue") + 
  geom_point(aes(x=agey, y=xmode)) + 
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") + 
  labs(x="Known Age [years]", y="Predicted Age [years]") + 
  theme_bw()

### Figure 11D
for(i in 1:338){
  
  q1[i] <- quantile(cop[,i], probs = c(0.025, 0.975), na.rm=T)[1]
  q2[i] <- quantile(cop[,i], probs = c(0.025, 0.975), na.rm=T)[2]
  
}

ggplot() + 
  geom_errorbar(aes(x=dat2$agey, ymin=q1, ymax=q2),col = "tomato") + 
  geom_point(aes(x=dat2$agey, y=colMeans(cop, na.rm=T))) + 
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") + 
  labs(x="Known Age [years]", y="Predicted Age [years]") + 
  theme_bw()

# Stop CLuster
stopImplicitCluster()

####################################END#########################################
