# Summary
The code and files associated with this repository may be used to recreate the results in *A Bayesian framework to infer ontogenetic relationships and predict associated parameters using human growth and development traits* by Christopher A. Wolfe and Kyra E. Stull. DOI from *Royal Society: Open Science* to follow upon final proof. All supplementary material from the text is stored here.

All data in the present analyses derives from the open-access repository associated with the **Subadult Virtual Anthropology Database**. The data file can be found at [Stull & Corron, 2022](https://zenodo.org/records/5193208). More information about the repository can be found at the [homepage](https://www.unr.edu/anthropology/research-and-facilities/subadult-database).
> Stull, K., & Corron, L. (2021). SVAD_US (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5193208.

> Stull, K., & Corron, L. (2022). Subadult Virtual Anthropology Database (SVAD) Data Collection Protocol: Epiphyseal Fusion, Diaphyseal Dimensions, Dental Development Stages, Vertebral Neural Canal Dimensions. Zenodo. https://doi.org/10.5281/zenodo.7293977.

# Reproducibility
Below we highlight the contents of each folder and describe the purpose of each file.

-  R Code
   -  `data_prep.R`
      - Step 1 of 6:
         - This file imports all necessary growth data, adjusts the ordinal scale to begin at 1, orders all factors (ordinal variables only) to ensure monotonicity, and removes the NA values in the continuous data only. Note, the NA values in the ordinal variables are transformed to the integer `99` as modified in files in `Stan_Models`.
         - All variables were chosen because they are the 'standard' growth and development variables collected as part of SVAD protocols. See Stull and Corron (2022).
   - `model_fit.R`
      - Step 2 of 6:
         - This file prepares the data into a list format for stan compilation and sampling.
         - All missing ordinal values are replaced with the integer `99`. This is for ease in the modeling step to account for missing threshold(s).
         - All continuous data is read in as a vector of present values (e.g., `FDL_complete`) and indices of missing and/or present values to piece together the complete response vector (e.g., `FDL_complete` or `FDL_missing`). This occurs because missing variables are treated as random variables and estimated as part of the modeling process.
         - The final part of this file completes the sampling step.
   - `correlation_plots.R`
      - Step 3 of 6:
         - This file must be run *after* sampling in Step 2. It imports all necessary variables to complete correlation matrices
         - Figure 3
         - Figure 5A, 6A, 7A
   - `eigen_analyses.R`
      - Step 4 of 6:
         - This file must be run *after* Step 3. It imports the correlation matrices and completes eigendecomposition.
         - Figure 4
         - Figure 5B, 6B, 7B
   - `imputation.R`
      - Step 5 of 6:
         - This file must be run *after* sampling in Step 2. It imports in the requisite draws and combines tthis information with 2 additional `.stan` files. This step completes data imputation using posterior draws from the model fit. 
         - `impute_1.stan` is used here for Figure 9A and 9B
         - `impute_2.stan` is used here for Figure 9C
   - `age_estim.R`
      - Step 6 of 6:
         - This file must be run *after* sampling in Step 2. It imports the requisite draws, combines this information with an additional `.stan` file, and runs the MCP model following steps laid out in this [vignette](https://rpubs.com/elainechu/mcp_vignette).
         - `cop_age_3var.stan` and associated code are used for Figure 11.
- Stan Models
   - `MixGaussCop_Growth.stan`
      - Samples the log probability density of a Gaussian copula. Necessary for all downstream analyses.
   - `impute1.stan` and `impute2.stan`
      - Takes the posterior samples from `MixGaussCop_Growth.stan` and uses them to impute missing data. The variables were selected at random. A user could modify these scripts for their own variable(s) and imputation task.
   - `cop_age_3var.stan`
      - Takes the posterior samples from `MixGaussCop_Growth.stan` and uses them to estimate age in a slice sampling routine. See the main text and supplementary material for specifics steps.
- `response_vars.csv` is a helper file to ensure a user does not always need to write out all variable names / abbreviations.      

\
\
*Note: The code here is designed to evaluate posterior samples of human growth and development data. That said, the basic mode fitting steps can be completed for any series of data types.*
