#################################################################################
# This is the 4th of 6 scripts to complete the analyses in Wolfe and Stull 2026.#
# The following script completes the eigendecomposition analyses. To do so, a   #
# user must used the mean poster matrix "cmat_mean" generated in                #
# "correlation_plots.R"                                                         #
#################################################################################

# Note - this script can be modified based on whether the user wants to get the
# loadings of all variables (Figure 4), skeletal growth (Figure 5B), skeletal
# development (Figure 6B), and dental development (Figure 7B).

## Load in the following from "correlation_plots.R":
## cmat_mean = All variables, Figure 4
## cmat_diaphyseal = Skeletal growth, Figure 5B
## cmat_skel_dev = Skeletal development, Figure 6B
## cmat_dental_dev = Dental development, Figure 7B

eig <- eigen(cmat_mean)

eig$values
eig$vectors


scree_data <- data.frame(
  Component = factor(1:length(eig$values)),
  Prop_Explained = eig$values/ sum(eig$values)
)
scree_data$cum_sum <- cumsum(scree_data$Prop_Explained)

ggplot(scree_data, aes(x = Component, y = Prop_Explained)) +
  geom_col() +
  labs(title = "Scree Plot", x = "Dimension", y = "% of Variance Explained") +
  scale_y_continuous(breaks = seq(0,1,0.1)) + theme_minimal()

var_contributions <- eig$vectors^2
var_contributions <- sweep(var_contributions, 2, eig$values, "*")
var_contributions <- sweep(var_contributions, 2, colSums(var_contributions), "/")
var_contributions %<>% as.data.frame %>% mutate(var = resp_vars)

## dim1
var_contributions %>% select(V1, V2, var) %>% ggplot(aes(x = fct_rev(fct_reorder(var, V1)), y = V1)) + geom_col() + geom_hline(yintercept = (1/54), col="tomato", lty = 2) + xlab("Variable") + ylab("Contributions (%)") + 
  theme(axis.text.x = element_text(angle = 90))

## dim2
var_contributions %>% select(V1, V2, var) %>% ggplot(aes(x = fct_rev(fct_reorder(var, V2)), y = V2)) + geom_col() + geom_hline(yintercept = (1/54), col="tomato", lty = 2) + xlab("Variable") + ylab("Contributions (%)") + 
  theme(axis.text.x = element_text(angle = 90))

## dim3
var_contributions %>% select(V1, V2,V3, var) %>% ggplot(aes(x = fct_rev(fct_reorder(var, V3)), y = V3)) + geom_col() + geom_hline(yintercept = (1/54), col="tomato", lty = 2) + xlab("Variable") + ylab("Contributions (%)") + 
  theme(axis.text.x = element_text(angle = 90))

var_contributions %>% select(V1, V2,V3,V4, var) %>% ggplot(aes(x = fct_rev(fct_reorder(var, V4)), y = V4)) + geom_col() + geom_hline(yintercept = (1/54), col="tomato", lty = 2) + xlab("Variable") + ylab("Contributions (%)") + 
  theme(axis.text.x = element_text(angle = 90))


# Prepare eigenvectors for arrows
loadings <- as.data.frame(eig$vectors)
loadings$Variables <- resp_vars

# Scaling the arrows for better visualization
loadings$PC1 <- loadings$V1 * sqrt(eig$values[1])
loadings$PC2 <- loadings$V2 * sqrt(eig$values[2])
loadings$PC3 <- loadings$V3 * sqrt(eig$values[3])

# Step 5: Create the biplot
ggplot() +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC3), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings, aes(x = PC1, y = PC3, label = Variables), 
                  hjust = 1.5, vjust =1.5, color = "red", max.overlaps = 30) +
  labs(title = "",
       x = paste("Dim 1 (", round(eig$values[1] / sum(eig$values) * 100, 1), "%)", sep = ""),
       y = paste("Dim 3 (", round(eig$values[3] / sum(eig$values) * 100, 1), "%)", sep = "")) +
  theme_minimal()
