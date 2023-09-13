library(tidyverse)
library(ggpubr)
library(doRNG)
library(doParallel)
library(LOCOM)
library(ANCOMBC)

# Simulation settings
data(throat.otu.table, package = "LOCOM")
# Use the URT data as the template to obtain mean vector and variance-covariance
# matrix. Discard OTUs that have less than 5% of prevalence across samples
prevalence = apply(t(throat.otu.table), 1, function(x)
  sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
tax_keep = which(prevalence >= 0.05)

set.seed(12345)
n = c(20, 40, 60, 100, 200)
d = length(tax_keep)
diff_prop = c(0.05, 0.1, 0.2, 0.5, 0.9)
iter_num = 100
seed = seq_len(iter_num)
df_sim_params = data.frame(expand.grid(n, diff_prop, seed)) %>%
  dplyr::rename(n = Var1, diff_prop = Var2, seed = Var3) %>%
  arrange(n, diff_prop, seed)
list_sim_params = apply(df_sim_params, 1, paste0, collapse = "_")

# Log-fold-changes for the binary exposure of DA taxa
lfc_value = c(-2, -1, 1, 2)
lfc_bin_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_bin_list[[i]] = sample(c(0, lfc_value), size = d, replace = TRUE,
                             prob = c(1 - diff_prop[i], 
                                      rep(diff_prop[i]/length(lfc_value), length(lfc_value))))
}
names(lfc_bin_list) = diff_prop

# Log-fold-changes for the continuous confounder
lfc_cont_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cont_list[[i]] = sample(c(0, 1), size = d, replace = TRUE,
                              prob = c(1 - diff_prop[i], diff_prop[i]))
}
names(lfc_cont_list) = diff_prop

cl = makeCluster(8)
registerDoParallel(cl)

res_sim = foreach(i = list_sim_params, .combine = rbind, .verbose = TRUE, 
                  .packages = c("ANCOMBC", "tidyverse")) %dorng% 
  {
    params = strsplit(i, "_")[[1]]
    n = as.numeric(params[1])
    diff_prop = as.numeric(params[2])
    seed = as.numeric(params[3])
    
    # Generate the true abundances
    set.seed(seed)
    abn_data = sim_plnm(abn_table = throat.otu.table, taxa_are_rows = FALSE, 
                        prv_cut = 0.05, n = n, lib_mean = 1e8, disp = 0.5)
    log_abn_data = log(abn_data + 1e-5)
    rownames(log_abn_data) = paste0("T", seq_len(d))
    colnames(log_abn_data) = paste0("S", seq_len(n))
    
    # Generate the sample and feature meta data
    # Sampling fractions are set to differ by the variable of interest
    smd = data.frame(sample = paste0("S", seq_len(n)),
                     samp_frac = log(c(runif(n/2, min = 1e-4, max = 1e-3),
                                       runif(n/2, min = 1e-3, max = 1e-2))),
                     cont_cov = rnorm(n),
                     bin_cov = as.factor(rep(seq_len(2), each = n/2)))
    
    d = nrow(abn_data) 
    lfc_bin = lfc_bin_list[[as.character(diff_prop)]]
    lfc_cont = lfc_cont_list[[as.character(diff_prop)]]
    fmd = data.frame(taxon = paste0("T", seq_len(d)),
                     seq_eff = log(runif(d, min = 0.1, max = 1)),
                     lfc_cont = lfc_cont,
                     lfc_bin = lfc_bin)
    
    # Add effect sizes of covariates to the true abundances
    smd_dmy = model.matrix(~ 0 + cont_cov + bin_cov, data = smd)
    log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
    log_abn_data = log_abn_data + outer(fmd$lfc_bin, smd_dmy[, "bin_cov2"])
    
    # Add sample- and taxon-specific biases
    log_otu_data = t(t(log_abn_data) + smd$samp_frac)
    log_otu_data = log_otu_data + fmd$seq_eff
    otu_data = round(exp(log_otu_data))
    
    # Check zero proportions
    overall_zero = sum(otu_data == 0)/(n * d)
    group1_zero = sum(otu_data[, seq_len(n/2)] == 0)/(n * d / 2)
    group2_zero = sum(otu_data[, seq(n/2, n)] == 0)/(n * d / 2)
    
    
    c(overall_zero, group1_zero, group2_zero)
  }

stopCluster(cl)

res_zero = data.frame(res_sim)
colnames(res_zero) = c("Overall", "Group 1", "Group 2")
simpattern = distinct(df_sim_params, n, diff_prop) %>%
  unite("setting", n:diff_prop, sep = ", ")
res_zero = res_zero %>%
  mutate(setting = rep(simpattern$setting, each = iter_num)) %>%
  separate(setting, c("n", "diff_prop"), ", ") %>%
  pivot_longer(cols = Overall:`Group 2`, names_to = "variable", values_to = "value")

fig = res_zero %>%
  ggboxplot(x = "n", y = "value", add = "jitter", color = "n", palette = "jco",
            xlab = "Sample Size", ylab = "Proportion of Zeros", 
            facet.by = "variable", nrow = 1, add.params = list(size = 0.1)) +
  guides(color = guide_legend(title = NULL))
fig

# This data set from Lahti et al. Nat. Comm. 5:4344, 2014 comes with 130 genus-like 
# taxonomic groups across 1006 western adults with no reported health complications
data(atlas1006, package = "microbiome")
feature_table = abundances(atlas1006)
sum(feature_table == 0)/(nrow(feature_table) * ncol(feature_table)) # 0.209

# This data set from O’Keefe et al. Nat. Comm. 6:6342, 2015 comes with 130 genus-like 
# taxonomic groups across 222 adults in a two-week diet swap study 
# between western (USA) and traditional (rural Africa) diets
data(dietswap, package = "microbiome")
feature_table = abundances(dietswap)
sum(feature_table == 0)/(nrow(feature_table) * ncol(feature_table)) # 0.206













