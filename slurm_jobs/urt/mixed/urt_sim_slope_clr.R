library(tidyverse)
library(mia)
library(ggpubr)
library(doRNG)
library(doParallel)
library(lmerTest)
library(compositions)

# Simulation settings
data(throat.otu.table, package = "LOCOM")
# Use the URT data as the template to obtain mean vector and variance-covariance
# matrix. Discard OTUs that have less than 5% of prevalence across samples
prevalence = apply(t(throat.otu.table), 1, function(x)
  sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
tax_keep = which(prevalence >= 0.05)

set.seed(12345)
n = c(30, 60, 90, 150, 300)
d = length(tax_keep)
diff_prop = c(0.05, 0.1, 0.2, 0.5, 0.9)
iter_num = 100
seed = seq_len(iter_num)
df_sim_params = data.frame(expand.grid(n, diff_prop, seed)) %>%
  dplyr::rename(n = Var1, diff_prop = Var2, seed = Var3) %>%
  arrange(n, diff_prop, seed)
list_sim_params = apply(df_sim_params, 1, paste0, collapse = "_")

# Log-fold-changes for the categorical exposure of DA taxa
lfc_value = c(-2, -1, 1, 2)
lfc_cat2_vs_1_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cat2_vs_1_list[[i]] = sample(c(0, lfc_value), size = d, replace = TRUE,
                                   prob = c(1 - diff_prop[i],
                                            rep(diff_prop[i]/length(lfc_value), length(lfc_value))))
}
names(lfc_cat2_vs_1_list) = diff_prop

lfc_cat3_vs_1_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cat3_vs_1_list[[i]] = sample(c(0, lfc_value), size = d, replace = TRUE,
                                   prob = c(1 - diff_prop[i],
                                            rep(diff_prop[i]/length(lfc_value), length(lfc_value))))
}
names(lfc_cat3_vs_1_list) = diff_prop

# Log-fold-changes for the continuous confounder
lfc_cont_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cont_list[[i]] = sample(c(0, 1), size = d, replace = TRUE,
                              prob = c(1 - diff_prop[i], diff_prop[i]))
}
names(lfc_cont_list) = diff_prop

cl = makeCluster(12)
registerDoParallel(cl)

res_sim = foreach(i = list_sim_params, .combine = rbind, .verbose = TRUE, 
                  .packages = c("ANCOMBC", "tidyverse", "compositions", "lmerTest")) %dorng% 
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
    
    # Generate the sample and feature meta data
    # Sampling fractions are set to differ by batches
    smd = data.frame(sample = paste0("S", seq_len(n)),
                     samp_frac = log(c(runif(n/3, min = 1e-4, max = 1e-3),
                                       runif(n/3, min = 1e-3, max = 1e-2),
                                       runif(n/3, min = 1e-2, max = 1e-1))),
                     cont_cov = rnorm(n),
                     cat_cov = as.factor(rep(seq_len(3), each = n/3)),
                     batch = rep(c("A", "B", "C"), n/3))
    
    d = nrow(abn_data)   
    lfc_cat2_vs_1 = lfc_cat2_vs_1_list[[as.character(diff_prop)]]
    lfc_cat3_vs_1 = lfc_cat3_vs_1_list[[as.character(diff_prop)]]
    lfc_cont = lfc_cont_list[[as.character(diff_prop)]]
    fmd = data.frame(taxon = paste0("T", seq_len(d)),
                     seq_eff = log(runif(d, min = 0.1, max = 1)),
                     lfc_cont = lfc_cont,
                     lfc_cat2_vs_1 = lfc_cat2_vs_1,
                     lfc_cat3_vs_1 = lfc_cat3_vs_1)
    
    # Add fixed effects
    smd_dmy = model.matrix(~ 0 + cont_cov + cat_cov, data = smd)
    log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
    log_abn_data = log_abn_data + outer(fmd$lfc_cat2_vs_1, smd_dmy[, "cat_cov2"])
    log_abn_data = log_abn_data + outer(fmd$lfc_cat3_vs_1, smd_dmy[, "cat_cov3"])
    
    # Add random effects
    covmat =  matrix(c(1, 0.5, 0.5, 2), 2, 2)
    rand_eff = MASS::mvrnorm(n = 3, mu = c(0, 0), Sigma = covmat)
    rand_int_mat = model.matrix(~ 0 + batch, data = smd)
    rand_slope_mat = model.matrix(~ 0 + batch, data = smd) * smd$cont_cov
    log_abn_data = t(apply(log_abn_data , 1, 
                           function(x) 
                             x + rand_int_mat %*% rand_eff[, 1] + rand_slope_mat %*% rand_eff[, 2]))
    rownames(log_abn_data) = paste0("T", seq_len(d))
    colnames(log_abn_data) = paste0("S", seq_len(n))
    
    # Add sample- and taxon-specific biases
    log_otu_data = t(t(log_abn_data) + smd$samp_frac)
    log_otu_data = log_otu_data + fmd$seq_eff
    otu_data = round(exp(log_otu_data))
    
    # Remove samples with low library sizes
    idx = which(colSums(otu_data) > 1000)
    otu_data = otu_data[, idx]
    smd = smd[idx, ]
    
    # Run linear mixed effects model on CLR-transformed data
    clr_data = t(clr(t(otu_data)))
    clr_data = as.matrix(clr_data)
    output = apply(clr_data, 1, function(y) {
      df = data.frame(y, smd)
      fit = lmer(y ~ cont_cov + cat_cov + (cont_cov | batch), data = df)
      summ = summary(fit)
      lfc = summ$coefficients[, "Estimate"]
      names(lfc) = paste0("lfc_", names(lfc))
      p_val = summ$coefficients[, "Pr(>|t|)"]
      names(p_val) = paste0("p_", names(p_val))
      return(c(lfc, p_val))
    })
    res = data.frame(t(output), check.names = FALSE)
    res = res %>%
      rownames_to_column("taxon") %>%
      mutate(across(starts_with("p_"), 
                    ~ p.adjust(., method = "holm"), 
                    .names = "q_{.col}")) %>%
      rename_with(~str_replace(., "q_p_", "q_"), starts_with("q_p_")) %>%
      mutate(diff_cat_cov2 = ifelse(q_cat_cov2 <= 0.05, TRUE, FALSE),
             diff_cat_cov3 = ifelse(q_cat_cov3 <= 0.05, TRUE, FALSE))
    
    
    res_merge = res %>%
      dplyr::transmute(taxon, 
                       lfc_est1 = lfc_cat_cov2 * diff_cat_cov2,
                       lfc_est2 = lfc_cat_cov3 * diff_cat_cov3) %>%
      dplyr::left_join(fmd %>%
                         dplyr::transmute(taxon, 
                                          lfc_true1 = lfc_cat2_vs_1,
                                          lfc_true2 = lfc_cat3_vs_1),
                       by = "taxon") %>%
      dplyr::transmute(taxon, 
                       lfc_est1 = case_when(lfc_est1 > 0 ~ 1,
                                            lfc_est1 < 0 ~ -1,
                                            TRUE ~ 0),
                       lfc_est2 = case_when(lfc_est2 > 0 ~ 1,
                                            lfc_est2 < 0 ~ -1,
                                            TRUE ~ 0),
                       lfc_true1 = case_when(lfc_true1 > 0 ~ 1,
                                             lfc_true1 < 0 ~ -1,
                                             TRUE ~ 0),
                       lfc_true2 = case_when(lfc_true2 > 0 ~ 1,
                                             lfc_true2 < 0 ~ -1,
                                             TRUE ~ 0))
    lfc_est1 = res_merge$lfc_est1
    lfc_true1 = res_merge$lfc_true1
    lfc_est2 = res_merge$lfc_est2
    lfc_true2 = res_merge$lfc_true2
    tp1 = sum(lfc_true1 == 1 & lfc_est1 == 1) +
      sum(lfc_true1 == -1 & lfc_est1 == -1)
    fp1 = sum(lfc_true1 == 0 & lfc_est1 != 0) +
      sum(lfc_true1 == 1 & lfc_est1 == -1) +
      sum(lfc_true1 == -1 & lfc_est1 == 1)
    fn1 = sum(lfc_true1 != 0 & lfc_est1 == 0)
    
    tp2 = sum(lfc_true2 == 1 & lfc_est2 == 1) +
      sum(lfc_true2 == -1 & lfc_est2 == -1)
    fp2 = sum(lfc_true2 == 0 & lfc_est2 != 0) +
      sum(lfc_true2 == 1 & lfc_est2 == -1) +
      sum(lfc_true2 == -1 & lfc_est2 == 1)
    fn2 = sum(lfc_true2 != 0 & lfc_est2 == 0)
    tp = tp1 + tp2
    fp = fp1 + fp2
    fn = fn1 + fn2
    power = tp/(tp + fn)
    fdr = fp/(tp + fp)
    
    c(power, fdr)
  }

stopCluster(cl)

write_csv(data.frame(res_sim), "urt_sim_slope_clr.csv")