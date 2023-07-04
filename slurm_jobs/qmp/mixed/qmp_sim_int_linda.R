library(tidyverse)
library(mia)
library(ggpubr)
library(doRNG)
library(doParallel)
library(ANCOMBC)
library(MicrobiomeStat)

# Simulation settings
set.seed(12345)
data(QMP, package = "ANCOMBC")
n = c(30, 60, 90, 150, 300)
d = ncol(QMP)
diff_prop = c(0.1, 0.2, 0.5, 0.9)
iter_num = 100
seed = seq_len(iter_num)
df_sim_params = data.frame(expand.grid(n, diff_prop, seed)) %>%
  dplyr::rename(n = Var1, diff_prop = Var2, seed = Var3) %>%
  arrange(n, diff_prop, seed)
list_sim_params = apply(df_sim_params, 1, paste0, collapse = "_")

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

lfc_cont_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cont_list[[i]] = sample(c(0, 1), size = d, replace = TRUE,
                              prob = c(1 - diff_prop[i], diff_prop[i]))
}
names(lfc_cont_list) = diff_prop

cl = makeCluster(4)
registerDoParallel(cl)

res_sim = foreach(i = list_sim_params, .combine = rbind, .verbose = TRUE, 
                  .packages = c("ANCOMBC", "tidyverse", "MicrobiomeStat")) %dorng% 
  {
    params = strsplit(i, "_")[[1]]
    n = as.numeric(params[1])
    diff_prop = as.numeric(params[2])
    seed = as.numeric(params[3])
    
    # Generate the true abundances
    set.seed(seed)
    abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05, 
                        n = n, lib_mean = 1e8, disp = 0.5)
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
    rand_eff = rnorm(3, mean = 0, sd = 1)
    rand_int_mat = model.matrix(~ 0 + batch, data = smd)
    log_abn_data = t(apply(log_abn_data , 1, function(x) x + rand_int_mat %*% rand_eff))
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
    
    # Run LinDA
    output = linda(feature.dat = otu_data, meta.dat = smd,
                   formula = "~ cont_cov + cat_cov + (1 | batch)",
                   alpha = 0.05, 
                   prev.filter = 0.10, 
                   mean.abund.filter = 0,
                   adaptive = TRUE,
                   max.abund.filter = 0,
                   p.adj.method = "holm",
                   n.cores = 1, 
                   verbose = FALSE)
    
    res = output$output
    res_merge = res$cat_cov2 %>%
      rownames_to_column("taxon") %>%
      dplyr::transmute(taxon, 
                       lfc_est1 = log2FoldChange * reject) %>%
      dplyr::left_join(fmd %>%
                         dplyr::transmute(taxon, 
                                          lfc_true1 = lfc_cat2_vs_1),
                       by = "taxon") %>%
      dplyr::transmute(taxon, 
                       lfc_est1 = case_when(lfc_est1 > 0 ~ 1,
                                            lfc_est1 < 0 ~ -1,
                                            TRUE ~ 0),
                       lfc_true1 = case_when(lfc_true1 > 0 ~ 1,
                                             lfc_true1 < 0 ~ -1,
                                             TRUE ~ 0)) %>%
      dplyr::left_join(
        res$cat_cov3 %>%
          rownames_to_column("taxon") %>%
          dplyr::transmute(taxon, 
                           lfc_est2 = log2FoldChange * reject) %>%
          dplyr::left_join(fmd %>%
                             dplyr::transmute(taxon, 
                                              lfc_true2 = lfc_cat3_vs_1),
                           by = "taxon") %>%
          dplyr::transmute(taxon, 
                           lfc_est2 = case_when(lfc_est2 > 0 ~ 1,
                                                lfc_est2 < 0 ~ -1,
                                                TRUE ~ 0),
                           lfc_true2 = case_when(lfc_true2 > 0 ~ 1,
                                                 lfc_true2 < 0 ~ -1,
                                                 TRUE ~ 0)),
        by = "taxon")
    
    lfc_est1 = res_merge$lfc_est1
    lfc_true1 = res_merge$lfc_true1
    lfc_est2 = res_merge$lfc_est2
    lfc_true2 = res_merge$lfc_true2
    
    tp1 = sum(lfc_true1 != 0 & lfc_est1 != 0)
    fp1 = sum(lfc_true1 == 0 & lfc_est1 != 0)
    fn1 = sum(lfc_true1 != 0 & lfc_est1 == 0)
    tp2 = sum(lfc_true2 != 0 & lfc_est2 != 0)
    fp2 = sum(lfc_true2 == 0 & lfc_est2 != 0)
    fn2 = sum(lfc_true2 != 0 & lfc_est2 == 0)
    tp = tp1 + tp2
    fp = fp1 + fp2
    fn = fn1 + fn2
    power = tp/(tp + fn)
    fdr = fp/(tp + fp)
    
    c(power, fdr)
  }

stopCluster(cl)

write_csv(data.frame(res_sim), "qmp_sim_int_linda.csv")