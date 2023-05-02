library(tidyverse)
library(mia)
library(ggpubr)
library(doRNG)
library(doParallel)
library(ANCOMBC)

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

lfc_value = seq(from = 0.5, to = 2, length.out = 4)
lfc_cat2_vs_1_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cat2_vs_1_list[[i]] = sample(c(0, lfc_value), size = d, replace = TRUE,
                                   prob = c(1 - diff_prop[i],
                                            rep(diff_prop[i]/length(lfc_value), length(lfc_value))))
}
names(lfc_cat2_vs_1_list) = diff_prop

lfc_cont_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cont_list[[i]] = sample(c(0, 1), size = d, replace = TRUE,
                              prob = c(1 - diff_prop[i], diff_prop[i]))
}
names(lfc_cont_list) = diff_prop

# detect all cpus per task
ncpus = Sys.getenv("SLURM_CPUS_PER_TASK")
# make cluster with ncpus
cl = makeCluster(strtoi(ncpus))
registerDoParallel(cl)

res_sim = foreach(i = list_sim_params, .combine = rbind, .verbose = TRUE, 
                  .packages = c("ANCOMBC", "tidyverse")) %dorng% 
  {
    Sys.sleep(4)
    
    params = strsplit(i, "_")[[1]]
    n = as.numeric(params[1])
    diff_prop = as.numeric(params[2])
    seed = as.numeric(params[3])
    
    # Generate the true abundances
    set.seed(seed)
    abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05, 
                        n = n, lib_mean = 1e8, disp = 0.5)
    log_abn_data = log(abn_data + 1e-5)
    rownames(log_abn_data) = paste0("T", seq_len(d))
    colnames(log_abn_data) = paste0("S", seq_len(n))
    
    # Generate the sample and feature meta data
    # Sampling fractions are set to differ by batches
    smd = data.frame(sample = paste0("S", seq_len(n)),
                     samp_frac = log(c(runif(n/3, min = 1e-4, max = 1e-3),
                                       runif(n/3, min = 1e-3, max = 1e-2),
                                       runif(n/3, min = 1e-2, max = 1e-1))),
                     cont_cov = rnorm(n),
                     cat_cov = as.factor(rep(seq_len(3), each = n/3)))
    
    d = nrow(abn_data)   
    lfc_cat2_vs_1 = lfc_cat2_vs_1_list[[as.character(diff_prop)]]
    lfc_cont = lfc_cont_list[[as.character(diff_prop)]]
    fmd = data.frame(taxon = paste0("T", seq_len(d)),
                     seq_eff = log(runif(d, min = 0.1, max = 1)),
                     lfc_cont = lfc_cont,
                     lfc_cat2_vs_1 = lfc_cat2_vs_1) %>%
      mutate(lfc_cat3_vs_1 = ifelse(lfc_cat2_vs_1 == 0, 0, lfc_cat2_vs_1 + 1))
    
    # Add effect sizes of covariates to the true abundances
    smd_dmy = model.matrix(~ 0 + cont_cov + cat_cov, data = smd)
    log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
    log_abn_data = log_abn_data + outer(fmd$lfc_cat2_vs_1, smd_dmy[, "cat_cov2"])
    log_abn_data = log_abn_data + outer(fmd$lfc_cat3_vs_1, smd_dmy[, "cat_cov3"])
    
    # Add sample- and taxon-specific biases
    log_otu_data = t(t(log_abn_data) + smd$samp_frac)
    log_otu_data = log_otu_data + fmd$seq_eff
    otu_data = round(exp(log_otu_data))
    
    # Create the tse object
    assays = S4Vectors::SimpleList(counts = otu_data)
    smd = S4Vectors::DataFrame(smd)
    tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
    
    # Run ANCOM-BC2
    set.seed(123)
    output = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                      fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
                      p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                      group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
                      alpha = 0.05, n_cl = 1, verbose = FALSE,
                      global = TRUE, pairwise = FALSE, 
                      dunnet = FALSE, trend = TRUE,
                      iter_control = list(tol = 1e-5, max_iter = 20, 
                                          verbose = FALSE),
                      em_control = list(tol = 1e-5, max_iter = 100),
                      lme_control = NULL, mdfdr_control = NULL, 
                      trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                  nrow = 2, byrow = TRUE)),
                                           node = list(2),
                                           solver = "ECOS",
                                           B = 100))
    
    res_trend = output$res_trend
    tab_sens = output$pseudo_sens_tab
    sens_cut = 1
    res_merge = res_trend %>%
      dplyr::transmute(taxon, 
                       diff_est = diff_abn) %>%
      dplyr::left_join(tab_sens %>%
                         dplyr::transmute(taxon, 
                                          sens = global),
                       by = "taxon") %>%
      dplyr::left_join(fmd %>%
                         dplyr::transmute(taxon, 
                                          diff_true = ifelse(lfc_cat2_vs_1 == 0, 0, 1)),
                       by = "taxon") %>%
      dplyr::transmute(taxon,
                       diff_est = diff_est * (sens < sens_cut),
                       diff_true)
    diff_est = res_merge$diff_est
    diff_true = res_merge$diff_true
    
    tp = sum(diff_true != 0 & diff_est != 0)
    fp = sum(diff_true == 0 & diff_est != 0)
    fn = sum(diff_true != 0 & diff_est == 0)
    power = tp/(tp + fn)
    fdr = fp/(tp + fp)
    
    c(power, fdr)
  }

write_csv(data.frame(res_sim), "sim_trend_ancombc2.csv")

stopCluster(cl)