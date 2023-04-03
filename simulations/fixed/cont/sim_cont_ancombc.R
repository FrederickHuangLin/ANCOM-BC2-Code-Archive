library(tidyverse)
library(mia)
library(ggpubr)
library(doRNG)
library(doParallel)
library(ANCOMBC)
library(MicrobiomeStat)
library(LOCOM)

logsumexp = function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax = function (x) {
  exp(x - logsumexp(x))
}

# Simulation settings
set.seed(12345)
data(QMP, package = "ANCOMBC")
n = c(10, 20, 30, 50, 100)
d = ncol(QMP)
diff_prop = c(0.1, 0.2, 0.5, 0.9)
iter_num = 100
seed = seq_len(iter_num)
df_sim_params = data.frame(expand.grid(n, diff_prop, seed)) %>%
  dplyr::rename(n = Var1, diff_prop = Var2, seed = Var3) %>%
  arrange(n, diff_prop, seed)
list_sim_params = apply(df_sim_params, 1, paste0, collapse = "_")

lfc_value = c(-2, -1, 1, 2)
lfc_cont_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_cont_list[[i]] = sample(c(0, lfc_value), size = d, replace = TRUE,
                              prob = c(1 - diff_prop[i], 
                                       rep(diff_prop[i]/length(lfc_value), length(lfc_value))))
}
names(lfc_cont_list) = diff_prop

lfc_bin_list = vector("list", length = length(diff_prop))
for (i in seq_along(diff_prop)) {
  lfc_bin_list[[i]] = sample(c(0, 1), size = d, replace = TRUE,
                             prob = c(1 - diff_prop[i], diff_prop[i]))
}
names(lfc_bin_list) = diff_prop

cl = makeCluster(4)
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
    abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05, 
                        n = n, lib_mean = 1e8, disp = 0.5)
    log_abn_data = log(abn_data + 1e-5)
    rownames(log_abn_data) = paste0("T", seq_len(d))
    colnames(log_abn_data) = paste0("S", seq_len(n))
    
    # Generate the sample and feature meta data
    # Sampling fractions are set to differ by the variable of interest
    smd = data.frame(sample = paste0("S", seq_len(n)),
                     cont_cov = rnorm(n),
                     bin_cov = as.factor(rep(seq_len(2), each = n/2))) %>%
      dplyr::mutate(samp_frac = log(softmax(cont_cov)/10))
    
    d = nrow(abn_data) 
    lfc_cont = lfc_cont_list[[as.character(diff_prop)]]
    lfc_bin = lfc_bin_list[[as.character(diff_prop)]]
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
    
    # Create the tse object
    assays = S4Vectors::SimpleList(counts = otu_data)
    smd = S4Vectors::DataFrame(smd)
    tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)
    
    # Run ANCOM-BC
    set.seed(123)
    output = ancombc(data = tse, assay_name = "counts", 
                     tax_level = NULL, phyloseq = NULL, 
                     formula = "cont_cov + bin_cov", 
                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                     group = "bin_cov", struc_zero = FALSE, neg_lb = FALSE, tol = 1e-5, 
                     max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE,
                     n_cl = 1, verbose = FALSE)
    
    res_prim = output$res
    res_merge = res_prim$lfc %>%
      dplyr::transmute(taxon, lfc_est = cont_cov) %>%
      dplyr::left_join(res_prim$diff_abn %>%
                         dplyr::transmute(taxon, diff = cont_cov),
                       by = "taxon") %>%
      dplyr::left_join(fmd %>%
                         dplyr::transmute(taxon, lfc_true = lfc_cont),
                       by = "taxon") %>%
      dplyr::transmute(taxon, 
                       lfc_est = lfc_est * diff,
                       lfc_est = case_when(lfc_est > 0 ~ 1,
                                           lfc_est < 0 ~ -1,
                                           TRUE ~ 0),
                       lfc_true = case_when(lfc_true > 0 ~ 1,
                                            lfc_true < 0 ~ -1,
                                            TRUE ~ 0))
    lfc_est = res_merge$lfc_est
    lfc_true = res_merge$lfc_true
    tp = sum(lfc_true != 0 & lfc_est != 0)
    fp = sum(lfc_true == 0 & lfc_est != 0)
    fn = sum(lfc_true != 0 & lfc_est == 0)
    power = tp/(tp + fn)
    fdr = fp/(tp + fp)
    
    c(power, fdr)
  }

stopCluster(cl)

write_csv(data.frame(res_sim), "sim_cont_ancombc.csv")