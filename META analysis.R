library(dplyr)
library(data.table)
library(readxl)
library(tibble)
library(tidyr)
library(metafor)

saliva_202 <- read.csv("saliva_202_without caries.csv", row.names = 1, check.names = FALSE)
group_saliva <- read_excel("group_saliva_cariesout.xlsx")

# sub_123 <- read.csv("sub_123.csv", row.names = 1, check.names = FALSE)
# group_sub <- read_excel("D:/Hp data/R_PROJECT OUTCOME/virus-site/input/group_sub.xlsx")

# take saliva for example below

profile <- saliva_202
min_value <- min(profile[profile > 0], na.rm = TRUE)
trans_profile <- log10(profile + min_value)

trans_df <- as.data.frame(trans_profile) %>%
  tibble::rownames_to_column("sample_ID")

data <- trans_df %>%
  inner_join(group_saliva %>% dplyr::select(sample_ID, condition, project), 
             by = "sample_ID") %>%
  setDT()  

feature_cols <- setdiff(colnames(data), c("sample_ID", "condition", "project"))

data <- data[, (feature_cols) := lapply(.SD, as.numeric), .SDcols = feature_cols]

meta_metafor <- function(dt, 
                         group = "group", 
                         group_pair = c("Disease", "Control"), 
                         proj = "proj", 
                         sample_id = "sample_ID",
                         measure = "SMD", 
                         method = "REML") {
  
  dt <- as.data.frame(dt)
  
  colnames(dt)[colnames(dt) == proj] <- "proj"
  colnames(dt)[colnames(dt) == group] <- "group"
  
  feature <- setdiff(colnames(dt), c(sample_id, "group", "proj"))
  
  results_list <- list()
  x_counter <- 1
  nfeature <- length(feature)
  
  for (i in feature) {
    cat("\rProcessing feature", x_counter, "/", nfeature)
    x_counter <- x_counter + 1
    
    idx1 <- dt$group %in% group_pair[1]
    if (sum(idx1) == 0) next 
    sub1 <- dt[idx1, c(i, "proj"), drop = FALSE]
    colnames(sub1)[1] <- "index"  
    
    tib <- sub1 %>%
      group_by(proj) %>%
      summarise(d_Mean = mean(index, na.rm = TRUE),
                d_Sd   = sd(index, na.rm = TRUE),
                d_N    = n(), .groups = "drop")
    
    idx2 <- dt$group %in% group_pair[2]
    if (sum(idx2) == 0) next
    sub2 <- dt[idx2, c(i, "proj"), drop = FALSE]
    colnames(sub2)[1] <- "index"
    
    tib2 <- sub2 %>%
      group_by(proj) %>%
      summarise(c_Mean = mean(index, na.rm = TRUE),
                c_Sd   = sd(index, na.rm = TRUE),
                c_N    = n(), .groups = "drop")
    
    meta_in <- merge(tib, tib2, by = "proj")
    if (nrow(meta_in) == 0) next

    smd_meta <- escalc(measure = measure, data = meta_in, append = TRUE,
                       m1i = d_Mean, m2i = c_Mean,
                       sd1i = d_Sd, sd2i = c_Sd,
                       n1i = d_N, n2i = c_N)
    
    smd_meta <- smd_meta %>% filter(!is.na(yi))
    if (nrow(smd_meta) < 2) next  
    
    smd_rma <- tryCatch(rma(yi, vi, method = method, data = smd_meta),
                        error = function(e) NULL)
    if (is.null(smd_rma)) next
    
    result <- smd_rma$data %>%
      mutate(measure     = measure,
             model       = "RM",
             method_tau2 = method,
             val_tau2    = as.numeric(smd_rma$tau2),
             I2          = paste0(round(smd_rma$I2, 2), "%"),
             Q           = smd_rma$QE,
             Q_pval      = smd_rma$QEp,
             feature     = i,
             estimate    = as.numeric(smd_rma$beta),
             ci_lb       = smd_rma$ci.lb,
             ci_ub       = smd_rma$ci.ub,
             pval        = smd_rma$pval)
    
    results_list[[i]] <- result
  }
  
  final_result <- bind_rows(results_list)
  
  cat("\n------------------------------------------------\n")
  cat("Interpretation:\n")
  cat("  estimate > 0  => feature is", group_pair[1], "enriched in\n")
  cat("  estimate < 0  => feature is", group_pair[2], "enriched in\n")
  cat("------------------------------------------------\n")
  
  return(final_result)
}

data_df <- as.data.frame(data)
x <- meta_metafor(data_df,
                  group       = "condition",
                  group_pair  = c("periodontitis", "healthy controls"),
                  proj        = "project",
                  sample_id   = "sample_ID",
                  measure     = "SMD",
                  method      = "REML")

write.table(x, "meta-analysis_saliva.csv", sep = ",", row.names = FALSE, quote = FALSE)