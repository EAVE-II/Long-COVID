##############################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 05_balance_checks - Assesses balance of covariates before and after matching by estimating
## standardised mean differences for treated and untreated. Exports tables and plots. 
##############################################################################################################

# 0. Set-up ------

# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(sandwich) # for estimating robust standard errors
library(Hmisc) # for wtd.mean
library(spatstat.geom) # for weighted variance (used to estimate standard mean differences = balance check)

# 1. Balance checks ---- 

## Function estimates standardized mean differences (SMD) in each covariate for treated (i.e. positive PCR test) versus untreated individuals in matched and
## unmatched samples, then outputs a table of results
### SMD: (weighted mean of the covariate for treated individuals minus weighted mean of the covariate for untreated individuals)/pooled SD of the covariate (all individuals)
### SMDs close to zero indicate good balance: .1 = an acceptable SMD threshold, or .05 for prognostically important covariates

balance_checks <- function(df_all_path,      # df_all_path = path of dataframe including matched and unmatched individuals 
                           df_matched_path){ # df_matched_path = path of dataframe including matched individuals only
                    
  # Read in data
  df_all <- readRDS(df_all_path)
  df_matched <- readRDS(df_matched_path)
  
  # Create treated vs untreated identifier that's consistent with matched dfs (pos_or_not = name included in matched dfs)
  df_all <- df_all %>% mutate(pos_or_not = ifelse(df_all$tes_res == "POSITIVE", 1, 0))
  
  # Identify dataframes to estimate SMDs for
  dfs <- list(df_all, df_matched)
  
  # Create empty dataframe to collect results in
  table <- NULL
  
  # Iterate over dataframes
  for (df in dfs) {
    
    # Prepare binary variables
    df <- df %>%
      mutate(mal = ifelse(sex == "M", 1, 0),
             fem = ifelse(sex == "F", 1, 0),
             age1 = ifelse(age_10 == "[18,28)", 1, 0),# 18 - 27
             age2 = ifelse(age_10 == "[28,38)", 1, 0),# 28 - 37
             age3 = ifelse(age_10 == "[38,48)", 1, 0),# 38 - 47
             age4 = ifelse(age_10 == "[48,58)", 1, 0),# 48 - 57
             age5 = ifelse(age_10 == "[58,68)", 1, 0),# 58 - 67
             age6 = ifelse(age_10 == "[68,78)", 1, 0),# 68 - 77
             age7 = ifelse(age_10 == "[78,88)", 1, 0),# 78 - 87
             age8 = ifelse(age_10 == "[88,98)" | age_10 == "[99,150]", 1, 0),# 88+
             SIMD1 = ifelse(simd == "1 - High", 1, 0), # Quintile 1 (most deprived)
             SIMD2 = ifelse(simd == 2, 1, 0), # Quintile 2 (most deprived)
             SIMD3 = ifelse(simd == 3, 1, 0), # Quintile 3 (most deprived)
             SIMD4 = ifelse(simd == 4, 1, 0), # Quintile 4 (most deprived)
             SIMD5 = ifelse(simd == "5-Low", 1, 0), # Quintile 5 (least deprived)
             hh1 = ifelse(hhold_bands == "1", 1, 0), # 1 individual in the household
             hh2 = ifelse(hhold_bands == "2", 1, 0), # 2 individuals in the household
             hh3 = ifelse(hhold_bands == "3-5", 1, 0), # 3-5 individuals in the household
             hh4 = ifelse(hhold_bands == "6-10", 1, 0), # 6-10 individuals in the household
             hh5 = ifelse(hhold_bands == "11-30" |
                          hhold_bands == "31-100" |
                          hhold_bands == "100+", 1, 0), # >10 individuals in the household
             ur1 = ifelse(ur6_lab == "1 Large Urban Areas", 1, 0), 
             ur2 = ifelse(ur6_lab == "2 Other Urban Areas", 1, 0),
             ur3 = ifelse(ur6_lab == "3 Accessible Small Towns", 1, 0),
             ur4 = ifelse(ur6_lab == "4 Remote Small Towns", 1, 0),
             ur5 = ifelse(ur6_lab == "5 Accessible Rural", 1, 0),
             ur6 = ifelse(ur6_lab == "6 Remote Rural", 1, 0),
             vac0 = ifelse(vac_n_clean == 0, 1, 0), # 0 vaccination doses received (up to 14 days before test)
             vac1 = ifelse(vac_n_clean == 1, 1, 0), # 1
             vac2 = ifelse(vac_n_clean == 2, 1, 0), # 2
             vac3 = ifelse(vac_n_clean == 3, 1, 0), # 3
             vac4 = ifelse(vac_n_clean == 4, 1, 0), # 4
             bmi_imp_bands = cut(bmi_imp, breaks = c(0, 18.5, 24.9, 29.9, 80), labels=c("Underweight","Normal weight","Overweight","Obese")),
             bmi_imp_bands = ifelse(is.na(bmi_imp_bands), "Unknown", bmi_imp_bands),
             bmi1 = ifelse(bmi_imp_bands == 1, 1, 0), # Underweight (< 18.5)
             bmi2 = ifelse(bmi_imp_bands == 2, 1, 0), # Normal weight (18.5 - 24.9)
             bmi3 = ifelse(bmi_imp_bands == 3, 1, 0), # Overweight (25 - 29.9)
             bmi4 = ifelse(bmi_imp_bands == 4, 1, 0), # Obese (>29.9) 
             n_risk0 = ifelse(n_risk_gps == "0", 1, 0), # 0 existing health conditions 
             n_risk1 = ifelse(n_risk_gps == "1", 1, 0), # 1 existing health conditions
             n_risk2 = ifelse(n_risk_gps == "2", 1, 0), # 2 existing health conditions
             n_risk3 = ifelse(n_risk_gps > 2, 1, 0) # 3+ existing health conditions - n is too small to split out 3-4 & 5+
             
      )
    
    # Create empty objects to collect values in
    wtd_mean_dif_vals <- NULL
    sd_pooled_vals <- NULL
    smd_vals <- NULL
    
    # Identify covariates to include in table and plots
    covariates <- list(#propensity score
                       df$prop_score,
                       # sex
                       df$mal, df$fem,
                       # age
                       df$age1, df$age2, df$age3, df$age4, df$age5, df$age6, df$age7, df$age8, 
                       # simd
                       df$SIMD1, df$SIMD2, df$SIMD3, df$SIMD4, df$SIMD5,
                       # household size
                       df$hh1, df$hh2, df$hh3, df$hh4, df$hh5,
                       # urban-rural settlement
                       df$ur1, df$ur2, df$ur3, df$ur4, df$ur5, df$ur6,
                       # vaccine doses (up to 14 days before test date)
                       df$vac0, df$vac1, df$vac2, df$vac3, df$vac4,
                       # bmi
                       df$bmi1, df$bmi2, df$bmi3, df$bmi4,
                       # n risk groups
                       df$n_risk0, df$n_risk1, df$n_risk2, df$n_risk3,  
                       # immunosuppressed, shielding
                       df$immuno_supp,
                       df$shielding,
                       # hospitalised in 365 days before test
                       df$hos_prior)
    
    # Iterate over covariates
    for (cov in covariates) {
      
      # Estimate standardised mean difference (SMD): (weighted mean treated - weighted mean untreated)/pooled standard deviation
      ## Difference in weighted means
      wtd_mean_pos <- wtd.mean(cov[df$pos_or_not == 1], df$eave_weight[df$pos_or_not == 1], na.rm = T)
      wtd_mean_not <- wtd.mean(cov[df$pos_or_not == 0], df$eave_weight[df$pos_or_not == 0], na.rm = T)
      wtd_mean_dif <- wtd_mean_pos - wtd_mean_not

      ## Weighted variance (used to estimated standard deviation)
      wtd_var_pos <- Hmisc::wtd.var(cov[df$pos_or_not == 1], df$eave_weight[df$pos_or_not == 1], na.rm = T)
      wtd_var_neg <- Hmisc::wtd.var(cov[df$pos_or_not == 0], df$eave_weight[df$pos_or_not == 0], na.rm = T)

      ## Pooled standard deviation
      sd_pooled <- sqrt((wtd_var_pos + wtd_var_neg) / 2)

      ## Standardised mean difference (SMD)
      smd <- wtd_mean_dif / sd_pooled

      ## Collect weightd mean, pooled standard deviation, SMD for each coefficient
      wtd_mean_dif_vals <- c(wtd_mean_dif_vals, wtd_mean_dif)
      sd_pooled_vals <- c(sd_pooled_vals, sd_pooled)
      smd_vals <- c(smd_vals, smd)
      vals <- cbind(wtd_mean_dif_vals, sd_pooled_vals, smd_vals)

    }

    # Collect results for all coefficients in each df
    table <- cbind(table, vals)
    table

  }
  
  colnames(table) <- c("Weighted mean differences: Full sample", 
                       "Standard deviation: Full sample",
                       "Standardised mean differences: Full sample", 
                       "Weighted mean differences: Matched sample", 
                       "Standard deviation: Matched sample",
                       "Standardised mean differences: Matched sample")
  
  table <- round(table,2)
  
  row.names(table) <- c("Propensity score",
                        "Male", "Female",
                        "18 - 27", "28 - 37", "38 - 47", "48 - 57", "58 - 67", "68 - 77", "78 - 87", "88 - 100",  
                        "SIMD 1 (Most deprived)", "SIMD 2", "SIMD 3", "SIMD 4", "SIMD 5 (Least deprived)", 
                        "Household size: 1", "Household size: 2", "Household size: 3-5", "Household size: 6-10", "Household size: >10",
                        "Large urban areas", "Other urban areas", "Accessible small towns", "Remote small towns", 
                        "Accessible rural", "Remote rural", 
                        "0 vaccine doses", "1 vaccine dose", "2 vaccine doses", "3 vaccine doses", "4 vaccine doses",
                        "Underweight (BMI <18.5)", "Normal weight (BMI 18.5 - 24.9)", "Overweight (BMI 25 - 29.9)", "Obese (BMI >29.9)", 
                        "0 risk groups", "1 risk group", "2 risk groups", "3+ risk groups", 
                        "Immunosuppressed",
                        "Shielding",
                        "Hospitalised during 12 months pre-test")
  
  table <- as.data.frame(table)
  
  return(table)
}


# 2. Plotting function: Generates and exports balance plots ----
balance_plots <- function(tab_neg,    # tab_neg = table of balance checks results - positive vs negative e.g. tab_neg, tab_neg_wild
                          tab_nt,     # tab_nt = table of balance checks results - positive vs never tested e.g. tab_nt, tab_nt_wild
                          plot_name){ # plot_name = name to save the plot with e.g. "Balance tests - full sample.png"
  
  # Format data for plotting
  bal_neg <- tab_neg %>%
    mutate(exposure_control = "Positive vs negative",
           "Covariate" = row.names(.))
  row.names(bal_neg) <- NULL         
  
  bal_nt <- tab_nt %>% 
    mutate(exposure_control = "Positive vs not yet tested",
           "Covariate" = row.names(.))
  row.names(bal_nt) <- NULL  
  
  plot_data <- rbind(bal_neg, bal_nt)
  
  plot_data <- plot_data %>% 
    rename("Full sample" = "Standardised mean differences: Full sample", 
           "Matched sample" = "Standardised mean differences: Matched sample") %>% 
    pivot_longer(cols = c("Full sample", "Matched sample"),
                 names_to = "Sample",
                 values_to = "SMD") %>% 
    select(Covariate, Sample, SMD, exposure_control)
  
  ## Set order of variables
  plot_data$Covariate <- factor(plot_data$Covariate, levels = rev(unique(plot_data$Covariate)))
  
  ## Plot
  balance_plots <- ggplot(plot_data, aes(x = SMD, y = Covariate, shape = Sample, colour = Sample), size=3) +
    geom_vline(xintercept = 0) + 
    geom_vline(xintercept = c(-0.1, 0.1), linetype=2) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype=3) +
    geom_point(size=2, position = position_jitter(w = 0, h = 0)) +
    theme_bw() +
    facet_grid(~exposure_control) +
    theme(strip.background = element_blank(),
          strip.text = element_text(colour = "black", hjust = 0, size=10), 
          legend.position = "bottom", 
          legend.title=element_blank(), 
          panel.spacing = unit(2, "lines")) + # adjust space between plots
    labs(x="Standardised mean differences", y="", title = "",
         shape = "Sample", color = "Sample")
  
  ## Export plots
  ggsave(plot_name,
         plot = balance_plots,
         width = 8, height = 10,
         units = "in", 
         dpi = 300)
}

# 3. Generate and export tables and plots ----

## Full dataset ----
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes")

### Generate tables
tab_neg <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs neg.rds",
                          df_matched_path = "Matched pairs - pos vs neg.rds") 

tab_nt <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs nt.rds",
                          df_matched_path = "Matched pairs - pos vs nt.rds") 

### Export tables
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/3. Balance checks")
write.csv(tab_neg, "Balance tests - pos vs neg.csv", row.names = T)
write.csv(tab_nt, "Balance tests - pos vs nt.csv", row.names = T)

### Generate and export plots
balance_plots(tab_neg = tab_neg, 
              tab_nt = tab_nt, 
              plot_name = "Balance tests - full sample.png")

### Remove files
rm(tab_neg, tab_nt)


## Wild ----
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes")

### Generate tables
tab_neg_wild <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs neg  - wild.rds",
                               df_matched_path = "Matched pairs - pos vs neg - wild.rds") 

tab_nt_wild <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs nt  - wild.rds",
                              df_matched_path = "Matched pairs - pos vs nt - wild.rds") 

### Export tables
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/3. Balance checks")
write_csv(tab_neg_wild, "Balance tests - pos vs neg - wild.csv")
write_csv(tab_nt_wild, "Balance tests - pos vs nt - wild.csv")

### Generate and export plots
balance_plots(tab_neg = tab_neg_wild, 
              tab_nt = tab_nt_wild, 
              plot_name = "Balance tests - wild.png")

### Remove files
rm(tab_neg_wild, tab_nt_wild)


## Alpha ----
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes")

### Generate tables
tab_neg_alpha <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs neg  - alpha.rds",
                               df_matched_path = "Matched pairs - pos vs neg - alpha.rds") 

tab_nt_alpha <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs nt  - alpha.rds",
                              df_matched_path = "Matched pairs - pos vs nt - alpha.rds") 

### Export tables
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/3. Balance checks")
write_csv(tab_neg_alpha, "Balance tests - pos vs neg - alpha.csv")
write_csv(tab_nt_alpha, "Balance tests - pos vs nt - alpha.csv")

### Generate and export plots
balance_plots(tab_neg = tab_neg_alpha, 
              tab_nt = tab_nt_alpha, 
              plot_name = "Balance tests - alpha.png")

### Remove files
rm(tab_neg_alpha, tab_nt_alpha)


## Delta ----
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes")

### Generate tables
tab_neg_delta <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs neg  - delta.rds",
                               df_matched_path = "Matched pairs - pos vs neg - delta.rds") 

tab_nt_delta <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs nt  - delta.rds",
                              df_matched_path = "Matched pairs - pos vs nt - delta.rds") 

### Export tables
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/3. Balance checks")
write_csv(tab_neg_delta, "Balance tests - pos vs neg - delta.csv")
write_csv(tab_nt_delta, "Balance tests - pos vs nt - delta.csv")

### Generate and export plots
balance_plots(tab_neg = tab_neg_delta, 
              tab_nt = tab_nt_delta, 
              plot_name = "Balance tests - delta.png")

### Remove files
rm(tab_neg_delta, tab_nt_delta)


## Omicron ----
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes")

### Generate tables
tab_neg_omicron <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs neg  - omicron.rds",
                               df_matched_path = "Matched pairs - pos vs neg - omicron.rds") 

tab_nt_omicron <- balance_checks(df_all_path = "Linked df with propensity scores - pos vs nt  - omicron.rds",
                              df_matched_path = "Matched pairs - pos vs nt - omicron.rds") 

### Export tables
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/3. Balance checks")
write_csv(tab_neg_omicron, "Balance tests - pos vs neg - omicron.csv")
write_csv(tab_nt_omicron, "Balance tests - pos vs nt - omicron.csv")

### Generate and export plots
balance_plots(tab_neg = tab_neg_omicron, 
              tab_nt = tab_nt_omicron, 
              plot_name = "Balance tests - omicron.png")

### Remove files
rm(tab_neg_omicron, tab_nt_omicron)


