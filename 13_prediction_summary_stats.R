#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 12_prediction_summary_stats - Summary stats showing breakdown of cohort used in prediction modelling.
#########################################################################################################################

# Set-up ------
# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(readxl)
library(zoo) # for dates
library(scales) # for 'label = comma' in ggplot
library(questionr) # for wtd.table

# Code for negation of %in%
`%!in%` = Negate(`%in%`) 

# Set location to export to
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling/Summary stats")

## Read in data prepared for prediction modelling
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/df_training_cleaned.rds") 
df_testing <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/df_testing_with_probs.rds")

## Set dependent variable
df$depvar = ifelse(df$opdef == 1 | df$code_txt_fitnote == 1, 1, 0)
df_testing$depvar = ifelse(df_testing$opdef == 1 | df_testing$code_txt_fitnote == 1, 1, 0)

## Check for missing values
missing_vals <- colSums(is.na(df))
print("Check for missing values")
missing_vals[missing_vals>0]

# 1. Descriptive stats - table ----
## Calculate descriptives
calculate_descriptives <- function(df, 
                                   mask){ # mask = the subset of interest e.g. df$`Operational definition` > 0
                                   
                                   
  ## Get population total
  pop_n <- sum(df$eave_weight) # NB all weights are set to 1
  
  ## Subset data of interest
  df <- df %>% filter(mask)
  
  ## Get N
  ### Sex ----
  sex_tab <- round(wtd.table(df$sex, exclude = NULL, weights = df$eave_weight), 0)
  sex_n <- ""
  fem_n <- sex_tab[1]
  mal_n <- sex_tab[2]
  
  ### Age ----
  age_tab <- round(wtd.table(df$age_10, exclude = NULL, weights = df$eave_weight), 0)
  age_n <- " "
  age1_n<- age_tab[1]
  age2_n<- age_tab[2]
  age3_n<- age_tab[3]
  age4_n<- age_tab[4]
  age5_n<- age_tab[5]
  age6_n<- age_tab[6]
  age7_n<- age_tab[7]
  age8_n<- sum(age_tab[8:9])
  
  ### SIMD ----
  simd_tab <- round(wtd.table(df$simd, exclude = NULL, weights = df$eave_weight), 0)
  simd_n <- ""
  simd1_n <- simd_tab[1]
  simd2_n <- simd_tab[2]
  simd3_n <- simd_tab[3]
  simd4_n <- simd_tab[4]
  simd5_n <- simd_tab[5]
  #simd6_n <- simd_tab[6]
  
  ### Household size ----
  hh_tab <- round(wtd.table(df$hhold_bands, weights = df$eave_weight, useNA = "always"), 0)
  hh_n <- ""
  hh1_n <- hh_tab[1] # 1
  hh2_n <- hh_tab[2] # 2
  hh3_n <- hh_tab[3] # 3-5
  hh4_n <- hh_tab[4] # 6-10
  hh5_n <- sum(hh_tab[5]) # 11+
  
  ### Urban-rural classification ----
  ur_tab <- round(wtd.table(df$ur6_lab, exclude = NULL, weights = df$eave_weight), 0)
  ur_n <- ""
  ur1_n <- ur_tab[1]
  ur2_n <- ur_tab[2]
  ur3_n <- ur_tab[3]
  ur4_n <- ur_tab[4]
  ur5_n <- ur_tab[5]
  ur6_n <- ur_tab[6]
  #ur7_n <- ur_tab[7]
  
  #### Tests ----
  #tes_tab <- round(wtd.table(df$tes_n, exclude = NULL, weights = df$eave_weight), 0)
  #tes_n <- ""
  #tes1_n <- ur_tab[1] # 1
  #tes2_n <- ur_tab[2] # 2
  #tes3_n <- ur_tab[3] # 3-4
  #tes4_n <- ur_tab[4] # 5-10
  #tes5_n <- ur_tab[5] # 11+
  
  ### Variant (based on test date) ----
  var_tab <- round(wtd.table(df$var_per, exclude = NULL, weights = df$eave_weight), 0)
  var_n <- ""
  var1_n <- var_tab[5] # Wild
  var2_n <- var_tab[1] # Alpha
  var3_n <- var_tab[2] # Delta
  var4_n <- var_tab[3] # Omicron
  var5_n <- var_tab[4] # Unknown
  
  ### Vaccine doses received by 14 days before test date ----
  vac_tab <- round(wtd.table(df$vac_dos, exclude = NULL, weights = df$eave_weight), 0)
  vac_n <- ""
  vac1_n <- vac_tab[1] # 0
  vac2_n <- vac_tab[2] # 1
  vac3_n <- vac_tab[3] # 2
  vac4_n <- vac_tab[4] # 3+
  
  ### Shielding ----
  she_tab <- round(wtd.table(df$shielding, exclude = NULL, weights = df$eave_weight), 0)
  she_n <- ""
  she1_n <- she_tab[2] # sheilding
  she2_n <- she_tab[1] # not shielding
  
  # Immunosuppressed ----
  imm_tab <- round(wtd.table(df$immuno_supp, exclude = NULL, weights = df$eave_weight), 0)
  imm_n <- ""
  imm1_n <- imm_tab[2] # immunosuppressed
  imm2_n <- imm_tab[1] # not immunosuppressed
  
  # Care home resident ----
  car_tab <- round(wtd.table(df$carehome, exclude = NULL, weights = df$eave_weight), 0)
  car_n <- ""
  car1_n <- car_tab[2] # carehome resident
  car2_n <- car_tab[1] # not carehome resident
  
  ### BMI (including imputed values) ----
  df <- df %>% 
    mutate(bmi_bands = cut(bmi_imp, breaks = c(0, 18.5, 24.9, 29.9, Inf), labels=c("Underweight","Normal weight","Overweight","Obese")),
           bmi_bands = ifelse(is.na(bmi_bands), "Unknown", bmi_bands))
  
  bmi_tab <- round(wtd.table(df$bmi_bands, exclude = NULL, weights = df$eave_weight), 0)
  bmi_n <- ""
  bmi1_n <- bmi_tab[1] # Underweight
  bmi2_n <- bmi_tab[2] # Normal weight
  bmi3_n <- bmi_tab[3] # Overweight
  bmi4_n <- bmi_tab[4] # Obese
  
  ### Comorbidities/Risk groups ----
  com_tab <- round(wtd.table(df$n_risk_gps, weights = df$eave_weight), 0)
  com_n <- ""
  com1_n <- com_tab[1] # 0
  com2_n <- com_tab[2] # 1
  com3_n <- com_tab[3]# 2
  com4_n <- sum(com_tab[-c(1:3)]) # 3+
  
 ### Severity of acute infection ----
 sev_tab <- round(wtd.table(df$severe, exclude = NULL, weights = df$eave_weight), 0)
 sev_n <- ""
 sev1_n <- sev_tab[2] # Hospitalised or admitted to ICU within 28 days of positive test
 sev2_n <- sev_tab[1] # Not hospitalised or admitted to ICU within 28 days of positive test
 
 ### Total ----
 tot_n <- round(sum(df$eave_weight), 0)
 
 N <- rbind(tot_n,
            sex_n, fem_n, mal_n,
            age_n, age1_n, age2_n, age3_n, age4_n, age5_n, age6_n, age7_n, age8_n,
            simd_n, simd1_n, simd2_n, simd3_n, simd4_n, simd5_n, #simd6_n, 
            hh_n, hh1_n, hh2_n, hh3_n, hh4_n, hh5_n,
            ur_n, ur1_n, ur2_n, ur3_n, ur4_n, ur5_n, ur6_n, #ur7_n,
            #tes_n, tes1_n, tes2_n, tes3_n, tes4_n, tes5_n,
            var_n, var1_n, var2_n, var3_n, var4_n, var5_n,
            vac_n, vac1_n, vac2_n, vac3_n, vac4_n, 
            she_n, she1_n, she2_n,
            imm_n, imm1_n, imm2_n,
            car_n, car1_n, car2_n,
            bmi_n, bmi1_n, bmi2_n, bmi3_n, bmi4_n,
            com_n, com1_n, com2_n, com3_n, com4_n, 
            sev_n, sev1_n, sev2_n)
 
  ## Suppress small numbers
  N <- as.numeric(N)
  N[!is.na(N) & N<5] <- -5
  
  ## Get %
  PC <- as_tibble(format(round(as.numeric(N)/tot_n*100, 1), nsmall=1L))
  colnames(PC) <- "PC"
  
  ## For total row, make % of total population
  PC$PC[1] <- round(as.numeric(N[1])/pop_n*100, 1)
  
  ## Labels
  labels <- rbind("Total (% of dataset)",
                  "Sex", "Female", "Male",
                  "Age", "18 - 27", "28 - 37", "38 - 47", "48 - 57", "58 - 67", "68 - 77", "78 - 87", "88 - 100",  
                  "SIMD quintiles", "1 - Most deprived", "2", "3", "4", "5 - Least deprived", #"Unknown", 
                  "Household size", "1", "2", "3-5", "6-10", "11+",
                  "Urban-Rural", "Large urban areas", "Other urban areas", "Accessible small towns", "Remote small towns", "Accessible rural", "Remote rural", #"Unknown", 
                  #"PCR tests by first positive PCR", "1", "2", "3-4", "5-10", "11+",
                  "Variant period", "Wild (up to 10/01/2021)", "Alpha (11/01/2021 - 09/05/2021)", "Delta (24/05/2021 - 28/11/2021)", "Omicron (20/12/2021 onwards)", "No dominant variant",
                  "Vaccination doses (up to 14 days before positive test/outcome)" , "0", "1", "2", "3+",
                  "Shielding", "Shielding", "Not shielding", 
                  "Immunosuppressed", "Immunosuppressed", "Not immunosuppressed",
                  "Carehome resident", "Carehome resident", "Not carehome resident",
                  "BMI", "Underweight (BMI < 18.5)", "Normal weight (BMI 18.5 - 24.9)", "Overweight (BMI 25 - 29.9)", "Obese (BMI >29.9)", 
                  "Comorbidities", "0", "1", "2", "3+",
                  "Severity of acute infection (positive cases only)", "Hospitalised within 28 days", "Not hospitalised within 28 days")
  
  ## Combine
  results <- cbind(labels, N, PC)
  results$N <- ifelse(is.na(results$N), "", results$N)
  results$PC <- gsub("NA", "", results$PC)
  
  return(results)
}

## Get descriptives
training_noLC <- data.frame(calculate_descriptives(df, mask = df$depvar == 0))
training_LC <- data.frame(calculate_descriptives(df, mask = df$depvar == 1))
holdout_noLC <- data.frame(calculate_descriptives(df_testing, mask = df_testing$depvar == 0))
holdout_LC <- data.frame(calculate_descriptives(df_testing, mask = df_testing$depvar == 1))


descr_stats <- cbind(training_noLC, 
                     training_LC[2:3],
                     holdout_noLC[2:3],
                     holdout_LC[2:3])

colnames(descr_stats) <- c("", "Training \nNo long COVID (N)", "Training \nNo long COVID  (%)", 
                          "Training \nLong COVID  (N)", "Training \nLong COVID  (%)",
                          "Hold out \nNo long COVID (N)", "Hold out \nNo long COVID  (%)", 
                          "Hold out \nLong COVID  (N)", "Hold out \nLong COVID  (%)")


View(descr_stats)

# Export
write_csv(descr_stats, "Summary stats.csv")



# 2. Q stats - table ----
## Calculate descriptives
Q_descriptives <- function(df, mask){ # mask = the subset of interest e.g. df$`Operational definition` > 1 
                           
                                   
  ## Get population total
  pop_n <- sum(df$eave_weight)
  
  ## Subset data of interest 
  df <- df %>% filter(mask) 
  
  ## Get N for each Q-covid variable
  q_vars <- names(df)[grep("^Q_", names(df))]
  
  N <- NULL
  
  for(var in q_vars){
    tab <- round(wtd.table(df[[var]], exclude = NULL, weights = df$eave_weight), 0)
    tab_n <- tab[2]
    N <- rbind(N, tab_n)
  }
  
  ### Total
  tot_n <- round(sum(df$eave_weight), 0)
  
  N <- rbind(tot_n, N)
  
  ## Suppress small numbers
  N <- as.numeric(N)
  N[!is.na(N) & N<5] <- -5
  
  ## Get %
  PC <- as_tibble(format(round(as.numeric(N)/tot_n*100, 1), nsmall=1L))
  colnames(PC) <- "PC"
  
  ## For total row, make % of total population
  PC$PC[1] <- round(as.numeric(N[1])/pop_n*100, 1)
  
  ## Labels
  labels <- "Total (% of dataset)"
  
  q_vars[q_vars == "Q_DIAG_DIABETES_1"] <- "Diabetes Type I"
  q_vars[q_vars == "Q_DIAG_DIABETES_2"] <- "Diabetes Type II"
  q_vars[q_vars == "Q_DIAG_AF"] <- "Atrial fibrillation"
  q_vars[q_vars == "Q_DIAG_ASTHMA"] <- "Asthma"
  q_vars[q_vars == "Q_DIAG_BLOOD_CANCER"] <- "Haematological cancer"
  q_vars[q_vars == "Q_DIAG_CCF"] <- "Heart failure"
  q_vars[q_vars == "Q_DIAG_CHD"] <- "Coronary heart disease"
  q_vars[q_vars == "Q_DIAG_COPD"] <- "Chronic obstructive pulmonary disease (COPD)"
  q_vars[q_vars == "Q_DIAG_DEMENTIA"] <- "Dementia"
  q_vars[q_vars == "Q_DIAG_EPILEPSY"] <- "Epilepsy"
  q_vars[q_vars == "Q_DIAG_FRACTURE"] <- "Fracture"
  q_vars[q_vars == "Q_DIAG_NEURO"] <- "Neurological disorder"
  q_vars[q_vars == "Q_DIAG_PARKINSONS"] <- "Parkinson’s disease"
  q_vars[q_vars == "Q_DIAG_PULM_HYPER"] <- "Pulmonary hypertension"
  q_vars[q_vars == "Q_DIAG_PULM_RARE"] <- "Rare pulmonary disease"
  q_vars[q_vars == "Q_DIAG_PVD"] <- "Peripheral vascular disease"
  q_vars[q_vars == "Q_DIAG_RA_SLE"] <- "Rheumatoid arthritis or systemic lupus erythematous (SLE)"
  q_vars[q_vars == "Q_DIAG_RESP_CANCER"] <- "Respiratory cancer"
  q_vars[q_vars == "Q_DIAG_SEV_MENT_ILL"] <- "Severe mental illness"
  q_vars[q_vars == "Q_DIAG_STROKE"] <- "Stroke/Transient Ischaemic Attack (TIA)"
  q_vars[q_vars == "Q_DIAG_VTE"] <- "Thrombosis or pulmonary embolus"
  q_vars[q_vars == "Q_DIAG_CKD_LEVEL"] <- "Chronic Kidney disease (level 3+)"
  
  labels <- c(labels, q_vars)
  
  ## Combine
  results <- cbind(labels, N, PC)
  results$PC <- gsub("NA", "", results$PC)
  results$N <- ifelse(is.na(results$N), "", results$N)
  
  return(results)
}

## Get Q descriptives
Q_training_noLC <- data.frame(Q_descriptives(df, mask = df$depvar == 0))
Q_training_LC <- data.frame(Q_descriptives(df, mask = df$depvar == 1))
Q_holdout_noLC <- data.frame(Q_descriptives(df_testing, mask = df_testing$depvar == 0))
Q_holdout_LC <- data.frame(Q_descriptives(df_testing, mask = df_testing$depvar == 1))


Q_stats <- cbind(Q_training_noLC, 
                 Q_training_LC[2:3],
                 Q_holdout_noLC[2:3],
                 Q_holdout_LC[2:3])


## Sort table so that rows are in alphabetical order with Total row at the top
colnames(Q_stats) <- c("label", colnames(descr_stats)[2:length(colnames(descr_stats))])

Q_stats_sorted <- Q_stats %>% 
  filter(label != "Total (% of dataset)") %>% 
  arrange(label)

Q_stats <- rbind(Q_stats[1,], Q_stats_sorted)

colnames(Q_stats) <- colnames(descr_stats)

# Export
write_csv(Q_stats, "Summary stats - Q.csv")


# 3.Prescriptions stats - table ----
## Calculate descriptives
PIS_descriptives <- function(df, mask){ # mask = the subset of interest e.g. df$`Operational definition` > 1 
  
  
  ## Get population total
  pop_n <- sum(df$eave_weight)
  
  ## Subset data of interest 
  df <- df %>% filter(mask) 
  
  ## Get N for each Q-covid variable
  pis_vars <- names(df)[grep("^PIS_", names(df))]
  
  N <- NULL
  
  for(var in pis_vars){
    tab <- round(wtd.table(df[[var]], exclude = NULL, weights = df$eave_weight), 0)
    tab_n <- tab[2]
    N <- rbind(N, tab_n)
  }
  
  ### Total
  tot_n <- round(sum(df$eave_weight), 0)
  
  N <- rbind(tot_n, N)
  
  ## Suppress small numbers
  N <- as.numeric(N)
  N[!is.na(N) & N<5] <- -5
  
  ## Get %
  PC <- as_tibble(format(round(as.numeric(N)/tot_n*100, 1), nsmall=1L))
  colnames(PC) <- "PC"
  
  ## For total row, make % of total population
  PC$PC[1] <- round(as.numeric(N[1])/pop_n*100, 1)
  
  ## Labels
  labels <- "Total (% of dataset)"
  
  pis_vars[pis_vars == "PIS_Lipid_regulating_drugs"] <- "Lipid-regulating drugs"
  pis_vars[pis_vars == "PIS_Angiotensin_converting_enzyme_inhibitors"] <- "Angiotensin-converting enzyme inhibitors"
  pis_vars[pis_vars == "PIS_Beta_adrenoceptor_blocking_drugs"] <- "Beta-adrenoceptor blocking drugs"
  pis_vars[pis_vars == "PIS_Oral_iron"] <- "Oral iron"
  pis_vars[pis_vars == "PIS_Selective_serotonin_re_uptake_inhibitors"] <- "Selective serotonin re-uptake inhibitors"
  pis_vars[pis_vars == "PIS_Loratadine"] <- "Loratadine (antihistamine)"
  pis_vars[pis_vars == "PIS_Direct_oral_anticoagulants"] <- "Direct oral anticoagulants"
  pis_vars[pis_vars == "PIS_Colchicine"] <- "Colchicine (anti-inflammatory)"
  pis_vars[pis_vars == "PIS_Antiplatelet_drugs"] <- "Antiplatelet drugs"
  pis_vars[pis_vars == "PIS_Alpha_adrenoceptor_blocking_drugs"] <- "Alpha-adrenoceptor blocking drugs"
  pis_vars[pis_vars == "PIS_Macrolides"] <- "Macrolides (antibacterial)"
  pis_vars[pis_vars == "PIS_Benzylpenicillin_and_phenoxymethylpenicillin"] <- "Benzylpenicillin and phenoxymethylpenicillin"
  pis_vars[pis_vars == "PIS_Leukotriene_receptor_antagonists"] <- "Leukotriene receptor antagonists"
  pis_vars[pis_vars == "PIS_Herpes_simplex_and_varicella_zoster"] <- "Herpes simplex and varicella-zoster (antiviral)"
  pis_vars[pis_vars == "PIS_Replacement_therapy"] <- "Corticosteroid replacement therapy"
  pis_vars[pis_vars == "PIS_Famotidine"] <- "Famotidine (histamine H2 receptor antagonist)"
  pis_vars[pis_vars == "PIS_Ursodeoxycholic_acid"] <- "Ursodeoxycholic acid"
  pis_vars[pis_vars == "PIS_Systemic_nasal_decongestants"] <- "Systemic nasal decongestants"
  pis_vars[pis_vars == "PIS_Warfarin_sodium"] <- "Warfarin sodium"
  pis_vars[pis_vars == "PIS_Parenteral_anticoagulants"] <- "Parenteral anticoagulants"
  pis_vars[pis_vars == "PIS_Compound_bronchodilator_preparations"] <- "Compound bronchodilator preparations"
  pis_vars[pis_vars == "PIS_Ranitidine_hydrochloride"] <- "Ranitidine hydrochloride"
  
  labels <- c(labels, pis_vars)
  
  ## Combine
  results <- cbind(labels, N, PC)
  results$PC <- gsub("NA", "", results$PC)
  results$N <- ifelse(is.na(results$N), "", results$N)
  
  return(results)
}

## Get PIS descriptives
PIS_training_noLC <- data.frame(PIS_descriptives(df, mask = df$depvar == 0))
PIS_training_LC <- data.frame(PIS_descriptives(df, mask = df$depvar == 1))
PIS_holdout_noLC <- data.frame(PIS_descriptives(df_testing, mask = df_testing$depvar == 0))
PIS_holdout_LC <- data.frame(PIS_descriptives(df_testing, mask = df_testing$depvar == 1))


PIS_stats <- cbind(PIS_training_noLC, 
                   PIS_training_LC[2:3],
                   PIS_holdout_noLC[2:3],
                   PIS_holdout_LC[2:3])

## Sort table so that rows are in alphabetical order with Total row at the top
colnames(PIS_stats) <- c("label", colnames(descr_stats)[2:length(colnames(descr_stats))])

PIS_stats_sorted <- PIS_stats %>% 
  filter(label != "Total (% of dataset)") %>% 
  arrange(label)

PIS_stats <- rbind(PIS_stats[1,], PIS_stats_sorted)

colnames(PIS_stats) <- colnames(descr_stats)

# Export
write_csv(PIS_stats, "Summary stats - PIS.csv")
