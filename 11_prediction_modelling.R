##########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 11_prediction_modelling - Checks for multicollinearity. Selects the optimal subset of predictors based
## on AIC/BIC and relaxed LASSO. Trains model using 10-fold cross validation. Get predicted probabilities using 
## (i) logistic regression, (ii) Naive Bayes Classfication, (iii) Gradient boosting decision trees.
## Sensitivity: runs main model (logistic regression) using alternative dependent variables. 
## Exports df with predicted probabilities from each version of the model.
##########################################################################################################################

# 0. Set-up ------
# Clear environment 
rm(list=ls())
setwd("/conf/EAVE/GPanalysis/analyses/long_covid")

# Libraries
library(tidyverse)
library(lubridate)
library(corrplot) # for correlograms
library(splines)
library(visreg) # for plotting splines
library(MASS) # for stepwise predictor selection using AIC and BIC ****Masks dplyr::select****
library(lmtest) # for likelihood ratio test
library(fmsb) # for Nagelkerke's R-squared
library(pROC) # for c-statistic/AUROC
library(ggplot2)
library(car) # for assessing collinearity between continuous and categorical variables using Variance Inflation Factor
library(pmsampsize) # to check for overfitting
library(DescTools) # to get Cox Snell Rsquared (for use in overfitting estimations)
library(stringr)
library(reshape2) # for melt function to reshape correlation matrix
library(caret) # for generating confusion matrix and k-fold cross-validation
library(yardstick) # for precision recall curves
library(glmnet) # for LASSO (use version 4.1-2 as C++17 not available)
library(e1071) # for Naive Bayes classifier
library(questionr) # for weighted tables
library(MLmetrics) # for F1_Score
library(naivebayes) # for Naive Bayes Classifier
library(xgboost) # for Gradient boosting decision trees
library(randomForest)


## Read in prepared data
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_predict_prepared.rds") # for main analysis

## Set dependent variable to use # e.g. df$opdef, df$opdef_no_blood, df$code_txt_fit
df$depvar <- ifelse(df$opdef == 1 | df$code_txt_fitnote == 1, 1, 0)
depvar_name <- "op def code txt sick" # for saving plots 

## Check for missing values
missing_vals <- colSums(is.na(df))
print("Check for missing values")
missing_vals[missing_vals>0]


## 0.1 Cleaning ----
df <- df %>% 
  # Remove Q_COVID variables that aren't clinically interesting with respect to long COVID)
  dplyr::select(-c("Q_DIAG_CEREBRALPALSY", "Q_DIAG_CIRRHOSIS", "Q_DIAG_CONGEN_HD", "Q_DIAG_SICKLE_CELL")) %>% 
  # Group 4+ vaccine doses (few observations have 4, 5, 6)
  dplyr::mutate(vac_dos = as.character(vac_dos),
                vac_dos = ifelse(vac_dos == "3" | vac_dos == "4" | vac_dos == "5" | vac_dos == "6", "3+", vac_dos),
                vac_dos = as.factor(vac_dos))%>% 
  # Drop 'unknown' SIMD and ur
  filter(ur != "Unknown")# %>% 

df$ur <- droplevels(df$ur)
df$simd <- droplevels(df$simd)


### Rename variables with spaces in their name (otherwise causes issues for caret::predict)
df <- df %>% 
  rename(PIS_Lipid_regulating_drugs = "PIS_Lipid-regulating drugs",
         PIS_Angiotensin_converting_enzyme_inhibitors = "PIS_Angiotensin-converting enzyme inhibitors",
         PIS_Beta_adrenoceptor_blocking_drugs = "PIS_Beta-adrenoceptor blocking drugs",
         PIS_Oral_iron = "PIS_Oral iron",
         PIS_Selective_serotonin_re_uptake_inhibitors = "PIS_Selective serotonin re-uptake inhibitors",
         PIS_Metformin_hydrochloride = "PIS_Metformin hydrochloride",
         PIS_Loratadine = "PIS_Loratadine",
         PIS_Direct_oral_anticoagulants = "PIS_Direct oral anticoagulants",
         PIS_Colchicine = "PIS_Colchicine",
         PIS_Antiplatelet_drugs = "PIS_Antiplatelet drugs",
         PIS_Alpha_adrenoceptor_blocking_drugs = "PIS_Alpha-adrenoceptor blocking drugs",
         PIS_Macrolides = "PIS_Macrolides",
         PIS_Benzylpenicillin_and_phenoxymethylpenicillin = "PIS_Benzylpenicillin and phenoxymethylpenicillin",
         PIS_Leukotriene_receptor_antagonists = "PIS_Leukotriene receptor antagonists",
         PIS_Herpes_simplex_and_varicella_zoster = "PIS_Herpes simplex and varicella-zoster",
         PIS_Replacement_therapy  = "PIS_Replacement therapy",
         PIS_Famotidine = "PIS_Famotidine",
         PIS_Ursodeoxycholic_acid = "PIS_Ursodeoxycholic acid",
         PIS_Systemic_nasal_decongestants = "PIS_Systemic nasal decongestants",
         PIS_Warfarin_sodium  = "PIS_Warfarin sodium",
         PIS_Parenteral_anticoagulants = "PIS_Parenteral anticoagulants",
         PIS_Compound_bronchodilator_preparations = "PIS_Compound bronchodilator preparations",
         PIS_Ranitidine_hydrochloride = "PIS_Ranitidine hydrochloride")


## 0.2 Split into training and testing 80:20 ----
set.seed(1234)

### 20% for testing
df_testing_id <- sample(df$EAVE_LINKNO, size = round(nrow(df)/5, 0), replace = F) 
df_testing <- df %>% filter(EAVE_LINKNO %in% df_testing_id)

### 80% for training
df <- df %>% filter(!(EAVE_LINKNO %in% df_testing_id))


## 0.3 Spline transformations ----
## Spline transform age and bmi variables (to avoid problems with spaces in ns(splines) variable names when using caret::predict)
### Specify variables to be spline-transformed
spline_vars <- c("age", "bmi_imp")

### Training: Apply spline transformation to the specified variables
for (var in spline_vars) {
  var_i <- which(colnames(df)==var) # get column number for indexing
  # 3 knots for age
  if(var == "age"){ df <- cbind(df, ns(df[, var_i], df = 3)) # knots at 34, 52
  start <- ncol(df)-2
  }
  # 2 knots for bmi
  if(var == "bmi_imp"){ df <- cbind(df, ns(df[, var_i], df = 2)) # knot at 28
  start <- ncol(df)-1} 
  # Rename new columns
  end <- ncol(df)
  colnames(df)[start:end] <- paste0(var, colnames(df)[start:end])
}


### Testing: Apply spline transformation to the specified variables
for (var in spline_vars) {
  var_i <- which(colnames(df_testing)==var) # get column number for indexing
  # 3 knots for age
  if(var == "age"){ df_testing <- cbind(df_testing, ns(df_testing[, var_i], df = 3))
  start <- ncol(df_testing)-2
  }
  # 2 knots for bmi
  if(var == "bmi_imp"){ df_testing <- cbind(df_testing, ns(df_testing[, var_i], df = 2))
  start <- ncol(df_testing)-1} 
  # Rename new columns
  end <- ncol(df_testing)
  colnames(df_testing)[start:end] <- paste0(var, colnames(df_testing)[start:end])
}


## 0.4 Weights ----
df$class_weights <- 1

## 0.5 Functions ----
## Function to run full model using different datasets 
run_full_mod <- function(data, # data = e.g. df, df[df$var_per == "Wild",] etc.
                         weight, # e.g. 1 or df$class_weight (for cost function customisation)
                         family){ # family = "binomial" or family = "quasibinomial" (if using weights)
  
  mod <- glm(depvar
             
             # Socio-demographic
             ~ relevel(sex, ref = "M")  # Force male to be reference category
             + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
             + relevel(simd, ref = 5) # Force least deprived to be reference category
             + hhold_bands
             
             # Geographic
             + relevel(ur, ref = 1) # Force remote rural to be reference category 
             # ur causes an error when using vif because coefficient on urUnknown = NA (not defined because of singularities)
             # + relevel(hb, ref = "NHS Lothian") # Force Lothian to be reference category -- exclude hb so model is more widely applicable
             
             # Testing
             #+ tes_n
             + relevel(var_per, ref = "Wild") # force Wild to be the reference category
             
             # Vaccinations
             + vac_dos
             #+ vac_elaps
             
             # Risk factors
             + shielding
             + immuno_supp
             + carehome 
             + bmi_imp1 + bmi_imp2 # splines in bmi 
             + Q_DIAG_AF 
             + Q_DIAG_ASTHMA 
             + Q_DIAG_BLOOD_CANCER 
             + Q_DIAG_CCF  
             + Q_DIAG_CHD  
             + Q_DIAG_CKD_LEVEL  
             + Q_DIAG_COPD 
             + Q_DIAG_DEMENTIA 
             + Q_DIAG_DIABETES_1 
             + Q_DIAG_DIABETES_2 
             + Q_DIAG_EPILEPSY 
             + Q_DIAG_FRACTURE 
             + Q_DIAG_NEURO  
             + Q_DIAG_PARKINSONS 
             + Q_DIAG_PULM_HYPER 
             + Q_DIAG_PULM_RARE  
             + Q_DIAG_PVD  
             + Q_DIAG_RA_SLE 
             + Q_DIAG_RESP_CANCER  
             + Q_DIAG_SEV_MENT_ILL 
             + Q_DIAG_STROKE 
             + Q_DIAG_VTE  
             
             # Prescriptions
             + PIS_Lipid_regulating_drugs
             + PIS_Angiotensin_converting_enzyme_inhibitors
             + PIS_Beta_adrenoceptor_blocking_drugs
             + PIS_Oral_iron
             + PIS_Selective_serotonin_re_uptake_inhibitors
             #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
             + PIS_Loratadine
             + PIS_Direct_oral_anticoagulants
             + PIS_Colchicine
             + PIS_Antiplatelet_drugs
             + PIS_Alpha_adrenoceptor_blocking_drugs
             + PIS_Macrolides
             + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
             + PIS_Leukotriene_receptor_antagonists
             + PIS_Herpes_simplex_and_varicella_zoster
             + PIS_Replacement_therapy
             + PIS_Famotidine
             + PIS_Ursodeoxycholic_acid
             + PIS_Systemic_nasal_decongestants
             + PIS_Warfarin_sodium
             + PIS_Parenteral_anticoagulants
             + PIS_Compound_bronchodilator_preparations
             + PIS_Ranitidine_hydrochloride
             
             # Severity of acute infection
             + severe, 
             
             # Model features 
             data = data,
             family = family, 
             weights = weight)
  
  return(mod)
}


coef_plot <- function(mod, # mod = the name of the model e.g. full_mod
                      plot_name, # the name to save the plot under e.g. "Full plot"
                      type = "standard", # type of model e.g. "standard" or "cv" (cross-validated)
                      add_to_x = 10){ # the amount of space to add to the x axis so that the labels fit
  # Plot
  ## Extract model coefficients, odds ratios, and confidence intervals (excluding NAs)
  var = names(coef(mod)[!is.na(coef(mod))])
  coef = exp(coef(mod)[!is.na(coef(mod))])
  ci = data.frame(exp(confint.default(mod)))
  ci <- ci %>% dplyr::filter(!is.na(X2.5..))
  pval = summary(mod)$coefficients[,4]
  
  ## Format pvals
  pval <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", format(round(pval, digits=3), nsmall = 2)))
  
  coef_df <- data.frame(var, coef, ci, pval)
  
  ## Remove intercept
  coef_df <- coef_df[-1,]
  
  ## Remove 'no dominant variant'
  coef_df <- coef_df[!(rownames(coef_df) == "`relevel(var_per, ref = \"Wild\")Unknown`"),]
  
  if(type == "standard"){
    ## Label categories
    coef_df <- coef_df %>% 
      mutate(subtitle = case_when(
        var == "relevel(sex, ref = \"M\")M"  ~	"Sex",
        var == "relevel(sex, ref = \"M\")F" ~	"Sex",
        var == "age1" ~ "Age",
        var == "age2" ~ "Age",
        var == "age3" ~ "Age",
        var == "relevel(simd, ref = 5)1" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "relevel(simd, ref = 5)2" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "relevel(simd, ref = 5)3" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "relevel(simd, ref = 5)4" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "relevel(simd, ref = 5)5" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "simdUnknown" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "hhold_bands.L" ~ "Household size (Ref: 1)",
        var == "hhold_bands.Q" ~ "Household size (Ref: 1)",
        var == "hhold_bands.C" ~ "Household size (Ref: 1)",
        var == "hhold_bands^4" ~ "Household size (Ref: 1)",
        var == "relevel(ur, ref = 1)2" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "relevel(ur, ref = 1)3" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "relevel(ur, ref = 1)4" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "relevel(ur, ref = 1)5" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "relevel(ur, ref = 1)6" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "relevel(ur, ref = 1)Unknown" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Ayrshire and Arran" ~	"Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Borders" ~	"Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Dumfries and Galloway" ~	     "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Fife" ~	"Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Forth Valley" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Grampian" ~	"Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Greater Glasgow and Clyde" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Highland" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Lanarkshire" ~	"Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Lothian" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Orkney" ~	"Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Shetland" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Tayside" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Western Isles" ~ "Health board (Ref: NHS Lothian)",
        var == "relevel(hb, ref = \"NHS Lothian\")Unknown" ~	"Health board (Ref: NHS Lothian)",
        var == "tes_n.L" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "tes_n.Q" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "tes_n.C" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "\`tes_n^4\`" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "relevel(var_per, ref = \"Wild\")Alpha"   ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "relevel(var_per, ref = \"Wild\")Delta"   ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "relevel(var_per, ref = \"Wild\")Omicron" ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "relevel(var_per, ref = \"Wild\")Unknown" ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "vac_dos1" ~ "Vaccine doses by 14 days before positive PCR test (Ref: 0)",
        var == "vac_dos2" ~ "Vaccine doses by 14 days before positive PCR test (Ref: 0)",
        var == "vac_dos3+" ~"Vaccine doses by 14 days before positive PCR test (Ref: 0)",
        var == "shielding" ~ "Risk factors",
        var == "immuno_supp" ~ "Risk factors",
        var == "carehome" ~ "Carehome resident",
        var == "bmi_imp1" ~ "Risk factors",
        var == "bmi_imp2" ~ "Risk factors",
        grepl("Q_", var) ~ "Risk factors",
        grepl("PIS_", var) ~ "Prescriptions dispensed in the 3 months prior to positive PCR test",
        var == "severe" ~ "Severity of acute infection"))
    
    ## Label variables
    coef_df <- coef_df %>% 
      mutate(label = case_when(
        var == "relevel(sex, ref = \"M\")F" ~	"Female",
        var == "relevel(sex, ref = \"M\")M" ~	"Male",
        var == "age1" ~ "Age 18 - 33 (spline 1)",
        var == "age2" ~ "Age 34 - 51 (spline 2)",
        var == "age3" ~ "Age 52+ (spline 3)",
        var == "relevel(simd, ref = 5)1" ~	"SIMD: Quintile 1 (most deprived)",
        var == "relevel(simd, ref = 5)2" ~	"SIMD: Quintile 2",
        var == "relevel(simd, ref = 5)3" ~ "SIMD: Quintile 3",
        var == "relevel(simd, ref = 5)4" ~ "SIMD: Quintile 4",
        var == "relevel(simd, ref = 5)5" ~ "SIMD: Quintile 5 (least deprived)",
        var == "simdUnknown" ~ "Unknown SIMD",
        var == "hhold_bands.L" ~ "Household size: 2",
        var == "hhold_bands.Q" ~ "Household size: 3-5",
        var == "hhold_bands.C" ~ "Household size: 6-10",
        var == "hhold_bands^4" ~ "Household size: >10",
        var == "relevel(ur, ref = 1)2" ~ "Other Urban Areas",
        var == "relevel(ur, ref = 1)3" ~ "Accessible Small Towns",
        var == "relevel(ur, ref = 1)4" ~ "Remote Small Towns",
        var == "relevel(ur, ref = 1)5" ~ "Accessible Rural",
        var == "relevel(ur, ref = 1)6" ~ "Remote Rural",
        var == "relevel(ur, ref = 1)Unknown" ~ "Unknown",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Ayrshire and Arran" ~	"NHS Ayrshire and Arran",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Borders" ~	"NHS Borders",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Dumfries and Galloway" ~	"NHS Dumfries and Galloway",      
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Fife" ~	"NHS Fife",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Forth Valley" ~ "NHS Forth Valley",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Grampian" ~	"NHS Grampian",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Greater Glasgow and Clyde" ~ "NHS Greater Glasgow and Clyde",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Highland" ~ "NHS Highland",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Lanarkshire" ~	"NHS Lanarkshire",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Lothian" ~ "NHS Lothian",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Orkney" ~	"NHS Orkney",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Shetland" ~ "NHS Shetland",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Tayside" ~ "NHS Tayside",
        var == "relevel(hb, ref = \"NHS Lothian\")NHS Western Isles" ~ "NHS Western Isles",
        var == "relevel(hb, ref = \"NHS Lothian\")Unknown" ~	"Unknown healthboard",
        var == "tes_n.L" ~ "2 PCR tests",
        var == "tes_n.Q" ~ "3-4 PCR tests",
        var == "tes_n.C" ~ "5-10 PCR tests",
        var == "tes_n^4" ~ ">10 PCR tests",
        var == "relevel(var_per, ref = \"Wild\")Alpha"   ~ "Alpha (11/01/2021 - 09/05/2021)",
        var == "relevel(var_per, ref = \"Wild\")Delta"   ~ "Delta (24/05/2021 - 05/12/2021)",
        var == "relevel(var_per, ref = \"Wild\")Omicron" ~ "Omicron (27/12/2021 - 30/04/2022)",
        var == "relevel(var_per, ref = \"Wild\")Unknown" ~ "No dominant variant",
        var == "vac_dos1" ~ "Vaccine doses: 1",
        var == "vac_dos2" ~ "Vaccine doses: 2",
        var == "vac_dos3+" ~ "Vaccine doses: 3+",
        var == "shielding" ~ "Shielding",
        var == "immuno_supp" ~ "Immunosuppressed",
        var == "carehome" ~ "Carehome resident",
        var == "bmi_imp1" ~ "Body Mass Index < 28 (spline 1)",
        var == "bmi_imp2" ~ "Body Mass Index 28+ (spline 2)",
        var == "Q_DIAG_AF" ~ "Atrial fibrillation",
        var == "Q_DIAG_ASTHMA" ~ "Asthma",
        var == "Q_DIAG_BLOOD_CANCER" ~ "Blood cancer",
        var == "Q_DIAG_CCF" ~ "Heart failure",
        var == "Q_DIAG_CHD" ~ "Coronary heart disease",
        var == "Q_DIAG_CKD_LEVEL" ~ "Kidney disease (level 3+)",
        var == "Q_DIAG_COPD" ~ "Chronic obstructive pulmonary disease (COPD)",
        var == "Q_DIAG_DEMENTIA" ~ "Dementia",
        var == "Q_DIAG_DIABETES_1" ~ "Diabetes Type I",
        var == "Q_DIAG_DIABETES_2" ~ "Diabetes Type II",
        var == "Q_DIAG_EPILEPSY" ~ "Epilepsy",
        var == "Q_DIAG_FRACTURE" ~ "Fracture",
        var == "Q_DIAG_NEURO" ~ "Neurological disorder",
        var == "Q_DIAG_PARKINSONS" ~ "Parkinsons",
        var == "Q_DIAG_PULM_HYPER" ~ "Pulmonary hypertension",
        var == "Q_DIAG_PULM_RARE" ~ "Pulmonary rare",
        var == "Q_DIAG_PVD" ~ "Peripheral vascular disease",
        var == "Q_DIAG_RA_SLE" ~ "Rheumatoid arthritis or SLE",
        var == "Q_DIAG_RESP_CANCER" ~ "Respiratory cancer",
        var == "Q_DIAG_SEV_MENT_ILL" ~ "Severe mental illness",
        var == "Q_DIAG_STROKE" ~ "Stroke/TIA",
        var == "Q_DIAG_VTE" ~ "A thrombosis or pulmonary embolus",
        var == "PIS_Lipid_regulating_drugs" ~ "Lipid-regulating drugs",
        var == "PIS_Angiotensin_converting_enzyme_inhibitors" ~ "Angiotensin-converting enzyme inhibitors",
        var == "PIS_Beta_adrenoceptor_blocking_drugs" ~ "Beta-adrenoceptor blocking drugs",
        var == "PIS_Oral_iron" ~ "Oral iron",
        var == "PIS_Selective_serotonin_re_uptake_inhibitors" ~ "Selective serotonin re-uptake inhibitors",
        var == "PIS_Loratadine" ~ "Loratadine (antihistamine)",
        var == "PIS_Direct_oral_anticoagulants" ~ "Direct oral anticoagulants",
        var == "PIS_Colchicine" ~ "Colchicine (anti-inflammatory)",
        var == "PIS_Antiplatelet_drugs" ~ "Antiplatelet drugs",
        var == "PIS_Alpha_adrenoceptor_blocking_drugs" ~ "Alpha-adrenoceptor blocking drugs",
        var == "PIS_Macrolides" ~ "Macrolides (antibacterial)",
        var == "PIS_Benzylpenicillin_and_phenoxymethylpenicillin" ~ "Benzylpenicillin and phenoxymethylpenicillin",
        var == "PIS_Leukotriene_receptor_antagonists" ~ "Leukotriene receptor antagonists",
        var == "PIS_Herpes_simplex_and_varicella_zoster" ~ "Herpes simplex and varicella-zoster (antiviral)",
        var == "PIS_Replacement_therapy" ~ "Corticosteroid replacement therapy",
        var == "PIS_Famotidine" ~ "Famotidine (histamine H2 receptor antagonist)",
        var == "PIS_Ursodeoxycholic_acid" ~ "Ursodeoxycholic acid",
        var == "PIS_Systemic_nasal_decongestants" ~ "Systemic nasal decongestants",
        var == "PIS_Warfarin_sodium" ~ "Warfarin sodium",
        var == "PIS_Parenteral_anticoagulants" ~ "Parenteral anticoagulants",
        var == "PIS_Compound_bronchodilator_preparations" ~ "Compound bronchodilator preparations",
        var == "PIS_Ranitidine_hydrochloride" ~ "Ranitidine hydrochloride",
        var == "severe" ~ "Hospitalised within 28 days of PCR test"))
  }
  
  if(type == "cv"){
    ## Label categories
    coef_df <- coef_df %>% 
      mutate(subtitle = case_when(
        var == "\`relevel(sex, ref = \"M\")M\`"  ~	"Sex",
        var == "\`relevel(sex, ref = \"M\")F\`" ~	"Sex",
        var == "age1" ~ "Age",
        var == "age2" ~ "Age",
        var == "age3" ~ "Age",
        var == "\`relevel(simd, ref = 5)1\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "\`relevel(simd, ref = 5)2\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "\`relevel(simd, ref = 5)3\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "\`relevel(simd, ref = 5)4\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "\`relevel(simd, ref = 5)5\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "simdUnknown" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
        var == "hhold_bands.L" ~ "Household size (Ref: 1)",
        var == "hhold_bands.Q" ~ "Household size (Ref: 1)",
        var == "hhold_bands.C" ~ "Household size (Ref: 1)",
        var == "\`hhold_bands^4\`" ~ "Household size (Ref: 1)",
        var == "\`relevel(ur, ref = 1)2\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "\`relevel(ur, ref = 1)3\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "\`relevel(ur, ref = 1)4\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "\`relevel(ur, ref = 1)5\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "\`relevel(ur, ref = 1)6\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "\`relevel(ur, ref = 1)Unknown\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Ayrshire and Arran\`" ~	"Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Borders\`" ~	"Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Dumfries and Galloway\`" ~	     "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Fife\`" ~	"Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Forth Valley\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Grampian\`" ~	"Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Greater Glasgow and Clyde\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Highland\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lanarkshire\`" ~	"Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lothian\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Orkney\`" ~	"Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Shetland\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Tayside\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Western Isles\`" ~ "Health board (Ref: NHS Lothian)",
        var == "\`relevel(hb, ref = \"NHS Lothian\")Unknown\`" ~	"Health board (Ref: NHS Lothian)",
        var == "tes_n.L" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "tes_n.Q" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "tes_n.C" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "\`tes_n^4\`" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
        var == "\`relevel(var_per, ref = \"Wild\")Alpha\`"   ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "\`relevel(var_per, ref = \"Wild\")Delta\`"   ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "\`relevel(var_per, ref = \"Wild\")Omicron\`" ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "\`relevel(var_per, ref = \"Wild\")Unknown\`" ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
        var == "vac_dos1" ~ "Vaccine doses by 14 days before positive PCR test (Ref: 0)",
        var == "vac_dos2" ~ "Vaccine doses by 14 days before positive PCR test (Ref: 0)",
        var == "\`vac_dos3+\`" ~"Vaccine doses by 14 days before positive PCR test (Ref: 0)",
        var == "shielding" ~ "Risk factors",
        var == "immuno_supp" ~ "Risk factors",
        var == "carehome" ~ "Carehome resident",
        var == "bmi_imp1" ~ "Risk factors",
        var == "bmi_imp2" ~ "Risk factors",
        grepl("Q_", var) ~ "Risk factors",
        grepl("PIS_", var) ~ "Prescriptions dispensed in the 3 months prior to positive PCR test",
        var == "severe" ~ "Severity of acute infection"))
    
    ## Label variables
    coef_df <- coef_df %>% 
      mutate(label = case_when(
        var == "\`relevel(sex, ref = \"M\")F\`" ~	"Female",
        var == "\`relevel(sex, ref = \"M\")M\`" ~	"Male",
        var == "age1" ~ "Age 18 - 33 (spline 1)",
        var == "age2" ~ "Age 34 - 51 (spline 2)",
        var == "age3" ~ "Age 52+ (spline 3)",
        var == "\`relevel(simd, ref = 5)1\`" ~	"SIMD: Quintile 1 (most deprived)",
        var == "\`relevel(simd, ref = 5)2\`" ~	"SIMD: Quintile 2",
        var == "\`relevel(simd, ref = 5)3\`" ~ "SIMD: Quintile 3",
        var == "\`relevel(simd, ref = 5)4\`" ~ "SIMD: Quintile 4",
        var == "\`relevel(simd, ref = 5)5\`" ~ "SIMD: Quintile 5 (least deprived)",
        var == "simdUnknown" ~ "Unknown SIMD",
        var == "hhold_bands.L" ~ "Household size: 2",
        var == "hhold_bands.Q" ~ "Household size: 3-5",
        var == "hhold_bands.C" ~ "Household size: 6-10",
        var == "\`hhold_bands^4\`" ~ "Household size: >10",
        var == "\`relevel(ur, ref = 1)2\`" ~ "Other Urban Areas",
        var == "\`relevel(ur, ref = 1)3\`" ~ "Accessible Small Towns",
        var == "\`relevel(ur, ref = 1)4\`" ~ "Remote Small Towns",
        var == "\`relevel(ur, ref = 1)5\`" ~ "Accessible Rural",
        var == "\`relevel(ur, ref = 1)6\`" ~ "Remote Rural",
        var == "\`relevel(ur, ref = 1)Unknown\`" ~ "Unknown",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Ayrshire and Arran\`" ~	"NHS Ayrshire and Arran",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Borders\`" ~	"NHS Borders",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Dumfries and Galloway\`" ~	"NHS Dumfries and Galloway",      
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Fife\`" ~	"NHS Fife",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Forth Valley\`" ~ "NHS Forth Valley",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Grampian\`" ~	"NHS Grampian",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Greater Glasgow and Clyde\`" ~ "NHS Greater Glasgow and Clyde",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Highland\`" ~ "NHS Highland",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lanarkshire\`" ~	"NHS Lanarkshire",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lothian\`" ~ "NHS Lothian",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Orkney\`" ~	"NHS Orkney",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Shetland\`" ~ "NHS Shetland",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Tayside\`" ~ "NHS Tayside",
        var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Western Isles\`" ~ "NHS Western Isles",
        var == "\`relevel(hb, ref = \"NHS Lothian\")Unknown\`" ~	"Unknown healthboard",
        var == "tes_n.L" ~ "2 PCR tests",
        var == "tes_n.Q" ~ "3-4 PCR tests",
        var == "tes_n.C" ~ "5-10 PCR tests",
        var == "\`tes_n^4\`" ~ ">10 PCR tests",
        var == "\`relevel(var_per, ref = \"Wild\")Alpha\`"   ~ "Alpha (11/01/2021 - 09/05/2021)",
        var == "\`relevel(var_per, ref = \"Wild\")Delta\`"   ~ "Delta (24/05/2021 - 05/12/2021)",
        var == "\`relevel(var_per, ref = \"Wild\")Omicron\`" ~ "Omicron (27/12/2021 - 30/04/2022)",
        var == "\`relevel(var_per, ref = \"Wild\")Unknown\`" ~ "No dominant variant",
        var == "vac_dos1" ~ "Vaccine doses: 1",
        var == "vac_dos2" ~ "Vaccine doses: 2",
        var == "\`vac_dos3+\`" ~ "Vaccine doses: 3+",
        var == "shielding" ~ "Shielding",
        var == "immuno_supp" ~ "Immunosuppressed",
        var == "carehome" ~ "Carehome resident",
        var == "bmi_imp1" ~ "Body Mass Index < 28 (spline 1)",
        var == "bmi_imp2" ~ "Body Mass Index 28+ (spline 2)",
        var == "Q_DIAG_AF" ~ "Atrial fibrillation",
        var == "Q_DIAG_ASTHMA" ~ "Asthma",
        var == "Q_DIAG_BLOOD_CANCER" ~ "Blood cancer",
        var == "Q_DIAG_CCF" ~ "Heart failure",
        var == "Q_DIAG_CHD" ~ "Coronary heart disease",
        var == "Q_DIAG_CKD_LEVEL" ~ "Kidney disease (level 3+)",
        var == "Q_DIAG_COPD" ~ "Chronic obstructive pulmonary disease (COPD)",
        var == "Q_DIAG_DEMENTIA" ~ "Dementia",
        var == "Q_DIAG_DIABETES_1" ~ "Diabetes Type I",
        var == "Q_DIAG_DIABETES_2" ~ "Diabetes Type II",
        var == "Q_DIAG_EPILEPSY" ~ "Epilepsy",
        var == "Q_DIAG_FRACTURE" ~ "Fracture",
        var == "Q_DIAG_NEURO" ~ "Neurological disorder",
        var == "Q_DIAG_PARKINSONS" ~ "Parkinsons",
        var == "Q_DIAG_PULM_HYPER" ~ "Pulmonary hypertension",
        var == "Q_DIAG_PULM_RARE" ~ "Pulmonary rare",
        var == "Q_DIAG_PVD" ~ "Peripheral vascular disease",
        var == "Q_DIAG_RA_SLE" ~ "Rheumatoid arthritis or SLE",
        var == "Q_DIAG_RESP_CANCER" ~ "Respiratory cancer",
        var == "Q_DIAG_SEV_MENT_ILL" ~ "Severe mental illness",
        var == "Q_DIAG_STROKE" ~ "Stroke/TIA",
        var == "Q_DIAG_VTE" ~ "A thrombosis or pulmonary embolus",
        var == "PIS_Lipid_regulating_drugs" ~ "Lipid-regulating drugs",
        var == "PIS_Angiotensin_converting_enzyme_inhibitors" ~ "Angiotensin-converting enzyme inhibitors",
        var == "PIS_Beta_adrenoceptor_blocking_drugs" ~ "Beta-adrenoceptor blocking drugs",
        var == "PIS_Oral_iron" ~ "Oral iron",
        var == "PIS_Selective_serotonin_re_uptake_inhibitors" ~ "Selective serotonin re-uptake inhibitors",
        var == "PIS_Loratadine" ~ "Loratadine (antihistamine)",
        var == "PIS_Direct_oral_anticoagulants" ~ "Direct oral anticoagulants",
        var == "PIS_Colchicine" ~ "Colchicine (anti-inflammatory)",
        var == "PIS_Antiplatelet_drugs" ~ "Antiplatelet drugs",
        var == "PIS_Alpha_adrenoceptor_blocking_drugs" ~ "Alpha-adrenoceptor blocking drugs",
        var == "PIS_Macrolides" ~ "Macrolides (antibacterial)",
        var == "PIS_Benzylpenicillin_and_phenoxymethylpenicillin" ~ "Benzylpenicillin and phenoxymethylpenicillin",
        var == "PIS_Leukotriene_receptor_antagonists" ~ "Leukotriene receptor antagonists",
        var == "PIS_Herpes_simplex_and_varicella_zoster" ~ "Herpes simplex and varicella-zoster (antiviral)",
        var == "PIS_Replacement_therapy" ~ "Corticosteroid replacement therapy",
        var == "PIS_Famotidine" ~ "Famotidine (histamine H2 receptor antagonist)",
        var == "PIS_Ursodeoxycholic_acid" ~ "Ursodeoxycholic acid",
        var == "PIS_Systemic_nasal_decongestants" ~ "Systemic nasal decongestants",
        var == "PIS_Warfarin_sodium" ~ "Warfarin sodium",
        var == "PIS_Parenteral_anticoagulants" ~ "Parenteral anticoagulants",
        var == "PIS_Compound_bronchodilator_preparations" ~ "Compound bronchodilator preparations",
        var == "PIS_Ranitidine_hydrochloride" ~ "Ranitidine hydrochloride",
        var == "severe" ~ "Hospitalised within 28 days of PCR test"))
  }
  
  # Set column names
  colnames(coef_df) <- c("var", "coef", "lower", "upper", "pval", "subtitle", "label")
  
  # Set order of facets to appear in plot
  coef_df$subtitle <- factor(coef_df$subtitle, levels=unique(coef_df$subtitle))
  subtitle <- factor(coef_df$subtitle, levels=unique(coef_df$subtitle)) # Needed for adjusting height of each facet
  
  # Set order variables to appear in plot (with Q_ and PIS_ vars ordered by coefficient size)
  ## Subset start of dataframe
  start <- coef_df[!(grepl("Q_|PIS_|hb", coef_df$var)),]
  
  ## Subset of dataframe where 'var' begins with "hb"
  #hb_subset <- coef_df[grep("hb", coef_df$var), ]
  #hb_subset <- hb_subset[order(hb_subset$coef), ]
  
  ## Subset of dataframe where 'var' begins with "Q_"
  Q_subset <- coef_df[grep("Q_", coef_df$var), ]
  Q_subset <- Q_subset[order(Q_subset$coef), ]
  
  ## Subset of dataframe where 'var' begins with "PIS_"
  PIS_subset <- coef_df[grep("PIS_", coef_df$var), ]
  PIS_subset <- PIS_subset[order(PIS_subset$coef), ]
  
  ## Combine the subsets and the remaining rows of original dataframe
  coef_df <- rbind(start, Q_subset, PIS_subset)
  
  ## Set order
  coef_df$label <- factor(coef_df$label, levels=rev(unique(coef_df$label)))
  
  # Create a forest plot using ggplot2
  # Calculate the expansion factor for x-axis
  expand_factor <- max(coef_df$upper) * 1.5 ## Edit this value if labels don't show on plot
  
  p <- ggplot(coef_df, aes(x = coef, y = label)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10(limits = c(min(coef_df$lower) - 0.05, max(coef_df$upper) + expand_factor),
                  expand = c(0, 0)) +  # Set the expand parameter to remove white space
    geom_text(aes(x = max(coef_df$upper),
                  label = sprintf("%.2f (%.2f, %.2f); %s", coef, lower, upper, pval)),
              hjust = 0,  # Adjust the value here to position the label correctly
              size = 3, position = position_nudge(x = 0.05)) +
    theme_bw() +
    ylab("") +
    xlab("Adjusted Odds Ratio (95% CI)") +
    theme(axis.text.y = element_text(size = 11, color = "black"))  # Set the font size for y-axis labels
  
  
  # Make it faceted (facet_grid sets heights for different numbers of variables in different facets)
  p.grid <- p + facet_grid(subtitle ~ ., scales = "free_y", space = "free_y")
  p.wrap <- p + facet_wrap(~ subtitle, ncol = 1, scales = "free_y")
  
  gp.grid <- ggplotGrob(p.grid)
  gp.wrap <- ggplotGrob(p.wrap)
  
  gp.wrap$heights[gp.wrap$layout[grep("panel", gp.wrap$layout$name), "t"]] <- 
    gp.grid$heights[gp.grid$layout[grep("panel", gp.grid$layout$name), "t"]]
  
  grid::grid.draw(gp.wrap)
  
  setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")
  ggsave(paste0(plot_name, " - ", depvar_name, ".jpg"), gp.wrap, height = 14, width = 11, units = "in")
  
  return(gp.wrap)
}


# 1. Multicollinnearity checks----
## 1.1 Correlations ----
## https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html 
## Select binarised versions of variables to check for correlations within each level of a variable
cor_df <- df %>% 
  dplyr::select(
    # Socio-demographics
    age, sex_m, simd1, simd2, simd3, simd4, simd5, hhold1, hhold2, hhold3, hhold4, hhold5, 
    # Geographic
    ur1, ur2, ur3, ur4, ur5, ur6, 
    #colnames(df[startsWith(colnames(df), "hb_")]), -hb_lab, -hb,
    # Testing
    #tes1, tes2, tes3, tes4, tes5, 
    wild, alpha, delta, omicron,
    # Vaccinations
    vac1, vac2, vac3, vac4, vac5, vac6,
    #vac_elaps1, vac_elaps2, vac_elaps3, vac_elaps4, vac_elaps5, 
    # Risk factors
    shielding, immuno_supp, carehome, bmi_imp,
    colnames(df[startsWith(colnames(df), "Q_")]), 
    # Severity of acute infection
    severe,
    # Prescriptions
    colnames(df[startsWith(colnames(df), "PIS_")]))

## Get correlations
cors <- cor(cor_df, method = "p") 

## Plot correlation matrix
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling/")

filetag <- "Correlations.pdf"

pdf(filetag, height = 10, width = 10)

corrplot(cors, method = 'circle', type = 'lower', diag = F, tl.cex = 0.6, tl.col = "black", mar = c(0,0,0,0))

dev.off()

## Inspect larger correlations (<-0.7, >0.7 generally considered large) 
up_lim <- 0.5
low_lim <- -0.5

cors[lower.tri(cors)]=NA
diag(cors)=NA
cors_df <- na.omit(melt(cors))
cors_df <- cors_df %>% 
  filter(value > 0.5 | value < -0.5)
cors_df

### Q_DIAG_DIABETES_2 x PIS_Metformin_hydrochloride (diabetes medication)
### Create combined diabetes variable
df$Q_DIAG_DIABETES_2 <- ifelse(df$PIS_Metformin_hydrochloride == 1, 1, df$Q_DIAG_DIABETES_2 )
df_testing$Q_DIAG_DIABETES_2 <- ifelse(df_testing$PIS_Metformin_hydrochloride == 1, 1, df_testing$Q_DIAG_DIABETES_2 )

## Clean up
rm(cors_df, cor_df, cors)


## 1.2 VIF: get Generalised Variance Inflation Factors for each predictor ----
## https://cran.r-project.org/web/packages/olsrr/vignettes/regression_diagnostics.html#:~:text=Collinearity%20is%20spotted%20by%20finding,range%20of%2030%20or%20larger.
## VIFs measure the inflation in the variances/elipsoid size of the parameter estimates due to collinearities that exist among the predictors. 
mod_vif <- run_full_mod(data = df[df$ur != 'Unknown',], weight = df$class_weights, family = "quasibinomial") 
summary(mod_vif) 

## GVIF(1/(2×Df) < 2 is not a problem
## 1. https://stats.stackexchange.com/questions/430412/vif-for-categorical-variable-with-more-than-2-categories?rq=1 
## 2. https://stats.stackexchange.com/questions/70679/which-variance-inflation-factor-should-i-be-using-textgvif-or-textgvif#:~:text=For%20the%20two%20continuous%20variables%2C%20GVIF(1,to%20the%20level%20of%20collinearity.
VIF <- vif(mod_vif)  
VIF

## Clean up
rm(mod_vif, VIF)


# 2. Model ----
## 2.1 Run model  ----
### Run with 'family = "binomial" to get AIC (needed for calculation of quasi AIC)
full_mod <- run_full_mod(data = df, weight = df$class_weight, family = "binomial") # NB df$class_weight = inverse class proportion weights
summary(full_mod)

### Run with 'family = "quasibinomial" to adjust for overdispersion 
full_mod_q <- run_full_mod(data = df, weight = df$class_weight, family = "quasibinomial") # NB df$class_weight = inverse class proportion weights
summary(full_mod_q)

### Plot
coef_plot(mod = full_mod_q, plot_name = "Coefficient plot - full model")

### Get quasi AIC
full_overdispersion <- deviance(full_mod_q) / full_mod_q$df.residual
full_quasi_AIC <- AIC(full_mod) + 2 * full_overdispersion
print(full_quasi_AIC) 

## 2.2 R2 ----
## Cox-Snell R2 (use Cox Snell given binary outcome)
full_mod_R2_cs <- PseudoR2(full_mod, which = "CoxSnell") 
full_mod_R2_cs 

### Nagelkerke's R2: a modification of the Cox-Snell R-squared that is anchored to a maximum value of 1. 
### More conservative than Cox-Snell and may be more appropriate for assessing model performance.
full_mod_R2_n <- NagelkerkeR2(full_mod)
full_mod_R2_n 


# 3. Predictor selection ----
# Aikaike/Bayseian Information Criterion (AIC/BIC) balances the fit of a model with the number of predictors
# AIC penalizes for more predictors
# BIC penalizes for larger sample sizes (and may lead to underfitting in very large samples)
# https://rforpoliticalscience.com/2020/10/23/choose-model-variables-by-aic-in-a-stepwise-algorithm-with-the-mass-package-in-r/

### 3.1 AIC-based selection ----
### Need to re-run full_mod without using function - otherwise stepAIC won't work
full_mod <- glm(depvar
                
                # Socio-demographic
                ~ relevel(sex, ref = "M")  # Force male to be reference category
                + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
                + relevel(simd, ref = 5) # Force least deprived to be reference category
                + hhold_bands
                
                # Geographic
                + relevel(ur, ref = 1) # Force remote rural to be reference category 
                # ur causes an error when using vif because coefficient on urUnknown = NA (not defined because of singularities)
                # + relevel(hb, ref = "NHS Lothian") # Force Lothian to be reference category -- exclude hb so model is more widely applicable
                
                # Testing
                #+ tes_n # not using due to removal of mass testing
                + relevel(var_per, ref = "Wild") # force Wild to be the reference category
                
                # Vaccinations
                + vac_dos
                #+ vac_elaps
                
                # Risk factors
                + shielding
                + immuno_supp
                + carehome 
                + bmi_imp1 + bmi_imp2 # splines in bmi 
                + Q_DIAG_AF 
                + Q_DIAG_ASTHMA 
                + Q_DIAG_BLOOD_CANCER 
                + Q_DIAG_CCF  
                + Q_DIAG_CHD  
                + Q_DIAG_CKD_LEVEL  
                + Q_DIAG_COPD 
                + Q_DIAG_DEMENTIA 
                + Q_DIAG_DIABETES_1 
                + Q_DIAG_DIABETES_2 
                + Q_DIAG_EPILEPSY 
                + Q_DIAG_FRACTURE 
                + Q_DIAG_NEURO  
                + Q_DIAG_PARKINSONS 
                + Q_DIAG_PULM_HYPER 
                + Q_DIAG_PULM_RARE  
                + Q_DIAG_PVD  
                + Q_DIAG_RA_SLE 
                + Q_DIAG_RESP_CANCER  
                + Q_DIAG_SEV_MENT_ILL 
                + Q_DIAG_STROKE 
                + Q_DIAG_VTE  
                
                # Prescriptions
                + PIS_Lipid_regulating_drugs
                + PIS_Angiotensin_converting_enzyme_inhibitors
                + PIS_Beta_adrenoceptor_blocking_drugs
                + PIS_Oral_iron
                + PIS_Selective_serotonin_re_uptake_inhibitors
                #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
                + PIS_Loratadine
                + PIS_Direct_oral_anticoagulants
                + PIS_Colchicine
                + PIS_Antiplatelet_drugs
                + PIS_Alpha_adrenoceptor_blocking_drugs
                + PIS_Macrolides
                + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
                + PIS_Leukotriene_receptor_antagonists
                + PIS_Herpes_simplex_and_varicella_zoster
                + PIS_Replacement_therapy
                + PIS_Famotidine
                + PIS_Ursodeoxycholic_acid
                + PIS_Systemic_nasal_decongestants
                + PIS_Warfarin_sodium
                + PIS_Parenteral_anticoagulants
                + PIS_Compound_bronchodilator_preparations
                + PIS_Ranitidine_hydrochloride
                
                # Severity of acute infection
                + severe, 
                
                # Model features 
                data = df,
                family = "binomial", 
                weights = df$class_weight)

## PLEASE NOTE
## Run time is v long (48-60 hours) --> I've commented out the stepwise selection code and 
## specified the optimal model that's identified by the commented out code 

### Use AIC to select the version of the model that has the smallest amount of error by varying predictor selection
### sink() saves the outputs to a text file as the code runs - to avoid losing progress if the R session crashes
#sink(file = "stepwise_aic_output.txt")
#stepwise_aic <- stepAIC(full_mod, direction="backward", trace=T) # "backward" = start with all variables and remove one by one (can also use "both") 
#summary(stepwise_aic)
#sink(file = NULL)


# Optimal model identified by running the commented out code above (the code above takes 48-60 hours to run)
run_aic_mod <- function(data,
                        family, # e.g. family = "binomial" or family = "quasibinomial"
                        weight){ # e.g. weight = 1 or weight = df$class_weights
  
  mod <- glm(depvar
             
             # Socio-demographic
             ~ relevel(sex, ref = "M") # Force male to be reference category
             + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
             + relevel(simd, ref = 5) # Force least deprived to be reference category
             + hhold_bands
             
             # Geographic
             + relevel(ur, ref = 1) # Force remote rural to be reference category
             #+ relevel(hb, ref = "NHS Lothian") # Force Lothian to be reference category
             
             # Testing
             #+ tes_n # not using
             + relevel(var_per, ref = "Wild") # force Wild to be the reference category
             
             # Vaccinations
             + vac_dos
             
             # Risk factors
             + shielding
             + immuno_supp
             + carehome 
             + bmi_imp1 + bmi_imp2 # splines in bmi 
             + Q_DIAG_AF 
             + Q_DIAG_ASTHMA 
             + Q_DIAG_BLOOD_CANCER 
             ### + Q_DIAG_CCF 
             + Q_DIAG_CHD  
             + Q_DIAG_CKD_LEVEL  
             + Q_DIAG_COPD 
             + Q_DIAG_DEMENTIA 
             + Q_DIAG_DIABETES_1 
             + Q_DIAG_DIABETES_2 
             ###+ Q_DIAG_EPILEPSY  
             ###+ Q_DIAG_FRACTURE 
             + Q_DIAG_NEURO   
             ###+ Q_DIAG_PARKINSONS 
             ### + Q_DIAG_PULM_HYPER  
             + Q_DIAG_PULM_RARE  
             ###+ Q_DIAG_PVD   
             + Q_DIAG_RA_SLE 
             ###+ Q_DIAG_RESP_CANCER   
             + Q_DIAG_SEV_MENT_ILL 
             ###+ Q_DIAG_STROKE  
             + Q_DIAG_VTE  
             
             # Prescriptions
             + PIS_Lipid_regulating_drugs
             + PIS_Angiotensin_converting_enzyme_inhibitors
             + PIS_Beta_adrenoceptor_blocking_drugs
             + PIS_Oral_iron
             + PIS_Selective_serotonin_re_uptake_inhibitors
             #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
             + PIS_Loratadine
             ###+ PIS_Direct_oral_anticoagulants 
             + PIS_Colchicine
             + PIS_Antiplatelet_drugs
             ###+ PIS_Alpha_adrenoceptor_blocking_drugs
             + PIS_Macrolides
             + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
             + PIS_Leukotriene_receptor_antagonists
             + PIS_Herpes_simplex_and_varicella_zoster
             ###+ PIS_Replacement_therapy 
             + PIS_Famotidine
             + PIS_Ursodeoxycholic_acid
             + PIS_Systemic_nasal_decongestants
             ###+ PIS_Warfarin_sodium 
             + PIS_Parenteral_anticoagulants 
             + PIS_Compound_bronchodilator_preparations
             ###+ PIS_Ranitidine_hydrochloride 
             
             # Severity of acute infection
             + severe, 
             
             # Model features 
             data = data,
             family = family,
             weights = weight)
  
  return(mod)
}

### Run with 'family = "binomial" to get AIC
stepwise_aic <- run_aic_mod(data = df, family = "binomial", weight = df$class_weights)
summary(stepwise_aic) 

### Plot
AIC_plot <- coef_plot(mod = stepwise_aic, plot_name = "Coefficient plot - stepwise AIC model", add_to_x = 11)


### 3.2 BIC-based selection ----
## Use BIC to select the version of the model that has the smallest amount of error by varying predictor selection
#sink(file = "stepwise_bic_output.txt")
#stepwise_bic <- stepAIC(full_mod, direction="backward", trace=T, k=log(nrow(df))) # "backward" = start with all variables and remove one by one (can also use "both") 
#summary(stepwise_bic)
#sink(file = NULL)

run_bic_mod <- function(data,
                        family, # e.g. family = "binomial" or family = "quasibinomial"
                        weight){ # e.g. weight = 1 or weight = df$class_weights
  
  mod <- glm(depvar
             
             # Socio-demographic
             ~ relevel(sex, ref = "M") # Force male to be reference category
             + age1 + age2 #+ age3 # age with natural splines (3 degrees of freedom)
             + relevel(simd, ref = 5) # Force least deprived to be reference category
             ###+ hhold_bands 
             
             # Geographic
             + relevel(ur, ref = 1) # Force remote rural to be reference category 
             + relevel(hb, ref = "NHS Lothian") # Force Lothian to be reference category
             
             # Testing
             #+ tes_n # do not include due to end of testing + is not important to model
             + relevel(var_per, ref = "Wild") # force Wild to be the reference category
             
             # Vaccinations
             + vac_dos 
             
             # Risk factors
             + shielding
             + immuno_supp
             ###+ carehome 
             + bmi_imp1 + bmi_imp2 # splines in bmi 
             ###+ Q_DIAG_AF 
             + Q_DIAG_ASTHMA 
             ###+ Q_DIAG_BLOOD_CANCER 
             ###+ Q_DIAG_CCF  
             + Q_DIAG_CHD  
             ###+ Q_DIAG_CKD_LEVEL 
             + Q_DIAG_COPD 
             + Q_DIAG_DEMENTIA 
             + Q_DIAG_DIABETES_1 
             + Q_DIAG_DIABETES_2 
             ###+ Q_DIAG_EPILEPSY 
             ###+ Q_DIAG_FRACTURE 
             ###+ Q_DIAG_NEURO  
             ###+ Q_DIAG_PARKINSONS 
             ###+ Q_DIAG_PULM_HYPER 
             ###+ Q_DIAG_PULM_RARE  
             ###+ Q_DIAG_PVD  
             ###+ Q_DIAG_RA_SLE 
             ###+ Q_DIAG_RESP_CANCER  
             + Q_DIAG_SEV_MENT_ILL 
             ###+ Q_DIAG_STROKE 
             ###+ Q_DIAG_VTE  
             
             # Prescriptions
             + PIS_Lipid_regulating_drugs
             + PIS_Angiotensin_converting_enzyme_inhibitors
             + PIS_Beta_adrenoceptor_blocking_drugs
             + PIS_Oral_iron
             + PIS_Selective_serotonin_re_uptake_inhibitors
             #+ PIS_Metformin_hydrochloride # combined with T2 diabetes
             + PIS_Loratadine
             ###+ PIS_Direct_oral_anticoagulants
             ###+ PIS_Colchicine 
             ###+ PIS_Antiplatelet_drugs 
             ###+ PIS_Alpha_adrenoceptor_blocking_drugs 
             + PIS_Macrolides
             + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
             + PIS_Leukotriene_receptor_antagonists
             + PIS_Herpes_simplex_and_varicella_zoster
             ###+ PIS_Replacement_therapy 
             ###+ PIS_Famotidine 
             ###+ PIS_Ursodeoxycholic_acid 
             + PIS_Systemic_nasal_decongestants
             ###+ PIS_Warfarin_sodium
             ###+ PIS_Parenteral_anticoagulants
             + PIS_Compound_bronchodilator_preparations
             ###+ PIS_Ranitidine_hydrochloride
             
             # Severity of acute infection
             + severe, 
             
             # Model features 
             data = data,
             family = family, # because binary dependent variable is weighted (accounts for overdispersion)
             weights = weight)
  
  return(mod)
  
}

### Run with 'family = "binomial" to get AIC
stepwise_bic <- run_bic_mod(data = df, family = "binomial", weight = df$class_weights)
summary(stepwise_bic) 

### Plot
BIC_plot <- coef_plot(mod = stepwise_bic, plot_name = "Coefficient plot - stepwise BIC model")


### 3.3 AIC/BIC R2 ----
### Cox-Snell R2 (use Cox Snell given binary outcome)
aic_R2_cs <- PseudoR2(stepwise_aic, which = "CoxSnell") 
aic_R2_cs  

bic_R2_cs <- PseudoR2(stepwise_bic, which = "CoxSnell") 
bic_R2_cs 


### Nagelkerke's R2: a modification of the Cox-Snell R-squared that is anchored to a maximum value of 1. 
### More conservative than Cox-Snell and may be more appropriate for assessing model performance.
aic_R2_n <- NagelkerkeR2(stepwise_aic) 
aic_R2_n  

bic_R2_n <- NagelkerkeR2(stepwise_bic) 
bic_R2_n 


### 3.4 Maximum likelihood ratio test ----
# Assess the fit of the different models (containing more or less predictors) using the maximum likelihood ratio test (MLRT)
# > 0.05 indicate a good fit

# Obtain log-likelihood values
full_loglik <- logLik(full_mod)
full_loglik 

aic_loglik <- logLik(stepwise_aic)
aic_loglik 

bic_loglik <- logLik(stepwise_bic)
bic_loglik 

# MLRT statistic and p-value 
lrtest(full_mod, stepwise_aic) 
lrtest(full_mod, stepwise_bic) 

## Clean up
rm(full_mod, stepwise_bic)

# 3.5 LASSO ----
# https://www.statology.org/lasso-regression-in-r/
# https://www.r-bloggers.com/2020/05/quick-tutorial-on-lasso-regression-with-example/

### Function to prepare predictors as a matrix (including splines) for LASSO model
LASSO_mat <- function(data){ # data = the data to be used e.g. df or df_testing
  
  x <- data %>% 
    
    dplyr::select(
      # Socio-demographics
      sex, age1, age2, age3, simd, hhold_bands,
      
      # Geographic
      ur, #hb,
      
      # Testing
      #tes_n, 
      var_per,
      
      # Vaccinations
      vac_dos,
      
      # Risk factors
      shielding, immuno_supp, carehome, bmi_imp1, bmi_imp2, 
      Q_DIAG_AF, Q_DIAG_ASTHMA, Q_DIAG_BLOOD_CANCER, Q_DIAG_CCF, Q_DIAG_CHD, 
      Q_DIAG_CKD_LEVEL, Q_DIAG_COPD, Q_DIAG_DEMENTIA, Q_DIAG_DIABETES_1, Q_DIAG_DIABETES_2, 
      Q_DIAG_EPILEPSY, Q_DIAG_FRACTURE, Q_DIAG_NEURO , Q_DIAG_PARKINSONS, Q_DIAG_PULM_HYPER, 
      Q_DIAG_PULM_RARE, Q_DIAG_PVD, Q_DIAG_RA_SLE, Q_DIAG_RESP_CANCER, Q_DIAG_SEV_MENT_ILL, 
      Q_DIAG_STROKE, Q_DIAG_VTE, 
      
      # Prescriptions
      PIS_Lipid_regulating_drugs, PIS_Angiotensin_converting_enzyme_inhibitors, PIS_Beta_adrenoceptor_blocking_drugs, 
      PIS_Oral_iron, PIS_Selective_serotonin_re_uptake_inhibitors, PIS_Loratadine, PIS_Direct_oral_anticoagulants, 
      PIS_Colchicine, PIS_Antiplatelet_drugs, PIS_Alpha_adrenoceptor_blocking_drugs, PIS_Macrolides, 
      PIS_Benzylpenicillin_and_phenoxymethylpenicillin, PIS_Leukotriene_receptor_antagonists, 
      PIS_Herpes_simplex_and_varicella_zoster, PIS_Replacement_therapy, PIS_Famotidine, PIS_Ursodeoxycholic_acid, 
      PIS_Systemic_nasal_decongestants, PIS_Warfarin_sodium, PIS_Parenteral_anticoagulants, 
      PIS_Compound_bronchodilator_preparations, PIS_Ranitidine_hydrochloride,
      
      # Severity of acute infection
      severe)
  
  return(data.matrix(x))
}

### Set the response variable
y <- df$depvar

### Set predictor variables (matrix)
x <- LASSO_mat(df)

### Perform k-fold cross-validation to find optimal lambda value
cv_lasso <- cv.glmnet(x, y,
                      weights = df$class_weights,
                      alpha = 1, # alpha=1 is the lasso penalty, and alpha=0 the ridge penalty
                      nfolds = 10)

### Identify optimal lambda that minimizes mean squared error (MSE)
lambda_optimal <- cv_lasso$lambda.min
lambda_optimal 

### Plot MSE by lambda value
plot(cv_lasso)

### Fit multiple LASSO models on 50% subsets of training data
# Prepare variable labels
LASSO_labs <- as.data.frame(colnames(x))

# Set variable labels
LASSO_labs <- LASSO_labs %>% 
  rename(label = `colnames(x)`) %>% 
  mutate(label = case_when(
    label == "sex"  ~	"Sex",
    label == "age1" ~ "Age 18 - 33 (spline 1)",
    label == "age2" ~ "Age 34 - 51 (spline 2)",
    label == "age3" ~ "Age 52+ (spline 3)",
    label == "simd" ~ "Scottish Index of Multiple Deprivation",
    label == "hhold_bands" ~ "Household size",
    label == "ur" ~ "Urban-Rural Classification",
    label == "tes_n" ~ "PCR tests taken (up to first positive PCR test)",
    label == "var_per"   ~ "Variant period",
    label == "vac_dos" ~ "Vaccine doses by 14 days before positive PCR test",
    label == "shielding" ~ "Shielding",
    label == "immuno_supp" ~ "Immunosuppressed",
    label == "carehome" ~ "Carehome resident",
    label == "bmi_imp1" ~ "Body Mass Index < 28 (spline 1)",
    label == "bmi_imp2" ~ "Body Mass Index 28+ (spline 2)",
    label == "Q_DIAG_AF" ~ "Atrial fibrillation",
    label == "Q_DIAG_ASTHMA" ~ "Asthma",
    label == "Q_DIAG_BLOOD_CANCER" ~ "Blood cancer",
    label == "Q_DIAG_CCF" ~ "Heart failure",
    label == "Q_DIAG_CHD" ~ "Coronary heart disease",
    label == "Q_DIAG_CKD_LEVEL" ~ "Kidney disease (level 3+)",
    label == "Q_DIAG_COPD" ~ "Chronic obstructive pulmonary disease (COPD)",
    label == "Q_DIAG_DEMENTIA" ~ "Dementia",
    label == "Q_DIAG_DIABETES_1" ~ "Diabetes Type I",
    label == "Q_DIAG_DIABETES_2" ~ "Diabetes Type II",
    label == "Q_DIAG_EPILEPSY" ~ "Epilepsy",
    label == "Q_DIAG_FRACTURE" ~ "Fracture",
    label == "Q_DIAG_NEURO" ~ "Neurological disorder",
    label == "Q_DIAG_PARKINSONS" ~ "Parkinsons",
    label == "Q_DIAG_PULM_HYPER" ~ "Pulmonary hypertension",
    label == "Q_DIAG_PULM_RARE" ~ "Pulmonary rare",
    label == "Q_DIAG_PVD" ~ "Peripheral vascular disease",
    label == "Q_DIAG_RA_SLE" ~ "Rheumatoid arthritis or SLE",
    label == "Q_DIAG_RESP_CANCER" ~ "Respiratory cancer",
    label == "Q_DIAG_SEV_MENT_ILL" ~ "Severe mental illness",
    label == "Q_DIAG_STROKE" ~ "Stroke/TIA",
    label == "Q_DIAG_VTE" ~ "A thrombosis or pulmonary embolus",
    label == "PIS_Lipid_regulating_drugs" ~ "Lipid-regulating drugs",
    label == "PIS_Angiotensin_converting_enzyme_inhibitors" ~ "Angiotensin-converting enzyme inhibitors",
    label == "PIS_Beta_adrenoceptor_blocking_drugs" ~ "Beta-adrenoceptor blocking drugs",
    label == "PIS_Oral_iron" ~ "Oral iron",
    label == "PIS_Selective_serotonin_re_uptake_inhibitors" ~ "Selective serotonin re-uptake inhibitors",
    label == "PIS_Loratadine" ~ "Loratadine (antihistamine)",
    label == "PIS_Direct_oral_anticoagulants" ~ "Direct oral anticoagulants",
    label == "PIS_Colchicine" ~ "Colchicine (anti-inflammatory)",
    label == "PIS_Antiplatelet_drugs" ~ "Antiplatelet drugs",
    label == "PIS_Alpha_adrenoceptor_blocking_drugs" ~ "Alpha-adrenoceptor blocking drugs",
    label == "PIS_Macrolides" ~ "Macrolides (antibacterial)",
    label == "PIS_Benzylpenicillin_and_phenoxymethylpenicillin" ~ "Benzylpenicillin and phenoxymethylpenicillin",
    label == "PIS_Leukotriene_receptor_antagonists" ~ "Leukotriene receptor antagonists",
    label == "PIS_Herpes_simplex_and_varicella_zoster" ~ "Herpes simplex and varicella-zoster (antiviral)",
    label == "PIS_Replacement_therapy" ~ "Corticosteroid replacement therapy",
    label == "PIS_Famotidine" ~ "Famotidine (histamine H2 receptor antagonist)",
    label == "PIS_Ursodeoxycholic_acid" ~ "Ursodeoxycholic acid",
    label == "PIS_Systemic_nasal_decongestants" ~ "Systemic nasal decongestants",
    label == "PIS_Warfarin_sodium" ~ "Warfarin sodium",
    label == "PIS_Parenteral_anticoagulants" ~ "Parenteral anticoagulants",
    label == "PIS_Compound_bronchodilator_preparations" ~ "Compound bronchodilator preparations",
    label == "PIS_Ranitidine_hydrochloride" ~ "Ranitidine hydrochloride",
    label == "severe" ~ "Severity of acute infection"))

# Set order
LASSO_order <- c("Sex", "Age 18 - 33 (spline 1)", "Age 34 - 51 (spline 2)", "Age 52+ (spline 3)", LASSO_labs$label[!LASSO_labs$label %in% c("Sex", "Age 18 - 33 (spline 1)", "Age 34 - 51 (spline 2)", "Age 52+ (spline 3)")])

# Set number of models to run
n_models <- 1000

# Create matrix to collect results of predictor selection
selected_vars <- matrix(0, nrow = ncol(x), ncol = n_models) 

# Run n LASSO models and collect predictor selection and coefficients
for (i in 1:n_models) { 
  idx <- sample(nrow(x), nrow(x) * 0.5) # Get indices for a random 50% of the rows in x
  x_i <- x[idx, ] # Select predictors at those indices
  y_i <- y[idx] # Select dependent vars at those indices
  weights_i <- df$class_weights[idx] # Select weights at those indices
  LASSO_model_i <- glmnet(x_i, y_i, weights = weights_i, alpha = 1, lambda = lambda_optimal) # LASSO model
  selected_vars[, i] <- as.numeric(coef(LASSO_model_i) != 0)[-1] # Identify whether each variable has a coefficient > 0 i.e. is not dropped
}

# Compute proportion of models in which each variable was selected
selected_vars <- t(selected_vars)
selected_vars <- as.data.frame(selected_vars)
colnames(selected_vars) <- c(LASSO_labs$label)
prop_selected <- colMeans(selected_vars)
prop_selected <- as.data.frame(prop_selected)
prop_selected$label <- colnames(selected_vars)
prop_selected$label <- factor(prop_selected$label, levels = rev(LASSO_order))


# Identify variables removed during AIC-predictor selection
AIC_removed <- c("Heart failure", "Epilepsy", "Fracture", "Parkinsons", "Pulmonary hypertension", "Peripheral vascular disease", 
                 "Respiratory cancer", "Stroke/TIA", "Direct oral anticoagulants", "Alpha-adrenoceptor blocking drugs", "Corticosteroid replacement therapy",  
                 "Warfarin sodium", "Ranitidine hydrochloride")

prop_selected$AIC_selected <- ifelse(prop_selected$label %in% AIC_removed, "Removed from main model during predictor selection", "Included in main model")
prop_selected$AIC_selected <- factor(prop_selected$AIC_selected)
color_pal <- c("Included in main model" = "orange", 
               "Removed from main model\nduring predictor selection" = "grey")

# Plot proportion of models each variable was selected in
prop_selected_plot <- ggplot(prop_selected, aes(x=label, y=prop_selected)) +
  geom_segment(aes(x=label, xend=label, y=0.0, yend=prop_selected), color="grey") +
  geom_point(aes(color=AIC_selected), size=3) +
  scale_color_manual(values=color_pal) +
  xlab("") +
  ylab("Proportion of LASSO models with predictor selected") +
  ylim(0.7, 1) +
  geom_hline(yintercept=0.95, linetype="dotted", color = "black") +
  coord_flip() +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title=element_blank(),
    legend.position = "bottom")

prop_selected_plot

setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")
ggsave(paste0("LASSO_selection - ", depvar_name, ".jpg"), prop_selected_plot, height = 10, width = 8, units = "in", dpi = 300)


# 4. Run model using k-fold cross validation ----
#https://daviddalpiaz.github.io/r4sl/the-caret-package.html
#http://www.sthda.com/english/articles/38-regression-model-validation/157-cross-validation-essentials-in-r/

## Make dependent variable a factor ----
df$depvar <- as.factor(df$depvar)

## Set control parameters for cross-validation ----
cv_results <- trainControl(method = "cv", number = 10, savePredictions = T) # cross-validation with 10 folds

## 4.1 AIC-fit ----
## Use predictors suggested by AIC selection
## Fit model using cross-validation (training data) to obtain results for each fold of the cross-validation
## Include 'relevel' for plotting (re-run without below - necessary so that the caret::predict works)
cv_model_plot <- train(depvar
                       
                       # Socio-demographic
                       ~ relevel(sex, ref = "M") # Force male to be reference category
                       + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
                       + relevel(simd, ref = 5) # Force least deprived to be reference category
                       + hhold_bands
                       
                       # Geographic
                       + relevel(ur, ref = 1) # Force remote rural to be reference category
                       #+ relevel(hb, ref = "NHS Lothian") # Force Lothian to be reference category
                       
                       # Testing
                       #+ tes_n # not using
                       + relevel(var_per, ref = "Wild") # force Wild to be the reference category
                       
                       # Vaccinations
                       + vac_dos
                       
                       # Risk factors
                       + shielding
                       + immuno_supp
                       + carehome 
                       + bmi_imp1 + bmi_imp2 # splines in bmi 
                       + Q_DIAG_AF 
                       + Q_DIAG_ASTHMA 
                       + Q_DIAG_BLOOD_CANCER 
                       ### + Q_DIAG_CCF 
                       + Q_DIAG_CHD  
                       + Q_DIAG_CKD_LEVEL  
                       + Q_DIAG_COPD 
                       + Q_DIAG_DEMENTIA 
                       + Q_DIAG_DIABETES_1 
                       + Q_DIAG_DIABETES_2 
                       ###+ Q_DIAG_EPILEPSY  
                       ###+ Q_DIAG_FRACTURE 
                       + Q_DIAG_NEURO   
                       ###+ Q_DIAG_PARKINSONS 
                       ### + Q_DIAG_PULM_HYPER  
                       + Q_DIAG_PULM_RARE  
                       ###+ Q_DIAG_PVD   
                       + Q_DIAG_RA_SLE 
                       ###+ Q_DIAG_RESP_CANCER   
                       + Q_DIAG_SEV_MENT_ILL 
                       ###+ Q_DIAG_STROKE  
                       + Q_DIAG_VTE  
                       
                       # Prescriptions
                       + PIS_Lipid_regulating_drugs
                       + PIS_Angiotensin_converting_enzyme_inhibitors
                       + PIS_Beta_adrenoceptor_blocking_drugs
                       + PIS_Oral_iron
                       + PIS_Selective_serotonin_re_uptake_inhibitors
                       #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
                       + PIS_Loratadine
                       ###+ PIS_Direct_oral_anticoagulants 
                       + PIS_Colchicine
                       + PIS_Antiplatelet_drugs
                       ###+ PIS_Alpha_adrenoceptor_blocking_drugs
                       + PIS_Macrolides
                       + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
                       + PIS_Leukotriene_receptor_antagonists
                       + PIS_Herpes_simplex_and_varicella_zoster
                       ###+ PIS_Replacement_therapy 
                       + PIS_Famotidine
                       + PIS_Ursodeoxycholic_acid
                       + PIS_Systemic_nasal_decongestants
                       ###+ PIS_Warfarin_sodium 
                       + PIS_Parenteral_anticoagulants 
                       + PIS_Compound_bronchodilator_preparations
                       ###+ PIS_Ranitidine_hydrochloride 
                       
                       # Severity of acute infection
                       + severe, 
                       
                       # Model features 
                       data = df,
                       family = "binomial", 
                       method = "glm",
                       weights = df$class_weight,
                       trControl = cv_results)

# Plot
cv_plot <- coef_plot(mod = cv_model_plot$finalModel, plot_name = "Coefficient plot - 10-fold", type = "cv")


# Run model without 'relevel' (so that stats::predict works)
cv_model <- train(depvar
                  
                  # Socio-demographic
                  ~ sex # Force male to be reference category
                  + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
                  + simd # Force least deprived to be reference category
                  + hhold_bands
                  
                  # Geographic
                  + ur # Force remote rural to be reference category
                  #+ hb, ref = "NHS Lothian") # Force Lothian to be reference category
                  
                  # Testing
                  #+ tes_n
                  + var_per # force Wild to be the reference category
                  
                  # Vaccinations
                  + vac_dos
                  
                  # Risk factors
                  + shielding
                  + immuno_supp
                  + carehome 
                  + bmi_imp1 + bmi_imp2 # splines in bmi 
                  + Q_DIAG_AF 
                  + Q_DIAG_ASTHMA 
                  + Q_DIAG_BLOOD_CANCER 
                  ### + Q_DIAG_CCF 
                  + Q_DIAG_CHD  
                  + Q_DIAG_CKD_LEVEL  
                  + Q_DIAG_COPD 
                  + Q_DIAG_DEMENTIA 
                  + Q_DIAG_DIABETES_1 
                  + Q_DIAG_DIABETES_2 
                  ###+ Q_DIAG_EPILEPSY  
                  ###+ Q_DIAG_FRACTURE 
                  + Q_DIAG_NEURO   
                  ###+ Q_DIAG_PARKINSONS 
                  ### + Q_DIAG_PULM_HYPER  
                  + Q_DIAG_PULM_RARE  
                  ###+ Q_DIAG_PVD   
                  + Q_DIAG_RA_SLE 
                  ###+ Q_DIAG_RESP_CANCER   
                  + Q_DIAG_SEV_MENT_ILL 
                  ###+ Q_DIAG_STROKE  
                  + Q_DIAG_VTE  
                  
                  # Prescriptions
                  + PIS_Lipid_regulating_drugs
                  + PIS_Angiotensin_converting_enzyme_inhibitors
                  + PIS_Beta_adrenoceptor_blocking_drugs
                  + PIS_Oral_iron
                  + PIS_Selective_serotonin_re_uptake_inhibitors
                  #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
                  + PIS_Loratadine
                  ###+ PIS_Direct_oral_anticoagulants 
                  + PIS_Colchicine
                  + PIS_Antiplatelet_drugs
                  ###+ PIS_Alpha_adrenoceptor_blocking_drugs
                  + PIS_Macrolides
                  + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
                  + PIS_Leukotriene_receptor_antagonists
                  + PIS_Herpes_simplex_and_varicella_zoster
                  ###+ PIS_Replacement_therapy 
                  + PIS_Famotidine
                  + PIS_Ursodeoxycholic_acid
                  + PIS_Systemic_nasal_decongestants
                  ###+ PIS_Warfarin_sodium 
                  + PIS_Parenteral_anticoagulants 
                  + PIS_Compound_bronchodilator_preparations
                  ###+ PIS_Ranitidine_hydrochloride 
                  
                  # Severity of acute infection
                  + severe, 
                  
                  # Model features 
                  data = df,
                  family = "binomial", 
                  method = "glm",
                  weights = df$class_weight,
                  trControl = cv_results)

## Get predicted probabilities
df_testing$AIC_pred <- predict(cv_model, newdata = df_testing, type = "prob")[,2]

## 4.2 LASSO-fit ----
## Use predictors suggested by LASSO
## Fit model using cross-validation (training data) to obtain results for each fold of the cross-validation
## Include 'relevel' for plotting (re-run without below - necessary so that the caret::predict works)
lasso_mod_for_plot_fun <- function(data){
  
  mod <- train(depvar
               
               # Socio-demographic
               ~ relevel(sex, ref = "M") # Force male to be reference category
               + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
               + relevel(simd, ref = 5) # Force least deprived to be reference category
               + hhold_bands
               
               # Geographic
               + relevel(ur, ref = 1) # Force remote rural to be reference category
               #+ relevel(hb, ref = "NHS Lothian") # Force Lothian to be reference category
               
               # Testing
               #+ tes_n # not using
               + relevel(var_per, ref = "Wild") # force Wild to be the reference category
               
               # Vaccinations
               + vac_dos
               
               # Risk factors
               + shielding
               + immuno_supp
               + carehome 
               + bmi_imp1 + bmi_imp2 # splines in bmi 
               + Q_DIAG_AF 
               + Q_DIAG_ASTHMA 
               + Q_DIAG_BLOOD_CANCER 
               ### + Q_DIAG_CCF 
               + Q_DIAG_CHD  
               # LASSO removed+ Q_DIAG_CKD_LEVEL  
               + Q_DIAG_COPD 
               + Q_DIAG_DEMENTIA 
               + Q_DIAG_DIABETES_1 
               + Q_DIAG_DIABETES_2 
               ###+ Q_DIAG_EPILEPSY  
               ###+ Q_DIAG_FRACTURE 
               # LASSO removed+ Q_DIAG_NEURO   
               ###+ Q_DIAG_PARKINSONS 
               ### + Q_DIAG_PULM_HYPER  
               + Q_DIAG_PULM_RARE  
               ###+ Q_DIAG_PVD   
               + Q_DIAG_RA_SLE 
               ###+ Q_DIAG_RESP_CANCER   
               + Q_DIAG_SEV_MENT_ILL 
               ###+ Q_DIAG_STROKE  
               + Q_DIAG_VTE  
               
               # Prescriptions
               + PIS_Lipid_regulating_drugs
               + PIS_Angiotensin_converting_enzyme_inhibitors
               + PIS_Beta_adrenoceptor_blocking_drugs
               + PIS_Oral_iron
               + PIS_Selective_serotonin_re_uptake_inhibitors
               #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
               + PIS_Loratadine
               ###+ PIS_Direct_oral_anticoagulants 
               + PIS_Colchicine
               # LASSO removed+ PIS_Antiplatelet_drugs
               ###+ PIS_Alpha_adrenoceptor_blocking_drugs
               + PIS_Macrolides
               + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
               + PIS_Leukotriene_receptor_antagonists
               + PIS_Herpes_simplex_and_varicella_zoster
               ###+ PIS_Replacement_therapy 
               + PIS_Famotidine
               + PIS_Ursodeoxycholic_acid
               + PIS_Systemic_nasal_decongestants
               ###+ PIS_Warfarin_sodium 
               + PIS_Parenteral_anticoagulants 
               + PIS_Compound_bronchodilator_preparations
               ###+ PIS_Ranitidine_hydrochloride 
               
               # Severity of acute infection
               + severe,   
               
               # Model features 
               data = data,
               family = "binomial", 
               method = "glm",
               weights = class_weights,
               trControl = cv_results)
  
  return(mod)
}

lasso_mod_for_plot <- lasso_mod_for_plot_fun(data = df)

## Save coefficients and CI
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")
sink(file = paste0("Main model - ", depvar_name, ".txt"))
lasso_mod_for_plot$finalModel
sink(file = NULL)

## Extract model coefficients, odds ratios, and confidence intervals (excluding NAs)
mod <- lasso_mod_for_plot$finalModel
var = names(coef(mod)[!is.na(coef(mod))])
coef = exp(coef(mod)[!is.na(coef(mod))])
ci = data.frame(exp(confint.default(mod)))
ci <- ci %>% dplyr::filter(!is.na(X2.5..))
pval = summary(mod)$coefficients[,4]
depvar_plot <- rep("Main outcome measure", length(pval))

## Combine
coef_df_main_mod <- data.frame(var, coef, ci, pval, depvar_plot)

# Plot
lasso_mod_plot <- coef_plot(mod = lasso_mod_for_plot$finalModel, plot_name = "Coefficient plot - 10-fold - LASSO selection", type = "cv")
rm(lasso_mod_plot)

## Sensitivity - re-do plot including only (i) individuals identified as having long COVID OR (ii) with 6 months follow up  ----
complete_follow_up_df <- df[df$depvar == 1 | df$incomplete_follow_up == 0,] # N = 827,750 (of which 6.6% long COVID (54,750))

lasso_mod_for_plot_complete <- lasso_mod_for_plot_fun(data = complete_follow_up_df)

## Save coefficients and CI
sink(file = paste0("Main model restricted to complete follow up - ", depvar_name, ".txt"))
lasso_mod_for_plot_complete$finalModel
sink(file = NULL)

## Extract model coefficients, odds ratios, and confidence intervals (excluding NAs)
mod_complete <- lasso_mod_for_plot_complete$finalModel
var_complete = names(coef(mod_complete)[!is.na(coef(mod_complete))])
coef_complete = exp(coef(mod_complete)[!is.na(coef(mod_complete))])
ci_complete = data.frame(exp(confint.default(mod_complete)))
ci_complete <- ci_complete %>% dplyr::filter(!is.na(X2.5..))
pval_complete = summary(mod_complete)$coefficients[,4]
depvar_plot_complete <- rep("Main outcome measure restricted to complete follow up", length(pval_complete))

## Combine
coef_df_main_mod_complete <- data.frame(var_complete, coef_complete, ci_complete, pval_complete, depvar_plot_complete)


# Plot
lasso_mod_plot_complete <- coef_plot(mod = lasso_mod_for_plot_complete$finalModel, plot_name = "Coefficient plot - restricted to complete follow up - 10-fold - LASSO selection", type = "cv")
rm(lasso_mod_plot_complete)

# Run model without 'relevel' (so that predict() works)
lasso_model_fun <- function(){
  mod <- train(depvar
               
               # Socio-demographic
               ~ sex # Force male to be reference category
               + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
               + simd # Force least deprived to be reference category
               + hhold_bands 
               
               # Geographic
               + ur 
               #+ hb, ref = "NHS Lothian") # Force Lothian to be reference category
               
               # Testing
               #+ tes_n
               + var_per # force Wild to be the reference category
               
               # Vaccinations
               + vac_dos
               
               # Risk factors
               + shielding
               + immuno_supp
               + carehome 
               + bmi_imp1 + bmi_imp2 # splines in bmi 
               + Q_DIAG_AF 
               + Q_DIAG_ASTHMA 
               + Q_DIAG_BLOOD_CANCER 
               ### + Q_DIAG_CCF 
               + Q_DIAG_CHD  
               # LASSO removed+ Q_DIAG_CKD_LEVEL  
               + Q_DIAG_COPD 
               + Q_DIAG_DEMENTIA 
               + Q_DIAG_DIABETES_1 
               + Q_DIAG_DIABETES_2 
               ###+ Q_DIAG_EPILEPSY  
               ###+ Q_DIAG_FRACTURE 
               # LASSO removed+ Q_DIAG_NEURO   
               ###+ Q_DIAG_PARKINSONS 
               ### + Q_DIAG_PULM_HYPER  
               + Q_DIAG_PULM_RARE  
               ###+ Q_DIAG_PVD   
               + Q_DIAG_RA_SLE 
               ###+ Q_DIAG_RESP_CANCER   
               + Q_DIAG_SEV_MENT_ILL 
               ###+ Q_DIAG_STROKE  
               + Q_DIAG_VTE  
               
               # Prescriptions
               + PIS_Lipid_regulating_drugs
               + PIS_Angiotensin_converting_enzyme_inhibitors
               + PIS_Beta_adrenoceptor_blocking_drugs
               + PIS_Oral_iron
               + PIS_Selective_serotonin_re_uptake_inhibitors
               #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
               + PIS_Loratadine
               ###+ PIS_Direct_oral_anticoagulants 
               + PIS_Colchicine
               # LASSO removed+ PIS_Antiplatelet_drugs
               ###+ PIS_Alpha_adrenoceptor_blocking_drugs
               + PIS_Macrolides
               + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
               + PIS_Leukotriene_receptor_antagonists
               + PIS_Herpes_simplex_and_varicella_zoster
               ###+ PIS_Replacement_therapy 
               + PIS_Famotidine
               + PIS_Ursodeoxycholic_acid
               + PIS_Systemic_nasal_decongestants
               ###+ PIS_Warfarin_sodium 
               + PIS_Parenteral_anticoagulants 
               + PIS_Compound_bronchodilator_preparations
               ###+ PIS_Ranitidine_hydrochloride 
               
               # Severity of acute infection
               + severe,   
               
               # Model features 
               data = df,
               family = "quasibinomial", 
               method = "glm",
               weights = df$class_weight,
               trControl = cv_results)
  return(mod)
}

lasso_model <- lasso_model_fun()

## Get predicted probabilities
df$LASSO_pred <- predict(lasso_model, newdata = df, type = "prob")[,2]
df_testing$LASSO_pred <- predict(lasso_model, newdata = df_testing, type = "prob")[,2]

# 5. Naive Bayes Classification ----
## Use predictors suggested by LASSO
## Not using cross validation
## Includes Laplace Correction so that feature value that appear in the test data but not the training data will not default to 0
## this reduces overfitting and results in more stable and robust prediction
df$depvar <- as.factor(df$depvar)

nb_model <- naive_bayes(depvar
                        
                        # Socio-demographic
                        ~ sex # Force male to be reference category
                        + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
                        + simd # Force least deprived to be reference category
                        + hhold_bands 
                        
                        # Geographic
                        + ur 
                        #+ hb, ref = "NHS Lothian") # Force Lothian to be reference category
                        
                        # Testing
                        #+ tes_n
                        + var_per # force Wild to be the reference category
                        
                        # Vaccinations
                        + vac_dos
                        
                        # Risk factors
                        + shielding
                        + immuno_supp
                        + carehome 
                        + bmi_imp1 + bmi_imp2 # splines in bmi 
                        + Q_DIAG_AF 
                        + Q_DIAG_ASTHMA 
                        + Q_DIAG_BLOOD_CANCER 
                        ### + Q_DIAG_CCF 
                        + Q_DIAG_CHD  
                        # LASSO removed+ Q_DIAG_CKD_LEVEL  
                        + Q_DIAG_COPD 
                        + Q_DIAG_DEMENTIA 
                        + Q_DIAG_DIABETES_1 
                        + Q_DIAG_DIABETES_2 
                        ###+ Q_DIAG_EPILEPSY  
                        ###+ Q_DIAG_FRACTURE 
                        # LASSO removed+ Q_DIAG_NEURO   
                        ###+ Q_DIAG_PARKINSONS 
                        ### + Q_DIAG_PULM_HYPER  
                        + Q_DIAG_PULM_RARE  
                        ###+ Q_DIAG_PVD   
                        + Q_DIAG_RA_SLE 
                        ###+ Q_DIAG_RESP_CANCER   
                        + Q_DIAG_SEV_MENT_ILL 
                        ###+ Q_DIAG_STROKE  
                        + Q_DIAG_VTE  
                        
                        # Prescriptions
                        + PIS_Lipid_regulating_drugs
                        + PIS_Angiotensin_converting_enzyme_inhibitors
                        + PIS_Beta_adrenoceptor_blocking_drugs
                        + PIS_Oral_iron
                        + PIS_Selective_serotonin_re_uptake_inhibitors
                        #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
                        + PIS_Loratadine
                        ###+ PIS_Direct_oral_anticoagulants 
                        + PIS_Colchicine
                        # LASSO removed+ PIS_Antiplatelet_drugs
                        ###+ PIS_Alpha_adrenoceptor_blocking_drugs
                        + PIS_Macrolides
                        + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
                        + PIS_Leukotriene_receptor_antagonists
                        + PIS_Herpes_simplex_and_varicella_zoster
                        ###+ PIS_Replacement_therapy 
                        + PIS_Famotidine
                        + PIS_Ursodeoxycholic_acid
                        + PIS_Systemic_nasal_decongestants
                        ###+ PIS_Warfarin_sodium 
                        + PIS_Parenteral_anticoagulants 
                        + PIS_Compound_bronchodilator_preparations
                        ###+ PIS_Ranitidine_hydrochloride 
                        
                        # Severity of acute infection
                        + severe,   
                        
                        # Model features 
                        data = df,
                        family = "binomial", 
                        weights = df$class_weight,
                        laplace = 1)

nb_coefs <- coef(nb_model)

## NB this predicts the classification of each individual - it doesn't produce coefficients for predictors

## Get predicted probabilities
df$NB_pred <- predict(nb_model, newdata = df, type = "prob")[,2]
df_testing$NB_pred <- predict(nb_model, newdata = df_testing, type = "prob")[,2]


# 6. Gradient boosting decision trees (xgboost) ----
## Convert data into a DMatrix object (LASSO_mat function helps with this)
df_matrix <- LASSO_mat(df) # omits df$depvar
dtrain <- xgb.DMatrix(data = df_matrix, label = as.numeric(as.character(df$depvar)))

## Specify the model parameters 
params <- list(objective = "binary:logistic", eval_metric = "logloss")

## Train the GBDT model with 10-fold cross validaiton to identify the optimal number of iterations 
gbdt_cv_model <- xgb.cv(data = dtrain, params = params, nfold = 10, nrounds = 100, metrics = "logloss", verbose = TRUE)
print(gbdt_cv_model)

## Get the optimal number of iterations
optimal_iterations <- which.min(gbdt_cv_model$evaluation_log$test_logloss_mean) # 51

## Train the final GBDT model with the optimal number of iterations (no cross validation i.e. use full dataset)
gbdt_model <- xgb.train(data = dtrain, params = params, nrounds = optimal_iterations)

## 6.1 Prepare data and plot feature importance scores ----
importance_scores <- xgb.importance(model = gbdt_model)
importance_scores <- as.data.frame(xgb.plot.importance(importance_matrix = importance_scores))

importance_scores_df <- importance_scores %>% 
  dplyr::select(Feature, Importance) %>% 
  mutate(Feature = case_when(
    Feature == "sex"  ~	"Sex",
    Feature == "age1" ~ "Age 18 - 33 (spline 1)",
    Feature == "age2" ~ "Age 34 - 51 (spline 2)",
    Feature == "age3" ~ "Age 52+ (spline 3)",
    Feature == "simd" ~ "Scottish Index of Multiple Deprivation",
    Feature == "hhold_bands" ~ "Household size",
    Feature == "ur" ~ "Urban-Rural Classification",
    Feature == "var_per"   ~ "Variant period",
    Feature == "vac_dos" ~ "Vaccine doses",
    Feature == "shielding" ~ "Shielding",
    Feature == "immuno_supp" ~ "Immunosuppressed",
    Feature == "carehome" ~ "Carehome resident",
    Feature == "bmi_imp1" ~ "BMI < 28 (spline 1)",
    Feature == "bmi_imp2" ~ "BMI 28+ (spline 2)",
    Feature == "Q_DIAG_AF" ~ "Atrial fibrillation",
    Feature == "Q_DIAG_ASTHMA" ~ "Asthma",
    Feature == "Q_DIAG_BLOOD_CANCER" ~ "Blood cancer",
    Feature == "Q_DIAG_CCF" ~ "Heart failure",
    Feature == "Q_DIAG_CHD" ~ "Coronary heart disease",
    Feature == "Q_DIAG_CKD_LEVEL" ~ "Kidney disease (level 3+)",
    Feature == "Q_DIAG_COPD" ~ "Chronic obstructive pulmonary disease (COPD)",
    Feature == "Q_DIAG_DEMENTIA" ~ "Dementia",
    Feature == "Q_DIAG_DIABETES_1" ~ "Diabetes Type I",
    Feature == "Q_DIAG_DIABETES_2" ~ "Diabetes Type II",
    Feature == "Q_DIAG_EPILEPSY" ~ "Epilepsy",
    Feature == "Q_DIAG_FRACTURE" ~ "Fracture",
    Feature == "Q_DIAG_NEURO" ~ "Neurological disorder",
    Feature == "Q_DIAG_PARKINSONS" ~ "Parkinsons",
    Feature == "Q_DIAG_PULM_HYPER" ~ "Pulmonary hypertension",
    Feature == "Q_DIAG_PULM_RARE" ~ "Pulmonary rare",
    Feature == "Q_DIAG_PVD" ~ "Peripheral vascular disease",
    Feature == "Q_DIAG_RA_SLE" ~ "Rheumatoid arthritis or SLE",
    Feature == "Q_DIAG_RESP_CANCER" ~ "Respiratory cancer",
    Feature == "Q_DIAG_SEV_MENT_ILL" ~ "Severe mental illness",
    Feature == "Q_DIAG_STROKE" ~ "Stroke/TIA",
    Feature == "Q_DIAG_VTE" ~ "A thrombosis or pulmonary embolus",
    Feature == "PIS_Lipid_regulating_drugs" ~ "Lipid-regulating drugs",
    Feature == "PIS_Angiotensin_converting_enzyme_inhibitors" ~ "Angiotensin-converting enzyme inhibitors",
    Feature == "PIS_Beta_adrenoceptor_blocking_drugs" ~ "Beta-adrenoceptor blocking drugs",
    Feature == "PIS_Oral_iron" ~ "Oral iron",
    Feature == "PIS_Selective_serotonin_re_uptake_inhibitors" ~ "Selective serotonin re-uptake inhibitors",
    Feature == "PIS_Loratadine" ~ "Loratadine (antihistamine)",
    Feature == "PIS_Direct_oral_anticoagulants" ~ "Direct oral anticoagulants",
    Feature == "PIS_Colchicine" ~ "Colchicine (anti-inflammatory)",
    Feature == "PIS_Antiplatelet_drugs" ~ "Antiplatelet drugs",
    Feature == "PIS_Alpha_adrenoceptor_blocking_drugs" ~ "Alpha-adrenoceptor blocking drugs",
    Feature == "PIS_Macrolides" ~ "Macrolides (antibacterial)",
    Feature == "PIS_Benzylpenicillin_and_phenoxymethylpenicillin" ~ "Benzylpenicillin and phenoxymethylpenicillin",
    Feature == "PIS_Leukotriene_receptor_antagonists" ~ "Leukotriene receptor antagonists",
    Feature == "PIS_Herpes_simplex_and_varicella_zoster" ~ "Herpes simplex and varicella-zoster (antiviral)",
    Feature == "PIS_Replacement_therapy" ~ "Corticosteroid replacement therapy",
    Feature == "PIS_Famotidine" ~ "Famotidine (histamine H2 receptor antagonist)",
    Feature == "PIS_Ursodeoxycholic_acid" ~ "Ursodeoxycholic acid",
    Feature == "PIS_Systemic_nasal_decongestants" ~ "Systemic nasal decongestants",
    Feature == "PIS_Warfarin_sodium" ~ "Warfarin sodium",
    Feature == "PIS_Parenteral_anticoagulants" ~ "Parenteral anticoagulants",
    Feature == "PIS_Compound_bronchodilator_preparations" ~ "Compound bronchodilator preparations",
    Feature == "PIS_Ranitidine_hydrochloride" ~ "Ranitidine hydrochloride",
    Feature == "severe" ~ "Severe acute infection"))

GBDT_plot <- ggplot(importance_scores_df, aes(x = Importance, y = reorder(Feature, Importance))) +
  geom_bar(stat = "identity", fill = "#00BFC4", width = 0.5) +  # Set width to control bar thickness
  labs(x = "Importance Score", y = NULL) +
  theme(axis.text.y = element_text(hjust = 1)) +
  theme_bw() +
  coord_cartesian(xlim = c(0, max(importance_scores_df$Importance) + 0.1))  # Set x-axis limits

print(GBDT_plot)

setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")
ggsave("Feature importance scores.jpg", GBDT_plot, height = 9, width = 8, units = "in")

## Convert testing data to matrix form and obtain predicted probabilities from final iteration
df_testing_matrix <- LASSO_mat(df_testing)
gbdt_pred <- predict(gbdt_model, newdata = df_testing_matrix, type = "prob")
df_testing$GBDT_pred <- gbdt_pred

## NB this predicts feature importance, but it doesn't produce coefficients for predictors

# 7. Sensitivity - comparing dependent variables ----
## 7.1 Prepare data for plotting and get predicted probabilities ----
## For each version of the dependent variable, run the main model and collect coef and CI, get predicted probabilities
depvar_variations <- function(plotname){ #plotname = label to identify the dependent variable used e.g. depvar = "Operational definition"
  
  ## Format dependent variable
  df$depvar <- as.factor(df$depvar)
  
  ## Run main model for plotting (LASSO with 10-fold CV, training data) 
  lasso_mod_for_plot <- lasso_mod_for_plot_fun()
  
  ## Extract model coefficients, odds ratios, and confidence intervals (excluding NAs)
  mod <- lasso_mod_for_plot$finalModel
  var = names(coef(mod)[!is.na(coef(mod))])
  coef = exp(coef(mod)[!is.na(coef(mod))])
  ci = data.frame(exp(confint.default(mod)))
  ci <- ci %>% dplyr::filter(!is.na(X2.5..))
  pval = summary(mod)$coefficients[,4]
  depvar_plot <- rep(plotname, length(pval)) 
  
  ## Combine
  coef_df <- data.frame(var, coef, ci, pval, depvar_plot)
  
  ## Remove model (to save memory)
  rm(lasso_mod_for_plot)
  
  ## Run main model for predicted probabilities (LASSO with 10-fold CV, training data)
  lasso_model_for_probs <- lasso_model_fun()
  
  ## Get predicted probabilities (testing data)
  probs <- predict(lasso_model_fun, newdata = df_testing, type = "prob")[,2]
  
  ## Remove model
  rm(lasso_model_for_probs)
  
  ## Return predicted probabilities
  return(list(probs = probs)) #, coef_df = coef_df)) 
}
  
## Depvar = clinical code, free text, long covid sick note
df$depvar <- df$code_txt_fitnote 
codes <- depvar_variations("Clinical code, free text, or sick note")
coef_df_codes <- codes$coef_df
df_testing$LASSO_pred_codes <- codes$probs

## Depvar = opdef
df$depvar <- ifelse(df$opdef == 1, 1, 0) 
opdef <- depvar_variations("Operational definition")
coef_df_opdef <- opdef$coef_df
df_testing$LASSO_pred_opdef <- opdef$probs

## Depvar = opdef without blood tests
df$depvar <- df$opdef_no_blood 
opdef_nb <- depvar_variations("Operational definition (excluding blood tests)")
coef_df_opdef_nb <- opdef_nb$coef_df
df_testing$LASSO_pred_opdef_nb <- opdef_nb$probs

## Depvar = opdef without blood tests OR clinical code, free text, long covid sick note
df$depvar <- ifelse(df$opdef_no_blood ==1 | df$code_txt_fitnote == 1, 1, 0) 
opdef_nb_codes <- depvar_variations("Operational definition (excluding blood tests)\nor clinical code, free text, or sick note")
coef_df_opdef_nb_codes <- opdef_nb_codes$coef_df
df_testing$LASSO_pred_opdef_nb_codes <- opdef_nb_codes$probs


### 7.2 Plot ---- 
## Combine data for plots
coef_df <- rbind(coef_df_main_mod, coef_df_codes, coef_df_opdef, coef_df_opdef_nb, coef_df_opdef_nb_codes)
row.names(coef_df) <- NULL

## Prepare data for plotting
coef_df <- coef_df %>% 
  
  # Remove intercept
  filter(var != "(Intercept)") %>% 
  
  # Remove 'no dominant variant'
  filter(var != ("`relevel(var_per, ref = \"Wild\")Unknown`")) %>% 
  
  # Label categories
  mutate(subtitle = case_when(
      var == "\`relevel(sex, ref = \"M\")M\`"  ~	"Sex",
      var == "\`relevel(sex, ref = \"M\")F\`" ~	"Sex",
      var == "age1" ~ "Age",
      var == "age2" ~ "Age",
      var == "age3" ~ "Age",
      var == "\`relevel(simd, ref = 5)1\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
      var == "\`relevel(simd, ref = 5)2\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
      var == "\`relevel(simd, ref = 5)3\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
      var == "\`relevel(simd, ref = 5)4\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
      var == "\`relevel(simd, ref = 5)5\`" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
      var == "simdUnknown" ~ "Scottish Index of Multiple Deprivation (Ref: Quintile 5 - least deprived)",
      var == "hhold_bands.L" ~ "Household size (Ref: 1)",
      var == "hhold_bands.Q" ~ "Household size (Ref: 1)",
      var == "hhold_bands.C" ~ "Household size (Ref: 1)",
      var == "\`hhold_bands^4\`" ~ "Household size (Ref: 1)",
      var == "\`relevel(ur, ref = 1)2\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
      var == "\`relevel(ur, ref = 1)3\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
      var == "\`relevel(ur, ref = 1)4\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
      var == "\`relevel(ur, ref = 1)5\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
      var == "\`relevel(ur, ref = 1)6\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
      var == "\`relevel(ur, ref = 1)Unknown\`" ~ "Urban-Rural Classification (Ref: Large Urban Area)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Ayrshire and Arran\`" ~	"Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Borders\`" ~	"Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Dumfries and Galloway\`" ~	     "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Fife\`" ~	"Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Forth Valley\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Grampian\`" ~	"Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Greater Glasgow and Clyde\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Highland\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lanarkshire\`" ~	"Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lothian\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Orkney\`" ~	"Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Shetland\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Tayside\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Western Isles\`" ~ "Health board (Ref: NHS Lothian)",
      var == "\`relevel(hb, ref = \"NHS Lothian\")Unknown\`" ~	"Health board (Ref: NHS Lothian)",
      var == "tes_n.L" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
      var == "tes_n.Q" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
      var == "tes_n.C" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
      var == "\`tes_n^4\`" ~ "PCR tests taken (up to first positive PCR test) (Ref: 1)",
      var == "\`relevel(var_per, ref = \"Wild\")Alpha\`"   ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
      var == "\`relevel(var_per, ref = \"Wild\")Delta\`"   ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
      var == "\`relevel(var_per, ref = \"Wild\")Omicron\`" ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
      var == "\`relevel(var_per, ref = \"Wild\")Unknown\`" ~ "Variant period (Ref: Wild-type (01/03/2020 - 10/01/2021))",
      var == "vac_dos1" ~ "Vaccine doses by 14 days before positive PCR test (Ref: 0)",
      var == "vac_dos2" ~ "Vaccine doses by 14 days before positive PCR test (Ref: 0)",
      var == "\`vac_dos3+\`" ~"Vaccine doses by 14 days before positive PCR test (Ref: 0)",
      var == "shielding" ~ "Risk factors",
      var == "immuno_supp" ~ "Risk factors",
      var == "carehome" ~ "Carehome resident",
      var == "bmi_imp1" ~ "Risk factors",
      var == "bmi_imp2" ~ "Risk factors",
      grepl("Q_", var) ~ "Risk factors",
      grepl("PIS_", var) ~ "Prescriptions dispensed in the 3 months prior to positive PCR test",
      var == "severe" ~ "Severity of acute infection")) %>% 
  
  ## Label variables
  mutate(label = case_when(
      var == "\`relevel(sex, ref = \"M\")F\`" ~	"Female",
      var == "\`relevel(sex, ref = \"M\")M\`" ~	"Male",
      var == "age1" ~ "Age 18 - 33 (spline 1)",
      var == "age2" ~ "Age 34 - 51 (spline 2)",
      var == "age3" ~ "Age 52+ (spline 3)",
      var == "\`relevel(simd, ref = 5)1\`" ~	"SIMD: Quintile 1 (most deprived)",
      var == "\`relevel(simd, ref = 5)2\`" ~	"SIMD: Quintile 2",
      var == "\`relevel(simd, ref = 5)3\`" ~ "SIMD: Quintile 3",
      var == "\`relevel(simd, ref = 5)4\`" ~ "SIMD: Quintile 4",
      var == "\`relevel(simd, ref = 5)5\`" ~ "SIMD: Quintile 5 (least deprived)",
      var == "simdUnknown" ~ "Unknown SIMD",
      var == "hhold_bands.L" ~ "Household size: 2",
      var == "hhold_bands.Q" ~ "Household size: 3-5",
      var == "hhold_bands.C" ~ "Household size: 6-10",
      var == "\`hhold_bands^4\`" ~ "Household size: >10",
      var == "\`relevel(ur, ref = 1)2\`" ~ "Other Urban Areas",
      var == "\`relevel(ur, ref = 1)3\`" ~ "Accessible Small Towns",
      var == "\`relevel(ur, ref = 1)4\`" ~ "Remote Small Towns",
      var == "\`relevel(ur, ref = 1)5\`" ~ "Accessible Rural",
      var == "\`relevel(ur, ref = 1)6\`" ~ "Remote Rural",
      var == "\`relevel(ur, ref = 1)Unknown\`" ~ "Unknown",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Ayrshire and Arran\`" ~	"NHS Ayrshire and Arran",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Borders\`" ~	"NHS Borders",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Dumfries and Galloway\`" ~	"NHS Dumfries and Galloway",      
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Fife\`" ~	"NHS Fife",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Forth Valley\`" ~ "NHS Forth Valley",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Grampian\`" ~	"NHS Grampian",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Greater Glasgow and Clyde\`" ~ "NHS Greater Glasgow and Clyde",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Highland\`" ~ "NHS Highland",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lanarkshire\`" ~	"NHS Lanarkshire",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Lothian\`" ~ "NHS Lothian",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Orkney\`" ~	"NHS Orkney",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Shetland\`" ~ "NHS Shetland",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Tayside\`" ~ "NHS Tayside",
      var == "\`relevel(hb, ref = \"NHS Lothian\")NHS Western Isles\`" ~ "NHS Western Isles",
      var == "\`relevel(hb, ref = \"NHS Lothian\")Unknown\`" ~	"Unknown healthboard",
      var == "tes_n.L" ~ "2 PCR tests",
      var == "tes_n.Q" ~ "3-4 PCR tests",
      var == "tes_n.C" ~ "5-10 PCR tests",
      var == "\`tes_n^4\`" ~ ">10 PCR tests",
      var == "\`relevel(var_per, ref = \"Wild\")Alpha\`"   ~ "Alpha (11/01/2021 - 09/05/2021)",
      var == "\`relevel(var_per, ref = \"Wild\")Delta\`"   ~ "Delta (24/05/2021 - 05/12/2021)",
      var == "\`relevel(var_per, ref = \"Wild\")Omicron\`" ~ "Omicron (27/12/2021 - 30/04/2022)",
      var == "\`relevel(var_per, ref = \"Wild\")Unknown\`" ~ "No dominant variant",
      var == "vac_dos1" ~ "Vaccine doses: 1",
      var == "vac_dos2" ~ "Vaccine doses: 2",
      var == "\`vac_dos3+\`" ~ "Vaccine doses: 3+",
      var == "shielding" ~ "Shielding",
      var == "immuno_supp" ~ "Immunosuppressed",
      var == "carehome" ~ "Carehome resident",
      var == "bmi_imp1" ~ "Body Mass Index < 28 (spline 1)",
      var == "bmi_imp2" ~ "Body Mass Index 28+ (spline 2)",
      var == "Q_DIAG_AF" ~ "Atrial fibrillation",
      var == "Q_DIAG_ASTHMA" ~ "Asthma",
      var == "Q_DIAG_BLOOD_CANCER" ~ "Blood cancer",
      var == "Q_DIAG_CCF" ~ "Heart failure",
      var == "Q_DIAG_CHD" ~ "Coronary heart disease",
      var == "Q_DIAG_CKD_LEVEL" ~ "Kidney disease (level 3+)",
      var == "Q_DIAG_COPD" ~ "Chronic obstructive pulmonary disease (COPD)",
      var == "Q_DIAG_DEMENTIA" ~ "Dementia",
      var == "Q_DIAG_DIABETES_1" ~ "Diabetes Type I",
      var == "Q_DIAG_DIABETES_2" ~ "Diabetes Type II",
      var == "Q_DIAG_EPILEPSY" ~ "Epilepsy",
      var == "Q_DIAG_FRACTURE" ~ "Fracture",
      var == "Q_DIAG_NEURO" ~ "Neurological disorder",
      var == "Q_DIAG_PARKINSONS" ~ "Parkinsons",
      var == "Q_DIAG_PULM_HYPER" ~ "Pulmonary hypertension",
      var == "Q_DIAG_PULM_RARE" ~ "Pulmonary rare",
      var == "Q_DIAG_PVD" ~ "Peripheral vascular disease",
      var == "Q_DIAG_RA_SLE" ~ "Rheumatoid arthritis or SLE",
      var == "Q_DIAG_RESP_CANCER" ~ "Respiratory cancer",
      var == "Q_DIAG_SEV_MENT_ILL" ~ "Severe mental illness",
      var == "Q_DIAG_STROKE" ~ "Stroke/TIA",
      var == "Q_DIAG_VTE" ~ "A thrombosis or pulmonary embolus",
      var == "PIS_Lipid_regulating_drugs" ~ "Lipid-regulating drugs",
      var == "PIS_Angiotensin_converting_enzyme_inhibitors" ~ "Angiotensin-converting enzyme inhibitors",
      var == "PIS_Beta_adrenoceptor_blocking_drugs" ~ "Beta-adrenoceptor blocking drugs",
      var == "PIS_Oral_iron" ~ "Oral iron",
      var == "PIS_Selective_serotonin_re_uptake_inhibitors" ~ "Selective serotonin re-uptake inhibitors",
      var == "PIS_Loratadine" ~ "Loratadine (antihistamine)",
      var == "PIS_Direct_oral_anticoagulants" ~ "Direct oral anticoagulants",
      var == "PIS_Colchicine" ~ "Colchicine (anti-inflammatory)",
      var == "PIS_Antiplatelet_drugs" ~ "Antiplatelet drugs",
      var == "PIS_Alpha_adrenoceptor_blocking_drugs" ~ "Alpha-adrenoceptor blocking drugs",
      var == "PIS_Macrolides" ~ "Macrolides (antibacterial)",
      var == "PIS_Benzylpenicillin_and_phenoxymethylpenicillin" ~ "Benzylpenicillin and phenoxymethylpenicillin",
      var == "PIS_Leukotriene_receptor_antagonists" ~ "Leukotriene receptor antagonists",
      var == "PIS_Herpes_simplex_and_varicella_zoster" ~ "Herpes simplex and varicella-zoster (antiviral)",
      var == "PIS_Replacement_therapy" ~ "Corticosteroid replacement therapy",
      var == "PIS_Famotidine" ~ "Famotidine (histamine H2 receptor antagonist)",
      var == "PIS_Ursodeoxycholic_acid" ~ "Ursodeoxycholic acid",
      var == "PIS_Systemic_nasal_decongestants" ~ "Systemic nasal decongestants",
      var == "PIS_Warfarin_sodium" ~ "Warfarin sodium",
      var == "PIS_Parenteral_anticoagulants" ~ "Parenteral anticoagulants",
      var == "PIS_Compound_bronchodilator_preparations" ~ "Compound bronchodilator preparations",
      var == "PIS_Ranitidine_hydrochloride" ~ "Ranitidine hydrochloride",
      var == "severe" ~ "Hospitalised within 28 days of PCR test"))

# Set column names
colnames(coef_df) <- c("var", "coef", "lower", "upper", "pval", "depvar", "subtitle", "label")

# Set order of facets to appear in plot
coef_df$depvar <- factor(coef_df$depvar, levels = c("Main outcome measure", 
                         "Clinical code, free text, or sick note", 
                         "Operational definition",
                         "Operational definition (excluding blood tests)\nor clinical code, free text, or sick note",
                         "Operational definition (excluding blood tests)"))
coef_df$subtitle <- factor(coef_df$subtitle, levels=unique(coef_df$subtitle))
subtitle <- factor(coef_df$subtitle, levels=unique(coef_df$subtitle)) # Needed for adjusting height of each facet

# Set order variables to appear in plot (with Q_ and PIS_ vars ordered by coefficient size)
## Subset start of dataframe
start <- coef_df[!(grepl("Q_|PIS_|hb", coef_df$var)),]

## Subset of dataframe where 'var' begins with "Q_"
Q_subset <- coef_df[grep("Q_", coef_df$var), ]
Q_subset <- Q_subset[order(Q_subset$coef), ]

## Subset of dataframe where 'var' begins with "PIS_"
PIS_subset <- coef_df[grep("PIS_", coef_df$var), ]
PIS_subset <- PIS_subset[order(PIS_subset$coef), ]

## Combine the subsets and the remaining rows of original dataframe
coef_df <- rbind(start, Q_subset, PIS_subset)

## Set order
coef_df$label <- factor(coef_df$label, levels=rev(unique(coef_df$label)))

## Set custom colours
cust_cols <- c("black", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")

# Create a forest plot using ggplot2
p <- ggplot(coef_df, aes(x = coef, y = label, color = depvar)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(limits = c(min(coef_df$lower) - 0.05, max(coef_df$upper))) +
  scale_color_manual(values = cust_cols) +  # Set the custom colors
  theme_bw() +
  ylab("") +
  xlab("Adjusted Odds Ratio (95% CI)") +
  theme(axis.text.y = element_text(size = 11, color = "black"), # Set the font size for y-axis labels
        legend.position = "none") +  # Remove the legend) 
  facet_wrap(~depvar, nrow = 1)

p

setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")
ggsave("Coefficient plot - comparing dependent variables.jpg", p, height = 11, width = 19, units = "in")


# Export data ----
saveRDS(df, "/conf/EAVE/GPanalysis/analyses/long_covid/data/df_training_cleaned.rds")
saveRDS(df_testing, "/conf/EAVE/GPanalysis/analyses/long_covid/data/df_testing_with_probs.rds")


