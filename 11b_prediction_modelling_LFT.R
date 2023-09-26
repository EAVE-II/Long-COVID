##########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 11b_prediction_modelling_LFT - Uses data on all individuals with a positive PCR OR LFT test.
## Trains the main model identified in 11 (LASSO-selected predictors with 10-fold cross validation). Gets predicted 
## probabilities. Sensitivity: re-runs main model using alterantive dependent variables (main version uses operational 
## definition): (i) opdef without blood tests and with clinical codes, free text, sick notes; (ii) clinical codes, free 
## text, sick notes; (iii) opef without blood tests. Exports df with predicted probabilities.
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


## Read in data prepared for prediction modelling
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_predict_prepared_LFT.rds") # for sensitivty including LFTs

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
  if(var == "age"){ df <- cbind(df, ns(df[, var_i], df = 3))
  start <- ncol(df)-2
  }
  # 2 knots for bmi
  if(var == "bmi_imp"){ df <- cbind(df, ns(df[, var_i], df = 2))
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
## Create inverse class proportion weights in the training data to allow cost function customisation to deal with case imbalance (training data)
## See Budiarto et al (2023) Handling Class Imbalance in Machine Learning-based Prediction Models
### Calculate class proportions
class_proportions <- table(df$depvar) / length(df$depvar)

### Calculate inverse class proportion weights
inverse_weights <- 1 / class_proportions

### Assign the inverse class proportion weights to a new variable in the dataframe
#df$class_weights <- ifelse(df$depvar == 0, inverse_weights[1], inverse_weights[2])
df$class_weights <- 1

## 0.4 Functions ----
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
              size = 3, position = position_nudge(x = 0.01)) +
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
  
  setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling/LFT")
  ggsave(paste0(plot_name, " - ", depvar_name, ".jpg"), gp.wrap, height = 14, width = 11, units = "in")
  
  return(gp.wrap)
}

### Create combined diabetes variable and remove PIS_Metformin_hydrochloride from model (implemented during Set-up)----
df$Q_DIAG_DIABETES_2 <- ifelse(df$PIS_Metformin_hydrochloride == 1, 1, df$Q_DIAG_DIABETES_2 )
df_testing$Q_DIAG_DIABETES_2 <- ifelse(df_testing$PIS_Metformin_hydrochloride == 1, 1, df_testing$Q_DIAG_DIABETES_2 )


# 1. Main model: LASSO-selected predictors with cross validation ----
## Use predictors suggested by LASSO
## Fit model using cross-validation (training data) to obtain results for each fold of the cross-validation
df$depvar <- as.factor(df$depvar)

## Set control parameters for cross-validation
cv_results <- trainControl(method = "cv", number = 10, savePredictions = T) # cross-validation with 10 folds

## Include 'relevel' for plotting (re-run without below - necessary so that the caret::predict works)
lasso_mod_for_plot_fun <- function(){
  
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
               data = df,
               family = "binomial", 
               method = "glm",
               weights = df$class_weight,
               trControl = cv_results)
  
  return(mod)
}


lasso_mod_for_plot <- lasso_mod_for_plot_fun()

## Save coefficients and CI
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

## Run model without 'relevel' (so that predict() works)
#lasso_model_fun <- function(){
#  mod <- train(depvar
#               
#               # Socio-demographic
#               ~ sex # Force male to be reference category
#               + age1 + age2 + age3 # age with natural splines (3 degrees of freedom)
#               + simd # Force least deprived to be reference category
#               + hhold_bands 
#               
#               # Geographic
#               + ur 
#               #+ hb, ref = "NHS Lothian") # Force Lothian to be reference category
#               
#               # Testing
#               #+ tes_n
#               + var_per # force Wild to be the reference category
#               
#               # Vaccinations
#               + vac_dos
#               
#               # Risk factors
#               + shielding
#               + immuno_supp
#               + carehome 
#               + bmi_imp1 + bmi_imp2 # splines in bmi 
#               + Q_DIAG_AF 
#               + Q_DIAG_ASTHMA 
#               + Q_DIAG_BLOOD_CANCER 
#               ### + Q_DIAG_CCF 
#               + Q_DIAG_CHD  
#               # LASSO removed+ Q_DIAG_CKD_LEVEL  
#               + Q_DIAG_COPD 
#               + Q_DIAG_DEMENTIA 
#               + Q_DIAG_DIABETES_1 
#               + Q_DIAG_DIABETES_2 
#               ###+ Q_DIAG_EPILEPSY  
#               ###+ Q_DIAG_FRACTURE 
#               # LASSO removed+ Q_DIAG_NEURO   
#               ###+ Q_DIAG_PARKINSONS 
#               ### + Q_DIAG_PULM_HYPER  
#               + Q_DIAG_PULM_RARE  
#               ###+ Q_DIAG_PVD   
#               + Q_DIAG_RA_SLE 
#               ###+ Q_DIAG_RESP_CANCER   
#               + Q_DIAG_SEV_MENT_ILL 
#               ###+ Q_DIAG_STROKE  
#               + Q_DIAG_VTE  
#               
#               # Prescriptions
#               + PIS_Lipid_regulating_drugs
#               + PIS_Angiotensin_converting_enzyme_inhibitors
#               + PIS_Beta_adrenoceptor_blocking_drugs
#               + PIS_Oral_iron
#               + PIS_Selective_serotonin_re_uptake_inhibitors
#               #+ PIS_Metformin_hydrochloride # combined with Q_DIAG_DIABETES_2
#               + PIS_Loratadine
#               ###+ PIS_Direct_oral_anticoagulants 
#               + PIS_Colchicine
#               # LASSO removed+ PIS_Antiplatelet_drugs
#               ###+ PIS_Alpha_adrenoceptor_blocking_drugs
#               + PIS_Macrolides
#               + PIS_Benzylpenicillin_and_phenoxymethylpenicillin
#               + PIS_Leukotriene_receptor_antagonists
#               + PIS_Herpes_simplex_and_varicella_zoster
#               ###+ PIS_Replacement_therapy 
#               + PIS_Famotidine
#               + PIS_Ursodeoxycholic_acid
#               + PIS_Systemic_nasal_decongestants
#               ###+ PIS_Warfarin_sodium 
#               + PIS_Parenteral_anticoagulants 
#               + PIS_Compound_bronchodilator_preparations
#               ###+ PIS_Ranitidine_hydrochloride 
#               
#               # Severity of acute infection
#               + severe,   
#               
#               # Model features 
#               data = df,
#               family = "quasibinomial", 
#               method = "glm",
#               weights = df$class_weight,
#               trControl = cv_results)
#  return(mod)
#}
#
#lasso_model <- lasso_model_fun()
#
## 2. Get predicted probabilities  ----
#df_testing$LASSO_pred <- predict(lasso_model, newdata = df_testing, type = "prob")[,2]

# 3. Sensitivity ----
## 8.1 Prepare data for plotting and get predicted probabilities ----
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
  
  # * Uncomment code below ----
  ## Run main model for predicted probabilities (LASSO with 10-fold CV, training data)
  #lasso_model_for_probs <- lasso_model_fun()
  
  ## Get predicted probabilities (testing data)
  #probs <- predict(lasso_model_0, newdata = df_testing, type = "prob")[,2]
  
  ## Remove model
  #rm(lasso_model_for_probs)
  
  ## Return predicted probabilities
  return(list(coef_df = coef_df))#, probs = probs))
}

## Depvar = clinical code, free text, long covid sick note
df$depvar <- df$code_txt_fitnote 
codes <- depvar_variations("Clinical code, free text, or sick note")
coef_df_codes <- codes$coef_df
#df_testing$LASSO_pred_codes <- codes$probs

## Depvar = opdef
df$depvar <- ifelse(df$opdef == 1, 1, 0) 
opdef <- depvar_variations("Operational definition")
coef_df_opdef <- opdef$coef_df
#df_testing$LASSO_pred_opdef <- opdef$probs

## Depvar = opdef without blood tests
df$depvar <- df$opdef_no_blood 
opdef_nb <- depvar_variations("Operational definition (excluding blood tests)")
coef_df_opdef_nb <- opdef_nb$coef_df
#df_testing$LASSO_pred_opdef_nb <- opdef_nb$probs

## Depvar = opdef without blood tests OR clinical code, free text, long covid sick note
df$depvar <- ifelse(df$opdef_no_blood ==1 | df$code_txt_fitnote == 1, 1, 0) 
opdef_nb_codes <- depvar_variations("Operational definition (excluding blood tests)\nor clinical code, free text, or sick note")
coef_df_opdef_nb_codes <- opdef_nb_codes$coef_df
#df_testing$LASSO_pred_opdef_nb_codes <- opdef_nb_codes$probs


## 8.2 Plot ---- 
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

setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling/LFT")
ggsave("Coefficient plot - comparing dependent variables.jpg", p, height = 11, width = 19, units = "in")


## 4. Export ----
### Analysis including those with positive PCR or LFT tests
#saveRDS(df, "/conf/EAVE/GPanalysis/analyses/long_covid/data/df_training_cleaned_LFT.rds")
#saveRDS(df_testing, "/conf/EAVE/GPanalysis/analyses/long_covid/data/df_testing_with_probs_LFT.rds")
