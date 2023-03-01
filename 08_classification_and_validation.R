#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 08_classification_and_validation - Classifies individuals as having long-COVID or not according to our 
## operatinal definition and checks for agreement with unambiguous indicators of long-COVID (Read code, free text, free 
## text on sick note). For sensitivity, repeats this analysis based on symptoms that were significantly higher at 12-26
## weeks only, and identifies cases as having long COVID based on symptoms recorded at 12-26 weeks only.
#########################################################################################################################

# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(ggplot2)
library(corrplot)

# 0. Set-up ------
# Read in df with counts of outcomes in each period for positive cases  
# NOTE: the raw data contains positive and negative cases only i.e. does not include 'not yet tested'
df_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/5. Cluster analysis/df_cluster_LFT.rds") %>% # counts of vars for PCR & LFT positive, prepared in 06_matched_analysis
  filter(tes_res == "POSITIVE") 

# Read in all variables that are significantly higher for positive cases
all_vars <- read.csv("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/2. Matched analysis/Significantly higher outcomes.csv") %>% 
  dplyr::select(x)

all_vars <- as.vector(all_vars$x)


# Sensitivity: Read in all variables that are significantly higher for positive cases @ >12-26 weeks only
vars_12 <- read.csv("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/2. Matched analysis/Significantly higher outcomes - 12-26.csv") %>% 
  dplyr::select(x)

vars_12 <- as.vector(vars_12$x)

# 2. Prepare df ----
## Update var names in df to match var labels used in results tables
df <- df_raw %>% 
  rename('Echocardiogram_4' = Ecg_4,
         'GP Interactions_4' = gp_n_4,
         'Inpatient Visits_4' = hos_4,
         'ICU Admissions_4' = icu_4,
         'Out of Hours_4' = ooh_4,
         'A&E Visits_4' = ae_4,
         'NHS 24 Calls_4' = n24_4,
         'Outpatient Visits_4' = op_4,
         'Echocardiogram_12' = Ecg_12,
         'GP Interactions_12' = gp_n_12,
         'Inpatient Visits_12' = hos_12,
         'ICU Admissions_12' = icu_12,
         'Out of Hours_12' = ooh_12,
         'A&E Visits_12' = ae_12,
         'NHS 24 Calls_12' = n24_12,
         'Outpatient Visits_12' = op_12)

colnames(df) <- gsub('Bloodtests', 'Blood Tests', colnames(df)) 
colnames(df) <- gsub('Consult Sickline', 'Sick Note', colnames(df)) 

## Get counts for each variable in 'all_vars' for 4-26 weeks (by summing _4 and _12 counts)
for(var in all_vars){
  df[[var]] <- df[[paste0(var,"_4")]] + df[[paste0(var,"_12")]]
}

## Sensitivity: Keep the 12-26 weeks counts + date variables 
df12 <- df %>% 
  dplyr::select(EAVE_LINKNO, tes_dat, all_of(paste0(vars_12, "_12")))

### Remove '_12' from colnames
colnames(df12) <- c("EAVE_LINKNO", "tes_dat", vars_12)


## Keep the 4-26 weeks counts + date variables (date is needed for plotting data in 09)
df <- df %>% 
  dplyr::select(EAVE_LINKNO, tes_dat, all_of(all_vars))


# 3. Classify individuals as having long-COVID or not according to our operational definition ----
## Classify variables according to our 3 categories
symptoms <- c("Fatigue", "Breathless", "Taste and Smell", "Consult Mental Health") 
investigations <- c(all_vars[grepl("Blood", all_vars)], "Chest Xray", "Echocardiogram") 
pres_sick <- c(setdiff(all_vars, c(symptoms, investigations)))

symptoms12 <- symptoms[symptoms != "Consult Mental Health"]
investigations12 <- investigations[investigations != "Blood Tests Biochem" & investigations != "Echocardiogram"] 
pres_sick12 <- pres_sick[pres_sick != "Corticosteroids (respiratory)" & pres_sick != "Tetracyclines" & pres_sick != "Coronavirus"]

## Check
length(symptoms) + length(investigations) + length(pres_sick) == length(all_vars)

length(symptoms12) + length(investigations12) + length(pres_sick12) == length(vars_12)


## Identify individuals who have 2 out of the following 3:
## - a symptom 
## - an investigation (blood test or other)
## - a prescription or sick note

# Create binary variables to signal whether each individual has each type of outcome (symptom, investigation, pres_sick) 
count_outcome_types <- function(df, outcome_type){ # e.g. outcome_type = symptoms
  df$n <- 0 # create empty column
  for(outcome in outcome_type){
    df$n <- ifelse(df[[outcome]] > 0, df$n + 1, df$n) # +1 for each outcome that is > 0
    df$n <- ifelse(df$n > 0, 1, 0) # make binary
  }
  
  print("% of individuals with >0 outcomes")
  print(round(prop.table(table(df$n))*100,1))
  return(df$n)
}

df$symptoms <- count_outcome_types(df, symptoms)
df$investigations <- count_outcome_types(df, investigations)
df$pres_sick <- count_outcome_types(df, pres_sick)

df12$symptoms <- count_outcome_types(df12, symptoms12)
df12$investigations <- count_outcome_types(df12, investigations12)
df12$pres_sick <- count_outcome_types(df12, pres_sick12)

## Identify individuals with 2 out of 3 of: symptoms, investigations, prescription/sick note
df$indicators <- df$symptoms + df$investigations + df$pres_sick
df$op_def <- ifelse(df$indicators >1, 1, 0)
table(df$op_def)
round(prop.table(table(df$op_def))*100,1) # 5.1% of positive cases

df12$indicators <- df12$symptoms + df$investigations + df$pres_sick
df12$op_def12 <- ifelse(df12$indicators >1, 1, 0)
table(df12$op_def12)
round(prop.table(table(df12$op_def12))*100,1) # 3.9% of positive cases

## Join op_def12 to df
df <- df %>% 
  left_join(dplyr::select(df12, EAVE_LINKNO, op_def12), by = "EAVE_LINKNO")

## Check overlap: 55210 of (55210+18465) = 74.9% of cases identified using only those symptoms that are sig higher at 12-26 weeks AND recorded at 12-26 weeks
xtab <- table(df$op_def, df$op_def12)
xtab
xtab[2,2]/(xtab[2,1] + xtab[2,2])

# 4. Export data ----
op_def_df <- df %>% 
  filter(op_def == 1) %>% 
  dplyr::select(EAVE_LINKNO, tes_dat, op_def) 

write_rds(op_def_df, "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/6. Long COVID indicators/opdef_classification.rds")
