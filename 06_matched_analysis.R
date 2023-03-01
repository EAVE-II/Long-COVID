#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 06_matched_analysis - For each period, incorporates counts of Read codes and other health outcomes for 
## matched pairs of individuals (pos vs neg and pos vs not yet tested). Estimates differences in health outcomes of those 
## who have and have not tested positive for COVID-19 using matched data.
#########################################################################################################################

# 0. Set-up ------

# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(sandwich) 
library(splines)


# Set chemical substances of interest -----
chemicals_of_interest <- c("Edoxaban") # We decided not to include analysis of chemicals in the 1st paper, but may want to run analysis on 
                                       # chemicals in the next stage of research, so leaving in a minimal version for now


# FUNCTIONS ----

# 1. Prepare GP data  ----
## Gets counts of Read codes and other health outcomes for each individual and time period
## NOTE: This step is performed after matching because stu_12w/stu_end is updated after matching
prepare_GP_data <- function(df, cluster){ # e.g. df = df_neg or df = df_nt
  # e.g. cluster = T if creating a df to be used for cluster analysis
  
  ## Read in GP data supplied by Albasoft
  GP_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/longcovid_261022.rds") 
  
  ## Clean up GP data
  # Format date variable
  GP_raw$EventDate <- ymd(GP_raw$EventDate)
  
  GP <- GP_raw %>% 
    # Keep EAVE_LINKNOs in df
    filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
    # Remove Read codes that are NA
    filter(!is.na(group)) %>% 
    # Select columns of interest
    select(EAVE_LINKNO, gp_read = group, gp_dat = EventDate) %>%
    # Incorporate stu_sta, stu_12w, stu_end for each matched pair 
    left_join(select(df, EAVE_LINKNO, stu_sta, stu_12w, stu_end), by = "EAVE_LINKNO") %>%
    # Remove Read codes recorded outside of study period
    filter(gp_dat >= stu_sta & gp_dat <= stu_end) 
  
  
  ## Count the number of times each Read code is recorded at 4-12 weeks ('_4') and 12-26 weeks ('_12') after testing
  GP_counts <- GP %>%
    # Identify Read codes at 4-12 weeks '_4' and >12 weeks '_12'
    mutate(gp_read_phase = ifelse(gp_dat <= stu_12w, paste(gp_read, "_4", sep = ""), # 4-12 weeks
                                  ifelse(gp_dat > stu_12w, paste(gp_read, "_12", sep =""), NA))) %>% # >12-26 weeks
    # Add counts of Read codes by patient and period
    group_by(EAVE_LINKNO, gp_read_phase) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(n = as.numeric(n)) %>%
    ## Transpose
    pivot_wider(
      names_from = gp_read_phase,
      values_from = n,
      values_fill = list(n = 0))
  
  ## Create a sick note variable for sick notes that don't mention long-COVID
  GP_counts <- GP_counts %>% 
    mutate(`Sick Note (except long-COVID)_4` = `Consult Sickline_4` - `Free text - fitnote_4`,
           `Sick Note (except long-COVID)_12` = `Consult Sickline_12` - `Free text - fitnote_12`)
  
  
  # Remove the free text variables when preparing data for regressions (but keep for cluster analysis)
  if(cluster == F){ 
    GP_counts <- GP_counts %>% 
      select(-c("Free text_4", "Free text - fitnote_4", "Free text_12", "Free text - fitnote_12"))
  }
  
  ## Join GP data to covariate data 
  df <- df %>% 
    left_join(GP_counts, by = "EAVE_LINKNO") %>% 
    mutate_if(is.numeric,coalesce,0) %>% # fill NAs with 0s
    mutate(n_risk_gps = ifelse(is.na(n_risk_gps), 0, n_risk_gps)) # fill NAs with 0s (non-numeric variable)
  
  ## Get gp_n (counts of GP interactions on unique days)
  ## NB does not capture GP interactions where Read codes not on our list are recorded. 
  ### 4-12 weeks
  gp_n_4 <- GP %>%
    filter(gp_dat <= stu_12w) %>%
    distinct(EAVE_LINKNO, gp_dat, .keep_all = T) %>% 
    add_count(EAVE_LINKNO, name = "gp_n_4") %>% 
    select(EAVE_LINKNO, gp_n_4) %>% 
    filter(!duplicated(EAVE_LINKNO))
  
  ### >12 - 26 weeks
  gp_n_12 <- GP %>%
    filter(gp_dat > stu_12w) %>%
    distinct(EAVE_LINKNO, gp_dat, .keep_all = T) %>% 
    add_count(EAVE_LINKNO, name = "gp_n_12") %>% 
    select(EAVE_LINKNO, gp_n_12) %>% 
    filter(!duplicated(EAVE_LINKNO))
  
  ### Join gp_n_4 and gp_n_12 to df
  df  <- df %>% 
    left_join(gp_n_4, by="EAVE_LINKNO") %>% 
    mutate(gp_n_4 = ifelse(is.na(gp_n_4), 0L, gp_n_4)) %>% 
    left_join(gp_n_12, by="EAVE_LINKNO") %>% 
    mutate(gp_n_12 = ifelse(is.na(gp_n_12), 0L, gp_n_12))
  
  return(df)
}


# 2. Prepare counts of other healthcare interactions & follow up periods  ----
prepare_interactions_data <- function(df, cluster){ # e.g. df = df_neg or df = df_nt
                                                    # cluster = T or F depending on whether or not a df for use in the cluster analysis is being prepared
  
  ## Capture test dates from df to use with inside function
  tes_dates <- df %>% select(EAVE_LINKNO, tes_dat, stu_sta, stu_12w, stu_end)
  
  ## Function to count number of times patients interact with different parts of the healthcare system ----
  count_interactions <- function(raw_data,       # the datafile e.g. hos
                                 date_variable) {# the variable containing the date of the interaction e.g. hos$ADMISSION_DATE
    
    ## Create variables
    dat <- raw_data %>%
      filter(!is.na(EAVE_LINKNO)) %>%
      select(EAVE_LINKNO, {{date_variable}})%>%
      # Incorporate test dates
      left_join(tes_dates, by = "EAVE_LINKNO") %>%
      # Remove interactions after stu_end
      filter({{date_variable}} <= stu_end) %>% 
      # Interactions within 4-12 and 12-26 weeks following PCR test (binary - used to create count)
      mutate(count_4 = ifelse({{date_variable}} >= stu_sta & {{date_variable}} <= stu_12w, 1, 0),
             count_12 = ifelse({{date_variable}} > stu_12w & {{date_variable}} <= stu_end, 1, 0))
    
    ## Translate binary to count - 4-12 weeks
    count_4 <- dat %>%
      filter(count_4 == 1) %>%
      add_count(EAVE_LINKNO, name = "var_4") %>%
      filter(!duplicated(EAVE_LINKNO)) %>%
      select(EAVE_LINKNO, var_4)
    
    ## Translate binary to count - >12-26 weeks
    count_12 <- dat %>%
      filter(count_12 == 1) %>%
      add_count(EAVE_LINKNO, name = "var_12") %>%
      filter(!duplicated(EAVE_LINKNO)) %>%
      select(EAVE_LINKNO, var_12)
    
    ## Link to df
    df <- df %>%
      left_join(count_4, by = "EAVE_LINKNO") %>%
      left_join(count_12, by = "EAVE_LINKNO") %>%
      mutate(var_4 = ifelse(is.na(var_4), 0L, var_4),
             var_12 = ifelse(is.na(var_12), 0L, var_12))
    
    return (df)
  }
  
  
  # 1.2.1 Hospitalisations - SMR01 data (in-patients) ----
  df <- count_interactions(raw_data = read_rds("/conf/EAVE/GPanalysis/analyses/long_covid/data/SMR01 - 2019-03-01 to 2022-09-21.rds"), 
                           date_variable = hos_dat)
  
  df <- df %>% rename(hos_4 = var_4, hos_12 = var_12)
  
  
  # 1.2.2 Intensive care admissions - SICSAG data ----
  df <- count_interactions(raw_data = read_rds("/conf/EAVE/GPanalysis/data/SICSAG_episode_level_.rds"), 
                           date_variable = AdmitUnit)
  
  df <- df %>% rename(icu_4 = var_4, icu_12 = var_12)
  
  
  # 1.2.3 Out of hours (OOH data mart) ---------
  
  ## Prepare OOH data
  ooh <- read_rds("/conf/EAVE/GPanalysis/analyses/long_covid/data/OOH_long_covid_20221209.rds") %>% 
    mutate(ooh_dat = date(contact_start_date_time)) %>% 
    # Remove multiple entries per individual per day i.e. remove duplicated (EAVE_LINKNO AND ooh_dat)
    distinct(EAVE_LINKNO, ooh_dat, .keep_all = TRUE)
  
  df <- count_interactions(raw_data = ooh, 
                           date_variable = ooh_dat)
  
  df <- df %>% rename(ooh_4 = var_4, ooh_12 = var_12)
  
  rm(ooh)
  
  # 1.2.4 Accident and emergency  ---------
  ## Prepare A&E data
  ae <- read_rds("/conf/EAVE/GPanalysis/analyses/long_covid/data/A&E_long_covid_20221209.rds") %>% 
    mutate(ae_dat = date(arrival_date_time)) %>% 
    # Remove multiple entries per individual per day i.e. remove duplicated (EAVE_LINKNO AND ooh_dat)
    distinct(EAVE_LINKNO, ae_dat, .keep_all = TRUE)
  
  df <- count_interactions(raw_data = ae, date_variable = ae_dat)
  
  df <- df %>% rename(ae_4 = var_4, ae_12 = var_12)
  
  rm(ae)
  
  
  # 1.2.5 NHS 24 -----
  ## Prepare NHS 24 data
  nhs24 <- read_rds("/conf/EAVE/GPanalysis/analyses/long_covid/data/NHS24_long_covid_20221209.rds") %>% 
    mutate(nhs24_dat = date(call_received_date_time)) %>% 
    # Remove multiple entries per individual per day i.e. remove duplicated (EAVE_LINKNO AND ooh_dat)
    distinct(EAVE_LINKNO, nhs24_dat, .keep_all = TRUE)
  
  df <- count_interactions(raw_data = nhs24, date_variable = nhs24_dat)
  
  df <- df %>% rename(n24_4 = var_4, n24_12 = var_12)
  
  rm(nhs24)
  
  
  # 1.2.6 Outpatient (SMR00) -----
  ## Prepare outpatient data
  op <- read_rds("/conf/EAVE/GPanalysis/analyses/long_covid/data/SMR00_long_covid_20221209.rds") %>% 
    mutate(op_dat = date(clinic_date)) %>% 
    # Remove multiple entries per individual per day i.e. remove duplicated (EAVE_LINKNO AND ooh_dat)
    distinct(EAVE_LINKNO, op_dat, .keep_all = TRUE)
  
  df <- count_interactions(raw_data = op, date_variable = op_dat)
  
  df <- df %>% rename(op_4 = var_4, op_12 = var_12)
  
  rm(op)
  
  rm(tes_dates)
  
  
  # 1.2.7 Follow up period ----
  # Add days of follow up for each period (to account for incomplete follow up within each matched pair)
  df <- df %>%
    ## Number of days of follow-up available at the end of 4-12 weeks (to account for incomplete follow up)
    mutate(follow_up_4 = ifelse(stu_12w < stu_end, stu_12w - stu_sta, stu_end - stu_sta),
           ## Number of days of follow-up available at the end of >12-26 weeks (to account for incomplete follow up)
           follow_up_12 = stu_end - stu_sta,
           ## Formatting
           follow_up_4 = as.numeric(follow_up_4),
           follow_up_12 = as.numeric(follow_up_12))
  
  # Check proportion where days of follow up for the pair is < 1
  prop.table(table(df$follow_up_4 < 0 | df$follow_up_12 < 1))[2] # 0%
  
  # Keep pairs with > 0 days follow up (i.e. data for 1+ days beyond the 4 weeks following their test date)
  df <- df %>% 
    filter(follow_up_4 > 0)
  
  
  # 1.2.8 Update vac_n_clean, hos_prior, icu_prior to match new tes_dat (covariates in regression) ----
  ## For controls, these were calculated (during matching) based on a tes_date selected at random within the month propensity scores 
  ## were estimated, but now controls' tes_dat has been updated to match their positive matches - so re-calculate these values
  
  ## Skip this step when prpearing datasets for cluster analysis (which don't have updated tes_dats, because they are prepared from the unmatched dfs)
  if(cluster == F){
    # Read in and clean up data  
    # Hospitalisations (1/3/2019 to latest test date)
    hos <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/SMR01 - 2019-03-01 to 2022-09-21.rds")
    
    # ICU admissions (1/3/19 to latest test date)
    icu <- readRDS("/conf/EAVE/GPanalysis/data/SICSAG_episode_level_.rds") %>% 
      filter(!is.na(EAVE_LINKNO)) %>% 
      select(EAVE_LINKNO, icu_dat = AdmitUnit) %>% 
      filter(icu_dat <= max(df$tes_dat))
    
    # Vaccinations
    vac_raw <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned.rds") 
    
    ## Clean up vaccinations data
    vac <- vac_raw %>% 
      # Select variables of interest
      select(EAVE_LINKNO, 
             vac_dos = vacc_dose_number,
             vac_dat = vacc_occurence_date, 
             vac_typ = vacc_product_name) %>%
      # Keep only vaccinations administered after 1st December 2020 (when rollout began) and 14 days before max test date
      filter(vac_dat >= "2020-12-01" & vac_dat <= (max(df$tes_dat - days(14)))) %>% 
      # Rename vaccine types
      mutate(vac_typ = case_when(vac_typ == "Covid-19 Vaccine AstraZeneca" ~ "AZ",
                                 vac_typ == "Covid-19 mRNA Vaccine Pfizer" ~ "PB",
                                 vac_typ == "Covid-19 mRNA Vaccine Moderna" ~ "Mo",
                                 TRUE ~ "Unknown")) %>% 
      # Remove duplicates (by EAVE_LINKNO --AND-- dose number)
      group_by(vac_dos) %>% 
      filter(!duplicated(EAVE_LINKNO)) %>% 
      ungroup %>% 
      # Remove rows with missing info
      filter(!is.na(vac_dos) & !is.na(vac_dat) & !is.na(vac_typ)) %>% 
      # Transpose (one row per individual)
      pivot_wider(id_cols = EAVE_LINKNO,
                  names_from = vac_dos,
                  values_from = c(vac_dat, vac_typ)) %>% 
      mutate(vac_incon = 0)
    
    if("vac_dat_3" %in% colnames(vac)){
      ## Flag vaccinations 1 & 2 or 2 & 3 <18 days apart
      vac <- vac %>% 
        mutate(vac_gap_1 = as.numeric(date(vac_dat_2) - date(vac_dat_1)), 
             vac_gap_2 = as.numeric(date(vac_dat_3) - date(vac_dat_2)),
             vac_incon = ifelse(vac_gap_1 <19 & !is.na(vac_gap_1), 1, 0))
    }
    
    if("vac_typ_4" %in% colnames(vac)){  
      vac <- vac %>% 
    ## Flag first dose ==PB or Mo, but second dose == AZ
      mutate(vac_incon = ifelse((vac_typ_1 == "PB" | vac_typ_1 == "Mo") & vac_typ_2 == "AZ", 1, 0),
             ## Flag no record of earlier doses
             vac_incon = ifelse((is.na(vac_typ_1) & (!is.na(vac_typ_2) | !is.na(vac_typ_3) | !is.na(vac_typ_4))) |    # No first dose, but subsequent doses
                                  ((is.na(vac_typ_1) | is.na(vac_typ_2)) & (!is.na(vac_typ_3) | !is.na(vac_typ_4))) |    # No first or second dose, but subsequent doses
                                  ((is.na(vac_typ_1) | is.na(vac_typ_2) | is.na(vac_typ_3)) & !is.na(vac_typ_4)), 1, 0))   # No first, second, or third dose, but subsequent doses
    }
      
    vac <- vac %>% 
      ## Keep data with no inconsistencies
      filter(vac_incon == 0) 
    
    ## Get max vac_dat column
    max_vac_dat_col <- colnames(vac)[which(colnames(vac) == "vac_typ_1")-1]
    
    vac <- vac %>% 
      ## Transpose to long form 
      pivot_longer(cols = c(vac_dat_1:max_vac_dat_col), # columns that should pivot from wide to long (unquoted)
                   names_to = "vac_dos", # name of the new category column as a quoted string
                   values_to = "vac_dat") # name of the new value column as a quoted string)
    
    vac <- vac %>% 
      filter(!is.na(vac_dat)) %>% 
      select(EAVE_LINKNO, vac_dat)
    
    rm(vac_raw)
    
    
    ## Prepare counts
    ## Hospitalisations in the year before tes_dat
    hos_n <- hos %>% 
      filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
      left_join(select(df, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
      filter(hos_dat <= tes_dat & hos_dat > (tes_dat - days(365))) %>% 
      group_by(EAVE_LINKNO) %>% 
      summarise(hos_prior = n())
    
    ## ICU admissions in the year before tes_dat
    icu_n <- icu %>% 
      filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
      left_join(select(df, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
      filter(icu_dat <= tes_dat & icu_dat > (tes_dat - days(365))) %>% 
      group_by(EAVE_LINKNO) %>% 
      summarise(icu_prior = n())
    
    # Vaccine doses by 14 days before tes_dat
    vac_n <- vac %>% 
      filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
      left_join(select(df, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
      filter(vac_dat <= (tes_dat - days(14))) %>% 
      group_by(EAVE_LINKNO) %>% 
      summarise(vac_n_clean = n())
    
    # Replace in df
    df <- df %>% 
      select(-c(hos_prior, icu_prior, vac_n_clean)) %>% 
      left_join(hos_n, by = "EAVE_LINKNO") %>%
      left_join(icu_n, by = "EAVE_LINKNO") %>% 
      left_join(vac_n, by = "EAVE_LINKNO") %>% 
      mutate(hos_prior = ifelse(is.na(hos_prior), 0L, hos_prior),
             icu_prior = ifelse(is.na(icu_prior), 0L, icu_prior),
             vac_n_clean= ifelse(is.na(vac_n_clean), 0L, vac_n_clean))
  }
  
  
  return(df)
  
}


# 3. Prepare prescriptions data (PIS) ----
prepare_pis <- function(df){ # e.g. df = df_neg or df = df_nt
  
  ## Read in PIS data (prescriptions)
  pis_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/PIS_07122022.rds")
  
  ## Format variable names to match BNF nomenclature 
  pis <- pis_raw %>% 
    rename(dispensed_full_date = `Disp Date`) %>% 
    rename(bnf_code = `PI BNF Item Code`) %>% 
    rename(bnf_sub_paragraph_code = `PI BNF Paragraph Code`) %>% 
    select(EAVE_LINKNO, dispensed_full_date, bnf_sub_paragraph_code, bnf_code)
  
  ## Get chemical substances codes
  pis$bnf_chem_code <- substr(pis$bnf_code,1,9)
  
  ## Lookup sub_paragraph and chemical substances names
  ## Read in lookup file and format to match pis
  bnf <- read_csv("/conf/EAVE/GPanalysis/analyses/long_covid/data/BNF lookup.csv") %>% 
    select(`BNF Subparagraph`,	`BNF Subparagraph Code`,	`BNF Chemical Substance`,	`BNF Chemical Substance Code`) %>% 
    # Rename cols to match pis
    rename(bnf_sub_paragraph_name = `BNF Subparagraph`) %>% 
    rename(bnf_sub_paragraph_code = `BNF Subparagraph Code`) %>% 
    rename(bnf_chem_name = `BNF Chemical Substance`) %>% 
    rename(bnf_chem_code = `BNF Chemical Substance Code`) %>% 
    # Add leading 0s to codes
    mutate(bnf_sub_paragraph_code = paste0(0, bnf_sub_paragraph_code))
  
  ## Integrate into pis
  bnf_sub <- bnf %>% 
    filter(!duplicated(bnf_sub_paragraph_code)) %>% 
    filter(bnf_sub_paragraph_code %in% pis$bnf_sub_paragraph_code) %>% 
    select(bnf_sub_paragraph_code, bnf_sub_paragraph_name)
  
  bnf_chem <- bnf %>%
    filter(!duplicated(bnf_chem_code)) %>% 
    filter(bnf_chem_code %in% pis$bnf_chem_code) %>% 
    select(bnf_chem_code, bnf_chem_name)
  
  rm(bnf)
  
  pis <- pis %>% 
    left_join(bnf_sub, by = "bnf_sub_paragraph_code") %>% 
    left_join(bnf_chem, by = "bnf_chem_code")
  
  ### Prepare test dates, study start, 12weeks, study end to incorporate
  dates <- df %>% 
    select(EAVE_LINKNO, tes_dat, stu_sta, stu_12w, stu_end)
  
  ## Rename main pis variable (sub-paragraph name)
  pis <- pis %>% 
    rename(bnf = bnf_sub_paragraph_name)
  
  # Identify prescriptions/chemical substances issued 12 months prior to test ('_prior'), 4-12 weeks after ('_4') and >12 weeks after ('_12')
  pis <- pis %>% 
    # Keep individuals in df only
    filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
    ## Rename date variable
    mutate(bnf_dat = dispensed_full_date) %>% 
    ## Incorporate dates
    left_join(dates, by = "EAVE_LINKNO") %>% 
    ## Keep only prescriptions issued 12 months before testing up to stu_end
    filter(bnf_dat > (tes_dat-365) & bnf_dat <= stu_end) %>%  
    # Identify when prescriptions were dispensed, relative to test date
    mutate(bnf_phase = ifelse(bnf_dat < tes_dat, paste(bnf, "_prior", sep=""), # 12 months prior
                              ifelse(bnf_dat >= stu_sta & bnf_dat <= stu_12w, paste(bnf, "_4", sep=""), # 4-12 weeks
                                     paste(bnf, "_12", sep=""))),  # >12-26 weeks
           # Identify when chemical substances were dispensed, relative to test date
           bnf_chem_phase = ifelse(bnf_dat < tes_dat, paste(bnf_chem_name, "_prior", sep=""), # 12 months prior
                                   ifelse(bnf_dat >= stu_sta & bnf_dat <= stu_12w, paste(bnf_chem_name, "_4", sep=""), # 4-12 weeks
                                          paste(bnf_chem_name, "_12", sep="")))) # >12-26 weeks
  
  
  ## Get chemical substances (processed later)
  pis_chems <- pis %>% 
    select(EAVE_LINKNO, bnf_chem_phase)
  
  
  ## Transpose data to wide form (one row per individual with counts by period)
  pis <- pis %>% 
    select(EAVE_LINKNO, bnf_phase) %>% 
    # Count prescriptions in each period
    add_count(EAVE_LINKNO, bnf_phase, name = "bnf_n") %>% 
    distinct(EAVE_LINKNO, bnf_phase, .keep_all = T) %>% 
    filter(!is.na(bnf_phase)) %>% 
    # Transpose (one column per bnf paragraph for each period, one row per EAVE_LINKNO)
    pivot_wider(names_from = bnf_phase, values_from = bnf_n, values_fill = list(bnf_n = 0)) %>% 
    # Rename "Replacement therapy" (give it a more instructive name)
    rename(`Fludrocortisone acetate_prior` = `Replacement therapy_prior`,
           `Fludrocortisone acetate_4` = `Replacement therapy_4`,
           `Fludrocortisone acetate_12` = `Replacement therapy_12`)
  
  
  ## Filter out prescriptions that are not new (i.e. remove prescriptions if they were also dispensed in the 12 months prior to testing)
  ### Get column names
  pis_cols <- c(colnames(pis))[-1]
  
  ### Get lists of prescriptions for each time period
  pis_4 <- pis_cols[grepl("_4", pis_cols)]
  pis_12 <- pis_cols[grepl("_12", pis_cols)]
  pis_prior <- pis_cols[grepl("_prior", pis_cols)]
  
  ### Filter out non-new prescriptions at 4-12 weeks 
  for(i in 1:length(pis_4)){
    var_4 <- pis_4[i]
    var_prior <- paste0(gsub('.{2}$', '', var_4), "_prior")
    
    if(var_prior %in% pis_prior){
      # If the prescription was dispensed in the previous 12 months (i.e. prescription_prior > 0), 
      # then don't count prescriptions issued in the follow up periods (i.e. make prescription_4 = 0 and prescription_12 = 0)  
      pis <- pis %>% 
        mutate(!!var_4 := ifelse(!!as.name(var_prior) > 0, 0, !!as.name(var_4)))
    }
  }
  
  
  ### Filter out non-new prescriptions at 12-26 weeks 
  for(i in 1:length(pis_12)){
    var_12 <- pis_12[i]
    var_prior <- paste0(gsub('.{3}$', '', var_12), "_prior")
    
    if(var_prior %in% pis_prior){
      # If the prescription was dispensed in the previous 12 months (i.e. prescription_prior > 0), 
      # then don't count prescriptions issued in the follow up periods (i.e. make prescription_4 = 0 and prescription_12 = 0)  
      pis <- pis %>% 
        mutate(!!var_12 := ifelse(!!as.name(var_prior) > 0, 0, !!as.name(var_12)))
    }
  }
  
  ## Join to df
  df <- df %>% 
    left_join(pis, by = "EAVE_LINKNO")
  
  ## Fill NAs with 0s
  df[, pis_cols][is.na(df[, pis_cols])] <- 0
  
  
  # Process chemical substances code
  ## Only keep chemical substances we're interested in
  chems <- c(paste0(chemicals_of_interest, "_prior"), paste0(chemicals_of_interest, "_4"), paste0(chemicals_of_interest, "_12"))
  
  pis_chems <- pis_chems %>% 
    mutate(bnf_chem_phase = ifelse(bnf_chem_phase %in% chems, bnf_chem_phase, NA)) %>% 
    filter(!is.na(bnf_chem_phase))
  
  ## Transpose chemical substances data to wide form (one row per individual with counts by period)
  pis_chems <- pis_chems %>% 
    select(EAVE_LINKNO, bnf_chem_phase) %>% 
    # Count prescriptions in each period
    add_count(EAVE_LINKNO, bnf_chem_phase, name = "bnf_chem_n") %>% 
    distinct(EAVE_LINKNO, bnf_chem_phase, .keep_all = T) %>% 
    filter(!is.na(bnf_chem_phase)) %>% 
    # Transpose (one column per bnf chemical substance for each period, one row per EAVE_LINKNO)
    pivot_wider(names_from = bnf_chem_phase, values_from = bnf_chem_n, values_fill = list(bnf_chem_n = 0))
  
  ## Filter out prescriptions that are not new (i.e. remove prescriptions if they were also dispensed in the 12 months prior to testing)
  ### Get column names
  pis_chem_cols <- c(colnames(pis_chems))[-1]
  
  ### Get lists of chemical substances for each time period
  chems_4 <- pis_chem_cols[grepl("_4", pis_chem_cols)]
  chems_12 <- pis_chem_cols[grepl("_12", pis_chem_cols)]
  chems_prior <- pis_chem_cols[grepl("_prior", pis_chem_cols)]
  
  ### Filter out non-new prescriptions at 4-12 weeks 
  for(i in 1:length(chems_4)){
    
    var_4 <- chems_4[i]
    var_prior <- paste0(gsub('.{2}$', '', var_4), "_prior")
    
    if(var_prior %in% chems_prior){
      # If the prescription was dispensed in the previous 12 months (i.e. prescription_prior > 0), 
      # then don't count prescriptions issued in the follow up periods (i.e. make prescription_4 = 0 and prescription_12 = 0)  
      pis_chems <- pis_chems %>% 
        mutate(!!var_4 := ifelse(!!as.name(var_prior) > 0, 0, !!as.name(var_4)))
    }
  }
  
  ### Filter out non-new prescriptions at >12-26 weeks 
  for(i in 1:length(chems_12)){
    var_12 <- chems_12[i]
    var_prior <- paste0(gsub('.{3}$', '', var_12), "_prior")
    
    if(var_prior %in% chems_prior){
      
      # If the prescription was dispensed in the previous 12 months (i.e. prescription_prior > 0), 
      # then don't count prescriptions issued in the follow up periods (i.e. make prescription_4 = 0 and prescription_12 = 0)  
      pis_chems <- pis_chems %>% 
        mutate(!!var_12 := ifelse(!!as.name(var_prior) > 0, 0, !!as.name(var_12)))
    }
  }
  
  ## Join to df
  df <- df %>% 
    left_join(pis_chems, by = "EAVE_LINKNO")
  
  ## Fill NAs with 0s
  df[, pis_chem_cols][is.na(df[, pis_chem_cols])] <- 0
  
  
  ## Return df with counts of prescriptions incorporated
  return(list("GP" = df, "pis" = pis, "chems" = pis_chems))
}


# 4. Dependent variables: Function to set dependent variables to be included in matched analysis for a given df ----- 
# Removes any variables where there are <5 non-zero values in the treated OR control group
set_dep_vars <- function(gp_df, pis_df, chems_df, comparison) {# e.g. gp_df = df$GP
  #      pis_df = df$pis
  #      chems_df = df$chems
  #      comparison = "pos vs neg"
  
  ## Get n of non-zero cases among treated (pos)
  nonzero_pos <- as.data.frame(colSums(gp_df != 0 & gp_df$pos_or_not == 1))
  colnames(nonzero_pos) <- "non-zeros"
  nonzero_pos <- tibble::rownames_to_column(nonzero_pos, "variable")
  
  ## Get n of non-zero cases among controls (neg / nt)
  nonzero_con <- as.data.frame(colSums(gp_df != 0 & gp_df$pos_or_not == 0))
  colnames(nonzero_con) <- "non-zeros"
  nonzero_con <- tibble::rownames_to_column(nonzero_con, "variable")
  
  ## Export counts of each variable
  var_counts <- merge(nonzero_pos, nonzero_con, by = "variable", all = TRUE) # full join
  
  var_counts <- var_counts[grep("_4|_12", var_counts$variable),] # keep only those rows with _4 or _12
  
  var_counts <- var_counts %>% 
    rename(Positive = `non-zeros.x`,
           Control = `non-zeros.y`) %>% 
    mutate(Total = Positive + Control)
  
  write_csv(var_counts, paste0("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/2. Matched analysis/Variable counts - ", comparison, ".csv"))
  
  ## Keep only dependent variables with >5 non-zero values in (i) treated & (ii) control groups 
  five_plus_pos <- nonzero_pos[nonzero_pos$`non-zeros`> 5,]
  five_plus_con <- nonzero_con[nonzero_con$`non-zeros`> 5,]
  
  
  # Get column names relating to Read codes and healthcare interactions 
  ## Select relevant variables
  start <- which(names(gp_df) == "weights") + 2
  end <- which(names(gp_df) == "follow_up_4") - 1
  
  vars <- colnames(gp_df[, start:end])
  vars_4 <- vars[grep("_4", vars)]
  vars_12 <- vars[grep("_12", vars)]
  
  ## Keep only variables with >5 non-zero values in treated OR control group
  vars_4 <- vars_4[vars_4 %in% five_plus_pos$variable | vars_4 %in% five_plus_con$variable]
  vars_12 <- vars_12[vars_12 %in% five_plus_pos$variable | vars_12 %in% five_plus_con$variable]
  
  
  # Chemical substances
  chems_vars <- colnames(chems_df[, 2:ncol(chems_df)])
  chems_4  <- chems_vars[grep("_4", chems_vars)]
  chems_12 <- chems_vars[grep("_12", chems_vars)]
  chems_prior <- chems_vars[grep("_prior", chems_vars)]
  
  ## Keep only variables with >5 non-zero values in treated OR control group
  chems_4 <- chems_4[chems_4 %in% five_plus_pos$variable | chems_4 %in% five_plus_con$variable]
  chems_12 <- chems_12[chems_12 %in% five_plus_pos$variable | chems_12 %in% five_plus_con$variable]
  chems_prior <- chems_prior[chems_prior %in% five_plus_pos$variable | chems_prior %in% five_plus_con$variable]
  
  # Prescriptions
  pis_vars <- colnames(pis_df[, 2:ncol(pis_df)])
  pis_4  <- pis_vars[grep("_4", pis_vars)]
  pis_12 <- pis_vars[grep("_12", pis_vars)]
  pis_prior <- pis_vars[grep("_prior", pis_vars)]
  
  ## Keep only variables with >5 non-zero values in treated OR control group
  pis_4 <- pis_4[pis_4 %in% five_plus_pos$variable | pis_4 %in% five_plus_con$variable]
  pis_12 <- pis_12[pis_12 %in% five_plus_pos$variable | pis_12 %in% five_plus_con$variable]
  pis_prior <- pis_prior[pis_prior %in% five_plus_pos$variable | pis_prior %in% five_plus_con$variable]
  
  
  # Create labels----
  labs_4 <- gsub('_4', '', vars_4) # variable names with "_4" removed
  labs_4 <- gsub('_', ' ', labs_4) # variable names with "_" removed
  labs_4 <- gsub('Ecg', 'Echocardiogram', labs_4) 
  labs_4 <- gsub('gp n', 'GP Interactions', labs_4)
  labs_4 <- gsub('hos', 'Inpatient Visits', labs_4)
  labs_4 <- gsub('icu', 'ICU Admissions', labs_4)
  labs_4 <- gsub('ooh', 'Out of Hours', labs_4)
  labs_4 <- gsub('ae', 'A&E Visits', labs_4)
  labs_4 <- gsub('n24', 'NHS 24 Calls', labs_4)
  labs_4 <- gsub('op', 'Outpatient Visits', labs_4)
  labs_4 <- gsub('Bloodtests HA&E Visitsmatology', 'Bloodtests Haematology', labs_4)
  labs_4 <- gsub('Bloodtests', 'Blood Tests', labs_4)
  labs_4 <- gsub('Consult Sickline', 'Sick Note', labs_4)
  
  labs_12 <- gsub('_12', '', vars_12)
  labs_12 <- gsub('_', ' ', labs_12)
  labs_12 <- gsub('Ecg', 'Echocardiogram', labs_12) 
  labs_12 <- gsub('gp n', 'GP Interactions', labs_12)
  labs_12 <- gsub('hos', 'Inpatient Visits', labs_12)
  labs_12 <- gsub('icu', 'ICU Admissions', labs_12)
  labs_12 <- gsub('ooh', 'Out of Hours', labs_12)
  labs_12 <- gsub('ae', 'A&E Visits', labs_12)
  labs_12 <- gsub('n24', 'NHS 24 Calls', labs_12)
  labs_12 <- gsub('op', 'Outpatient Visits', labs_12)
  labs_12 <- gsub('Bloodtests HA&E Visitsmatology', 'Bloodtests Haematology', labs_12)
  labs_12 <- gsub('Bloodtests', 'Blood Tests', labs_12)
  labs_12 <- gsub('Consult Sickline', 'Sick Note', labs_12)
  
  labs_chems_4 <- gsub('_4', '', chems_4)
  labs_chems_12 <- gsub('_12', '', chems_12) 
  labs_chems_prior <- gsub('_prior', '', chems_prior)
  
  labs_pis_4 <- gsub('_4', '', pis_4)
  labs_pis_12 <- gsub('_12', '', pis_12)
  labs_pis_prior <- gsub('_prior', '', pis_prior)
  
  
  return(list("vars_4" = vars_4, "vars_12" = vars_12, 
              "pis_4" = pis_4, "pis_12" = pis_12, "pis_prior" = pis_prior, 
              "chems_4" = chems_4, "chems_12" = chems_12, "chems_prior" = chems_prior,
              "labs_4" = labs_4, "labs_12" = labs_12, 
              "labs_pis_4" = labs_pis_4, "labs_pis_12" = labs_pis_12, "labs_pis_prior" = labs_pis_prior,
              "labs_chems_4" = labs_chems_4, "labs_chems_12" = labs_chems_12, "labs_chems_prior" = labs_chems_prior))
  
}


# 5. Estimate models and plot results  -----
plot_model_results <- function(dep_vars, follow_up, labels, include_extremes, df){ 
  # dep_vars = a list of dependent variables e.g. list(df$bloodteststaken_4, df$long_covid_4) 
  # follow_up = days of follow up each individual completed by endpoint e.g. df$follow_up_4
  # labels = Grouped Read codes'/healthcare interactions names e.g. c("bloodteststaken", "longcovid")
  # include_extremes = include extreme values (< -100, > 100) e.g. TRUE/FALSE
  # df = dataframe to be used e.g. df_neg_m or df_nt_m
  
  
  # Specfy model: Function to specify model to be used to perform matched analysis -----
  ## Poisson regression 
  ## https://bookdown.org/roback/bookdown-BeyondMLR/ch-poissonreg.html 
  ## Interpret coefficients as rate ratios (for binary dependent variables, coefficients would be interpreted as risk ratios )
  
  # Health outcome ~ exposure status + covariates
  run_model <- function(y, follow_up, df){ # y = dependent variable e.g. bloodteststaken_4; 
    # follow_up = days of follow up achieved by 4-12 weeks and >12-26 weeks  
    ## follow_up is used to account for individuals who have an incomplete follow-up period 
    ## (early stu_end date due to death, control testing positive, or testing > 26 weeks before 
    ## the data extract was done)
    # df = dataframe to be used e.g. df_neg_m or df_nt_m
    
    
    ## As long as vac_n_clean has more than one level, run the model as ususal... 
    ## deals with the issue that arises during sub-sample analysis by variant e.g. no-one was vaccinated in 'Wild' period
    if(dim(table(df$vac_n_clean))>1){
      
      mod <- glm(get(y) ~ pos_or_not
                 # Covariates
                 + ns(age, df = 3) # age with natural splines (3 degrees of freedom)
                 + sex
                 + simd
                 + ur6_lab 
                 + council     
                 + tes_n # number of tests taken by tes_dat
                 + Q_DIAG_CKD_LEVEL
                 + Q_DIAG_AF
                 + Q_DIAG_ASTHMA
                 + Q_DIAG_BLOOD_CANCER
                 + Q_DIAG_CCF
                 + Q_DIAG_CEREBRALPALSY
                 + Q_DIAG_CHD
                 + Q_DIAG_CIRRHOSIS
                 + Q_DIAG_CONGEN_HD
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
                 + Q_DIAG_SICKLE_CELL
                 + Q_DIAG_STROKE
                 + Q_DIAG_VTE
                 + ns(bmi_imp, df = 3) # imputes missing bmi data for 54%, includes natural splines (3 degrees of freedom)
                 + immuno_supp
                 + shielding
                 + hhold_bands # no. individuals in the household
                 + vac_n_clean # number of vaccine doses received
                 + hos_prior # hospitalisations in the year prior to testing 
                 + icu_prior # icu admissions in the year prior to testing 
                 + tes_mon # controls for seasonality
                 ,weights = weights # matching weights incorporate sampling weights where appropriate
                 ,family = "quasipoisson" # use quasi to adjust for overdispersion
                 ,offset = log(follow_up)             
                 ,data = df)
      #print(summary(mod))
      return(mod)
    }
    
    ## If vac_n_clean only has one level (e.g. if no one was vaccinated in that period), run the model without it
    if(dim(table(df$vac_n_clean))<=1){
      mod <- glm(get(y) ~ pos_or_not
                 # Covariates
                 + ns(age, df = 3) # age with natural splines (3 degrees of freedom)
                 + sex
                 + simd
                 + ur6_lab 
                 + council     
                 + tes_n # number of tests taken by tes_dat
                 + Q_DIAG_CKD_LEVEL
                 + Q_DIAG_AF
                 + Q_DIAG_ASTHMA
                 + Q_DIAG_BLOOD_CANCER
                 + Q_DIAG_CCF
                 + Q_DIAG_CEREBRALPALSY
                 + Q_DIAG_CHD
                 + Q_DIAG_CIRRHOSIS
                 + Q_DIAG_CONGEN_HD
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
                 + Q_DIAG_SICKLE_CELL
                 + Q_DIAG_STROKE
                 + Q_DIAG_VTE
                 + ns(bmi_imp, df = 3) # imputes missing bmi data for 54%, includes natural splines (3 degrees of freedom)
                 + immuno_supp
                 + shielding
                 + hhold_bands # no. individuals in the household
                 #+ vac_n_clean # number of vaccine doses received - intentionally omitted (see note above)
                 + hos_prior # hospitalised in the year prior to testing 
                 + icu_prior 
                 + tes_mon # controls for seasonality
                 ,weights = weights # matching weights incorporate sampling weights where appropriate
                 ,family = "quasipoisson" # use quasi to adjust for overdispersion
                 ,offset = log(follow_up)             
                 ,data = df)
      #print(summary(mod))
      return(mod)
    }
    
    
  }
  
  # Run model for each dependent variable, and:
  ## Check for overdispersion
  ## Goodness of fit test
  ## Collect results in 'mods' 
  
  # Create empty data frames to store results in
  #goodness_of_fit <- data.frame()
  #overdispersion <- data.frame() # Not sure if this is needed, given overdispersion adjustment in model e.g. family = "quasipoisson"
  mods <- data.frame()
  
  for(var in dep_vars){
    
    # Run model
    mod <- run_model(y = var, follow_up, df)
    
    # Estimate SE
    # https://stats.oarc.ucla.edu/r/dae/poisson-regression/ 
    cov.mod <- vcovHC(mod, type="HC0")
    std.err <- sqrt(diag(cov.mod))
    Estimate <- coef(mod)[!is.na(coef(mod))]
    robust <- data.frame("Estimate" = Estimate, 
                         "Robust SE" = std.err,
                         "Pr(>|z|)" = 2 * pnorm(abs(Estimate/std.err), lower.tail=FALSE),
                         LL = Estimate - 1.96 * std.err,
                         UL = Estimate + 1.96 * std.err)
    
    # Collect coefficient and confidence interval
    coef <- exp(robust$Estimate[2])
    conf <- exp(c(robust$LL[2], robust$UL[2]))
    pval <- robust$Pr...z..[2]
    mod_output <- c(coef, conf, pval)
    
    # Save in mods
    mods <- rbind(mods, mod_output)
    
    # Goodness of fit test
    # https://stats.oarc.ucla.edu/r/dae/poisson-regression/
    gf <- cbind(res.deviance = mod$deviance, 
                deg_free = mod$df.residual,
                p = pchisq(mod$deviance, mod$df.residual, lower.tail=FALSE))
    goodness_of_fit <- rbind(goodness_of_fit, gf)
    
    
    ## Estimate overdispersion
    ### Guidance: https://towardsdatascience.com/adjust-for-overdispersion-in-poisson-regression-4b1f52baa2f1 
    od_val <- sum(residuals(mod, type ="pearson")^2)/mod$df.residual
    od_sig <- pchisq(mod$deviance, mod$df.residual, lower.tail=F)
    od <- c(od_val, od_sig)
    overdispersion <- rbind(overdispersion, od)
  }
  
  mods$labs <- labels
  colnames(mods) <- c("RR", "lower", "upper", "pval", "grp")
  
  # Adjust p-values to minimise false discovery rate
  mods$pval <- p.adjust(mods$pval, method = "BH")
  
  # Arrange by Rate Ratio  
  mods <- arrange(mods, RR)
  order <- c(mods$grp)
  mods$grp <- factor(mods$grp, levels = order)
  
  # Remove extreme values if specified
  extremes <- include_extremes
  
  if(extremes == F){
    mods <- mods %>%  filter(upper < 1000) %>% filter(lower > 0.01)
    print("Extreme values removed from plot")}
  
  # Plot figure
  fig <- ggplot(data=mods, aes(x=grp, y=RR, ymin=lower, ymax=upper)) +
    geom_pointrange() + 
    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    scale_y_log10() + # make axis logarithmic
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Health outcome") + ylab("Adjusted rate ratio (95% CI), log scale") +
    theme_bw() # use a white background
  
  # Format table (has to come after plot is made because of renaming grp to Outcome)
  mods <- mods %>% 
    mutate(Outcome = grp,
           `Sig` = ifelse(lower > 1 & pval <=.05, "Significantly higher", 
                          ifelse(upper < 1 & pval <=.05, "Significantly lower", "Not significant"))) %>% 
    select(Outcome, RR, lower, upper, pval, Sig)
  
  
  print("Coefficients and confidence intervals (robust SEs)")
  print(mods)
  
  print("Goodness of fit")
  print("Statistically significant chi-square (p < .05) indicates that the data do not fit the model well")
  print(goodness_of_fit)
  
  
  print("Checks for overdispersion (value, significance)")
  print("Values significantly > 1 suggest overdispersion")
  print(overdispersion)
  
  return(list("table" = mods, "plot" = fig))
  
}


# 6. Prepare data and perform analyses: Function runs all previous functions on the specified dataframe----
analysis <- function(df_path, # e.g. "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs neg.rds"
                     comparison, # string for plot and table titles e.g. "pos vs neg"
                     period){ # string for plot and table titles e.g. "alpha"
  
  # Read in data
  df <- readRDS(df_path)
  
  # Remove unnecessary variables
  df <- df %>% 
    select(-c(subclass)) 
  
  # Prepare count data
  print("Preparing GP data")
  df <- prepare_GP_data(df,  cluster = F)
  print("Preparing interactions data")
  df <- prepare_interactions_data(df, cluster = F)
  print("Preparing PIS data")
  df <- prepare_pis(df) 
  print("Preparing variables")
  df_vars <- set_dep_vars(df$GP, df$pis, df$chems, comparison) 
  
  # Run matched analysis, create and export plots
  ## 4-12 weeks
  print("Preparing plot and table - 4-12 weeks")
  res_4 <- plot_model_results(dep_vars = df_vars$vars_4,
                              follow_up = df$GP$follow_up_4, 
                              labels = df_vars$labs_4, 
                              include_extremes = F,
                              df = df$GP)
  
  ggsave(paste0("4-12 weeks - ", comparison, " - matched sample - ", period, ".pdf"),
         plot = res_4$plot, width = 6, height = 8, units = "in", dpi = 300)
  
  print("Writing to file")
  write.csv(res_4$table, paste0("4-12 weeks - ", comparison, " - matched sample - ", period, ".csv"), row.names = F)
  
  ## >12-26 weeks
  print("Preparing plot and table - 12-26 weeks")
  res_12 <- plot_model_results(dep_vars = df_vars$vars_12,
                               follow_up = df$GP$follow_up_12, 
                               labels = df_vars$labs_12, 
                               include_extremes = F,
                               df = df$GP)
  
  ggsave(paste0("12-26 weeks - ",  comparison, " - matched sample - ", period, ".pdf"),
         plot = res_12$plot, width = 6, height = 8, units = "in", dpi = 300)
  
  print("Writing to file")
  write.csv(res_12$table, paste0("12-26 weeks - ",  comparison, " - matched sample - ", period, ".csv"), row.names = F)
  
  
  ### Prescriptions 4-12 weeks
  print("Preparing prescriptions plot and table - 4-12 weeks")
  pis_4  <- plot_model_results(dep_vars = df_vars$pis_4,
                               follow_up = df$GP$follow_up_4, 
                               labels = df_vars$labs_pis_4, 
                               include_extremes = F,
                               df = df$GP)
  
  ggsave(paste0("PIS 4-12 weeks - ", comparison, " - matched sample - ", period, ".pdf"),
         plot = pis_4$plot, width = 6, height = 8, units = "in", dpi = 300)
  
  print("Writing to file")
  write.csv(pis_4$table, paste0("PIS 4-12 weeks - ", comparison, " - matched sample - ", period, ".csv"), row.names = F)
  
  
  ## Prescriptions >12-26 weeks
  print("Preparing prescriptions plot and table - 12-26 weeks")
  pis_12  <- plot_model_results(dep_vars = df_vars$pis_12,
                                follow_up = df$GP$follow_up_12, 
                                labels = df_vars$labs_pis_12, 
                                include_extremes = F,
                                df = df$GP)
  
  ggsave(paste0("PIS 12-26 weeks - ", comparison, " - matched sample - ", period, ".pdf"),
         plot = pis_12$plot, width = 6, height = 8, units = "in", dpi = 300)
  
  print("Writing to file")
  write.csv(pis_12$table, paste0("PIS 12-26 weeks - ", comparison, " - matched sample - ", period, ".csv"), row.names = F)
  
  
  ### Chemical substances 4-12 weeks
  print("Preparing chemical substances plot and table - 4-12 weeks")
  chems_4  <- plot_model_results(dep_vars = df_vars$chems_4,
                                 follow_up = df$GP$follow_up_4, 
                                 labels = df_vars$labs_chems_4, 
                                 include_extremes = F,
                                 df = df$GP)
  
  ggsave(paste0("Chems 4-12 weeks - ", comparison, " - matched sample - ", period, ".pdf"),
         plot = chems_4$plot, width = 6, height = 8, units = "in", dpi = 300)
  
  print("Writing to file")
  write.csv(chems_4$table, paste0("Chems 4-12 weeks - ", comparison, " - matched sample - ", period, ".csv"), row.names = F)
  
  
  ### Chemical substances >12-26 weeks
  print("Preparing chemical substances plot and table - 12-26 weeks")
  chems_12  <- plot_model_results(dep_vars = df_vars$chems_12,
                                  follow_up = df$GP$follow_up_12, 
                                  labels = df_vars$labs_chems_12, 
                                  include_extremes = F,
                                  df = df$GP)
  
  ggsave(paste0("Chems 12-26 weeks - ", comparison, " - matched sample - ", period, ".pdf"),
         plot = chems_12$plot, width = 6, height = 8, units = "in", dpi = 300)
  
  print("Writing to file")
  write.csv(chems_12$table, paste0("Chems 12-26 weeks - ", comparison, " - matched sample - ", period, ".csv"), row.names = F)
  
}


# EXPORT DATA FOR CLUSTER ANALYSIS/OPERATIONAL DEFINITION ----
# Prepares and exports a df that attaches the counts of dependent variables to the linked data for the unmatched sample (i.e. the full cohort after censoring)
# Does this once each for dfs derived using: (i) PCR test results, and (ii) PCR and LFT results
cluster_prep <- function(df){

  ## Prepare counts of dependent variables
  df <- prepare_GP_data(df, cluster = T)
  df <- prepare_interactions_data(df, cluster = T)
  df <- prepare_pis(df) 
  
  ## Prepare df to export
  df <- df$GP
  
  ## Remove "_prior" columns
  col_names <- colnames(df)
  valid_names <- col_names[!grepl("_prior", col_names)]
  df <- df[, valid_names]
  
  return(df)
}


## (i) df derived using PCR test results (remove positive cases duplicated for use as controls in matched analysis)
## NB the analysis below removes 'not yet tested' form the dataset - so only positive and negative cases remain (but that's ok)
df_pcr_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_linked.rds") %>% 
  filter(pos_control == 0) # N = 5,056,071
df_pcr <- cluster_prep(df_pcr_raw)
saveRDS(df_pcr, "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/5. Cluster analysis/df_cluster_PCR.rds")


## (ii) df derived using PCR & LFT results
df_lft_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_linked_LFT.rds")
df_lft <- cluster_prep(df_lft_raw)
saveRDS(df_lft, "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/5. Cluster analysis/df_cluster_LFT.rds")


# RUN ANALYSES ----
## Matched analysis for each df - results are exported to:
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/2. Matched analysis")

### All ----
## Pos vs neg
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs neg.rds",
         comparison = "pos vs neg", period = "all")

## Pos vs not yet tested
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs nt.rds",
         comparison = "pos vs nt", period = "all")


### Wild ----
## Pos vs neg
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs neg - wild.rds",
         comparison = "pos vs neg", period = "wild")

## Pos vs not yet tested
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs nt - wild.rds",
         comparison = "pos vs nt", period = "wild")



### Alpha ----
## Pos vs neg
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs neg - alpha.rds",
         comparison = "pos vs neg", period = "alpha")

## Pos vs not yet tested
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs nt - alpha.rds",
         comparison = "pos vs nt", period = "alpha")



### Delta ----
## Pos vs neg
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs neg - delta.rds",
         comparison = "pos vs neg", period = "delta")

## Pos vs not yet tested
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs nt - delta.rds",
         comparison = "pos vs nt", period = "delta")



### Omicron ----
## Pos vs neg
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs neg - omicron.rds",
         comparison = "pos vs neg", period = "omicron")

## Pos vs not yet tested
analysis(df_path = "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes/Matched pairs - pos vs nt - omicron.rds",
         comparison = "pos vs nt", period = "omicron")

# 7. Plot by variant periods (pos vs neg only) ----
var_per_plot <- function(wild_path, alpha_path, delta_path, omicron_path){ # filepaths for the relevant data
  
  # Read in dataa and identify the Period
  wild <- read_csv(wild_path) %>% mutate(Period = "Wild")
  alpha <- read_csv(alpha_path) %>% mutate(Period = "Alpha")
  delta <- read_csv(delta_path) %>% mutate(Period = "Delta")
  omicron <- read_csv(omicron_path) %>% mutate(Period = "Omicron")
  
  # Combine
  df <- rbind(wild, alpha, delta, omicron)
  
  # Remove unwanted rows and columns
  df <- df %>% 
    select(-c(pval, Sig)) %>% 
    filter(Outcome != "Sick Note (except long-COVID)") %>% 
    mutate(Outcome = ifelse(Outcome == "Coronavirus", "Antivirals for Coronavirus", Outcome))
  
  # Arrange by RR
  df <- arrange(df, RR)
  order <- unique(c(df$Outcome))
  df$Outcome <- factor(df$Outcome, levels = order)
  
  # Set order of legend labels
  df$Period <- factor(df$Period, levels = c("Wild", "Alpha", "Delta", "Omicron"))

  ## Staggered points on a single plot
  #fig <- ggplot(data=df, aes(x=Outcome, y=RR, ymin=lower, ymax=upper, color = Period, shape = Period)) +
  #  geom_pointrange(position = position_dodge(width = 0.9)) +
  #  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  #  scale_y_log10() + # make axis logarithmic
  #  coord_flip() +  # flip coordinates (puts labels on y axis)
  #  xlab("Health outcome") + ylab("Adjusted rate ratio (95% CI), log scale") +
  #  theme_bw() #+# use a white background
  #  #scale_colour_viridis_d(option = "inferno") + scale_shape_discrete(solid = T)
  
  # Faceted plot
  fig <- ggplot(data = df, aes(x = RR, y=Outcome, color = Period)) +
    geom_point() + 
    geom_pointrange(aes(xmin=lower, xmax=upper)) +
    geom_vline(xintercept=1, lty=2) +  # add a dotted line at x=1 after flip
    scale_x_log10() + # make axis logarithmic
    ylab("Health outcome") + xlab("Adjusted rate ratio (95% CI), log scale") +
    theme_bw() +# use a white background
    facet_wrap(~Period, nrow = 1)
    
  return(fig)
}

## 4-12 weeks
var_per_4 <- var_per_plot(wild_path = "4-12 weeks - pos vs neg - matched sample - wild.csv",
                          alpha_path = "4-12 weeks - pos vs neg - matched sample - alpha.csv",
                          delta_path = "4-12 weeks - pos vs neg - matched sample - delta.csv",
                          omicron_path = "4-12 weeks - pos vs neg - matched sample - omicron.csv")

ggsave("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/1. Variant periods/Variant periods - 4-12 weeks - pos vs neg.jpg",
       plot = var_per_4, width = 15, height = 9, units = "in", dpi = 300)

## 4-12 weeks: PIS
pis_per_4 <- var_per_plot(wild_path = "PIS 4-12 weeks - pos vs neg - matched sample - wild.csv",
                          alpha_path = "PIS 4-12 weeks - pos vs neg - matched sample - alpha.csv",
                          delta_path = "PIS 4-12 weeks - pos vs neg - matched sample - delta.csv",
                          omicron_path = "PIS 4-12 weeks - pos vs neg - matched sample - omicron.csv")

ggsave("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/1. Variant periods/Variant periods - PIS 4-12 weeks - pos vs neg.jpg",
       plot = pis_per_4, width = 15, height = 9, units = "in", dpi = 300)


## 12-26 weeks
var_per_12 <- var_per_plot(wild_path = "12-26 weeks - pos vs neg - matched sample - wild.csv",
                          alpha_path = "12-26 weeks - pos vs neg - matched sample - alpha.csv",
                          delta_path = "12-26 weeks - pos vs neg - matched sample - delta.csv",
                          omicron_path = "12-26 weeks - pos vs neg - matched sample - omicron.csv")

ggsave("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/1. Variant periods/Variant periods - 12-26 weeks - pos vs neg.jpg",
       plot = var_per_12, width = 15, height = 9, units = "in", dpi = 300)


## 12-26 weeks: PIS
pis_per_12 <- var_per_plot(wild_path = "PIS 12-26 weeks - pos vs neg - matched sample - wild.csv",
                          alpha_path = "PIS 12-26 weeks - pos vs neg - matched sample - alpha.csv",
                          delta_path = "PIS 12-26 weeks - pos vs neg - matched sample - delta.csv",
                          omicron_path = "PIS 12-26 weeks - pos vs neg - matched sample - omicron.csv")

ggsave("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/1. Variant periods/Variant periods - PIS 12-26 weeks - pos vs neg.jpg",
       plot = pis_per_12, width = 15, height = 9, units = "in", dpi = 300)

