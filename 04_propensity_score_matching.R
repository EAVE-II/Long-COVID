#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 04_propensity_score_matching - Estimates time-varying (by month) propensity to test positive for 
## COVID-19 and matches individuals who have and have not tested positive for COVID-19 on propensity scores and other covariates.
## Exports df of matched pairs (positive versus negative and positive versus not yet tested) to be used to identify 
## differences in health outcomes between groups 4-12 weeks and >12-26 weeks following testing.
#########################################################################################################################

# 0. Set-up ------
# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(sandwich) # for estimating robust standard errors
library(MatchIt) # for propensity-score matching https://sejdemyr.github.io/r-tutorials/statistics/tutorial8.html
library(splines) # for modelling non-linear age and BMI

# Read in linked data
# NB this file contains duplicate EAVE_LINKNOs where an individual can be used as a negative case and later as a positive case
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_linked.rds") 

# Split df into positive vs negative and positive vs not yet tested dataframes
df_neg <- df %>% filter(tes_res != "NEVER TESTED") # Positive vs negative
df_nt <- df # Positive vs not yet tested - positive cases will be matched to controls that have not tested by the date of the positive case's test.
            # Those controls may go on to test positive or negative later --> don't filter out any data here.
rm(df)

## Set periods for sub-sample analysis (where the variant represents >60% of positive cases, based on analysis of whole genome sequencing data in 02_variant_periods)
wild_start <- ymd("2020-03-01")
wild_end <-   ymd("2021-01-10")

alpha_start <-  ymd("2021-01-11")
alpha_end <-  ymd("2021-05-09")

delta_start <-  ymd("2021-05-24")
delta_end <-  ymd("2021-12-05")

omicron_start <-  ymd("2021-12-27")
omicron_end <-  ymd("2022-04-30")# date when PCR testing ended in Scotland 


# 1. Read in & prepare data used for estimating propensity scores ----
## Do this to get counts of (i) vaccinations up to 14 days before tes_dat and (ii) hospitalisations and (iii) icu admissions in the year before tes_dat 
## (used in propensity score model). Do this inside ps_match (rather than 03_data_setup) so that a random date in the relevant month can been used 
## as tes_dat for cases that have 'not yet tested'

# Hospitalisations (1/3/2019 to latest test date)
hos <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/SMR01 - 2019-03-01 to 2022-09-21.rds")

# ICU admissions (1/3/19 to latest test date)
icu <- readRDS("/conf/EAVE/GPanalysis/data/SICSAG_episode_level_.rds") %>% 
  filter(!is.na(EAVE_LINKNO)) %>% 
  select(EAVE_LINKNO, icu_dat = AdmitUnit) %>% 
  filter(icu_dat <= max(df_neg$tes_dat))

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
  filter(vac_dat >= "2020-12-01" & vac_dat <= (max(df_neg$tes_dat - days(14)))) %>% 
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
  ## Flag vaccinations 1 & 2 or 2 & 3 <18 days apart
  mutate(vac_gap_1 = as.numeric(date(vac_dat_2) - date(vac_dat_1)), 
         vac_gap_2 = as.numeric(date(vac_dat_3) - date(vac_dat_2)),
         vac_incon = ifelse(vac_gap_1 <19 & !is.na(vac_gap_1), 1, 0)) %>% 
  ## Flag first dose ==PB or Mo, but second dose == AZ
  mutate(vac_incon = ifelse((vac_typ_1 == "PB" | vac_typ_1 == "Mo") & vac_typ_2 == "AZ", 1, 0),
         ## Flag no record of earlier doses
         vac_incon = ifelse((is.na(vac_typ_1) & (!is.na(vac_typ_2) | !is.na(vac_typ_3) | !is.na(vac_typ_4))) |    # No first dose, but subsequent doses
                              ((is.na(vac_typ_1) | is.na(vac_typ_2)) & (!is.na(vac_typ_3) | !is.na(vac_typ_4))) |    # No first or second dose, but subsequent doses
                              ((is.na(vac_typ_1) | is.na(vac_typ_2) | is.na(vac_typ_3)) & !is.na(vac_typ_4)), 1, 0)) %>%  # No first, second, or third dose, but subsequent doses
  ## Keep data with no inconsistencies
  filter(vac_incon == 0) %>% 
  ## Transpose to long form 
  pivot_longer(cols = c(vac_dat_1:vac_dat_7), # columns that should pivot from wide to long (unquoted)
               names_to = "vac_dos", # name of the new category column as a quoted string
               values_to = "vac_dat") %>% # name of the new value column as a quoted string)
  filter(!is.na(vac_dat)) %>% 
  select(EAVE_LINKNO, vac_dat)

rm(vac_raw)


# 2. Function - Time varying propensity score matching ----

## Performs time-varying (by month) propensity score matching on a specified dataframe
ps_match <- function(df, control_type, start, end){ # df = the dataframe to use e.g. df_neg
                                                    # control_type = whether controls are negative or not yet tested (to allow propensity scores with/without 
                                                    # exact match on week of test) e.g. "neg" or "nt"  
                                                    # start = the start of the period of interest e.g. alpha_start
                                                    # end = the end of the period of itnerest e.g. omicron_end
  
  
  ## For each month in which there is a tes_dat in the dataset:
  ## (i) estimate individuals' propensity to test positive for COVID 
  ## (ii) collect pairs of individuals who did and did not test positive, matched on propensity scores and additional covariates
  ## (iii) identify the individual matched to
  ## Then merge matched pairs from each month into a single df 
  
  ## Prepare time periods (months)
  sta_date <- start
  end_date <- end
  months <- c(seq(sta_date, end_date, by='months'), end_date)
  
  ## Create empty object to capture matched and unmatched individuals' results for each month 
  all_patients <- NULL
  
  ## Create empty object to capture matched pairs for each month 
  all_matches <- NULL
  
  ## Loop over each month in the dataset and: 
  ## (i) estimate propensity to test positive for COVID 
  ## (ii) perform nearest neighbour matching on propensity scores and exact matching on additional covariates
  for(month in 1:(length(months)-1)){ # starting point for each month
 
    # Extract month start and end date 
    month_start <- months[month]
    month_end <- months[month + 1]
    
    # Print progress
    print(paste0("Time-period: ", month, " of ", length(months)-1))
    
    # Select individuals to include in analysis for the month and identify as exposed or control (pos_or_not)
    df_month <- df %>% 
      # Keep every 'never tested' case (for use as controls) and... 
      filter(tes_res == "NEVER TESTED" |
               # Keep every positive case that tested positive during the month (for use as exposed cases)
               (tes_res == "POSITIVE" & tes_dat >= month_start & tes_dat < month_end) |
               # Keep every positive case that tested positive after month_end (for use as neg/nt controls) *but* 
               # only keep the subset whose first test (pos or neg) was also after month end (because positive cases 
               # that had a negative test earlier were duplicated and marked as 'negative' - so these will be picked up with the negative cases)
               (tes_res == "POSITIVE" & tes_dat >= month_end & first_tes >= month_end) |
               # Keep every negative case that tested negative during the month (for use as negative controls) -- this picks up individuals who 
               # go on to test positive after month_end
               (tes_res == "NEGATIVE" & tes_dat >= month_start & tes_dat < month_end) | 
               # Keep every negative case that tested negative after month end (for use as nt controls)
               (tes_res == "NEGATIVE" & tes_dat >= month_end)) %>% 
      # Identify individuals as being in the exposed or control group
      mutate(pos_or_not = ifelse(tes_res == "POSITIVE" & !is.na(tes_dat) & tes_dat >= month_start & tes_dat < month_end, 1, 0)) %>% 
      # Remove cases used as controls used in earlier months
      filter(!(EAVE_LINKNO %in% all_matches$EAVE_LINKNO)) 
    
    
    # Skip any months where there are not 2 values of tes_res (i.e. there must be positive AND non-positive test results)
    if(dim(table(df_month$pos_or_not)) == 2) { 
      
      # Set whether controls used are individuals with a negative PCR test or never tested
      control = control_type
      
      ## Set tes_dat to a random date within the month for cases that have 'not yet tested' to give a date for counting 
      ## vac_n_clean, hos_prio, icu_prior & update tes_n (excluding it from prop scores, but necessary as a control in regressions later)
      if(control == "nt"){
        
        # Define range of dates to select from
        start_date <- as.Date(month_start)
        end_date <- as.Date(month_end - days(1))
        
        # Generate random dates in range
        random_dates <- sample(seq(start_date, end_date, by="day"), nrow(df_month[df_month$pos_or_not==0,]), replace = T)
        
        #assign random dates to each row in dataframe that has not yet tested
        df_month$tes_dat[df_month$pos_or_not == 0] <- random_dates
        
        # Update existing tes_n (which contains tes_n by first negative or first positive test) so that it is correct for controls (i.e. all controls 
        # have not yet tested, so tes_n = 0)
        df_month <- df_month %>% 
          mutate(tes_n = ifelse(pos_or_not == 1, tes_n, 0L))
        
      }
      
      # Prepare variables needed for propensity score estimation 
      ## Do this here rather than during 03_data_setup so that counts for 'not yet tested' cases are in relation to month_start (given absence of a relevant tes_dat)
      
      ## Hospitalisations in the year before tes_dat
      hos_n <- hos %>% 
        filter(EAVE_LINKNO %in% df_month$EAVE_LINKNO) %>% 
        left_join(select(df_month, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
        filter(hos_dat <= tes_dat & hos_dat > (tes_dat - days(365))) %>% 
        group_by(EAVE_LINKNO) %>% 
        summarise(hos_prior = n())
      
      ## ICU admissions in the year before tes_dat
      icu_n <- icu %>% 
        filter(EAVE_LINKNO %in% df_month$EAVE_LINKNO) %>% 
        left_join(select(df_month, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
        filter(icu_dat <= tes_dat & icu_dat > (tes_dat - days(365))) %>% 
        group_by(EAVE_LINKNO) %>% 
        summarise(icu_prior = n())
      
      # Vaccine doses by 14 days before tes_dat
      vac_n <- vac %>% 
        filter(EAVE_LINKNO %in% df_month$EAVE_LINKNO) %>% 
        left_join(select(df_month, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
        filter(vac_dat <= (tes_dat - days(14))) %>% 
        group_by(EAVE_LINKNO) %>% 
        summarise(vac_n_clean = n())
      
      # Join to df_month
      df_month <- df_month %>% 
        left_join(hos_n, by = "EAVE_LINKNO") %>%
        left_join(icu_n, by = "EAVE_LINKNO") %>% 
        left_join(vac_n, by = "EAVE_LINKNO") %>% 
        mutate(hos_prior = ifelse(is.na(hos_prior), 0L, hos_prior),
               icu_prior = ifelse(is.na(icu_prior), 0L, icu_prior),
               vac_n_clean= ifelse(is.na(vac_n_clean), 0L, vac_n_clean))
      

      # For negative controls, include exact matching on week of test
      if(control == "neg") {
        
        # Estimate propensity scores 
        # MatchIt documentation: https://cran.r-project.org/web/packages/MatchIt/vignettes/
        # Guidance: https://stats.stackexchange.com/questions/562403/propensity-score-match-with-exact-matching-on-one-variable
        matches <- matchit(
          
          
          ## Estimate propensity to test positive for COVID and match pos & neg on propensity scores (nearest neighbour matching) 
          pos_or_not ~ ns(age, df = 3) # age with natural splines (3 degrees of freedom)
          + sex 
          + simd
          + ur6_lab 
          + council
          + tes_n # n tests by tes_dat - NB impossible to achieve balance on (neg always = 1, nt always = 0) so not included in balance checks
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
          + ns(bmi_imp, df = 3) # imputes missing bmi data includes natural splines (3 degrees of freedom)
          + immuno_supp
          + shielding
          + hhold_bands # no. individuals in the household
          + vac_n_clean # number of vaccine doses received
          + hos_prior # hospitalised in the year prior to testing
          + icu_prior # admitted to ICU in the year prior to testing
          
          , method = "nearest" # matching method: nearest neighbour matching (on propensity score) 
          , distance = "glm" # method for estimating the propensity score: 'glm' = logit
          , caliper = 0.8 # only allow matches within x standard deviations of one another
          
          # Also perform exact matching on additional variables
          , exact = ~ tes_wk # week of testing (1:52) 
          + tes_yr # year of testing 
          + age # truncated older age groups
          
          , reuse.max = 1 # use each control once only (per month) - extra code (further down) ensures controls are not used more than once across months
          , data = df_month
          , s.weights = ~ eave_weight # include sampling weights - not essential, but recommended if it improves matching
          
        )
      }
      
      
      # For not yet tested controls, do not include exact maching on week of test (because the test date in the df is assigned at random)
      # NB although there's no matching on week of test, later code ensures we follow up matched pairs over identical dates 
      # e.g. a matched pair including an individual who tested positive on 1st February 2021 will be followed-up from 4 - 26 weeks following 1st February 2021 
      if(control == "nt") {
        
        # Estimate propensity scores 
        # MatchIt documentation: https://cran.r-project.org/web/packages/MatchIt/vignettes/
        # Guidance: https://stats.stackexchange.com/questions/562403/propensity-score-match-with-exact-matching-on-one-variable
        matches <- matchit(
          
          ## Estimate propensity to test positive for COVID and match of pos & not pos on propensity scores (nearest neighbour matching) 
          pos_or_not ~ ns(age, df = 3) # age with natural splines (3 degrees of freedom)
          + sex 
          + simd
          + ur6_lab 
          + council
          #+ tes_n # n tests by tes_dat - NB impossible to match positive (always >1) and not yet tested (always 0) so omit 
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
          + ns(bmi_imp, df = 3) # imputes missing bmi data includes natural splines (3 degrees of freedom)
          + immuno_supp
          + shielding
          + hhold_bands # no. individuals in the household
          + vac_n_clean # number of vaccine doses received
          + hos_prior # hospitalised in the year prior to testing
          + icu_prior # admitted to ICU in the year prior to testing
          
          , method = "nearest" # matching method: nearest neighbour matching (on propensity score) 
          , distance = "glm" # method for estimating the propensity score: 'glm' = logit
          , caliper = 0.8 # only allow matches within x standard deviations of one another
          
          # Also perform exact matching on additional variables
          , exact = ~ age # truncated older age groups
          , reuse.max = 1 # use each control once only (per month) - extra code ensures controls are not used more than once across months
          , data = df_month
          , s.weights = ~ eave_weight # include sampling weights - not essential, but recommended if it improves matching
          
        )
      }
      
      
      ## Get matched pairs only and rename 'distance' to 'prop_score'
      month_matches <- match.data(matches, distance = "prop_score")
      
      ## Incorporate EAVE_LINKNO of individuals matched to 
      ### Extract df_month row numbers for matched pairs from matches$match.matrix
      EAVE_match <- as.data.frame(na.omit(matches$match.matrix))
      
      ### Make row numbers containing location of matched individuals into a column
      EAVE_match <- tibble::rownames_to_column(EAVE_match, "V0")
      
      ### Get corresponding EAVE_LINKNOs
      EAVE_match <- EAVE_match %>% 
        mutate(V0 = as.numeric(V0),
               V1 = as.numeric(levels(V1))[V1],
               V0 = df_month$EAVE_LINKNO[V0],
               V1 = df_month$EAVE_LINKNO[V1])
      
      ### Integrate matched individuals' EAVE_LINKNO into 'month_matches'
      month_matches <- month_matches %>% 
        left_join(EAVE_match, by = c("EAVE_LINKNO" = "V0")) %>% 
        left_join(EAVE_match, by = c("EAVE_LINKNO" = "V1")) %>% 
        mutate(m_EAVE_LINKNO = ifelse(is.na(V0), V1, V0)) %>% 
        select(-c(V0, V1))  
      
      
      ## Integrate stu_sta, stu_12w, stu_end, dea_dat, tes_res, tes_dat of matched individuals  
      ## this allows updating of stu_sta, stu_12w, stu_end so that exposed and controls are followed up over identical dates (exposed match's dates)
      ## but pairs are censored if one individual tests positive later, dies, or the end of the follow up is reached (due to data extract before stu_end)
      month_matches <- month_matches %>% 
        left_join(select(df_month, EAVE_LINKNO, stu_sta, stu_12w, stu_end, dea_dat, tes_dat, tes_res, pos_dat, pos_control), by = c("m_EAVE_LINKNO" = "EAVE_LINKNO")) %>% 
        rename(stu_sta = stu_sta.x,
               stu_12w = stu_12w.x,
               stu_end = stu_end.x,
               dea_dat = dea_dat.x,
               tes_dat = tes_dat.x,
               tes_res = tes_res.x,
               pos_dat = pos_dat.x,
               pos_control = pos_control.x,
               m_stu_sta = stu_sta.y,
               m_stu_12w = stu_12w.y,
               m_stu_end = stu_end.y,
               m_dea_dat = dea_dat.y,
               m_tes_dat = tes_dat.y,
               m_tes_res = tes_res.y,
               m_pos_dat = pos_dat.y,
               m_pos_control = pos_control.y) 
      
      ## Update tes_dat, stu_sta, stu_12w, stu_end so that controls are followed for the same dates as their positive matches
      ## and so that stu_12w and stu_end is censored to be the earliest of all pairs' stu_12w and stu_end 
      ## NB: 01_data_setup censors end dates for deaths/incomplete follow up/subsequent positive test for those with pos/neg tests
      ## but the code below updates censoring for 'not yet tested'
      month_matches <- month_matches %>%  
        
        mutate(tes_dat = ifelse(pos_or_not == 0, m_tes_dat, tes_dat), # update controls' dates to be the same as the positive case they are matched to      
               tes_dat = as.Date(tes_dat, origin = "1970-01-01"), 
               stu_sta = ifelse(pos_or_not == 0, m_stu_sta, stu_sta),
               stu_sta = as.Date(stu_sta, origin = "1970-01-01"),
               stu_12w = ifelse(pos_or_not == 0, m_stu_12w, stu_12w),
               stu_12w = as.Date(stu_12w, origin = "1970-01-01"),
               stu_end = ifelse(pos_or_not == 0, m_stu_end, stu_end),
               stu_end = as.Date(stu_end, origin = "1970-01-01")) %>% 
        
        ## Censor pairs for controls' death
        filter(is.na(dea_dat) | dea_dat > stu_sta) %>% # removes controls that die before their stu_sta (which has been updated to match the positive case's stu_sta)
        filter(is.na(m_dea_dat) | m_dea_dat > stu_sta) %>% # removes positive cases where the matched control dies before the positive case's stu_sta 
        mutate(stu_12w = ifelse(pos_or_not == 0 & !is.na(dea_dat) & dea_dat < stu_12w, dea_dat, stu_12w), # note: if pos_or_not == 0, dea_dat = control's dea_dat
               stu_12w = ifelse(pos_or_not == 1 & !is.na(m_dea_dat) & m_dea_dat < stu_12w, m_dea_dat, stu_12w), # note: if pos_or_not == 1, m_dea_dat = control's dea_dat
               stu_12w = as.Date(stu_12w, origin = "1970-01-01"),
               stu_end = ifelse(pos_or_not == 0 & !is.na(dea_dat) & dea_dat < stu_end, dea_dat, stu_end),
               stu_end = ifelse(pos_or_not == 1 & !is.na(m_dea_dat) & m_dea_dat < stu_end, m_dea_dat, stu_end),
               stu_end = as.Date(stu_end, origin = "1970-01-01")) %>% 
        
        # Censor pairs for controls' subsequent positive tests
        filter(!(pos_control == 1 & pos_dat < stu_sta)) %>% # don't keep controls that tested positive before their (updated) study start
        filter(!(m_pos_control == 1 & m_pos_dat < stu_sta)) %>% #don't keep positive cases matched to a control that tested pos before study start
        mutate(stu_12w = ifelse(pos_control == 1 & pos_dat < stu_12w, pos_dat, stu_12w), # note: if pos_control == 1, pos_dat = control's positive test date
               stu_12w = ifelse(m_pos_control == 1 & m_pos_dat < stu_12w, m_pos_dat, stu_12w), # note: if m_pos_control == 1, m_tes_dat = control's positive test date
               stu_12w = as.Date(stu_12w, origin = "1970-01-01"),
               stu_end = ifelse(pos_control == 1 & pos_dat < stu_end, pos_dat, stu_end),
               stu_end = ifelse(m_pos_control == 1 & m_pos_dat < stu_end, m_pos_dat, stu_end),
               stu_end = as.Date(stu_end, origin = "1970-01-01")) %>%  
        
        # remove cols
        select(-c(m_stu_sta, m_stu_12w, m_stu_end, m_dea_dat, m_tes_dat, m_tes_res, m_pos_dat, m_pos_control))
      
      
      ## Merge matched pairs for month into df with previous months' pairs
      all_matches <- rbind(all_matches, month_matches)
      
      ## Get output for ALL individuals (not just matched) and rename 'distance' to 'prop_score' (for balance checks)
      month_patients <- match.data(matches, distance = "prop_score", drop.unmatched = FALSE, data = df_month)
      
      ## Keep only those whose tes_dat is in the current month
      month_patients <- month_patients %>% filter(tes_dat >= month_start & tes_dat < month_end)
      
      ## Merge all patients for month into df with previous months' patients
      all_patients <- rbind(all_patients, month_patients)
      
      ## Remove files
      rm(df_month, month_patients, month_matches, EAVE_match)
      
    }
  }
  
  # Print % of positive cases preserved after matching
  pos_cases_matched <- round(length(unique(all_matches$EAVE_LINKNO[all_matches$pos_or_not == 1]))/length(unique(all_patients$EAVE_LINKNO[all_patients$tes_res == "POSITIVE"]))*100, 1)
  
  print(paste0(pos_cases_matched, "% of positive cases matched to a control"))
  
  return(list("all_matches" = all_matches, "all_patients" = all_patients))
}



# 3. MATCHING (Estimate propensity scores then perform time-varying propensity score matching) ----
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/2. Matching/2. Matched dataframes")

## Full dataset ----
### Pos vs neg
### Match positive and negative cases
df_neg_matching <- ps_match(df_neg, 
                            control_type = "neg",
                            start = wild_start, 
                            end = omicron_end) # 54.4% of positive cases matched to a control

### Select matched pairs and export
df_neg_m <- df_neg_matching$all_matches
saveRDS(df_neg_m, "Matched pairs - pos vs neg.rds")

### Select all patients (with propensity scores, weights etc.) and export
df_neg_all <- df_neg_matching$all_patients
saveRDS(df_neg_all, "Linked df with propensity scores - pos vs neg.rds")

### Remove files
rm(df_neg_matching, df_neg_m, df_neg_all)


### Pos vs nt
### Match positive and not yet tested cases
df_nt_matching <- ps_match(df_nt, 
                           control_type = "nt",
                           start = wild_start,
                           end = omicron_end) # 84.4% of positive cases matched to a control

### Select matched pairs
df_nt_m <- df_nt_matching$all_matches
saveRDS(df_nt_m, "Matched pairs - pos vs nt.rds")

### Select all patients (with propensity scores, weights etc.)
df_nt_all <- df_nt_matching$all_patients
saveRDS(df_nt_all, "Linked df with propensity scores - pos vs nt.rds")

## Remove files
rm(df_nt_matching, df_nt_m, df_nt_all)


## Wild ----
### Pos vs neg
### Matching
wild_matching_neg <- ps_match(df = df_neg,
                              control_type = "neg",
                              start = wild_start,
                              end = wild_end) 
### Matched df
df_neg_m_wild <- wild_matching_neg$all_matches
saveRDS(df_neg_m_wild, "Matched pairs - pos vs neg - wild.rds")

### Full df (with propensity scores, weights etc.)
df_neg_wild <- wild_matching_neg$all_patients
saveRDS(df_neg_wild, "Linked df with propensity scores - pos vs neg  - wild.rds")

## Remove files
rm(wild_matching_neg, df_neg_m_wild, df_neg_wild)


### Pos vs nt
### Matching
wild_matching_nt <- ps_match(df = df_nt,
                             control_type = "nt",
                             start = wild_start,
                             end = wild_end) 
### Matched df
df_nt_m_wild <- wild_matching_nt$all_matches
saveRDS(df_nt_m_wild, "Matched pairs - pos vs nt - wild.rds")

### Full df (with propensity scores, weights etc.)
df_nt_wild <- wild_matching_nt$all_patients
saveRDS(df_nt_wild, "Linked df with propensity scores - pos vs nt  - wild.rds")

## Remove files
rm(wild_matching_nt, df_nt_m_wild, df_nt_wild)


## Alpha ----
### Pos vs neg
### Matching
alpha_matching_neg <- ps_match(df = df_neg,
                               control_type = "neg",
                               start = alpha_start,
                               end = alpha_end) 
### Matched df
df_neg_m_alpha <- alpha_matching_neg$all_matches
saveRDS(df_neg_m_alpha, "Matched pairs - pos vs neg - alpha.rds")

### Full df (with propensity scores, weights etc.)
df_neg_alpha <- alpha_matching_neg$all_patients
saveRDS(df_neg_alpha, "Linked df with propensity scores - pos vs neg  - alpha.rds")

## Remove files
rm(alpha_matching_neg, df_neg_m_alpha, df_neg_alpha)


### Pos vs nt
### Matching
alpha_matching_nt <- ps_match(df = df_nt,
                              control_type = "nt",
                              start = alpha_start,
                              end = alpha_end) 
### Matched df
df_nt_m_alpha <- alpha_matching_nt$all_matches
saveRDS(df_nt_m_alpha, "Matched pairs - pos vs nt - alpha.rds")

### Full df (with propensity scores, weights etc.)
df_nt_alpha <- alpha_matching_nt$all_patients
saveRDS(df_nt_alpha, "Linked df with propensity scores - pos vs nt  - alpha.rds")

## Remove files
rm(alpha_matching_nt, df_nt_m_alpha, df_nt_alpha)


## Delta ----
### Pos vs neg
### Matching
delta_matching_neg <- ps_match(df = df_neg,
                               control_type = "neg",
                               start = delta_start,
                               end = delta_end) 
### Matched df
df_neg_m_delta <- delta_matching_neg$all_matches
saveRDS(df_neg_m_delta, "Matched pairs - pos vs neg - delta.rds")

### Full df (with propensity scores, weights etc.)
df_neg_delta <- delta_matching_neg$all_patients
saveRDS(df_neg_delta, "Linked df with propensity scores - pos vs neg  - delta.rds")

## Remove files
rm(delta_matching_neg, df_neg_m_delta, df_neg_delta)


### Pos vs nt
### Matching
delta_matching_nt <- ps_match(df = df_nt,
                              control_type = "nt",
                              start = delta_start,
                              end = delta_end) 
### Matched df
df_nt_m_delta <- delta_matching_nt$all_matches
saveRDS(df_nt_m_delta, "Matched pairs - pos vs nt - delta.rds")

### Full df (with propensity scores, weights etc.)
df_nt_delta <- delta_matching_nt$all_patients
saveRDS(df_nt_delta, "Linked df with propensity scores - pos vs nt  - delta.rds")

## Remove files
rm(delta_matching_nt, df_nt_m_delta, df_nt_delta)


## Omicron ----
### Pos vs neg
### Matching
omicron_matching_neg <- ps_match(df = df_neg,
                                 control_type = "neg",
                                 start = omicron_start,
                                 end = omicron_end) 
### Matched df
df_neg_m_omicron <- omicron_matching_neg$all_matches
saveRDS(df_neg_m_omicron, "Matched pairs - pos vs neg - omicron.rds")

### Full df (with propensity scores, weights etc.)
df_neg_omicron <- omicron_matching_neg$all_patients
saveRDS(df_neg_omicron, "Linked df with propensity scores - pos vs neg  - omicron.rds")

## Remove files
rm(omicron_matching_neg, df_neg_m_omicron, df_neg_omicron)


### Pos vs nt
### Matching
omicron_matching_nt <- ps_match(df = df_nt,
                                control_type = "nt",
                                start = omicron_start,
                                end = omicron_end)  
### Matched df
df_nt_m_omicron <- omicron_matching_nt$all_matches
saveRDS(df_nt_m_omicron, "Matched pairs - pos vs nt - omicron.rds")

### Full df (with propensity scores, weights etc.)
df_nt_omicron <- omicron_matching_nt$all_patients
saveRDS(df_nt_omicron, "Linked df with propensity scores - pos vs nt  - omicron.rds")

## Remove files
rm(omicron_matching_nt, df_nt_m_omicron, df_nt_omicron)
