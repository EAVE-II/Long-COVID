#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 09_long_covid_indicators - Descriptive analysis of individuals classified as having long COVID according to
## our operational definition, long-COVID Read codes, free text, and sick note free text.
## Plots counts by month for each outcome variable and prepares a table showing the demographic breakdown.
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

# Set location to export figures to
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/6. Long COVID indicators")


# 1. Prepare data for plot of new counts of each outcome variable by month ----
## Set study start
start <- ymd("2020-03-01")
end <- ymd("2022-10-20")

## Read in EAVE_LINKNO refresh file from PHS (use this to clean out erroneous EAVE_LINKNOs and update EAVE_LINKNOs)
EAVE_raw <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_LINKNO_refresh.rds")

## Keep individuals who:
EAVE <- EAVE_raw %>% # 6,988,571
  ## have a validated CHI number
  filter(unvalidatedCHI_flag != 1) %>% #6,978,344
  ## were alive at the start of the study
  filter(is.na(NRS.Date.Death) | NRS.Date.Death > start) %>% #6,972,704
  ## did not transfer in after the study start date
  filter(is.na(DATE_TRANSFER_IN) | DATE_TRANSFER_IN < start) %>%  #6,476,394
  ## did not transfer out before or during the study period
  filter(is.na(DATE_TRANSFER_OUT) | DATE_TRANSFER_OUT > end) %>% # 6,143,690
  ## Remove duplicates
  distinct(EAVE_LINKNO_old, .keep_all = T) %>% 
  distinct(EAVE_LINKNO, .keep_all = T) # 6,116,287

## Read in GP data supplied by Albasoft 
GP_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/longcovid_261022.rds") # full cohort, extracted 20/10/22
  
## Remove erroneous EAVE_LINKNOs, update gp to use new EAVE_LINKNOs, and remove under 18s
gp <- GP_raw %>%  
  filter(EAVE_LINKNO %in% EAVE$EAVE_LINKNO_old) %>% 
  left_join(dplyr::select(EAVE, "EAVE_LINKNO" = EAVE_LINKNO_old, "EAVE_LINKNO_new" = EAVE_LINKNO), by = "EAVE_LINKNO") %>% 
  mutate("old_EAVE_LINKNO" = EAVE_LINKNO, # keep old EAVE_LINKNO
         "EAVE_LINKNO" = EAVE_LINKNO_new) %>% 
  dplyr::select(-EAVE_LINKNO_new) %>% 
  mutate(chi_age = as.numeric(chi_age)) %>% 
  filter(chi_age >= 18 & chi_age <=110) 

denominator <- length(unique(gp$EAVE_LINKNO)) # 5,139,041

rm(EAVE_raw, EAVE)

# Clean up and keep info on variables of interest only
gp <- GP_raw %>% 
  mutate(gp_dat = ymd(EventDate)) %>% 
  filter(gp_dat >= ymd(start) & gp_dat <= ymd(end)) %>% 
  dplyr::select(EAVE_LINKNO, gp_dat, group) %>% 
  filter(group == "Free text" | group == "Free text - fitnote" | group == "Long Covid") %>% 
  distinct(EAVE_LINKNO, gp_dat, group) # remove duplicate entries

# Remove GP_raw 
rm(GP_raw)


# 2. Incorporate operational definition  ----
opdef <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/6. Long COVID indicators/opdef_classification.rds") # PCR or LFT


opdef_1 <- opdef %>%   # Format to match gp data
  rename("gp_dat" = tes_dat, "group" = op_def) %>% 
  dplyr::select(-op_def_no_blood) %>% 
  mutate(group = "Operational definition") %>% 
  ## Add 28 days to gp_dat for operational definition variables (since they are based on tes_dat) 
  ## Do this to remove some of the lag between opdef and other outcome vars when plotting
  mutate(gp_dat = gp_dat + days(28))

opdef_2 <- opdef %>%   # Format to match gp data
  rename("gp_dat" = tes_dat, "group" = op_def_no_blood) %>% 
  dplyr::select(-op_def) %>% 
  filter(group == 1) %>% 
  mutate(group = "Operational definition - no blood tests") %>% 
  ## Add 28 days to gp_dat for operational definition variables (since they are based on tes_dat) 
  ## Do this to remove some of the lag between opdef and other outcome vars when plotting
  mutate(gp_dat = gp_dat + days(28))

gp <- gp %>% bind_rows(opdef_1, opdef_2)

rm(opdef)

## Get counts of individuals' first instance of each indicator, by month 
first_indicator <- gp %>%
  # Keep only the first occurence of any outcome measure for each individual
  arrange(gp_dat) %>% 
  group_by(EAVE_LINKNO, group) %>% 
  slice(1) %>% 
  ungroup() 

## Export for use in prediction modelling
saveRDS(first_indicator, "Outcome measure dates.rds")

first_indicator <- first_indicator %>% 
  # Get counts per month
  dplyr::select(-EAVE_LINKNO) %>% 
  mutate(gp_dat = as.yearmon(gp_dat, "%m/%Y"))%>% 
  group_by(gp_dat, group) %>% 
  mutate(n= n()) %>% 
  arrange(gp_dat) %>% 
  distinct(gp_dat, group, .keep_all=T) %>% 
  rename("Date" = gp_dat) %>% 
  filter(Date != "Oct 2022") # remove final month with incomplete data


## Get counts of individuals' first validation measure ('Long-COVID code, free text or sick note'), by month 
first_valid <- gp %>%
  filter(group == "Free text" | group == "Free text - fitnote" | group == "Long Covid") %>% 
  arrange(gp_dat) %>% 
  group_by(EAVE_LINKNO) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(-EAVE_LINKNO) %>% 
  mutate(gp_dat = as.yearmon(gp_dat, "%m/%Y"))%>% 
  group_by(gp_dat) %>% 
  mutate(n= n()) %>% 
  arrange(gp_dat) %>% 
  distinct(gp_dat, .keep_all=T) %>% 
  mutate(group = 'Long-COVID code, free text or sick note') %>% 
  rename("Date" = gp_dat) %>% 
  filter(Date != "Oct 2022") # remove final month with incomplete data


## Get counts of individuals' first instance of any measure, by month 
first_indicator_or_valid <- gp %>%
  # Keep only the first occurence of any outcome measure for each individual
  arrange(gp_dat) %>% 
  group_by(EAVE_LINKNO) %>% 
  slice(1) %>% 
  ungroup() %>% 
  # Get counts per month
  dplyr::select(-EAVE_LINKNO) %>% 
  mutate(gp_dat = as.yearmon(gp_dat, "%m/%Y"))%>% 
  group_by(gp_dat) %>% 
  mutate(n= n()) %>% 
  arrange(gp_dat) %>% 
  distinct(gp_dat, .keep_all=T) %>% 
  mutate(group = "Any outcome measure") %>% 
  rename("Date" = gp_dat) %>% 
  filter(Date != "Oct 2022") # remove final month with incomplete data

individuals_data <- rbind(first_indicator, first_valid, first_indicator_or_valid)

individuals_data <- individuals_data %>% arrange(Date)

rm(first_indicator, first_valid, first_indicator_or_valid)

# Recode 'group' using labels that will show in plot legends
individuals_relabelled <- individuals_data %>% 
  filter(group != "Operational definition - no blood tests") %>% 
  mutate(group = recode(group, 
                        'Free text'='Long-COVID in free text', 
                        'Free text - fitnote'='Long-COVID on sick note',
                        'Long Covid' = 'Long-COVID clinical code'))

## Set the order that outcomes will appear in the plot legend
individuals_relabelled$group <- factor(individuals_relabelled$group, 
                                       levels = c("Any outcome measure",
                                                  'Operational definition',
                                                  'Operational definition (PCR test)',
                                                  'Long-COVID code, free text or sick note',
                                                  'Long-COVID on sick note',
                                                  'Long-COVID in free text', 
                                                  'Long-COVID clinical code'))

# Export plot data
write_csv(individuals_relabelled, "plot data.csv")

# Function to produce plots
make_plot <- function(plot_data, plot_name){
  
  plot <- ggplot(plot_data, aes(x=Date, y=n, group=group, color=group)) +
    geom_line()+
    theme_bw() +
    scale_x_yearmon(n = 42) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "") +
    theme(legend.title=element_blank()) +
    theme(legend.position = "right") +
    geom_vline(xintercept = as.yearmon("Apr 2022", "%b %Y"), linetype="dashed", color = "dark grey") + # add vertical dashed line
    labs(vline="Testing Withdrawn (April 2022)") +
    scale_y_continuous("New cases per month", label = comma, 
                       expand = c(0, 0), limits = c(0, max(plot_data$n) + 100)) # to ensure y-axis starts at 0
    #scale_y_log10() # make axis logarithmic +
    
  
  
  ggsave(paste0(plot_name, ".png"), plot = plot, width = 12, height = 6, units = "in", dpi = 300)
  
  return(plot)
}

## Plot all variables
### Select data
all_vars <- individuals_relabelled %>%
  filter(group == 'Long-COVID in free text' | 
         group == 'Long-COVID on sick note' |
         group == 'Long-COVID clinical code' |
         group == 'Operational definition' |
         #group == 'Operational definition (PCR test)' |
         group == 'Long-COVID code, free text or sick note' |
         group == "Any outcome measure")

### Generate and export plot
all_vars_plot <- make_plot(plot_data = all_vars,
                           plot_name = "All variables")
all_vars_plot


# 3. Descriptive stats: Prepare data ----

## 3.1 Prepare counts of outcome measures for the full cohort ----
## Collapse data on long-COVID indicators for each indvidual to counts, then join to linked data to allow analysis by demographic groups
gp_counts <- gp %>% 
  # Add counts of Read codes by patient and period
  group_by(EAVE_LINKNO, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(n = as.numeric(n)) %>%
  ## Transpose
  pivot_wider(
    names_from = group,
    values_from = n,
    values_fill = list(n = 0))

## Read in linked data 
## NB using the 'uncensored' linked data containing demographic variables for the full cohort (including never tested)
df_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_linked_LFT_uncensored.rds")

## There is a single duplicate - remove
df_raw <- df_raw %>% distinct(EAVE_LINKNO, .keep_all = T) # N = 5,094,522

## Join counts of long-COVID indicators to linked data and fill NAs with 0s and create an overall indicator
df <- df_raw %>% 
  left_join(gp_counts, by = "EAVE_LINKNO") %>% 
  mutate(`Long-COVID on sick note` = ifelse(is.na(`Free text - fitnote`), 0, `Free text - fitnote`),
         `Long-COVID in free text` = ifelse(is.na(`Free text`), 0, `Free text`),
         `Long-COVID clinical code` = ifelse(is.na(`Long Covid`), 0, `Long Covid`),
         `Long-COVID code, free text or sick note` = ifelse(`Long-COVID on sick note` > 0 | `Long-COVID in free text` > 0 | `Long-COVID clinical code` > 0, 1, 0),
         `Operational definition` = ifelse(is.na(`Operational definition`), 0, `Operational definition`),
         `Operational definition - no blood tests` = ifelse(is.na(`Operational definition - no blood tests`), 0, `Operational definition - no blood tests`),
         `Any outcome measure` = ifelse(`Operational definition` > 0 |  `Long-COVID code, free text or sick note` > 0, 1, 0)) %>% 
  dplyr::select(-c(`Free text - fitnote`, `Free text`, `Long Covid`))#, `Operational definition (PCR)`))


## 3.2 Identify positive cases by variant period ----
## These stats will be presented as % of all positive cases e.g. % of positive with a long COVID read code etc.

## Get periods where the variant represents >60% of positive cases, based on analysis of whole genome sequencing data in 02_variant_periods)
wild_start <- ymd("2020-03-01")
wild_end <-   ymd("2021-01-10")

alpha_start <-  ymd("2021-01-11")
alpha_end <-  ymd("2021-05-09")

delta_start <-  ymd("2021-05-24")
delta_end <-  ymd("2021-12-05")

omicron_start <-  ymd("2021-12-27")
omicron_end <-  ymd("2022-04-30") # widespread PCR testing ended

## For individuals with a positive test, identify variant period
df <- df %>% 
  mutate(var_per = ifelse(is.na(tes_dat) | tes_res == "NEGATIVE", NA,
                          ifelse(tes_dat >= wild_start & tes_dat <= wild_end, "Wild",
                              ifelse(tes_dat >= alpha_start & tes_dat <= alpha_end, "Alpha",
                                 ifelse(tes_dat >= delta_start & tes_dat <= delta_end, "Delta",
                                        ifelse(tes_dat >= omicron_start & tes_dat <= omicron_end, "Omicron",
                                               NA))))))


## 3.3 Vaccination doses ----
## This will be presented as number of doses up to 14 days before each outcome measure was recorded (won't be presented for full population)
vac_raw <- readRDS("/conf/EAVE/GPanalysis/data/temp/vaccine_cleaned.rds")  # cleaned file prepared by Vera @ PHS

## Prepare vac data
vac <- vac_raw %>% 
  ## Select vars of interest
  dplyr::select(c(EAVE_LINKNO, starts_with("date_vacc")))

## Pivot longer
vac <- pivot_longer(vac,
                    cols = date_vacc_1:date_vacc_sb,
                    values_to = "vac_dat")

## Clean up
vac <- vac %>% 
  dplyr::select(-name) %>% 
  filter(!is.na(vac_dat)) %>% 
  filter(vac_dat <= "2022-10-20")

# Count vaccine doses by 14 days before tes_dat (if individual is identified by operational defintion), 
# or by 14 days before outcome date (if identified by clinical code/free text)

## Get date of first recorded outcome measure (opdef/free text etc.)
gp_dat <- gp %>% 
  arrange(gp_dat) %>% 
  distinct(EAVE_LINKNO, .keep_all = T)

vac_n <- vac %>% 
  filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
  left_join(dplyr::select(gp_dat, EAVE_LINKNO, gp_dat), by = "EAVE_LINKNO") %>% 
  filter(vac_dat < (gp_dat - days(14))) %>% 
  group_by(EAVE_LINKNO) %>% 
  summarise(vac_dos = n())

# Replace in df
df <- df %>% 
  left_join(vac_n, by = "EAVE_LINKNO") %>% 
  mutate(vac_dos = ifelse(is.na(vac_dos), 0L, vac_dos))  

rm(vac_raw, vac, vac_n)


## 3.4 Severity of acute infection (admitted to hospital or ICU within 28 days of testing positive) ----
### Hospitalisations (1/3/2019 to latest test date)
hos <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/SMR01 - 2019-03-01 to 2022-09-21.rds")

### ICU admissions (1/3/19 to latest test date)
icu <- readRDS("/conf/EAVE/GPanalysis/data/SICSAG_episode_level_.rds") %>% 
  filter(!is.na(EAVE_LINKNO)) %>% 
  dplyr::select(EAVE_LINKNO, hos_dat = AdmitUnit)

## Get EAVE_LINKNOs admitted to hospital or ICU within 28 days of positive test
severe <- hos %>% 
  bind_rows(icu) %>% 
  ## Get hospitalisations/ICU admissions within 28 days of positive test
  left_join(dplyr::select(df, EAVE_LINKNO, tes_res, tes_dat), by = "EAVE_LINKNO") %>% 
  filter(tes_res == "POSITIVE") %>% 
  filter(hos_dat >= tes_dat & hos_dat <= (tes_dat + days(28))) %>% 
  distinct(EAVE_LINKNO)

## Flag severe cases in df
df$severe <- ifelse(df$EAVE_LINKNO %in% severe$EAVE_LINKNO, 1, 0)


# 4. Descriptive stats - table ----
## Calculate descriptives
calculate_descriptives <- function(df, 
                                   mask, # mask = the subset of interest e.g. df$`Operational definition` > 0
                                   vac_dat, # vac_dat = logical - get stats for no. of vaccinations or not (don't estiamte these for the full population)
                                   sev_dat){ # sev_dat = logical - include severity of accute infection data or not (only possible for the subset with a +ve test)
                                   
  ## Get population total
  pop_n <- sum(df$eave_weight) # NB weighting is no longer necessary so these are all set to 1
  
  ## Subset data of interest
  df <- df %>% filter(mask)
  
  ## Get N
  ### Sex 
  sex_tab <- round(wtd.table(df$sex, exclude = NULL, weights = df$eave_weight), 0)
  sex_n <- ""
  fem_n <- sex_tab[1]
  mal_n <- sex_tab[2]
  
  ### Age
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
  
  ### Test result (PCR or LFT)
  tes_res_tab <- round(wtd.table(df$tes_res, exclude = NULL, weights = df$eave_weight), 0)
  tes_res_n <- ""
  pos_n <- tes_res_tab["POSITIVE"]
  neg_n <- tes_res_tab["NEGATIVE"]
  nt_n  <- tes_res_tab["NEVER TESTED"]
  
  ### SIMD
  simd_tab <- round(wtd.table(df$simd, exclude = NULL, weights = df$eave_weight), 0)
  simd_n <- ""
  simd1_n <- simd_tab[1]
  simd2_n <- simd_tab[2]
  simd3_n <- simd_tab[3]
  simd4_n <- simd_tab[4]
  simd5_n <- simd_tab[5]
  simd6_n <- simd_tab[6]
  
  ### Urban-rural classification
  ur_tab <- round(wtd.table(df$ur6_lab, exclude = NULL, weights = df$eave_weight), 0)
  ur_n <- ""
  ur1_n <- ur_tab[1]
  ur2_n <- ur_tab[2]
  ur3_n <- ur_tab[3]
  ur4_n <- ur_tab[4]
  ur5_n <- ur_tab[5]
  ur6_n <- ur_tab[6]
  ur7_n <- ur_tab[7]
  
  ### Household size
  hh_tab <- round(wtd.table(df$hhold_raw_bands, weights = df$eave_weight, useNA = "always"), 0)
  hh_n <- ""
  hh1_n <- hh_tab[1] # 1
  hh2_n <- hh_tab[2] # 2
  hh3_n <- hh_tab[3] # 3-5
  hh4_n <- hh_tab[4] # 6-10
  hh5_n <- sum(hh_tab[5:7]) # 11+
  hh6_n <- hh_tab[8] # Unknown
  
  ### BMI - not using imputed values (38%) 
  df <- df %>% 
    mutate(bmi_bands = cut(bmi, breaks = c(0, 18.5, 24.9, 29.9, 80), labels=c("Underweight","Normal weight","Overweight","Obese")),
           bmi_bands = ifelse(is.na(bmi_bands), "Unknown", bmi_bands))
  
  bmi_tab <- round(wtd.table(df$bmi_bands, exclude = NULL, weights = df$eave_weight), 0)
  bmi_n <- ""
  bmi1_n <- bmi_tab[1] # Underweight
  bmi2_n <- bmi_tab[2] # Normal weight
  bmi3_n <- bmi_tab[3] # Overweight
  bmi4_n <- bmi_tab[4] # Obese
  bmi5_n <- bmi_tab[5] # Unknown
  
  ### Comorbidities/Risk groups
  com_tab <- round(wtd.table(df$n_risk_gps, exclude = NULL, weights = df$eave_weight), 0)
  com_n <- ""
  com1_n <- com_tab[1] # 0
  com2_n <- com_tab[2] # 1
  com3_n <- com_tab[3]# 2
  com4_n <- sum(com_tab[4:5]) # 3+
  
  ### Sheilding 
  she_tab <- round(wtd.table(df$shielding, exclude = NULL, weights = df$eave_weight), 0)
  she_n <- ""
  she1_n <- she_tab[2] # sheilding
  she2_n <- she_tab[1] # not shielding
  
  # Immunosuppressed
  imm_tab <- round(wtd.table(df$immuno_supp, exclude = NULL, weights = df$eave_weight), 0)
  imm_n <- ""
  imm1_n <- imm_tab[2] # immunosuppressed
  imm2_n <- imm_tab[1] # not immunosuppressed
  
  ### Vaccine doses received by 14 days before test date -- can't be included for the full sample because we don't have test dates for nt
  if(vac_dat == T){
    
    vac_tab <- round(wtd.table(df$vac_dos, exclude = NULL, weights = df$eave_weight), 0)
    
    vac_n <- ""
    vac1_n <- vac_tab[1] # 0
    vac2_n <- vac_tab[2] # 1
    vac3_n <- vac_tab[3] # 2
    vac4_n <- vac_tab[4] # 3
    vac5_n <- sum(vac_tab[5:6]) # 4+
    vac6_n <- vac_tab[7] # Unknown
  }
  
  if(vac_dat == F){
    vac_n <- ""
    vac1_n <- ""
    vac2_n <- ""
    vac3_n <- ""
    vac4_n <- ""
    vac5_n <- ""
    vac6_n <- ""
  }
  
  
  ### Variant (based on test date)
  if(sev_dat == T){
    
    # Restrict to the subset of individuals with a positive test result
    var_tab <- round(wtd.table(df$var_per[df$tes_res=="POSITIVE"], exclude = NULL, weights = df$eave_weight[df$tes_res=="POSITIVE"]), 0)
    
    var_n <- ""
    var1_n <- var_tab[4] # Wild
    var2_n <- var_tab[1] # Alpha
    var3_n <- var_tab[2] # Delta
    var4_n <- var_tab[3] # Omicron
    
  }
  
  if(sev_dat == F){
    var_n <- ""
    var1_n <- ""
    var2_n <- ""
    var3_n <- ""
    var4_n <- ""
    
  }
  
  ### Severity of acute infection
  if(sev_dat == T){
    
    # Restrict to the subset of individuals with a positive test result
    sev_tab <- round(wtd.table(df$severe[df$tes_res=="POSITIVE"], exclude = NULL, weights = df$eave_weight[df$tes_res=="POSITIVE"]), 0)
    
    sev_n <- ""
    sev1_n <- sev_tab[2] # Hospitalised or admitted to ICU within 28 days of positive test
    sev2_n <- sev_tab[1] # Not hospitalised or admitted to ICU within 28 days of positive test
    
  }
  
  if(sev_dat == F){
    sev_n <- ""
    sev1_n <- ""
    sev2_n <- ""
    
  }
  
  ### Total
  tot_n <- round(sum(df$eave_weight), 0)
  
  N <- rbind(tot_n,
             sex_n, fem_n, mal_n,
             age_n, age1_n, age2_n, age3_n, age4_n, age5_n, age6_n, age7_n, age8_n,
             tes_res_n, pos_n, neg_n, nt_n,
             simd_n, simd1_n, simd2_n, simd3_n, simd4_n, simd5_n, simd6_n, 
             ur_n, ur1_n, ur2_n, ur3_n, ur4_n, ur5_n, ur6_n, ur7_n, 
             hh_n, hh1_n, hh2_n, hh3_n, hh4_n, hh5_n, hh6_n,
             bmi_n, bmi1_n, bmi2_n, bmi3_n, bmi4_n, bmi5_n,
             com_n, com1_n, com2_n, com3_n, com4_n, 
             she_n, she1_n, she2_n,
             imm_n, imm1_n, imm2_n,
             vac_n, vac1_n, vac2_n, vac3_n, vac4_n, vac5_n, vac6_n,
             var_n, var1_n, var2_n, var3_n, var4_n,
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
  labels <- rbind("Total (% of population)",
                  "Sex", "Female", "Male",
                  "Age", "18 - 27", "28 - 37", "38 - 47", "48 - 57", "58 - 67", "68 - 77", "78 - 87", "88 - 100",  
                  "Testing (PCR or LFT)", "Positive", "Negative (and not positive)", "No tests recorded",
                  "SIMD quintiles", "1 - Most deprived", "2", "3", "4", "5 - Least deprived", "Unknown", 
                  "Urban-Rural", "Large urban areas", "Other urban areas", "Accessible small towns", "Remote small towns", "Accessible rural", "Remote rural", "Unknown", 
                  "Household size", "1", "2", "3-5", "6-10", "11+", "Unknown",
                  "BMI", "Underweight (BMI < 18.5)", "Normal weight (BMI 18.5 - 24.9)", "Overweight (BMI 25 - 29.9)", "Obese (BMI >29.9)", "Unknown",
                  "Comorbidities", "0", "1", "2", "3+",
                  "Shielding", "Shielding", "Not shielding", 
                  "Immunosuppressed", "Immunosuppressed", "Not immunosuppressed",
                  "Vaccination doses (up to 14 days before positive test/outcome)" , "0", "1", "2", "3", "4+", "Unknown",
                  "Variant period (positive cases only)", "Wild (up to 10/01/2021)", "Alpha (11/01/2021 - 09/05/2021)", "Delta (24/05/2021 - 28/11/2021)", "Omicron (20/12/2021 onwards)",
                  "Severity of acute infection (positive cases only)", "Hospitalised within 28 days", "Not hospitalised within 28 days")
  
  ## Combine
  results <- cbind(labels, N, PC)
  results$N <- ifelse(is.na(results$N), "", results$N)
  results$PC <- gsub("NA", "", results$PC)
  
  ## Use positive cases denominator for variant and severity rows
  denom_pos = round(sum(df$eave_weight[df$tes_res=="POSITIVE"]), 0)
  
  for(row in 1:nrow(results)){
    results$PC[row] <- ifelse(results$labels[row] %in% c("Wild (up to 10/01/2021)", "Alpha (11/01/2021 - 09/05/2021)", "Delta (24/05/2021 - 28/11/2021)", 
                                                         "Omicron (20/12/2021 onwards)", "Hospitalised within 28 days", "Not hospitalised within 28 days"), 
                              format(round(as.numeric(results$N[row])/denom_pos*100, 1), nsmall=1L), results$PC[row])
  }
  
  
  
  
  return(results)
}

## Get descriptives
pop_summary <- data.frame(calculate_descriptives(df, mask = !is.na(df$EAVE_LINKNO), vac_dat = F, sev_dat = T))
any_summary <- data.frame(calculate_descriptives(df, mask = df$`Any outcome measure` > 0, vac_dat = T, sev_dat = T))
cod_summary <- data.frame(calculate_descriptives(df, mask = df$`Long-COVID clinical code`> 0, vac_dat = T, sev_dat = T))
txt_summary <- data.frame(calculate_descriptives(df, mask = df$`Long-COVID in free text` > 0, vac_dat = T, sev_dat = T))
fit_summary <- data.frame(calculate_descriptives(df, mask = df$`Long-COVID on sick note` > 0, vac_dat = T, sev_dat = T))
def_summary <- data.frame(calculate_descriptives(df, mask = df$`Operational definition` > 0, vac_dat = T, sev_dat = T))


descr_stats <- cbind(pop_summary, 
                     any_summary[2:3],
                     cod_summary[2:3],
                     txt_summary[2:3],
                     fit_summary[2:3],
                     def_summary[2:3])

colnames(descr_stats) <- c("", "Full sample (N)", "Full sample (%)", 
                          "Any outcome measure (N)", "Any outcome measure (%)",
                          "Long-COVID clinical code (N)", "Long-COVID clinical code (%)",
                          "Long-COVID in free text (N)", "Long-COVID in free text (%)",
                          "Long-COVID on sick note (N)", "Long-COVID on sick note (%)",
                          "Operational definition (N)", "Operational definition (%)")


View(descr_stats)

# Export
write_csv(descr_stats, "Summary stats.csv")



## 5. Q stats - table ----
## Calculate descriptives
Q_descriptives <- function(df, mask){ # mask = the subset of interest e.g. df$`Operational definition` > 1 
                           
                                   
  ## Get population total
  pop_n <- sum(df$eave_weight)
  
  ## Subset data of interest (i.e. remove those not clinically interesting with respect to long COVID)
  df <- df %>% filter(mask) %>% 
    dplyr::select(-c("Q_DIAG_CEREBRALPALSY", "Q_DIAG_CIRRHOSIS", "Q_DIAG_CONGEN_HD", "Q_DIAG_SICKLE_CELL"))
  
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
  labels <- "Total (% of population)"
  
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
pop_Q <- data.frame(Q_descriptives(df, mask = !is.na(df$EAVE_LINKNO)))
any_Q <- data.frame(Q_descriptives(df, mask = df$`Any outcome measure` > 0))
cod_Q <- data.frame(Q_descriptives(df, mask = df$`Long-COVID clinical code`> 0))
txt_Q <- data.frame(Q_descriptives(df, mask = df$`Long-COVID in free text` > 0))
fit_Q <- data.frame(Q_descriptives(df, mask = df$`Long-COVID on sick note` > 0))
def_Q <- data.frame(Q_descriptives(df, mask = df$`Operational definition` > 0))


Q_stats <- cbind(pop_Q, 
                 any_Q[2:3],
                 cod_Q[2:3],
                 txt_Q[2:3],
                 fit_Q[2:3],
                 def_Q[2:3])

colnames(Q_stats) <- colnames(descr_stats)


View(Q_stats)

# Export
write_csv(Q_stats, "Summary stats - Q.csv")


# 6.1 Health board table ----
## Get hb labels
df <- df %>% 
  mutate(hb_lab = case_when(hb == "S08000015" ~ "NHS Ayrshire and Arran",
                        hb == "S08000016" ~ "NHS Borders",
                        hb == "S08000017" ~ "NHS Dumfries and Galloway",
                        hb == "S08000029" ~ "NHS Fife",
                        hb == "S08000019" ~ "NHS Forth Valley",
                        hb == "S08000020" ~ "NHS Grampian" ,
                        hb == "S08000031" ~ "NHS Greater Glasgow and Clyde",
                        hb == "S08000022" ~ "NHS Highland",
                        hb == "S08000032" ~ "NHS Lanarkshire",
                        hb == "S08000024" ~ "NHS Lothian",
                        hb == "S08000025" ~ "NHS Orkney",
                        hb == "S08000026" ~ "NHS Shetland",
                        hb == "S08000030" ~ "NHS Tayside",
                        hb == "S08000028" ~ "NHS Western Isles",
                        hb == "Unknown" ~ "Unknown"))

hbs <- c("NHS Ayrshire and Arran" , "NHS Borders", "NHS Dumfries and Galloway", "NHS Fife", "NHS Forth Valley", "NHS Grampian",
         "NHS Greater Glasgow and Clyde", "NHS Highland", "NHS Lanarkshire", "NHS Lothian", "NHS Orkney", "NHS Shetland", "NHS Tayside",
         "NHS Western Isles", "Unknown")

## Calculate descriptives
hb_descriptives <- function(df, mask){ # mask = the subset of interest e.g. df$`Operational definition` > 1 
  
  ## Get population total
  pop_n <- sum(df$eave_weight)
  
  ## Subset data of interest
  df <- df %>% filter(mask)
  
  ## Get N by hb
  hb_tab <- table(factor(df$hb_lab, levels = hbs))
  hb1_n<-  hb_tab[1]
  hb2_n<-  hb_tab[2]
  hb3_n<-  hb_tab[3]
  hb4_n<-  hb_tab[4]
  hb5_n<-  hb_tab[5]
  hb6_n<-  hb_tab[6]
  hb7_n<-  hb_tab[7]
  hb8_n<-  hb_tab[8]
  hb9_n<-  hb_tab[9]
  hb10_n<- hb_tab[10]
  hb11_n<- hb_tab[11]
  hb12_n<- hb_tab[12]
  hb13_n<- hb_tab[13]
  hb14_n<- hb_tab[14]
  hb15_n<- hb_tab[15]
  
  ### Total
  tot_n <- round(sum(df$eave_weight), 0)

  N <- rbind(tot_n,
             hb1_n, hb2_n, hb3_n, hb4_n, hb5_n, hb6_n, hb7_n, hb8_n, hb9_n, hb10_n,
             hb11_n, hb12_n, hb13_n, hb14_n, hb15_n)
  
  ## Suppress small numbers
  N <- as.numeric(N)
  N[!is.na(N) & N<5] <- -5
  
  ## Get %
  PC <- as_tibble(format(round(as.numeric(N)/tot_n*100, 1), nsmall=1L))
  colnames(PC) <- "PC"
  
  ## For total row, make % of total population
  PC$PC[1] <- round(as.numeric(N[1])/pop_n*100, 1)
  
  ## Labels
  labels <- c("Total (% of population)", row.names(hb_tab))
  
  ## Combine
  results <- cbind(labels, N, PC)
  results$N <- ifelse(is.na(results$N), "", results$N)
  results$PC <- gsub("NA", "", results$PC)
  
  return(results)
}

## Get descriptives
pop_hb <- data.frame(hb_descriptives(df, mask = !is.na(df$EAVE_LINKNO)))
any_hb <- data.frame(hb_descriptives(df, mask = df$`Any outcome measure` > 0))
cod_hb <- data.frame(hb_descriptives(df, mask = df$`Long-COVID clinical code`> 0))
txt_hb <- data.frame(hb_descriptives(df, mask = df$`Long-COVID in free text` > 0))
fit_hb <- data.frame(hb_descriptives(df, mask = df$`Long-COVID on sick note` > 0))
def_hb <- data.frame(hb_descriptives(df, mask = df$`Operational definition` > 0))


hb_stats <- cbind(pop_hb, 
                  any_hb[2:3],
                  cod_hb[2:3],
                  txt_hb[2:3],
                  fit_hb[2:3],
                  def_hb[2:3])

colnames(hb_stats) <- c("", "Full sample (N)", "Full sample (%)", 
                           "Any outcome measure (N)", "Any outcome measure (%)",
                           "Long-COVID clinical code (N)", "Long-COVID clinical code (%)",
                           "Long-COVID in free text (N)", "Long-COVID in free text (%)",
                           "Long-COVID on sick note (N)", "Long-COVID on sick note (%)",
                           "Operational definition (N)", "Operational definition (%)")


View(hb_stats)

# Export
write_csv(hb_stats, "Healthboard stats.csv")


# 7. Export df (for use in prediction modelling) ----
saveRDS(df, "/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_predict.rds", compress = TRUE)


# 8. Validation ----
## How well does our definition identify individuals who we're confident have/have had long-COVID? 
round(prop.table(table(df$`Operational definition`, df$`Long-COVID code, free text or sick note`))*100,2)
# True negatives: 98.3%
# True positives: 0.1%
# Op def, not validation: 1.2%
# Validation, not op def: 0.4%

## Check for opdef without blood tests
round(prop.table(table(df$`Operational definition - no blood tests` > 0, df$`Long-COVID code, free text or sick note`))*100,2)
# True negatives: 98.3%
# True positives: 0.1%
# Op def, not validation: 1.2%
# Validation, not op def: 0.4%


## As above, but for working age adults
round(prop.table(table(df$`Operational definition`[df$age<=66], df$`Long-COVID code, free text or sick note`[df$age<=66]))*100,2)
# True negatives: 98.1%
# True positives: 0.2%
# Op def, not validation: 1.3%
# Validation, not op def: 0.4%

## Validation against Long COVID read code OR Long COVID term in GP free text OR term in sick note free text
print("All validation measures")
df_valid <- df[df$`Long-COVID code, free text or sick note`> 0,]
round(prop.table(table(df_valid$`Operational definition` > 0))*100,1) # 27.0% are identified by the operational defintion

## Validation against Long COVID Read code
print("Long COVID Read code")
df_lc_read <- df[df$`Long-COVID clinical code` > 0,]
round(prop.table(table(df_lc_read$`Operational definition` > 0))*100,1) # 19.3% are identified by the operational defintion

## Validation against Long COVID term in GP free text
print("Long COVID in free text")
df_text <- df[df$`Long-COVID in free text` > 0,]
round(prop.table(table(df_text$`Operational definition` > 0))*100,1) # 21.6% are identified by the operational defintion

## Validation against Long COVID term in sick note free text
print("Long COVID in fitnote free text")
df_note <- df[df$`Long-COVID on sick note` > 0,]
round(prop.table(table(df_note$`Operational definition` > 0))*100,1) # 31.3% are identified by the operational defintion


# What's the overlap between our operational definition and our concrete measures of long-COVID?
print("Classified as having long COVID AND have a validation measure")
df_op_def <- df[df$`Operational definition` > 0,]
round(prop.table(table(df_op_def$`Long-COVID code, free text or sick note` > 0))*100,1) # 9.8% of those identified by the operational defintion are identified
                                                                                        # by another outcome measure

# 9. Validation plot (Sid's code) ----
## Prepare data for Sid
val_plot_df <- df %>% 
  dplyr::select(EAVE_LINKNO, old_EAVE_LINKNO, council, hb, `Long-COVID clinical code`, `Long-COVID in free text`, `Long-COVID on sick note`, `Operational definition`) %>% 
  mutate(`Any outcome measure` = ifelse(`Long-COVID clinical code` > 0 |
                                        `Long-COVID in free text`> 0 |
                                        `Long-COVID on sick note`> 0 |
                                        `Operational definition`> 0, 1, 0))


# Sid's code
library(networkD3)
library(tidygraph)
library(ggalluvial)
library(igraph)

# Read in df with counts of outcomes in each period  
# * note this read in positive and negative cases only
dfRAW <- val_plot_df 
rm(val_plot_df)
#df <- dfRAW %>% select(`Long-COVID clinical code`, `Long-COVID in free text`, `Long-COVID on sick note`, `Operational definition`, `Any outcome measure`)
df <- dfRAW
df$`Long-COVID clinical code`[df$`Long-COVID clinical code` > 1] <- 1
df$`Long-COVID in free text`[df$`Long-COVID in free text` > 1] <- 1
df$`Long-COVID on sick note`[df$`Long-COVID on sick note` > 1] <- 1
df$`Operational definition`[df$`Operational definition` > 1] <- 1
df$`Any outcome measure`[df$`Any outcome measure` > 1] <- 1

dfgraph <- df %>% 
  dplyr::select(`Long-COVID clinical code`, `Long-COVID in free text`, `Long-COVID on sick note`, `Operational definition`, `Any outcome measure`) #%>% 
#   count(`Long-COVID clinical code`, `Long-COVID in free text`, `Long-COVID on sick note`, `Operational definition`, `Any outcome measure`)
# colSums(df)

# create adjacency matrix and then plot with circlize library
nc = ncol(dfgraph)
nr = nrow(dfgraph)
m3b <- matrix(0, nrow=nc, ncol=nc)
for(i in seq(1,nc)){
  for(j in seq(1,nc)){
    t3 <- table(dfgraph[,i], dfgraph[,j])
    m3b[i,j] = t3[2,2]
  }
}
#m3b <- m3b[,-(1:2)]/1e03
m3df <- data.frame(
  order=1:5,
  classes = c("Long-COVID clinical code", "Long-COVID in free text", "Long-COVID on sick note", "Operational definition", "Any outcome measure"), 
  m3b,
  r = c(255,255,255,153, 51),
  g = c(51, 153, 255, 255, 255),
  b = c(51, 51, 51, 51, 153),
  stringsAsFactors = FALSE)

refdf <- m3df[, c(1,2, 8:10)]

dimnames(m3b) <- list(orig=m3df$classes, dest = m3df$classes)
#Sort order of data.frame and matrix for plotting in circos
m3df <- arrange(m3df, order)
m3df$classes <- factor(m3df$classes, levels = m3df$classes)
m3b <- m3b[levels(m3df$classes),levels(m3df$classes)]

### Define ranges of circos sectors and their colors (both of the sectors and the links)
m3df$xmin <- 0
m3df$xmax <- rowSums(m3b) + colSums(m3b)
n <- nrow(m3df)
m3df$rcol<-rgb(m3df$r, m3df$g, m3df$b, max = 255)
m3df$lcol<-rgb(m3df$r, m3df$g, m3df$b, alpha=200, max = 255)
border_mat <- matrix("red", nrow = 1, ncol = ncol(m3b))
rownames(border_mat) = rownames(m3b)[4]
colnames(border_mat) = colnames(m3b)

library(circlize)
png(
  "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/6. Long COVID indicators/symmetricmyplot.png",
  width     = 6,
  height    = 4,
  units     = "in",
  res       = 1200,
  pointsize = 5
)
circos.par(gap.degree = 8)
chordDiagram(m3b, grid.col = 1:5, symmetric = TRUE, directional = TRUE, annotationTrack = "grid",
             preAllocateTracks = list(list(track.height = 0.05),
                                      list(track.height = 0.05)),
             link.lwd = 2, link.lty = 2, link.border = border_mat)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
  dd <- ifelse(theta < 90 || theta > 270, "downward", "downward")
  aa = c(1, 0.5)
  if(theta < 90 || theta > 270)  aa = c(0, 0.5)
  circos.text(x=mean(xlim), y=1, labels=sector.index, facing = dd, cex=0.8,  adj = aa)
  #circos.text(mean(xlim), mean(ylim), sector.index, facing = "downward", cex = 1, niceFacing = TRUE)
}, bg.border = NA)
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  circos.axis("bottom", major.tick.length = 0.2, labels.cex = 0.4)
}, bg.border = NA)
circos.clear()
dev.off()



