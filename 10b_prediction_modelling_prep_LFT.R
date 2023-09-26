#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 10b_prediction_modelling_prep_LFT - Prepares variables for prediction modelling. Gets data for all individuals
## with a positive PCR or LFT test AND identified as having long COVID by the operational definition, clinical codes,
## free text, or long covid sick note. Incorporates prescriptions data. Formats variables.
#########################################################################################################################

# Set-up ------
# Clear environment 
rm(list=ls())
setwd("/conf/EAVE/GPanalysis/analyses/long_covid")

# Libraries
library(tidyverse)
library(lubridate)
library(stringr)
library(ggplot2)
library(splines)
library(MatchIt)

## Read in data: covariates + long COVID classification for all individuals with a positive PCR or LFT result
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_predict.rds") %>% 
  filter(tes_res == "POSITIVE" | `Long-COVID code, free text or sick note` > 0) # 1,462,676


# 1. Update tes_dat to equal gp_dat for those never tested
## For those who have never tested (but have a long covid code, free text or LC sick note) make tes_dat = the date their 
## long covid code, free text or LC sick note was recorded. This will be used to count tests taken by that date & number of vaccinations

## Read in dates long covid code, free text or LC sick note was recorded (preapred in 09)
## Don't keep operational definition, as that is already recorded in df
dates <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/6. Long COVID indicators/Outcome measure dates.rds") %>% 
  filter(group != "Operational definition" & group != "Operational definition - no blood tests") %>% 
  
  ## Keep earliest recorded indicator of long COVID
  arrange(gp_dat) %>% 
  distinct(EAVE_LINKNO,.keep_all = T) %>% 
  select(-group)

# Update tes_dat if it's missing (i.e. no tests take) OR if long covid code, free text or LC sick note is earlier than tes_dat
df <- df %>% 
  left_join(dates, by = "EAVE_LINKNO") %>% 
  mutate(tes_dat = ifelse(is.na(tes_dat), gp_dat, tes_dat),
         tes_dat = ifelse(!is.na(gp_dat) & gp_dat < tes_dat, gp_dat, tes_dat), 
         tes_dat = as.Date(tes_dat, origin = "1970-01-01")) %>% 
  select(-gp_dat)


# 2. Get number of vaccine doses by test date (PCR test date) or by validation measure date ----
## Read in vaccinations data for individuals in df
vac_raw <- readRDS("/conf/EAVE/GPanalysis/data/temp/vaccine_cleaned.rds") %>% # cleaned file prepared by Vera @ PHS
  filter(EAVE_LINKNO %in% df$EAVE_LINKNO)

## Prepare vac data
vac <- vac_raw %>% 
  ## Select vars of interest
  dplyr::select(c(EAVE_LINKNO, starts_with("date_vacc")))

## Pivot longer
vac <- pivot_longer(vac,
                    cols = date_vacc_1:date_vacc_sb,
                    names_to = "vac_dos",
                    values_to = "vac_dat")


## Keep last vaccination dose up to 14 days before test date/date LC code, free text, sick note recorded
vac <- vac %>% 
  filter(!is.na(vac_dat)) %>% 
  left_join(dplyr::select(df, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
  filter(vac_dat < date(tes_dat - days(14))) %>% 
  mutate(vac_elaps = as.numeric(tes_dat - vac_dat)) %>% 
  arrange(desc(vac_dat)) %>% 
  group_by(EAVE_LINKNO) %>% 
  summarise(vac_dos = n(), vac_elaps = first(vac_elaps)) %>% 
  dplyr::select(EAVE_LINKNO, vac_dos, vac_elaps) %>% 
  rename(vac_dos_clean = vac_dos)

## Join to df (keep unknowns - these are people cleaned out, as opposed to people not vaccinated)
df <- df %>% 
  left_join(vac, by = "EAVE_LINKNO") %>% 
  mutate(vac_dos = ifelse(vac_dos == "Unknown" & is.na(vac_dos_clean), "Unknown", vac_dos_clean),
         vac_dos = ifelse(vac_dos != "Unknown", vac_dos_clean, "Unknown")) %>% 
  dplyr::select(-vac_dos_clean) %>% 
  mutate(vac_dos = ifelse(is.na(vac_dos), 0, vac_dos),
         vac_dos = ifelse(vac_dos == "Unknown", NA, vac_dos)) 

## Remove files
rm(vac_raw, vac)


# 3. Format predictor vars for modelling and create binary variables for correlations ----
df <- df %>% 
  mutate(
    # Age: There are very few age 100+, combine with the 97-99 category to make 97+
    age = ifelse(age > 97, 97, age),
    # Sex
    sex = as.factor(sex),
    sex_m = ifelse(sex == "M", 1, 0),
    # SIMD
    simd = case_when(simd == "1 - High" ~ 1,
                     simd == "2" ~ 2,
                     simd == "3" ~ 3,
                     simd == "4" ~ 4,
                     simd == "5-Low" ~ 5),
    simd = ifelse(is.na(simd), "Unknown", simd),
    simd = as.factor(simd),
    simd1 = ifelse(simd == 1, 1, 0),
    simd2 = ifelse(simd == 2, 1, 0),
    simd3 = ifelse(simd == 3, 1, 0),
    simd4 = ifelse(simd == 4, 1, 0),
    simd5 = ifelse(simd == 5, 1, 0),
    # Urban Rural classification
    ur = str_sub(ur6_lab, 1, 1),
    ur = ifelse(ur == "U", "Unknown", ur),
    ur = as.factor(ur),
    ur1 = ifelse(ur == 1, 1, 0),
    ur2 = ifelse(ur == 2, 1, 0),
    ur3 = ifelse(ur == 3, 1, 0),
    ur4 = ifelse(ur == 4, 1, 0),
    ur5 = ifelse(ur == 5, 1, 0),
    ur6 = ifelse(ur == 6, 1, 0),
    # Healthboards
    hb = as.factor(hb_lab),
    hb_ayr = ifelse(hb_lab == "NHS Ayrshire and Arran", 1, 0),
    hb_bor = ifelse(hb_lab == "NHS Borders", 1, 0),
    hb_dum = ifelse(hb_lab == "NHS Dumfries and Galloway", 1, 0),
    hb_fif = ifelse(hb_lab == "NHS Fife", 1, 0),
    hb_for = ifelse(hb_lab == "NHS Forth Valley", 1, 0),
    hb_gra = ifelse(hb_lab == "NHS Grampian", 1, 0),
    hb_gla = ifelse(hb_lab == "NHS Greater Glasgow and Clyde", 1, 0),
    hb_hig = ifelse(hb_lab == "NHS Highland", 1, 0),
    hb_lan = ifelse(hb_lab == "NHS Lanarkshire", 1, 0),
    hb_lot = ifelse(hb_lab == "NHS Lothian", 1, 0),
    hb_ork = ifelse(hb_lab == "NHS Orkney", 1, 0),
    hb_she = ifelse(hb_lab == "NHS Shetland", 1, 0),
    hb_tay = ifelse(hb_lab == "NHS Tayside", 1, 0),
    hb_wes = ifelse(hb_lab == "NHS Western Isles", 1, 0),
    # Household size
    hhold_bands = case_when(hhold_n == 1 ~ "1",
                            hhold_n == 2 ~ "2",
                            hhold_n > 2 & hhold_n < 6 ~ "3-5",
                            hhold_n > 5 & hhold_n < 11 ~ "6-10",
                            hhold_n > 10 ~ ">10"),
    hhold_bands = ordered(hhold_bands, levels = c("1", "2", "3-5", "6-10", ">10")),
    hhold1 = ifelse(hhold_bands == "1", 1, 0),
    hhold2 = ifelse(hhold_bands == "2", 1, 0),
    hhold3 = ifelse(hhold_bands == "3-5", 1, 0),
    hhold4 = ifelse(hhold_bands == "6-10", 1, 0),
    hhold5 = ifelse(hhold_bands == ">10", 1, 0),
    # Home type
    carehome = ifelse(hom_typ == 1, 1, 0),
    # Vaccinations
    vac_dos = ifelse(is.na(vac_dos), "Unknown", vac_dos),
    vac_dos = ordered(vac_dos, levels = c(0:6)),
    vac1 = ifelse(vac_dos == 1, 1, 0),
    vac2 = ifelse(vac_dos == 2, 1, 0),
    vac3 = ifelse(vac_dos == 3, 1, 0),
    vac4 = ifelse(vac_dos == 4, 1, 0),
    vac5 = ifelse(vac_dos == 5, 1, 0),
    vac6 = ifelse(vac_dos == 6, 1, 0),
    # Variant period
    var_per = ifelse(is.na(var_per), "Unknown", var_per),
    var_per = as.factor(var_per),
    wild = ifelse(var_per == "Wild", 1, 0),
    wild = ifelse(is.na(wild), 0, wild),
    alpha = ifelse(var_per == "Alpha", 1, 0),
    alpha = ifelse(is.na(alpha), 0, alpha),
    delta = ifelse(var_per == "Delta", 1, 0),
    delta = ifelse(is.na(delta), 0, delta),
    omicron = ifelse(var_per == "Omicron", 1, 0),
    omicron = ifelse(is.na(omicron), 0, omicron)) %>% 
  # Rename outcome measures (for ease of coding)
   rename("opdef" = `Operational definition`,
          "opdef_no_blood" = `Operational definition - no blood tests`,
          "code_txt_fitnote" = `Long-COVID code, free text or sick note`) %>% 
  # Make outcome measures binary
  mutate(opdef = ifelse(opdef > 0, 1, 0),
         opdef_no_blood = ifelse(opdef_no_blood == 2, 1, 0),
         code_txt_fitnote = ifelse(code_txt_fitnote > 0, 1, 0))


# 4. Tidy up df ----
df <- df %>% 
  dplyr::select(EAVE_LINKNO, old_EAVE_LINKNO, 
                # Outcome measure
                opdef, opdef_no_blood, code_txt_fitnote,
                # Socio-demographics
                age, age_10, sex, sex_m, colnames(df[startsWith(colnames(df), "simd")]), hhold_bands, hhold1, hhold2, hhold3, hhold4, hhold5,
                # Geographic
                colnames(df[startsWith(colnames(df), "ur")]), colnames(df[startsWith(colnames(df), "hb")]),
                # Testing
                tes_dat, tes_res, 
                # Variant
                var_per, wild, alpha, delta, omicron,
                # Vaccinations
                colnames(df[startsWith(colnames(df), "vac")]),
                # Risk factors
                shielding, immuno_supp, carehome, bmi_imp, colnames(df[startsWith(colnames(df), "Q_")]), n_risk_gps,
                # Severity of acute infection
                severe,
                # Weights
                eave_weight)


# Checks for duplicates 
cat("The file contains no duplicate EAVE_LINKNOs: T/F?", length(unique(df$EAVE_LINKNO)) == length(df$EAVE_LINKNO))

# Check for missing values
missing_vals <- colSums(is.na(df))
print("Check for missing values")
missing_vals[missing_vals>0]

## Mean impute missing bmi_imp (38,497 = 3.5%). These are individuals not in the demographics file (validated cohort) but 
## who have interacted with the healthcare system in recent years
df <- df %>% 
  mutate(bmi_imp = ifelse(is.na(bmi_imp), mean(bmi_imp, na.rm = T), bmi_imp))

# Check values for each variable
# Set model variables
mod_vars <- c("opdef", "opdef_no_blood", "code_txt_fitnote","age", "sex", "simd", "hhold_bands", "ur", "hb", "var_per", "vac_dos", "vac_elaps", "shielding", "immuno_supp", "carehome", 
              "bmi_imp", colnames(df[startsWith(colnames(df), "Q_")]), "severe")

for(var in mod_vars){
  print(var)
  print(table(df[[var]], exclude = NULL))
}


# 5. Integrate PIS data ----
## Set number of days before test date to consider
n_days <- 52/4*7 # days (7) in a quarter of a year (52/4) e.g. days in 3 months

## Read in PIS data supplied Dec 22
pis1 <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/PIS_07122022.rds")

## Keep data for EAVE_LINKNOs in df, for the 365 days before tes_dat
pis1 <- pis1 %>% 
  filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
  left_join(dplyr::select(df, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
  filter(`Disp Date` < tes_dat) %>% 
  filter(`Disp Date` >= (tes_dat - days(n_days))) %>% 
  # Format
  dplyr::select(EAVE_LINKNO, pis_dat = `Disp Date`, item = `PI BNF Item Code`)


## Read in PIS data supplied Mar 23 (additional prescriptions)
pis2 <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/Long_COVID_PIS_extract_2023-03-17.rds")

## Keep data for EAVE_LINKNOs in df, for the 365 days before tes_dat
pis2 <- pis2 %>% 
  filter(EAVE_LINKNO %in% df$EAVE_LINKNO) %>% 
  left_join(dplyr::select(df, EAVE_LINKNO, tes_dat), by = "EAVE_LINKNO") %>% 
  filter(Disp.Date < tes_dat) %>% 
  filter(Disp.Date >= (tes_dat - days(n_days))) %>% 
  # Format
  dplyr::select(EAVE_LINKNO, pis_dat = Disp.Date, item = PI.BNF.Item.Code)

## Merge pis1 and pis2 and remove duplicates
pis <- rbind(pis1, pis2) %>% 
  distinct(EAVE_LINKNO, pis_dat, item)
rm(pis1, pis2)

## Read in lookup file (for labels) and format to match pis
bnf <- read_csv("/conf/EAVE/GPanalysis/analyses/long_covid/data/BNF lookup.csv") %>% 
  dplyr::select(subpar = `BNF Subparagraph Code`, subpar_lab = `BNF Subparagraph`, chem = `BNF Chemical Substance Code`, chem_lab = `BNF Chemical Substance`, item = `BNF Presentation Code`) %>% 
  mutate(subpar = as.character(subpar))

## Get subparagraph and chemical substance codes from item number
pis <- pis %>% 
  mutate(subpar = substr(item, 2, 7),
         chem = substr(item, 1, 9))

## Add sub-paragraph labels
subpar_labs <- bnf %>% 
  dplyr::select(subpar, subpar_lab) %>%
  distinct(subpar, .keep_all = T) %>% 
  filter(subpar %in% pis$subpar) %>% 
  mutate(subpar = as.character(subpar))

pis <- pis %>% 
  mutate(subpar = as.character(subpar)) %>% 
  left_join(subpar_labs, by = "subpar")

## Add chemical substance labels
chem_labs <- bnf %>% 
  dplyr::select(chem, chem_lab) %>% 
  distinct(chem, .keep_all = T) %>% 
  filter(chem %in% pis$chem) %>% 
  mutate(chem = as.character(chem))

pis <- pis %>%  
  mutate(chem = as.character(chem)) %>% 
  left_join(chem_labs, by = "chem")

## Remove files
rm(bnf, chem_labs, subpar_labs)

## Set prescriptions of interest
subpars <- c("204000", "205040", "205051", "208010", "209000", "212000", "301040",  "303020", 
              "310000", "501011", "501050", "503021", "603010", "901011", "403030")
### Don't include those used in the operational definition: "301011", "302000", "309010", "309020", "501030" 

chems <- c("0304010D0", "0601022B0", "0109010U0", "0103010H0", "0103010T0", "1001040G0", 
           # Direct oral anticoagulatns
           "0208020Z0", "0208020X0", "0208020AA", "0208020Y0", 
           # Warfarin
           "0208020V0")

## Capture either sub-paragraph or chemical substance of interest in a single variable
pis <- pis %>% 
  filter(subpar %in% subpars | chem %in% chems)

pis <- pis %>% 
  mutate(pis = ifelse(subpar %in% subpars, subpar,
                      ifelse(chem %in% chems, chem, NA)),
         pis_lab = ifelse(subpar %in% subpars, subpar_lab,
                          ifelse(chem %in% chems, chem_lab, NA)))

## Create custom categories for direct oral anticoagulants and other oral anticoagulants
pis <- pis %>% 
  mutate(
         # Direct oral anticoagulants
         pis = ifelse(chem %in% c("0208020Z0", "0208020X0", "0208020AA", "0208020Y0"), "208020a", pis),
         pis_lab = ifelse(pis == "208020a", "Direct oral anticoagulants", pis_lab),
         # Warfarin
         pis = ifelse(chem %in% c("0208020V0"), "208020c", pis),
         pis_lab = ifelse(pis == "208020b", "Warfarin sodium", pis_lab))

## Get counts of prescriptions dispensed in the period before testing (for subparagraph OR chemical substance of interest, as captured by pis)
pis <- pis %>% 
  group_by(EAVE_LINKNO, pis_lab) %>% 
  summarise(n = n(), .groups = "drop") %>%
  mutate(n = as.numeric(n)) %>%
  ## Transpose
  pivot_wider(
    names_from = pis_lab,
    values_from = n,
    values_fill = list(n = 0))

# Add PIS prefix to colnames 
colnames(pis) <- c("EAVE_LINKNO", paste0("PIS_", colnames(pis)[2:ncol(pis)]))

# Get distribution of number of prescriptions
pis_vars <- colnames(pis)[2:ncol(pis)]

plot_data <- pivot_longer(pis, pis_vars, names_to = "Prescription", values_to = "N")

plot_data <- plot_data %>% 
  filter(N > 0) %>% 
  dplyr::select(-EAVE_LINKNO) %>% 
  filter(Prescription != "PIS_Coronavirus") %>% 
  filter(N < 30) %>% 
  rename("N prescriptions" = N) 
  
ggplot(plot_data, aes(x=`N prescriptions`))+
  geom_histogram() +
  facet_wrap(~Prescription)

# Get counts of individuals in receipt of each prescription
summary_data <- plot_data %>% 
  group_by(Prescription) %>% 
  summarise(Individuals = n()) 


# Join to df, fill NAs with 0s, make binary
df <- df %>% 
  left_join(pis, by = "EAVE_LINKNO") %>% 
  mutate_at(vars(`PIS_Selective serotonin re-uptake inhibitors`:`PIS_Ranitidine hydrochloride`), ~replace_na(.,0)) 

for(var in pis_vars){
  df[[var]] <- ifelse(df[[var]] > 0, 1, 0)
}

rm(pis)


# 6. Export data ----
## Check for missing values
missing_vals <- colSums(is.na(df))
print("Check for missing values")
missing_vals[missing_vals>0]

saveRDS(df, "/conf/EAVE/GPanalysis/analyses/long_covid/data/LongCOVID_predict_prepared_LFT.rds")
