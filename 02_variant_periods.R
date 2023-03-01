##############################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk
## Description: 02_variant_periods - Use whole genome sequencing data to identify periods when different 
## variants were dominant (>60% of cases) - for use in sub-sample analysis by variant period
## Info on variants: https://www.who.int/activities/tracking-SARS-CoV-2-variants
##############################################################################################################

# 0. Set-up ------

# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(lubridate)

# Read in data (from University of Edinburgh, Centre of Genomics)
wgs_raw <- readRDS("/conf/EAVE/GPanalysis/data/WGS_latest.rds")


# 1. Get proportion of variants each week ----
## Keep one entry per individual per month (so that multiple tests per case don't skew proportions - allow a month for reinfection) + remove data beyond 20th Oct 22
wgs <- wgs_raw %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>%
  dplyr::rename(specimen_date = Collection_Date) %>%
  filter(specimen_date <= "2022-10-20") %>% 
  mutate(month = format(as.Date(specimen_date), "%Y-%m")) %>% 
  filter(!duplicated(EAVE_LINKNO, month)) %>% 
  select(EAVE_LINKNO, VariantShorthand, specimen_date)

## Check proportion for each variant
sort(round(prop.table(table(wgs$VariantShorthand))*100,1), decreasing = T)

## Identify main variants (>0.5% of cases)
wgs <- wgs %>% 
  mutate(variant = ifelse(VariantShorthand %in% c('BA.2', # omicron = 29.6% (at 20/10/22)
                                                  'B.1.617.2', # delta = 26.9%
                                                  'BA.1', # omicron = 20.9%
                                                  'VOC1', # alpha = 7.9% 
                                                  'BA.5', # omicron = 5.1% 
                                                  'AY.4.2', # delta = 3.8% 
                                                  'B.1.1.529', # omicron = 0.7% 
                                                  'BA.4'), # omicron = 0.7%
                          VariantShorthand, 
                          ifelse(VariantShorthand=="NA", "Unknown", "Other"))) 


## Estimate dominant variant by week (>=60% of cases that week) (use to impute missing variant)
wgs <- wgs %>%
  mutate(year = year(specimen_date),
         week = week(specimen_date),
         year_week = paste(year, week, sep="-")) %>% 
  # Instances of variant by week
  add_count(year_week, variant, name = "var_total") %>% 
  # Results by week
  add_count(year_week, name = "week_total") %>% 
  # % variant by week
  mutate(var_week_pc = var_total/week_total,
         # If variant is >60% of cases in a week, make it the dominant variant for that week
         dominant_var = ifelse(var_week_pc >= 0.6, variant, "Unknown")) %>% 
  select(EAVE_LINKNO, specimen_date, year_week, week_total, variant, dominant_var)

## Identify peak periods for different variants
peaks <- wgs %>% 
  select(year_week, week_total, dominant_var) %>% 
  filter(!duplicated(year_week)) %>% 
  arrange(year_week)
## NB not all tests had wgs performed --> 'week total' isn't the total number of cases - just the total number that were sequenced

## Export
write_csv(peaks, "/conf/EAVE/GPanalysis/analyses/long_covid/outputs/1. Misc/Dominant variant by week.csv")

rm(peaks, wgs, wgs_raw)
