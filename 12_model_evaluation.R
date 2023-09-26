##########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 12_model_evaluation - Gets AUROC and precision-recall curve. Identifies optimal discrimination threshold 
## and estimates stats for model evaluation. Prepares calibration plots. Plots predicted probabiites for sub-samples.
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
library(pmsampsize) # power analysis
library(DescTools) # to get Cox Snell Rsquared (for use in overfitting estimations)
library(stringr)
library(reshape2) # for melt function to reshape correlation matrix
library(caret) # for generating confusion matrix and k-fold cross-validation
library(yardstick) # for precision recall curves
library(glmnet) # for LASSO (use version 4.1-2 as C++17 not available)
library(e1071) # for Naive Bayes classifier
library(questionr) # for weighted tables
library(MLmetrics) # for F1_Score
library(PRROC)  # for pr_auc function
library(boot)   # for bootstrapping

## Set WD
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")

# 1. Set data, dependent variable, and predictors to use ----

# Set dataset to use
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/df_training_cleaned.rds") # Training
#df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/df_testing_with_probs.rds") # Holdout

# Set predicted probabilities to use
#df$pred <- df$AIC_pred
df$pred <- df$LASSO_pred  # opdef or code, txt, sick note (main analysis - all others are sensitivity tests)
#df$pred <- df$GBDT_pred
#df$pred <- df$NB_pred
#df$pred <- df$LASSO_pred_opdef # opdef 
#df$pred <- df$LASSO_pred_opdef_nb_codes # opdef no bloods or code, txt, sick note
#df$pred <- df$LASSO_pred_codes # code, txt, sick note
#df$pred <- df$LASSO_pred_opdef_nb # opdef no bloods

# Set model name (for labels)
#mod_name <- "AIC - opt" # e.g. "AIC opt", "LASSO opt", "Naive Bayes", "LASSO - no blood code txt sick", "LASSO code txt sick", "LASSO no blood"
mod_name <- "LASSO - opt" # (main analysis - all others are sensitivity tests)
#mod_name <- "GBDT"
#mod_name <- "Naive Bayes"
#mod_name <- "LASSO - op def"
#mod_name <- "LASSO - no bloods code txt sick"
#mod_name <- "LASSO - code txt sick"
#mod_name <- "LASSO - no bloods"

# Set dependent variable (to match predicted probabilities being used)
df$depvar <- ifelse(df$opdef == 1 | df$code_txt_fitnote == 1, 1, 0) # for LASSO_pred_opdef (6.2%)
#df$depvar <- df$opdef # for LASSO_pred_opdef (4.7%)
#df$depvar <- ifelse(df$opdef_no_blood == 1 | df$code_txt_fitnote == 1, 1, 0) # for LASSO_pred_opdef_nb_codes (2.8%)
#df$depvar <- df$code_txt_fitnote # for LASSO_pred_codes (2.1%)
#df$depvar <- df$opdef_no_blood # for LASSO_pred_opdef_nb (0.9%)


# 2. AUROC and c-stat (concordance) ----
auroc <- pROC::roc(df$depvar, df$pred, direction = "<", weights = df$class_weight)
cstat <- auroc$auc
cstat # AIC = 0.7135; LASSO = 0.7135; GBDT = 0.726; NB = 0.690
# LASSO_pred_opdef (opdef) = 0.707
# LASSO_pred_opdef_nb_codes (opdef no bloods, or LC indicators) = 0.718
# LASSO_pred_codes (LC indicators, no opdef) = 0.756
# LASSO_pred_opdef_nb (opdef no bloods) = 0.725

# Calculate the confidence interval using bootstrap resampling (2,000 samples)
ci <- pROC::ci(auroc, method = "bootstrap", of = "auc", level = 0.95)

# Extract the lower and upper bounds of the confidence interval
lower_bound <- ci[1]
lower_bound

upper_bound <- ci[3]
upper_bound

AUROC_plot <- ggroc(auroc) +
  theme_minimal() + 
  ggtitle("Receiver Operator Curve") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") 

AUROC_plot

ggsave(paste0("AUROC - ", mod_name, ".jpg"), AUROC_plot)


# 3. Precision-recall curve ----
# https://yardstick.tidymodels.org/reference/pr_curve.html
# Plots the share of cases identified as true positives (sensitivity) (y-axis) against the  
# precision rate (true positives as a share of all predicted positives i.e. the proportion of positive predictions that were correct). 
# Helpful because it does not include true negatives in its calculation --> is not affected by imbalance in the share of true positives and true negatives. 

## Area under the precision-recall curve
df_pr <- data.frame(obs = as.factor(df$depvar), pred = as.numeric(df$pred))
au_prc <- pr_auc(df_pr, obs, pred, event_level = "second")
au_prc # AIC = 0.141, LASSO = 0.141, GBDT = 0.155; NB = 0.152
# LASSO_pred_opdef (opdef) = 0.107
# LASSO_pred_opdef_nb_codes (opdef no bloods, or LC indicators) = 0.0673
# LASSO_pred_codes (LC indicators, no opdef) = 0.081
# LASSO_pred_opdef_nb (opdef no bloods) = 0.028

precision_recall_plot <- autoplot(pr_curve(df_pr, obs, pred, event_level = "second"))
precision_recall_plot
ggsave(paste0("Precision recall - ", mod_name, ".jpg"), precision_recall_plot)


# 4. Evaluation stats ----
## Set discrimination threshold
dt <- length(df$depvar[df$depvar == 1]) / length(df$depvar)

## Get predicted and actual values
predicted <- ifelse(df$pred >= dt, 1, 0)
actual <- df$depvar

## Define functions to get stats + 95% CI
## Sensitivity
calculateSensitivity <- function(predicted, actual) {
  true_positives <- sum(predicted[actual == 1] == 1)
  all_positives <- sum(actual == 1)
  if (all_positives == 0) return(list(Sensitivity = NA, CI = NA))
  sensitivity <- true_positives / all_positives
  
  # Calculate the 95% confidence interval 
  ci <- binom.test(true_positives, all_positives, conf.level = 0.95)
  
  return(data.frame(stat = sensitivity, lower = ci$conf.int[1], upper = ci$conf.int[2]))
}

# Specificity
calculateSpecificity <- function(predicted, actual) {
  true_negatives <- sum(predicted[actual == 0] == 0)
  all_negatives <- sum(actual == 0)
  if (all_negatives == 0) return(NA)
  specificity <- true_negatives / all_negatives
  
  # Calculate the 95% confidence interval
  ci <- binom.test(true_negatives, all_negatives, conf.level = 0.95)
  
  return(data.frame(stat = specificity, lower = ci$conf.int[1], upper = ci$conf.int[2]))
}

# Accuracy
calculateAccuracy <- function(predicted, actual){
  correct_pos_predictions <- length(predicted[predicted == 1 & actual == 1])
  correct_neg_predictions <- length(predicted[predicted == 0 & actual == 0])
  correct_predictions <- correct_pos_predictions + correct_neg_predictions
  accuracy <- correct_predictions/length(predicted)
  
  # Calculate the 95% confidence interval
  ci <- binom.test(correct_predictions, length(predicted), conf.level = 0.95)
  
  return(data.frame(stat = accuracy, lower = ci$conf.int[1], upper = ci$conf.int[2]))
  
}

# Positive predicted values 
calculatePPV<- function(predicted, actual){
  correct_pos_predictions <- length(predicted[predicted == 1 & actual == 1])
  all_pos_predictions <- sum(predicted[predicted == 1])
  PPV <- correct_pos_predictions/all_pos_predictions
  
  # Calculate the 95% confidence interval
  ci <- binom.test(correct_pos_predictions, all_pos_predictions, conf.level = 0.95)
  
  return(data.frame(stat = PPV, lower = ci$conf.int[1], upper = ci$conf.int[2]))
  
}

# Negative predicted values
calculateNPV<- function(predicted, actual){
  correct_neg_predictions <- length(predicted[predicted == 0 & actual == 0])
  all_neg_predictions <- length(predicted[predicted == 0])
  NPV <- correct_neg_predictions/all_neg_predictions
  
  # Calculate the 95% confidence interval
  ci <- binom.test(correct_neg_predictions, all_neg_predictions, conf.level = 0.95)
  
  return(data.frame(stat = NPV, lower = ci$conf.int[1], upper = ci$conf.int[2]))
  
}

# F1 Score
calculateF1 <- function(predicted, actual) {
  precision <- calculateSensitivity(predicted, actual)[1]
  recall <- calculateSensitivity(actual, predicted)[1]
  if (precision + recall == 0) return(NA)
  F1 = (2 * (precision * recall) / (precision + recall))
  
  # Calculate the 95% confidence interval 
  n <- length(predicted)
  z <- qnorm(0.975)  # Z-score for 95% confidence interval
  se <- sqrt((F1 * (1 - F1)) / n)
  margin_of_error <- z * se
  lower <- F1 - margin_of_error
  upper <- F1 + margin_of_error
  
  F1stats <- data.frame(stat = F1, lower = lower, upper = upper)
  
  colnames(F1stats) <- c("stat", "lower", "upper")
  
  return(F1stats)
}

# Brier Score
calculateBrierScore <- function(predicted, actual) {
  Brier <- mean((as.numeric(as.character(actual)) - as.numeric(predicted))^2)
  
  # Calculate the 95% confidence interval 
  n <- length(predicted)
  num_iterations <- 1000  # 1,000 is recommended
  
  # Initialize an empty vector to store bootstrapped Brier Scores
  bootstrapped_scores <- numeric(num_iterations)
  
  for (i in 1:num_iterations) {
    # Sample with replacement from the observed data
    sampled_indices <- sample(1:n, replace = TRUE)
    sampled_actual <- actual[sampled_indices]
    sampled_predicted <- predicted[sampled_indices]
    
    # Calculate Brier Score for the bootstrapped sample
    bootstrapped_scores[i] <- mean((as.numeric(as.character(sampled_actual)) - as.numeric(sampled_predicted))^2)
  }
  
  # Calculate the lower and upper bounds of the 95% confidence interval
  lower <- quantile(bootstrapped_scores, 0.025)
  upper <- quantile(bootstrapped_scores, 0.975)
  
  
  return(data.frame(stat = Brier, lower = lower, upper = upper))
}


## Return stats and CI
sensitivity <- data.frame(stat_name = "sensitivity", round(calculateSensitivity(predicted, actual), 2))
specificity <- data.frame(stat_name = "specificity", round(calculateSpecificity(predicted, actual),2))
accuracy <- data.frame(stat_name = "accuracy", round(calculateAccuracy(predicted, actual),2))
ppv <- data.frame(stat_name = "PPV", round(calculatePPV(predicted, actual),2))
npv <- data.frame(stat_name = "NPV", round(calculateNPV(predicted, actual),2))
F1 <- data.frame(stat_name = "F1 Score", round(calculateF1(predicted, actual),2))
brier <- data.frame(stat_name = "Brier", round(calculateBrierScore(predicted, actual),2))

eval_stats <- rbind(sensitivity, specificity, accuracy, ppv, npv, F1, brier)
eval_stats


# 6. Confusion matrices ----
## Function to generate confusion matrix and evaluation stats with 95% CI at a given discrimination threshold
conf_mat_fun <- function(dt){ # dt = discrimination threshold e.g. 0.062
  
  ### Create binary predictions
  predicted <- as.factor(ifelse(df$pred >= dt, 1, 0))
  
  ### Make depvar a factor
  df$depvar <- as.factor(df$depvar)
  
  ### Confusion matrix
  conf_mat <- confusionMatrix(predicted, df$depvar) 
  
  ### Get stats at dt
  predicted <- as.numeric(ifelse(df$pred >= dt, 1, 0))
  stats <- c(Sensitivity = sum(predicted[df$depvar == 1] == 1) / sum(df$depvar == 1), # true pos predictions/observed pos
             Specificity = sum(predicted[df$depvar == 0] == 0) / sum(df$depvar == 0), # true neg predictions/observed neg
             Accuracy    = sum(predicted == df$depvar) / length(df$depvar), # correct predictions
             `Positive predicted values`= ifelse(is.nan(sum(predicted[df$depvar == 1] == 1) / sum(predicted)), NA,
                                                 sum(predicted[df$depvar == 1] == 1) / sum(predicted)), # true pos predictions/all pos predictions
             `Negative predicted values` = ifelse(dt == 0, 0,
                                                  as.numeric(table(predicted[df$depvar == 0] == 0)["TRUE"] / table(predicted[predicted == 0]))),# true neg predictions/all neg predictions
             F1= 2/((1/(sum(predicted[df$depvar == 1] == 1) / sum(df$depvar == 1))) + 
                      (1/(sum(predicted[df$depvar == 1] == 1) / sum(predicted)))), # weighted avg of precision and sensitivity
             # Matthew's Correlation Coefficient
             # https://towardsdatascience.com/the-best-classification-metric-youve-never-heard-of-the-matthews-correlation-coefficient-3bf50a2f3e9a
             `Matthew's Correlation Coefficient` = ifelse(dim(table(predicted)) > 1, # Can only be estimated if 1s and 0s are predicted
                                                          mcc_vec(truth = as.factor(df$depvar), estimate = as.factor(predicted)), NA),
             `Brier Score` = mean((as.numeric(as.character(df$depvar)) - as.numeric(predicted))^2)
  )
  
  return(list(conf_mat = conf_mat, stats = round(stats,2)))
}
## Threshold == prevalence
conf_mat_prev <- conf_mat_fun(prop.table(table(df$depvar))[2]) # 0.062
conf_mat_prev$conf_mat 
conf_mat_prev$stats


# 7. Calibration plot and slope - expected vs observed ----
## Calculate the mean predicted probabilities and the observed frequencies of positive outcomes (testing data)

## Check distribution of predicted probabilities
hist_pp <- hist(df$pred) # 4x4 Histogram of predicted probabilities - AIC opt

## Plot observations by deciles of predicted probabilities e.g. first 10% of predicted probabilities, next 10%, ..., last 10%
calib_data <- df %>% 
  mutate(depvar = as.numeric(as.character(depvar))) %>% # change from factor to numeric
  group_by(decile = ntile(pred, 10)) %>% # group by deciles of predicted probabilities
  summarize(pred = mean(pred), # get mean predicted probability in each decile
            obs_pred = mean(depvar)) # get the rate of positive observations in each decile (mean = rate, given binary variable)

## Calibration plot 
calib_plot <- ggplot(calib_data, aes(x = pred, y = obs_pred)) +
  # Add smooth line between points
  geom_smooth(method = "loess", # Loess smoothing is a non-parametric form of regression that uses a 
              # weighted, sliding-window, average to calculate a line of best fit. 
              # Within each "window", a weighted average is calculated, and the sliding 
              # window passes along the x-axis.
              color = "grey", se = F) + # se = logical - whether or not to show 95% confidence interval
  # Add points over the top
  geom_point() +
  # Add 45 degree 'ideal' line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  # Label axes
  labs(x = "Predicted Probability", y = "Observed Frequency of Positive Outcomes") +
  theme_minimal() 

# Plot notes: Calibration plot of actual outcomes versus predictions for the model. N = xxx (testing data).
# Points indicate the observed frequencies by deciles of predicted probabilitys. The Loess smoother is shown in gray. 
# The dahsed line represents the ideal 45° line.

# Display the calibration plot
print(calib_plot)

# Export
ggsave(paste0("Calibration plot (deciles) - ", mod_name, ".jpg"), calib_plot, height = 4, width = 4, units = "in")

# Calculate the slope & intercept of the calibration line NB using individual points, rather than deciles (used in plot)
df$depvar <- as.numeric(as.character(df$depvar))

calib_mod <- summary(lm(depvar ~ pred, data = df))

calib_slope <- coef(calib_mod)[2,1]
calib_slope # 0.93

calib_intercept <- coef(calib_mod)[1,1]
calib_intercept # 0.00

# Calculate the confidence interval of calib_slope
ci <- confint(lm(depvar ~ pred, data = df), "pred", level = 0.95)

# Extract the lower and upper bounds of the confidence interval
lower_bound <- ci[1]
lower_bound # 0.954

upper_bound <- ci[2]
upper_bound # 1.140

# Print the calibration slope
cat("Calibration slope:", calib_slope, "\n") # AIC = 1.047 / LASSO = 1.047 / NB = 0.097 (??) / GBDT = 0.991
# LASSO_pred_opdef (opdef) = 0.78
# LASSO_pred_opdef_nb_codes (opdef no bloods, or LC indicators) = 0.51
# LASSO_pred_codes (LC indicators, no opdef) = 1.03
# LASSO_pred_opdef_nb (opdef no bloods) = 1.03



# 8. Observed and predicted probabilities, by vigentile of predicted probabilities ----
vigentile_plot <-
  function(mask, plot_name) {
    # e.g. mask = !is.na(df$EAVE_LINKNO)
    #      plot_name = "full sample"
    
    ## Get mean depvar and mean pred at each vigentile of pred
    df_20 <- df[mask,] %>%
      dplyr::select(pred, depvar) %>%
      arrange(pred) %>%
      mutate(vigentile = ntile(pred, 20)) %>%
      group_by(vigentile) %>%
      summarise(mean_depvar = mean(as.numeric(as.character(depvar))),
                mean_pred = mean(pred))
    
    ## Pivot long
    df_20 <- pivot_longer(
      df_20,
      cols = mean_depvar:mean_pred,
      names_to = "Category",
      values_to = "Probability"
    )
    
    plot <- ggplot(data = df_20,
                   aes(x = vigentile, y = Probability, colour = Category)) +
      geom_point() + labs(y = "Probability", x = "Vigentile Predicted Probability", colour = "") +
      theme_minimal() +
      scale_color_manual(
        values = c("mean_depvar" = "#F8766D", "mean_pred" = "#00BFC4"),
        labels = c("Observed", "Predicted")
      )
    
    setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/7. Prediction modelling")
    ggsave(
      paste0("Vigentile probabilities - ", plot_name, ".jpg"),
      plot,
      height = 4, width = 6, units = "in"
    )
    print(plot)
  }

## Full sample
vigentile_all <- vigentile_plot(mask = !is.na(df$EAVE_LINKNO), plot_name = "Full sample")

## Under 50s
vigentile_u50 <- vigentile_plot(mask = df$age < 50, plot_name = "Under 50s")

## Over 50s
vigentile_o50 <- vigentile_plot(mask = df$age >= 50, plot_name = "Over 50s")

## Wild
vigentile_wild <- vigentile_plot(mask = df$var_per == "Wild", plot_name = "Wild")

## Alpha
vigentile_alpha <- vigentile_plot(mask = df$var_per == "Alpha", plot_name = "Alpha")

## Delta
vigentile_delta <- vigentile_plot(mask = df$var_per == "Delta", plot_name = "Delta")

## Omicron
vigentile_omicron <- vigentile_plot(mask = df$var_per == "Omicron", plot_name = "Omicron")



# 9. Predicted probabilities (testing data) ----
## 9.1 By age and sex ----
### Summarise data by age and sex
male_data <- subset(df, sex_m == 1) %>% 
  group_by(age) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred)) 

female_data <- subset(df, sex_m == 0) %>% 
  group_by(age) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

### Create line chart
age_sex_plot <- ggplot() +
  # Add male line with CI
  geom_line(data = male_data, aes(x = age, y = pred, color = "Male")) +
  geom_ribbon(data = male_data, aes(x = age, ymin = ci_low, ymax = ci_high), fill = "#00BFC4", alpha = 0.2) +
  # Add female line with CI
  geom_line(data = female_data, aes(x = age, y = pred, color = "Female")) +
  geom_ribbon(data = female_data, aes(x = age, ymin = ci_low, ymax = ci_high), fill = "#F8766D", alpha = 0.2) +
  # Set axis labels and chart title
  xlab("Age") +
  ylab("Predicted probability of long COVID") +
  # Set x-axis tick marks for every fifth tick
  scale_x_continuous(breaks = seq(0, max(df$age), by = 5)) +
  # Remove legend title
  guides(color=guide_legend(title=NULL)) +
  theme_minimal()

age_sex_plot

ggsave(paste0("Prob by age and sex - ", mod_name,".jpg"), 
       age_sex_plot,
       width = 7, height = 7, units = "in")


## 9.2 By age and variant period ----
### Summarise data by age and variant period
wild_data <- subset(df, var_per == "Wild") %>% 
  group_by(age) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred)) 

alpha_data <- subset(df, var_per == "Alpha") %>% 
  group_by(age) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

delta_data <- subset(df, var_per == "Delta") %>% 
  group_by(age) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

omicron_data <- subset(df, var_per == "Omicron") %>% 
  group_by(age) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

### Create line chart
age_variant_plot <- ggplot() +
  # Add Wild line with CI
  geom_line(data = wild_data, aes(x = age, y = pred, color = "Wild")) +
  geom_ribbon(data = wild_data, aes(x = age, ymin = ci_low, ymax = ci_high), fill = "#C77CFF", alpha = 0.2) +
  # Add Alpha line with CI
  geom_line(data = alpha_data, aes(x = age, y = pred, color = "Alpha")) +
  geom_ribbon(data = alpha_data, aes(x = age, ymin = ci_low, ymax = ci_high), fill = "#F8766D", alpha = 0.2) +
  # Add Delta line with CI
  geom_line(data = delta_data, aes(x = age, y = pred, color = "Delta")) +
  geom_ribbon(data = delta_data, aes(x = age, ymin = ci_low, ymax = ci_high), fill = "#7CAE00", alpha = 0.2) +
  # Add Omicron line with CI
  geom_line(data = omicron_data, aes(x = age, y = pred, color = "Omicron")) +
  geom_ribbon(data = omicron_data, aes(x = age, ymin = ci_low, ymax = ci_high), fill = "#00BFC4", alpha = 0.2) +
  # Set axis labels and chart title
  xlab("Age") +
  ylab("Predicted probability of long COVID") +
  # Set x-axis tick marks for every fifth tick
  scale_x_continuous(breaks = seq(0, max(df$age), by = 5)) +
  # Remove legend title
  guides(color=guide_legend(title=NULL)) +
  theme_minimal()

age_variant_plot

ggsave(paste0("Prob by age and variant - ", mod_name,".jpg"), 
       age_variant_plot,
       width = 7, height = 7, units = "in")


## 9.3 By BMI and sex ----
### Summarise data by age and sex
male_data <- subset(df, sex_m == 1) %>% 
  group_by(bmi_imp) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred)) 

female_data <- subset(df, sex_m == 0) %>% 
  group_by(bmi_imp) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

### Create line chart
bmi_sex_plot <- ggplot() +
  # Add male line with CI
  geom_line(data = male_data, aes(x = bmi_imp, y = pred, color = "Male")) +
  geom_ribbon(data = male_data, aes(x = bmi_imp, ymin = ci_low, ymax = ci_high), fill = "#00BFC4", alpha = 0.2) +
  # Add female line with CI
  geom_line(data = female_data, aes(x = bmi_imp, y = pred, color = "Female")) +
  geom_ribbon(data = female_data, aes(x = bmi_imp, ymin = ci_low, ymax = ci_high), fill = "#F8766D", alpha = 0.2) +
  # Set axis labels and chart title
  xlab("Body Mass Index") +
  ylab("Predicted probability of long COVID") +
  # Set x-axis tick marks for every fifth tick
  scale_x_continuous(breaks = seq(0, max(df$bmi_imp), by = 5)) +
  # Remove legend title
  guides(color=guide_legend(title=NULL)) +
  theme_minimal()

bmi_sex_plot

ggsave(paste0("Prob by BMI and sex - ", mod_name,".jpg"), 
       bmi_sex_plot,
       width = 7, height = 7, units = "in")


## 9.4 By BMI and variant period ----
### Summarise data by age and variant period
wild_data <- subset(df, var_per == "Wild") %>% 
  group_by(bmi_imp) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred)) 

alpha_data <- subset(df, var_per == "Alpha") %>% 
  group_by(bmi_imp) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

delta_data <- subset(df, var_per == "Delta") %>% 
  group_by(bmi_imp) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

omicron_data <- subset(df, var_per == "Omicron") %>% 
  group_by(bmi_imp) %>% 
  summarise(ci_low = mean(pred) - 1.96*sd(pred)/sqrt(n()), # add lower CI
            ci_high = mean(pred) + 1.96*sd(pred)/sqrt(n()), # add upper CI
            pred = mean(pred))

### Create line chart
bmi_variant_plot <- ggplot() +
  # Add Wild line with CI
  geom_line(data = wild_data, aes(x = bmi_imp, y = pred, color = "Wild")) +
  geom_ribbon(data = wild_data, aes(x = bmi_imp, ymin = ci_low, ymax = ci_high), fill = "#C77CFF", alpha = 0.2) +
  # Add Alpha line with CI
  geom_line(data = alpha_data, aes(x = bmi_imp, y = pred, color = "Alpha")) +
  geom_ribbon(data = alpha_data, aes(x = bmi_imp, ymin = ci_low, ymax = ci_high), fill = "#F8766D", alpha = 0.2) +
  # Add Delta line with CI
  geom_line(data = delta_data, aes(x = bmi_imp, y = pred, color = "Delta")) +
  geom_ribbon(data = delta_data, aes(x = bmi_imp, ymin = ci_low, ymax = ci_high), fill = "#7CAE00", alpha = 0.2) +
  # Add Omicron line with CI
  geom_line(data = omicron_data, aes(x = bmi_imp, y = pred, color = "Omicron")) +
  geom_ribbon(data = omicron_data, aes(x = bmi_imp, ymin = ci_low, ymax = ci_high), fill = "#00BFC4", alpha = 0.2) +
  # Set axis labels and chart title
  xlab("Body Mass Index") +
  ylab("Predicted probability of long COVID") +
  # Set x-axis tick marks for every fifth tick
  scale_x_continuous(breaks = seq(0, max(df$bmi_imp), by = 5)) +
  # Remove legend title
  guides(color=guide_legend(title=NULL)) +
  theme_minimal()

bmi_variant_plot

ggsave(paste0("Prob by bmi and variant - ", mod_name,".jpg"), 
       bmi_variant_plot,
       width = 7, height = 7, units = "in")

