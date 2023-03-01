#########################################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk 
## Description: 07_cluster_analysis - Performs pam (aka k-medoids) and hierarchical clustering on outcomes
## that matched analysis suggests are significantly higher following a positive PCR test. Estimates average
## silhouette width and Dunn index to assess clusters' internal validity. 
#########################################################################################################################

# 0. Set-up ------

# Clear environment 
rm(list=ls())

# Libraries
library(tidyverse)
library(ggplot2)
library(stats) # for k-means clustering and hclust and pca
library(cluster) # for k-medoids clustering (pam) using clara
library(ggdendro) # allows creation of ggplot dendrograms (which have the benefit of being savable as objects)
library(clValid) # for checking the validity of clusters e.g. silhouette, Dunn index
library(poLCA) # for latent class analysis (masks dplyr::select)
library(psych) # for PCA
library(corrplot) # for correlograms

# Set seed
set.seed(1234)

# Read in outcome counts (NB reads in positive and negative cases only - according to PCR + LFT) 
df_raw <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/5. Cluster analysis/df_cluster_LFT.rds") # prepared in 06_matched_analysis

## Select positive cases only
df <- df_raw[df_raw$tes_res == "POSITIVE",]


# 1. Identify dependent variables to analyse ----
# Function that reads in tables of results from matched analysis and keeps dependent variables that are significantly higher among positive cases
get_depvars <- function(dep_m, pis_m){ # dep_m: file containing matched analysis of health outcomes e.g. "df.csv"
                                       # pis_m: file containing matched analysis of prescriptions data e.g. "pis_df.csv"
  dep_m_df <- read_csv(dep_m) # Read codes and healthcare intereactions
  pis_m_df <- read_csv(pis_m) # Newly dispensed prescriptions
  
  depvar <- dep_m_df %>% 
    
  filter(Sig == "Significantly higher") %>% 
    filter(Outcome != "Long Covid") %>% # We don't want to cluster on this
    filter(Outcome != "Free text") %>% # We don't want to cluster on this
    filter(Outcome != "Free text - fitnote") %>% # We don't want to cluster on this
    filter(Outcome != "Covid") %>% # We don't want to cluster on this
    filter(Outcome != "Sick Note (except long-COVID)") %>% 
    dplyr::select(Outcome)
  
  pisvar <- pis_m_df %>% 
    filter(Sig == "Significantly higher") %>% 
    dplyr::select(Outcome)
  
  vars <- c(depvar$Outcome, pisvar$Outcome)

  return(vars)
  
}

## Get dependent variables that are significantly higher (or lower) among positive cases
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/3. Regression outputs/2. Matched analysis")

# Update var names to match var labels used in tables read in below
df <- df %>% 
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

## 4-12 weeks: Positive vs negative
neg_4 <- get_depvars(dep_m = "4-12 weeks - pos vs neg - matched sample - all.csv"
                     ,pis_m = "PIS 4-12 weeks - pos vs neg - matched sample - all.csv")

## >12-26 weeks: Positive vs negative
neg_12 <- get_depvars(dep_m = "12-26 weeks - pos vs neg - matched sample - all.csv"
                      ,pis_m = "PIS 12-26 weeks - pos vs neg - matched sample - all.csv") 

## Get variables that are significantly higher in both periods and export for use in 08_classification_and_validation
all_vars <- union(neg_4, neg_12)
write.csv(all_vars, "Significantly higher outcomes.csv")

## Get variables that are significantly higher at >12-26 weeks for sensitivity analysis in 08_classification_and_validation
write.csv(neg_12, "Significantly higher outcomes - 12-26.csv")


# 2. Prepare dataframes ----
## Function selects specified variables for cluster analysis for a given follow-up period and transposes df in preparation for clustering
prep_df <- function(period, # period = the period under investigation (4-12 or >12-26 weeks) e.g. "_4" or "_12"
                    control, # control = control group e.g. "neg" or "nt"
                    depvars) { 
  
  ## Save EAVE_LINKNO + keep variables for the relevant period
  EAVE_LINKNO <- df$EAVE_LINKNO
  df <- df[, grepl(period, names(df))]
  
  ## Clean up column names
  colnames(df) <- sub(period, "", colnames(df)) # remove "_4" or "_12"
  colnames(df) <- sub("_", " ", colnames(df)) # remove "_"
  
  ## Select dependent variables of interest (those that are significantly higher for positive cases in matched analysis)
  df <- df %>% dplyr::select(c(all_of(depvars)))
  
  ## Re-attach EAVE_LINKNO
  df <- cbind(EAVE_LINKNO, df)
  
  # Check counts of non-zero values for each variable
  print(sort(colSums(df[2:ncol(df)]!= 0)))
  
  # Remove rows that are all 0s (makes df smaller to work with)
  df <- df[rowSums(df[2:ncol(df)])!=0,]
  
  # Standardize variables (subtract mean, divide by standard deviation to give mean = 0, sd = 1)
  df <- as.data.frame(scale(df[2:ncol(df)]))
  
  # Binarize variables # (Alternative to standardizing, but standardizing gives us more information)
  df_bin <- as.data.frame(ifelse(df[,2:ncol(df)] > 0, 1, 0))
  
 
  # Transpose
  df <- t(df)
  df_bin <- t(df_bin)
  
  # Return dataframe
  return(list("df" = as.data.frame(df), "df_bin" = as.data.frame(df_bin)))
}


# WARNING: Don't try to view the dfs produced below - after transposing, they are >1M columns View(df) is likely to crash R.
## 4-12 weeks: Positive vs negative
df_neg_4 <- prep_df(period = "_4", control = "neg", depvars = neg_4)

## >12-26 weeks: Positive vs negative
df_neg_12 <- prep_df(period = "_12", control = "neg", depvars = neg_12)

## Remove files
rm(df, neg_4, neg_12)


# 3 Correlograms ----
# To explore the relationships between variables
# NB data is very sparse so there correlations are extremely weak (except on blood tests)
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html 
cor_plot <- function(df,  # e.g. df = df_neg_4$df_bin (binary data)
                     lower){ # whether to include significantly lower outcomes e.g. lower = T
  
  # Transpose df
  df <- as.data.frame(t(df))
  # Exclude 'lower'
  if(lower == F) { df <- df[, !grepl("lower", colnames(df))] }
  # Shorten names (otherwise plot won't generate)
  names(df) <- substring(names(df), 1,25) #shorten names (otherwise plot won't generate)
  # Get correlations
  cors <- cor(df)
  # Plot
  plot <- corrplot(cors, method = 'circle', type = 'full', diag = F,    
                   tl.cex = 0.6, tl.col = "black", mar = c(0,0,0,0))
  # Print colsums
  print(colSums(df))
  
  return(plot)
}

# Generate plots
neg_4_cor <- cor_plot(df_neg_4$df_bin, lower = F)
neg_12_cor <- cor_plot(df_neg_12$df_bin, lower = F)


# 4. K-medoids clustering ----
## Function returns variables' cluster membership and average silhouette width for a given df, and for 2:K clusters
k_medoids <- function(df, K){ # df = e.g. df = df_neg_4
                              # K = number of clusters to check (min = 2) e.g. K = 15
  
  ## Create empty df with variables as row names to capture results in
  all_results <- data.frame(c(row.names(df), "Silhouette", "Dunn"))
  colnames(all_results) <- "Outcome"
  
  ## Create dissimilarity matrix
  diss <- daisy(df, metric="gower")
  
  for(k in 2:K){
    k_meds <- cluster::pam(x = diss, 
                           k = k, # number of clusters
                           diss = T, # indicates whether x is a dissimilarity matrix or not 
                           stand = F) # whether to standardise (normalise) data - already done

    # Get variables' cluster membership
    clusters <- as.data.frame(k_meds$clustering)
    colnames(clusters) <- paste0(k, " clusters")
    
    # Get silhouette width
    sil <- as.data.frame(k_meds$silinfo$avg.width)
    row.names(sil) <- "Silhouette"
    colnames(sil) <- paste0(k, " clusters")
    
    # Get Dunn index
    dunn <- as.data.frame(clValid::dunn(distance = diss, clusters = clusters))
    row.names(dunn) <- "Dunn Index"
    colnames(dunn) <- paste0(k, " clusters")
    
    # Collect results
    results <- rbind(clusters, sil, dunn)
    all_results <- cbind(all_results, results)
  }
  
  return(all_results)
}

## Get clusters, average silhouette width and Dunn index for 2:12 clusters for each df and export
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/5. Cluster analysis/K-medoids")

## Set parameters
k=12 # Max number of clusters to create

kmeds_neg_4 <- k_medoids(df_neg_4$df, k) 
write_csv(kmeds_neg_4, "k-medoids - pos vs neg - 4-12 weeks.csv") # optimal silhouette width  = 2 clusters; optimal Dunn index = 6 clusters

kmeds_neg_12 <- k_medoids(df_neg_12$df, k)
write_csv(kmeds_neg_12, "k-medoids - pos vs neg - 12-26 weeks.csv") # optimal silhouette width  = 2 clusters; optimal Dunn index = 4 clusters


# Helper functions for plotting dendrograms and saving as objects ----
# https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- 90 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}
plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 0.8,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(segment(hcdata)$y, n = 5)
  ymax      <- max(segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  "solid",
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  0.75)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}


# 5. Hierarchical clustering ----
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
# https://www.datacamp.com/tutorial/hierarchical-clustering-R 
hierarchical <- function(df, # e.g. df_neg_4
                         name, # detail included in the names of exported files e.g. "pos vs neg - 4-12 weeks"
                         dist_method, # method used for calculating dissimilarities between observations
                         cluster_method){  # method to be used for clustering e.g. "ward.D" - uses increases in the error sum of squares to determine which clusters should be fused
                         
  # Create distance matrix
  d <- daisy(df, metric = dist_method)
  
  # Hierarchical Clustering
  fit <- hclust(d, method = cluster_method) 
  
  # Create empty df to store silhouette widths in
  sil_means <- NULL
  
  # Get optimal silhouette width
  for(k in 2:12){
    
    # Get silhouette width of each cluster
    sil_width <- silhouette(cutree(fit, k = k), dist = d)[, 3]
    
    # Collect mean silhouette width
    sil_means <- c(sil_means, mean(sil_width))
    
  }
  
  sil_means_tab <- as.data.frame(cbind("K" = c(2:12), "Sil" = sil_means))
  print(sil_means_tab)
  
  # Get k that maximises silhouette width
  k <- sil_means_tab$K[sil_means_tab$Sil == max(sil_means_tab$Sil)]
  
  # Plot dendrogram coloured by k clusters 
  dend_dat <- dendro_data_k(fit, k)
  
  dend_k <- plot_ggdendro(dend_dat,
                          direction   = "tb", # plot orientation (alternative = "lr")
                          #fan = T, # creates a circular plot
                          label.size  = 4, # x-axis font size
                          nudge.label = 0.025, # gap between plot and x-axis labels
                          expand.y    = 2.0) + # space for x-axis labels
    ylab("Cluster distance") +
    xlab(" ") +
    theme_minimal() +
    theme(text = element_text(size=12)) # font size for y-axis title and labels
  
  ## Export dendrogram coloured by k clusters
  dend_k_filename <- paste0("Dendrogram (", k, " clusters) - ", name, " - ", dist_method, " - ", cluster_method, ".png")
  ggsave(dend_k_filename, plot = dend_k, width = 10, height = 10, units = "in", dpi = 300)
  
}


# Prepare and export hierarchical clusters: dendrograms, dendrograms for k clusters, cluster membership
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/5. Cluster analysis/Hierarchical")

## Set parameters
dist_method <- "gower" 
cluster_method <- "ward.D2" 

## Run analysis and export plots and tables (check objects' elements using $fit, $dend, $cluster_membership, $dend_k)
hier_neg_4 <- hierarchical(df_neg_4$df, name = "pos vs neg - 4-12 weeks", dist_method, cluster_method) # Positive vs negative, 4-12 weeks

hier_neg_12 <- hierarchical(df_neg_12$df, name = "pos vs neg - 12-26 weeks", dist_method, cluster_method) # Positive vs negative, >12-26 weeks


