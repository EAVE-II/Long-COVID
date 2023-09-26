##############################################################################################################
## Project: Long COVID
## Code author(s): karen.jeffrey@ed.ac.uk
## Description: 01_text_analysis - Identify frequently-used terms that occurr in GP free text alongside the 
## long-COVID Read code. Terms will be used by Albasoft to identify long-COVID in GP free text in the full dataset.
##############################################################################################################

# Install.packages
library(tidyverse)
library(stringi) # for cleaning data
library(quanteda) # for text analysis
library(wordcloud) # word-cloud generator 


## lemmatize = combine words with the same stems https://stackoverflow.com/questions/61257802/r-text-mining-grouping-similar-words-using-stemdocuments-in-tm-package 
## https://data.library.virginia.edu/a-beginners-guide-to-text-analysis-with-quanteda/ basics
## https://quanteda.io/articles/quickstart.html detail
## Corpus = all texts to be analysed
## Types = unique words
## Tokens = words

# Read in data ----------
df <- readRDS("/conf/EAVE/GPanalysis/analyses/long_covid/data/free_text.rds") # redacted sample of GP free text wint long COVID Read codes supplied by Albasoft
# df <- read.delim("/conf/EAVE/GPanalysis/analyses/long_covid/data/Sick notes - free text.txt", header = TRUE, sep = "\t", dec = ".") # sample of sick notes text
colnames(df) <- "Txt"
df$Txt <- as.character(df$Txt)

# Clean text -------
df$Txt <- stringi::stri_replace_all_regex(df$Txt, "\\d", "") # Remove numbers
df$Txt <- stringi::stri_replace_all_regex(df$Txt, "[\\p{p}\\p{S}]", "") # Remove punctuation and symbols
df$Txt <- stringi::stri_trans_tolower(df$Txt) # Make all lowercase

# Text analysis --------
## Create a corpus and identify text field 
corp <- corpus(df, text_field = "Txt")

## Identify words as tokens
doc.tokens <- tokens(corp, what = "word")

## Remove custom words and stopwords
doc.tokens <- tokens_remove(doc.tokens, pattern=c("c", stopwords('en')), padding = TRUE) 

## Collect n-grams
doc.tokens <- tokens_ngrams(doc.tokens, n=1:6, concatenator = " ") # collect 1,2,3,4,5,6-word phrases

## Create a document feature matrix (dfm)
doc.dfm <- dfm(doc.tokens)

## Explore
topfeatures(doc.dfm, 100) # Top 100 n-grams
top100 <- data.frame(ngram = labels(topfeatures(doc.dfm, 100)), count = topfeatures(doc.dfm, 100))

## Export
#setwd("/conf/EAVE/GPanalysis/analyses/long_covid/outputs/4. Free text")
setwd("/conf/EAVE/GPanalysis/analyses/long_covid/code_review/Review outputs/4. Free text")
write_csv(top100, "top_100_ngrams.csv")
write_csv(df, "free_text_cleaned.csv")

## Wordcloud --------
# http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know#step-5-generate-the-word-cloud
set.seed(1234)
wordcloud <- wordcloud(words = top100$ngram, freq = top100$count, min.freq = 1,
          max.words=800, random.order=FALSE, colors=brewer.pal(4, "RdYlBu"), scale=c(3.,0.25))
