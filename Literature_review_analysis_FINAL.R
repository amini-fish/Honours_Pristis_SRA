## Packages 
install.packages("ggpubr")
install.packages("gt")
install.packages("gtExtras")

if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("reshape2")) install.packages("reshape2")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(reshape2)
require(gridExtra)
require(cowplot)
require(gt)
require(gtExtras)

#-------------------------------------------------------------------------------
## Set working directory

# My desktop - or wherever you have the data stored

setwd("C:/Users/samue/Desktop")

#-------------------------------------------------------------------------------

## Data Loading 

# Dataset 1 - this is used for analysis where we don't need to split studies by category (i.e. duplicate the study in each research category)

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)

# Data_2 will be used for analysis when we need to group things/count by research category. This means we have duplicated the studies where there is > 1 category 

data_2 <- read.csv("relatedness_literature_review_working3.csv", stringsAsFactors = T)

#-------------------------------------------------------------------------------

## Data Cleaning & Filtering 

# Now we can delete the studies that didn't make the cut (e.g. grey lit)

data <- data %>% 
  filter(Index == "Keep")%>%
  droplevels

# estimator filtering 

data$Kinship.Method <- recode_factor(data$Kinship.Method, "Allele counts  " = "Allele counts", "GERUD & COLONY " = "GERUD + COLONY", "COLONY + CERVUS" = "COLONY + CERVUS", "Allele counts & GERUD 1" = "Allele counts + GERUD", "Allele counts & GERUD 2.0" = "Allele counts + GERUD", "GERUD 2.0" = "GERUD", "GERUD 2.0 & COLONY" = "GERUD + COLONY", "Allele counts, GERUD 2.0 & COLONY2" = "Allele counts + GERUD 2.0 + COLONY", "IR Values & KINSHIP 1.3" = "Kinship 1.3", "Kinship 1.3 + Cervus 2.0 " = "Kinship 1.3 + Cervus", "GERUD 1.0, COLONY, STORM & allele counts" = "Allele counts + COLONY + GERUD + STORM", "GERUD 2.0, CERVUS 3.0.7 and COLONY" = "CERVUS + COLONY + GERUD", "Allele counts, GERUD 2.0, COLONY" = "Allele counts + COLONY + GERUD", "GERUD & COLONY" = "COLONY + GERUD", "Allele counts, GERUD 2.0 & COLONY" = "Allele counts + COLONY + GERUD", "Sequoia, COLONY, dartR*" = "COLONY + dartR + Sequoia", "KINFERENCE" = "Kinference", "COLONY2" = "COLONY", "COLONY v.2." = "COLONY", "MLRELATE, COLONY v1.2, KINGROUP 1" = "ML-Relate + COLONY + KinGroup", "KinGroup  " = "KinGroup", "Allele Counts + COLONY + GERUD" = "Allele counts + COLONY + GERUD", "Coancestry + COLONY+ CERVUS " = "Coancestry + COLONY + CERVUS", "COLONY & CERVUS" = "COLONY + CERVUS", "Kinship 1.3 & Cervus" = "Kinship 1.3 + Cervus", "Allele Counts" = "Allele counts", "Allele Counts + ML-Relate" = "Allele counts + ML-Relate")

# Same again for data_2

data_2 <- data_2 %>% 
  filter(Index == "Keep")%>%
  droplevels

## Fixing estimator names

data_2$Kinship.Method <- recode_factor(data_2$Kinship.Method, "Allele counts  " = "Allele counts", "GERUD & COLONY " = "GERUD + COLONY", "COLONY + CERVUS" = "COLONY + CERVUS", "Allele counts & GERUD 1" = "Allele counts + GERUD", "Allele counts & GERUD 2.0" = "Allele counts + GERUD", "GERUD 2.0" = "GERUD", "GERUD 2.0 & COLONY" = "GERUD + COLONY", "Allele counts, GERUD 2.0 & COLONY2" = "Allele counts + GERUD 2.0 + COLONY", "IR Values & KINSHIP 1.3" = "Kinship 1.3", "Kinship 1.3 + Cervus 2.0 " = "Kinship 1.3 + Cervus", "GERUD 1.0, COLONY, STORM & allele counts" = "Allele counts + COLONY + GERUD + STORM", "GERUD 2.0, CERVUS 3.0.7 and COLONY" = "CERVUS + COLONY + GERUD", "Allele counts, GERUD 2.0, COLONY" = "Allele counts + COLONY + GERUD", "GERUD & COLONY" = "COLONY + GERUD", "Allele counts, GERUD 2.0 & COLONY" = "Allele counts + COLONY + GERUD", "Sequoia, COLONY, dartR*" = "COLONY + dartR + Sequoia", "KINFERENCE" = "Kinference", "COLONY2" = "COLONY", "COLONY v.2." = "COLONY", "MLRELATE, COLONY v1.2, KINGROUP 1" = "ML-Relate + COLONY + KinGroup", "KinGroup  " = "KinGroup", "Allele Counts + COLONY + GERUD" = "Allele counts + COLONY + GERUD", "Coancestry + COLONY+ CERVUS " = "Coancestry + COLONY + CERVUS", "COLONY & CERVUS" = "COLONY + CERVUS", "Kinship 1.3 & Cervus" = "Kinship 1.3 + Cervus", "Allele Counts" = "Allele counts", "Allele Counts + ML-Relate" = "Allele counts + ML-Relate")

#-------------------------------------------------------------------------------

#####################
## Section 1 - ALL ##
#####################

nlevels(data$Title) # Number of papers 

nrow(data) # Number of studies

table(data$Super_Order) # Superorder 

nlevels(data$Order) # Number of orders
 
nlevels(data$Family) # Number of families

nlevels(data$Species) # Number of species

min(data$Year) #Year min

max(data$Year) #Year max 

# Table of Order, Family, Species

# Number of studies by order

data %>%
  group_by(Order) %>%
  count(Order, sort = T)

# To family level 

data %>%
  group_by(Order, Family)%>%
  count(Family, sort = T)%>%
  print(n = 50)

# To species level 

data %>%
  group_by(Order, Family, Species)%>%
  count(Family, sort = T)%>%
  print(n = 55)

## Extinction risk - all

IUCN <- data.frame(data %>% 
                     group_by(IUCN.Status) %>%
                     count(IUCN.Status, sort = T))

IUCN <- IUCN %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100) 

IUCN$Percent <- round(IUCN$Percent, 2); IUCN

Threatened <- 18+14+22 ; Threatened # total threatened 

82 - Threatened #total non-threatened (exlcuding DD)

# Marker use - all

markers_all <- data.frame(data %>% 
  group_by(Markers) %>%
  count(Markers))

markers_all %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100)

# msat summary stats
msat_all <- data %>%
  filter(Markers == "mSat") %>%
  select(No..mSats)

summary(msat_all$No..mSats)
sd(msat_all$No..mSats)

# snp summary stats 

snp_all <- data %>%
  filter(Markers == "SNP") %>% 
  select(No..SNPs)

summary(snp_all$No..SNPs)
sd(snp_all$No..SNPs)

# Estimators :) 

# We also need the breakdown of which studies used one, two, or three estimators to get the number of studies and the % 

analyses_num <- data.frame(data %>% 
                             group_by(No..analyses) %>%
                             count(No..analyses))

analyses_num %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100)

# But which estimators? 

# Make a dataset for it 

estimator <- data %>% 
  group_by(Kinship.Method)%>%
  count(Kinship.Method, sort= T)
 
print(estimator, n = 45) # quick visual check

estimator$Kinship.Method <- as.character(estimator$Kinship.Method)

# Here we're essentially telling R to split anything with a + and place it into the column named Kinship.Method and then selecting only estiamtor names and their counts

estimator <- estimator %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n)

# But there is a catch 
# There will be some data cleaning to do now that we have split things...it's just a part of it but can be prone to errors.

estimator$Kinship.Method <- recode_factor(estimator$Kinship.Method, "CERVUS " = "CERVUS", "Allele counts  " = "Allele counts", "CKMRsim " = "CKMRsim", "KinGroup  " = "KinGroup", "Sequioa " = "Sequoia")

# Now, after some more visual checks we're confident we have removed any typos we can transform and count our estimator use properly

estimator$Kinship.Method <- as.character(estimator$Kinship.Method)

#Group it by the name, summarise it by calculating the sum of N for every row that each estimator occurs in (rather than counting how many rows it occurs in - easy mistake to make). Arrange it in descending order. 

estimator <- estimator %>%
  group_by(Kinship.Method) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE)) %>%
  arrange(desc(Total_Frequency))

# Print it :)

print(estimator, n = 100)

# The above should give us a broken down count of how many times each estimator was used, whether in isolation or in combination with another, this means the number wont add to 83 (because more than one estimator was used  57% of the time). 

#-------------------------------------------------------------------------------

## Section 2 - Categories (will get to this later as less urgent)

# Rinse and repeat

est_foc <- data_2 %>% 
  group_by(Kinship.Method, Focus)%>%
  count(Kinship.Method, sort= T)

est_foc$Kinship.Method <- recode_factor(est_foc$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS")

est_foc$Kinship.Method <- as.character(est_foc$Kinship.Method)
est_foc$Focus <- as.character(est_foc$Focus)

expanded_est <- est_foc %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n, Focus)

expanded_est$Kinship.Method <- recode_factor(expanded_est$Kinship.Method, "CERVUS " = "CERVUS", "Allele counts  " = "Allele counts", "CKMRsim " = "CKMRsim", "KinGroup  " = "KinGroup", "Sequioa " = "Sequoia")

expanded_est$Kinship.Method <- as.character(expanded_est$Kinship.Method)
expanded_est$Focus <- as.character(expanded_est$Focus)

expanded_est <- expanded_est %>%
  group_by(Kinship.Method, Focus) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

expanded_est

expanded_est$Focus <- recode_factor(expanded_est$Focus, "Popgen" = "Population Genetics", "Social" = "Sociality")

expanded_est$Focus <- factor(expanded_est$Focus, levels = c("Reproduction", "Population Genetics", "Demography", "Sociality"))

print(expanded_est, n = 50)


#-------------------------------------------------------------------------------

## Section 3 - Reproductive behavior - use data_2 

# first thing is first - we need to get all the studies focused on repro 

reproduction <- data_2 %>%
  filter(Focus == "Reproduction")

reproduction

# Rinse and repeat 

# Taxonomy

reproduction %>%
  group_by(Order) %>%
  count(Order, sort = T)

# To family level 

reproduction %>%
  group_by(Order, Family)%>%
  count(Family, sort = T)%>%
  print(n = 50)

# To species level 

reproduction %>%
  group_by(Order, Family, Species)%>%
  count(Family, sort = T)%>%
  print(n = 55)

## Extinction risk 

IUCN_repro <- data.frame(reproduction %>% 
                     group_by(IUCN.Status) %>%
                     count(IUCN.Status, sort = T))

IUCN_repro <- IUCN_repro %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100) 

IUCN_repro$Percent <- round(IUCN_repro$Percent, 2); IUCN_repro

10 + 10 + 6

20.41 + 20.41 + 12.24 #threatened %

7 + 15

30.61 + 14.29 + 2.04

# Alrighty - estimators...

# we just need to use our tibble named expanded_est here...filter it by the column "focus" and we have our info

expanded_est %>%
  filter(Focus == "Reproduction") %>%
  arrange(desc(Total_Frequency))

#BOOM!!

# Sample sizes

# 

#-------------------------------------------------------------------------------

## Section 4 - Population Genetics 

popgen <- data_2 %>%
  filter(Focus == "Popgen")

popgen

# Taxonomy 

  ## Order

popgen %>%
  group_by(Order) %>%
  count(Order, sort = T)

  ## Family

popgen %>%
  group_by(Order, Family)%>%
  count(Family, sort = T)%>%
  print(n = 50)

  ## Species 

popgen %>%
  group_by(Order, Family, Species)%>%
  count(Family, sort = F)%>%
  print(n = 55)

# Extinction risk 

IUCN_pg <- data.frame(popgen %>% 
                           group_by(IUCN.Status) %>%
                           count(IUCN.Status, sort = T))

IUCN_pg <- IUCN_pg%>% mutate(SUM=sum(n),
                                    Percent=n/SUM*100) 

IUCN_pg$Percent <- round(IUCN_pg$Percent, 2); IUCN_pg

  ## Quick calculations 

14 +10+3 # Threatened studies
43.75 + 31.25 + 9.38 #threatened % 
9.38 + 6.25 # rest % 

# Markers 

# Sample sizes 

# Estimators 

#-------------------------------------------------------------------------------

## Section 5 - Demography

# Subset 

demo <- data_2 %>%
  filter(Focus == "Demography")

# Taxonomy 

# IUCN 

IUCN_demo <- data.frame(demo %>% 
                           group_by(IUCN.Status) %>%
                           count(IUCN.Status, sort = T))

IUCN_demo <- IUCN_demo %>% mutate(SUM=sum(n),
                                    Percent=n/SUM*100) 

IUCN_demo$Percent <- round(IUCN_demo$Percent, 2); IUCN_demo

# Markers

# Samples Sizes 

# Estimators

#-------------------------------------------------------------------------------

## Section 6 - Sociality 

