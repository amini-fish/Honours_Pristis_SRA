#### Install Packages ####

install.packages("ggpubr")
install.packages("gt")
install.packages("gtExtras")

if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("reshape2")) install.packages("reshape2")

#### Load Packages into R ####

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(reshape2)
require(gridExtra)
require(cowplot)
require(gt)
require(gtExtras)
library(viridis)

#### Set working directory ####

# My desktop - or wherever you have the data stored

setwd("C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/Metadata")

#### Data Loading ####

# Dataset 1 - this is used for analysis where we don't need to split studies by category (i.e. duplicate the study in each research category)

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)

# Data_2 will be used for analysis when we need to group things/count by research category. This means we have duplicated the studies where there is > 1 category 

data_2 <- read.csv("relatedness_literature_review_working3.csv", stringsAsFactors = T)

#### Data Cleaning & Filtering ####

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
  count(Family, sort = F)%>%
  print(n = 50)

# To species level 

  data %>%
  group_by(Order, Family, Species)%>%
  count(Family, sort = F)%>%
  print(n = 70)

## Extinction risk - all

IUCN <- data.frame(data %>% 
                     group_by(IUCN.Status) %>%
                     count(IUCN.Status, sort = T))

IUCN <- IUCN %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100) 

IUCN$Percent <- round(IUCN$Percent, 2); IUCN

Threatend <- 23+19+29 ; Threatened # total threatened 

(Threatend/sum(IUCN$n))*100

106 - Threatened #total non-threatened (exlcuding DD)

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

summary(snp_all$No..SNPs, na.rm = T)
sd(snp_all$No..SNPs, na.rm = T)
base::sd(snp_all$No..SNPs)

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
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% 
  unnest(Kinship.Method) %>%
  group_by(Kinship.Method) %>%
  summarise(n = sum(n), .groups = "drop")%>%
  arrange(desc(n))

print(estimator, n = 70)
# But there is a catch 
# There will be some data cleaning to do now that we have split things...it's just a part of it but can be prone to errors.

estimator$Kinship.Method <- recode_factor(estimator$Kinship.Method, "CERVUS " = "CERVUS", "Allele counts  " = "Allele counts", "CKMRsim " = "CKMRsim", "KinGroup  " = "KinGroup", "Sequioa " = "Sequoia", "Cervus" = "CERVUS")

# Now, after some more visual checks we're confident we have removed any typos we can transform and count our estimator use properly

estimator$Kinship.Method <- as.character(estimator$Kinship.Method)

print(estimator, n = 50)

#Group it by the name, summarise it by calculating the sum of N for every row that each estimator occurs in (rather than counting how many rows it occurs in - easy mistake to make). Arrange it in descending order. 

estimator <- estimator %>%
  group_by(Kinship.Method) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE)) %>%
  arrange(desc(Total_Frequency))

estimator$Kinship.Method

Kinship.Method <- unique(estimator$Kinship.Method)

Estimators


Cat_Con <- c("Cat", "Cat", "Con", "Cat", "Cat/Con", 
                       "Cat", "Con", "Con","Cat", "Con", "Cat", 
                       "Cat", "Cat", "Con", "Cat", "Cat", "Con", "Con", "Con")


Cat_Con_frame <- data.frame(Kinship.Method, Cat_Con)

table(Cat_Con)

# Print it :)

print(estimator, n = 100)

# The above should give us a broken down count of how many times each estimator was used, whether in isolation or in combination with another, this means the number wont add to 83 (because more than one estimator was used  57% of the time). 


## Marker power data 

power_comp <- data %>%
  select(Power_comp, Markers)

summary(power_comp$Power_comp, na.rm = T)
sd(power_comp$Power_comp, na.rm = T)

# this works, but the range is so great that it jumps up a couple of orders of magnitude...I could log transfrom it? 

# it is also a good idea to do some kind of statistical test (t-test) 

####--------------------Section 2 Categories--------------------####

# Overall number of studies by research category

data_2 %>%
  group_by(Focus) %>%
  count(Focus, sort = T) %>%
  mutate(Focus, SUM = sum(n), 
         Percent = n/SUM*100)

# Taxonomy by focus - species

# Order 
data_2 %>%
  group_by(Focus, Order) %>%
  count(Order) %>%
  print(n = 100)

# Family 

data_2 %>%
  group_by(Focus, Family) %>%
  count(Family) %>%
  print( n = 100)

# Species

data_2 %>%
  group_by(Focus,Species) %>%
  count(Species) %>%
  print(n = 100)

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

expanded_est$Kinship.Method <- recode_factor(expanded_est$Kinship.Method, "CERVUS " = "CERVUS", "Allele counts  " = "Allele counts", "CKMRsim " = "CKMRsim", "KinGroup  " = "KinGroup", "Sequioa " = "Sequoia", "Cervus" = "CERVUS")

expanded_est$Kinship.Method <- as.character(expanded_est$Kinship.Method)
expanded_est$Focus <- as.character(expanded_est$Focus)

expanded_est <- expanded_est %>%
  group_by(Kinship.Method, Focus) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

expanded_est

expanded_est$Focus <- recode_factor(expanded_est$Focus, "Popgen" = "Population Genetics", "Social" = "Sociality")

expanded_est$Focus <- factor(expanded_est$Focus, levels = c("Reproduction", "Population Genetics", "Demography", "Sociality"))

print(expanded_est, n = 50)

####----------- Section 3 - Reproductive behavior ----------------####

# first thing is first - we need to get all the studies focused on repro 

reproduction <- data_2 %>%
  filter(Focus == "Reproduction") %>%
  droplevels()

reproduction

# Rinse and repeat 

nlevels(reproduction$Order) # Number of orders

nlevels(reproduction$Family) # Number of families

nlevels(reproduction$Species)

# Taxonomy

reproduction %>%
  group_by(Order) %>%
  count(Order, sort = T)

# To family level 

reproduction %>%
  group_by(Order, Family)%>%
  count(Family, sort = F)%>%
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

Rep_Th <- 6 + 10 + 11

Rep_Th/50

20.41 + 20.41 + 12.24 #threatened %

7 + 15

30.61 + 14.29 + 2.04

## Markers

  # Marker types 

reproduction %>%
  group_by(No..Marker.Types) %>%
  count(No..Marker.Types)

reproduction %>%
  group_by(Markers) %>%
  count(Markers)
  
  # Summary Stats

summary(reproduction$No..mSats, na.rm = T)
sd(reproduction$No..mSats, na.rm = T)

summary(reproduction$No..SNPs, na.rm = T)
sd(reproduction$No..SNPs, na.rm = T)

  # Marker Power


## Sample sizes

str(reproduction)

summary(as.numeric(as.character(reproduction$n_samples, na.rm = T)))

sd(as.numeric(as.character(reproduction$n_samples, na.rm = T)))

sd(reproduction$n_samples)

mode <- function(x, na.rm = FALSE) {
  
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}

  
mode(as.character(reproduction$n_samples, na.rm = T))

mode(reproduction$N_litters, na.rm = T)



# Estimators 

# we just need to use our tibble named expanded_est here...filter it by the column "focus" and we have our info

nest_rep <- data.frame(reproduction %>%
             group_by(No..analyses) %>%
             count(No..analyses))

nest_rep %>% mutate(SUM=sum(n),Percent=n/SUM*100)


expanded_est %>%
  filter(Focus == "Reproduction") %>%
  arrange(desc(Total_Frequency)) #BOOM!!



#-------------------------------------------------------------------------------

## Section 4 - Population Genetics 

popgen <- data_2 %>%
  filter(Focus == "Popgen")%>%
  droplevels()

popgen

# Taxonomy 


nlevels(popgen$Order) # Number of orders

nlevels(popgen$Family) # Number of families

nlevels(popgen$Species)

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

# Sub-category breakdown 

popgen %>%
  group_by(Popgen_focus) %>%
  count(Popgen_focus)


# Extinction risk 

IUCN_pg <- data.frame(popgen %>% 
                           group_by(IUCN.Status) %>%
                           count(IUCN.Status, sort = T))

IUCN_pg <- IUCN_pg%>% mutate(SUM=sum(n),
                                    Percent=n/SUM*100) 

IUCN_pg$Percent <- round(IUCN_pg$Percent, 2); IUCN_pg

  ## Quick calculations 

18 + 8 + 16 # Threatened studies

42/54 

# Markers 

popgen %>%
  group_by(No..Marker.Types) %>%
  count(No..Marker.Types)

popgen %>%
  group_by(Markers) %>%
  count(Markers)

# Summary Stats

summary(popgen$No..mSats, na.rm = T)
sd(popgen$No..mSats, na.rm = T)

summary(popgen$No..SNPs, na.rm = T)
sd(popgen$No..SNPs, na.rm = T)

## Sample sizes

str(reproduction)

summary(as.numeric(as.character(popgen$n_samples, na.rm = T)))

sd(as.numeric(as.character(popgen$n_samples, na.rm = T)))

sd(reproduction$n_samples)

# Sample sizes 

# Estimators 

## Number of estimators 
nest_pg <- data.frame(popgen %>%
                         group_by(No..analyses) %>%
                         count(No..analyses))

nest_pg %>% mutate(SUM=sum(n),Percent=n/SUM*100)

  ## Detailed breakdown

expanded_est %>%
  filter(Focus == "Population Genetics") %>%
  arrange(desc(Total_Frequency))

#-------------------------------------------------------------------------------

## Section 5 - Demography

# Subset 

demo <- data_2 %>%
  filter(Focus == "Demography")%>%
  droplevels()

nlevels(demo$Order) # Number of orders

nlevels(demo$Family) # Number of families

nlevels(demo$Species)

## Order

demo%>%
  group_by(Order) %>%
  count(Order, sort = T)

## Family

demo %>%
  group_by(Order, Family)%>%
  count(Family, sort = T)%>%
  print(n = 50)

## Species 

demo %>%
  group_by(Order, Family, Species)%>%
  count(Family, sort = F)%>%
  print(n = 55)

# IUCN 

IUCN_demo <- data.frame(demo %>% 
                           group_by(IUCN.Status) %>%
                           count(IUCN.Status, sort = T))

IUCN_demo <- IUCN_demo %>% mutate(SUM=sum(n),
                                    Percent=n/SUM*100) 

IUCN_demo$Percent <- round(IUCN_demo$Percent, 2); IUCN_demo

57.14 + 28.57 

# Markers

  ## SNP number 

summary(demo$No..SNPs, na.rm = T)
sd(demo$No..SNPs, na.rm = T)

# Samples Sizes 

summary(as.numeric(demo$n_samples), na.rm = T)
sd(as.numeric(demo$n_samples), na.rm = T)

# Estimators

  ## Detailed breakdown

expanded_est %>%
  filter(Focus == "Demography") %>%
  arrange(desc(Total_Frequency))

#-------------------------------------------------------------------------------

## Section 6 - Sociality 

# Subset 

soc <- data_2 %>%
  filter(Focus == "Social")

# Taxonomy 

## Order

soc %>%
  group_by(Order) %>%
  count(Order, sort = T)

## Family

soc %>%
  group_by(Order, Family)%>%
  count(Family, sort = T)%>%
  print(n = 50)

## Species 

soc %>%
  group_by(Order, Family, Species)%>%
  count(Family, sort = F)%>%
  print(n = 55)

# Sub-category breakdown 


# Extinction Risk

IUCN_soc <- data.frame(soc %>% 
                          group_by(IUCN.Status) %>%
                          count(IUCN.Status, sort = T))

IUCN_soc <- IUCN_soc %>% mutate(SUM=sum(n),
                                  Percent=n/SUM*100) 

IUCN_soc$Percent <- round(IUCN_soc$Percent, 2); IUCN_soc

42.86 + 28.57 + 14.29


# Markers 

soc %>%
  group_by(No..Marker.Types) %>%
  count(No..Marker.Types)

soc %>%
  group_by(Markers) %>%
  count(Markers)

summary(soc$No..mSats, na.rm = T)
sd(soc$No..mSats, na.rm = T)

summary(soc$No..SNPs, na.rm = T)
sd(soc$No..SNPs, na.rm = T)

soc %>%
  select(No..SNPs)

# Samples Sizes 

summary(as.numeric(as.character(soc$n_samples)))

sd(as.numeric(as.character(soc$n_samples)))

# Sample Sizes 

# Estimators

N_est_soc <- data.frame(soc %>%
                         group_by(No..analyses) %>%
                         count(No..analyses))

N_est_soc %>% mutate(SUM=sum(n),Percent=n/SUM*100)


expanded_est %>%
  filter(Focus == "Sociality") %>%
  arrange(desc(Total_Frequency)) #BOOM!!


(1361 + 3057)/2

sd(c(1361,3057))
#-------------------------------------------------------------------------------

## Section 7 - Miscellaneous 

# Species table 1 

Species_table <- data_2 %>%
  group_by(Order, Family, Species, Focus) %>%
  count(Species) %>%
  pivot_wider(names_from = Focus,values_from=n) %>%
  mutate(across(everything(), as.character)) %>%
  mutate_at(c("Reproduction", "Popgen", "Demography", "Social"), ~replace_na(.,""))

Species_table

write.csv(Species_table, file = "species_table.csv", row.names = F)

## Suplementary Figures

## Supp Fig 1 - IUCN status of all & each category...

n_est <- data_2 %>%
  select(No..analyses, Focus)%>%
  count(No..analyses, Focus)

n_est

n_est$Focus <- recode_factor(n_est$Focus, "Popgen" = "Population Genetics", "Social" = "Sociality")

n_est$Focus <- factor(n_est$Focus, levels = c("Reproduction", "Population Genetics", "Demography", "Sociality"))


supp_1 <- ggbarplot(n_est,
                   x = "No..analyses", 
                   y = "n", 
                   fill = "Focus", 
                   palette = "Spectral", 
                   label = F,
                   ylab = "No. Studies",
                   xlab = "No. Analyses",
                   position = position_dodge(0.9)) +
  theme_bw()

supp_1 <- supp_1 + theme(legend.position = "bottom", 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 12), 
                         axis.title = element_text(size = 12), 
                         axis.title.y = element_text(margin = margin(0,15,0,0)), 
                         plot.margin = unit(c(1,1,1,1), "cm")) 

print(supp_1)


ggsave("supp_1.tiff",
       plot = supp_1,
       width = 28,
       height = 20, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 1000
)


## Supp Fig 2 - Sample Sizes

sample_size <- data_2 %>%
  select(n_samples, Focus)

sample_size$Focus <- recode_factor(sample_size$Focus, "Popgen" = "Population Genetics", "Social" = "Sociality")

sample_size$Focus <- factor(sample_size$Focus, levels = c("Reproduction", "Population Genetics", "Demography", "Sociality"))

sample_size$log_n <- log(sample_size$n_samples)

supp2 <- ggboxplot(sample_size,
          x = "Focus", 
          y = "n_samples", 
          fill = "Focus", 
          ylab = "No. Samples",
          palette = "Spectral") +
  theme_bw()

supp2 <- supp2 + 
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        axis.title.y = element_text(margin = margin(0,15,0,0)), 
        plot.margin = unit(c(1,1,1,1), "cm")) +
  guides(size = FALSE, fill = FALSE)

print(supp2)

ggsave("supp_2.tiff",
       plot = supp2,
       width = 28,
       height = 20, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 1000
)

# Log

supp_2.2 <- ggboxplot(sample_size,
                   x = "Focus", 
                   y = "log_n", 
                   fill = "Focus", 
                   ylab = "log(No. Samples)",
                   palette = "Spectral") +
  theme_bw()

print(supp_2.2)

supp_2.2 <- supp_2.2 + 
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        axis.title.y = element_text(margin = margin(0,15,0,0)), 
        plot.margin = unit(c(1,1,1,1), "cm")) +
  guides(size = FALSE, fill = FALSE)

print(supp_2.2)

ggsave("supp_2.2.tiff",
       plot = supp_2.2,
       width = 28,
       height = 20, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 1000
)


## Supp Fig 3 - 

## Supp Fig 4 - 

## Supp Fig 5 - 

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)
data$Method_PedvsMark

df <- data %>%
  dplyr::group_by(Method_PedvsMark)%>%
  count(Method_PedvsMark)


df$Method_PedvsMark <- as.character(df$Method_PedvsMark)

df


## Seperate the values and get counts by Focus for each type but need to combine mutliples of same estomator type

df_2 <- df %>%
  mutate(Method_PedvsMark = strsplit(Method_PedvsMark, " \\+ ")) %>% # Split the combinations
  unnest(Method_PedvsMark) %>%
  select(Method_PedvsMark, n)

## Merge the same estimator types using the group_by function...then sum the n column for these rows to get the total count for each estimator type, in each research category.
df_2$Method_PedvsMark <- recode_factor(df_2$Method_PedvsMark, "Continuous " = "Continuous") 

df_3 <- df_2 %>%
  group_by(Method_PedvsMark) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

head(df_3)


# Extract unique items
unique_items <- unique(unlist(strsplit(paste(df$Method_PedvsMark, collapse = " + "), " \\+ ")))

unique_items


df_3$Focus <- recode_factor(df_3$Focus, "Popgen" = "Population Genetics", "Social" = "Social Behaviour", "Reproduction" = "Reproductive Behaviour")

df_3$Focus <- factor(df_3$Focus, levels = c("Reproductive Behaviour", "Population Genetics", "Demography", "Social Behaviour"))

df_3

##------------------------- Without focus but with year -------------------------##

data <- data %>%
  mutate(new_bin = cut(Year, breaks = seq(2000, 2025, by = 5),
                       labels = c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025")))

## Just to get how many of each estimator 

df <- data %>%
  dplyr::group_by(Method_PedvsMark, new_bin)%>%
  count(Method_PedvsMark)

df$Method_PedvsMark <- as.character(df$Method_PedvsMark)

df

## Seperate the values and get counts by Focus for each type but need to combine mutliples of same estomator type

df_2 <- df %>%
  mutate(Method_PedvsMark = strsplit(Method_PedvsMark, " \\+ ")) %>% # Split the combinations
  unnest(Method_PedvsMark) %>%
  select(new_bin, Method_PedvsMark, n)

df_2

## Merge the same estimator types using the group_by function...then sum the n column for these rows to get the total count for each estimator type, in each research category.
df_2$Method_PedvsMark <- recode_factor(df_2$Method_PedvsMark, "Continuous " = "Continuous") 

df_3 <- df_2 %>%
  group_by(new_bin, Method_PedvsMark) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

print(df_3)

palette_sa <- c ("grey","orange", "skyblue")

plot <- ggbarplot(df_3, x = "Method_PedvsMark", y = "Total_Frequency",
                    fill = "Method_PedvsMark",
                  width = 0.5,
                    position = position_dodge(preserve = "single"),
                    facet.by = "new_bin", 
                    palette = palette_sa,
                    ylab = "Frequency of use", 
                  xlab = "Estimator Type",
                    rotate = F,
                    dot.size = 10,
                    ggtheme = theme_bw())

plot <- plot + theme(legend.position = "none",
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 11), 
                         axis.title = element_text(size = 11), 
                         axis.text.x = element_text(angle = 0),
                         axis.title.y = element_text(margin = margin(0,15,0,0)), 
                         plot.margin = unit(c(1,1,1,1), "cm"), 
                         strip.text = element_text(size = 11)) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3))

print(plot)

ggsave("estimator_type_temporal.png",
       plot = plot,
       width = 20,
       height = 20,
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 2000
)

##------------------------- Estimators over time -------------------------##

estimator <- data %>% 
  group_by(Kinship.Method, new_bin)%>%
  count(Kinship.Method, sort= T)

print(estimator, n = 45) # quick visual check

estimator$Kinship.Method <- as.character(estimator$Kinship.Method)

# Here we're essentially telling R to split anything with a + and place it into the column named Kinship.Method and then selecting only estiamtor names and their counts

estimator <- estimator %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% 
  unnest(Kinship.Method) %>%
  group_by(Kinship.Method, new_bin) %>%
  summarise(n = sum(n), .groups = "drop")%>%
  arrange(desc(n))

print(estimator, n = 70)
# But there is a catch 
# There will be some data cleaning to do now that we have split things...it's just a part of it but can be prone to errors.

estimator$Kinship.Method <- recode_factor(estimator$Kinship.Method, "CERVUS " = "CERVUS", "Allele counts  " = "Allele counts", "CKMRsim " = "CKMRsim", "KinGroup  " = "KinGroup", "Sequioa " = "Sequoia", "Cervus" = "CERVUS")

# Now, after some more visual checks we're confident we have removed any typos we can transform and count our estimator use properly

estimator$Kinship.Method <- as.character(estimator$Kinship.Method)

estimator <- estimator %>%
  group_by(Kinship.Method, new_bin) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE)) %>%
  arrange(desc(Total_Frequency))

estimator <- estimator %>%
  left_join(Cat_Con_frame, by = "Kinship.Method") %>%
  arrange(Kinship.Method)

print(estimator)

palette_sa <- c ("skyblue", "orange", "grey")

plot <- ggplot(estimator, aes(x = Kinship.Method, y = Total_Frequency, fill = Cat_Con.y)) +
  geom_bar(stat = "identity", col = "black", width = 0.7, position = position_dodge(preserve = "single")) +
  facet_wrap(~ new_bin, scales = "free_x") +
  scale_fill_manual(values = palette_sa) +
  labs(y = "Frequency of use", x = "Estimator Type") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
    axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    strip.text = element_text(size = 11)
  ) + guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3))

print(plot)

?scale_fill_viridis
## Without facet

library(ggpubr)
