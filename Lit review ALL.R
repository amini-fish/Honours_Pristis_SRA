#install.packages("tidyverse")

library(ggplot2)
library(tidyverse)
library(dplyr)


## Custom function to extract certain values in an easier way ##

# Function to search for a keyword in a column and sum corresponding values from another column

sum_if_keyword <- function(data, keyword_column, value_column, keyword) {
  # Filter rows that contain the keyword in the specified column
  matching_rows <- data[grepl(keyword, data[[keyword_column]], ignore.case = TRUE), ]
  
  # Sum the values from the value_column for the matching rows
  total <- sum(matching_rows[[value_column]], na.rm = TRUE)
  
  return(total)
}


## -----------------------------------------------------------------------------

setwd("C:/Users/samue/Desktop")

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)

data <- data %>% 
  filter(Index == "Keep")%>%
  droplevels

## We now have all the data to assess general patterns in species, conservation and order etc...

## Number of studies

nrow(data)

## How many pubs...

nlevels(data$Title)

## How many species

nlevels(data$Species)

## Across how many orders

nlevels(data$Order)

nlevels(data$Family)

## Sharks and Rays 

table(data$Super_Order)

## Now lets look at what papers did what...

table(data$Order)

data %>%
  group_by(Order) %>%
  count(Order, sort = T)

data %>%
  group_by(Order, Family, Species)%>%
  count(Family)%>%
  print(n = 55)

tab_1 <- data %>%
  group_by(Order, Family)%>%
  count(Family)%>%
  print(n = 50)

write.csv(tab_1, file = "litrev_table1.csv")

## Plot 1 - meh 

ggplot(tab_1, 
       aes(forcats::fct_reorder(Family, -n, sum), n, fill = Order)) +
         geom_bar(stat = "identity") +
  theme_classic()


## Super order

table(data$Super_Order)

species_superorder <- data %>%
  group_by(Super_Order) %>%
  count(Species)

print(species_superorder, n = 49)

species_order <- data %>%
  group_by(Super_Order, Order, Species, Focus) %>%
  count(Order)

print(species_order, n = 49)

order_plot <- ggplot(species_order, 
       aes(forcats::fct_reorder(Order, -n, sum), n, fill = Focus)) +
  geom_col() + 
  xlab("Order") +
  ylab("Number of Studies") +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 65, by = 10)) + 
  ggtitle("Number of marker types used in relation to research focus") + 
  theme_classic()

order_plot + theme(axis.text.x = element_text(angle = 90),
                   axis.text = element_text(size = 12),
                   axis.title = element_text(size = 14),
                  plot.title = element_text(hjust = 0.5, size = 16))


species <- data %>%
  group_by(Order) %>%
  count(Species)

## ------------------------------------------------------------------------------

## Results section 1 (all) ###

## I have calculated these for the specific families that have had new studies added in to speed up the recalc process...

levels(data$Order)

## Carcharhiniformes 

Car_all <- data %>%
  filter(Order == c("Carcharhiniformes")) %>%
  group_by(Family, Species)%>%
  count(Family)%>%
  print(n = 50)

tapply(Car_all$n, Car_all$Family, sum)


## Squaliformes

Squ_all <- data %>%
  filter(Order == c("Squaliformes")) %>%
  group_by(Family, Species)%>%
  count(Family)%>%
  print(n = 50)

tapply(Squ_all$n, Squ_all$Family, sum)

## Myliobatiformes 

Myl_all <- data %>% 
  filter(Order == "Myliobatiformes") %>%
  group_by(Family, Species) %>%
  count(Family) %>%
  print(n = 50)


## ------------------------------------------------------------------------------

## Now lets redo the conservation data 

IUCN <- data.frame(data %>% 
  group_by(IUCN.Status) %>%
  count(IUCN.Status, sort = T))


## make a function to create a % value of data
IUCN <- IUCN %>% mutate(SUM=sum(n),
                Percent=n/SUM*100) 

IUCN$Percent <- round(IUCN$Percent, 2)

IUCN

CC.THR <- 18 +14 + 23 #num threatened CE + EN + VU

(CC.THR/83)*100 #61% threatened with extinction

## Prep for plot

IUCN$fraction <- IUCN$n / sum(IUCN$n)
IUCN$ymax <- cumsum(IUCN$fraction)
IUCN$ymin <- c(0, head(IUCN$ymax, n = -1))

## Make a donut plot


Cons_plot <- ggplot(IUCN, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill= IUCN.Status)) +
  geom_rect() +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) +
  ggtitle("Conservation Status of Species Examined") +
  theme_void()

Cons_plot + theme(
  title = element_text(hjust = 1)
)


## Rework some old code to make it nicer

Status <- c("DD", "LC", "NT", "VU", "EN", "CR")

## No species (DD to EX)
count <- c(1, 17,10,23,14,17)

## Merge into a dataframe

data <- data.frame(cbind(Status, count))

## Get the proportions 

str(data)

data$fraction <- as.numeric(data$count)/sum(as.numeric(data$count))

View(data)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

colo <- c("grey", "green", 'lightgreen', 'yellow', 'orange', 'red')


# Make the plot
status <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill = colo)) +
  geom_rect(fill = colo) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=12)) +
  ggtitle("Concervation Status of Elasmobranchs")

status

scale_color_manual(name='Status',
                   breaks=c('DD','LC','NT', 'VU', 'EN', 'CR', 'EX'),
                   values=c('DD' = 'grey', "LC" = 'green', "NT" = 'lightgreen', "VU"= 'yellow', "EN" = 'orange', "CR" = "red", "EX" = 'black'))

status 
## -------------------------------------------------------------------------------

## Marker types used

data %>%
  group_by(Markers)%>%
  count(Markers)

(62/83)*100 # mSat 
(21/83)*100

## Get the average number of each marker summed over the entire dataset 

msat_all <- data %>%
  filter(Markers == "mSat") %>%
  select(No..mSats)

summary(msat_all$No..mSats)
sd(msat_all$No..mSats)

snp_all <- data %>%
  filter(Markers == "SNP") %>% 
  select(No..SNPs)

summary(snp_all$No..SNPs)
sd(snp_all$No..SNPs)

## ------------------------------------------------------------------------------

## ESTIMATOR TIME - see other working code for an interesting outcome...

data$Kinship.Method <- recode_factor(data$Kinship.Method, "Allele counts  " = "Allele counts", "GERUD & COLONY " = "GERUD + COLONY", "COLONY + CERVUS" = "COLONY + CERVUS", "Allele counts & GERUD 1" = "Allele counts + GERUD", "Allele counts & GERUD 2.0" = "Allele counts + GERUD", "GERUD 2.0" = "GERUD", "GERUD 2.0 & COLONY" = "GERUD + COLONY", "Allele counts, GERUD 2.0 & COLONY2" = "Allele counts + GERUD 2.0 + COLONY", "IR Values & KINSHIP 1.3" = "Kinship 1.3", "Kinship 1.3 + Cervus 2.0 " = "Kinship 1.3 + Cervus", "GERUD 1.0, COLONY, STORM & allele counts" = "Allele counts + COLONY + GERUD + STORM", "GERUD 2.0, CERVUS 3.0.7 and COLONY" = "CERVUS + COLONY + GERUD", "Allele counts, GERUD 2.0, COLONY" = "Allele counts + COLONY + GERUD", "GERUD & COLONY" = "COLONY + GERUD", "Allele counts, GERUD 2.0 & COLONY" = "Allele counts + COLONY + GERUD", "Sequoia, COLONY, dartR*" = "COLONY + dartR + Sequoia", "KINFERENCE" = "Kinference", "COLONY2" = "COLONY", "COLONY v.2." = "COLONY", "MLRELATE, COLONY v1.2, KINGROUP 1" = "ML-Relate + COLONY + KinGroup", "KinGroup  " = "KinGroup", "Allele Counts + COLONY + GERUD" = "Allele counts + COLONY + GERUD", "Coancestry + COLONY+ CERVUS " = "Coancestry + COLONY + CERVUS", "Allele Counts" = "Allele counts", "COLONY & CERVUS" = "COLONY + CERVUS", "Kinship 1.3 & Cervus" = "Kinship 1.3 + Cervus")


## Can add in Focus to get specific groupings, or could try tapply and use as grouping factor


estimator <- data %>% 
  group_by(Kinship.Method)%>%
  count(Kinship.Method, sort= T)

print(estimator, n = 45)

order <- estimator$n
order

## Do it this shitty way until I remember the r solution

COLONY_all <- sum_if_keyword(estimator, "Kinship.Method", "n", "COLONY"); COLONY_all

GERUD_all <- sum_if_keyword(estimator, "Kinship.Method", "n", "GERUD"); GERUD_all

  
## -----------------------------------------------------------------------------

## Get COLONY combinations 

COL_all <- data %>%
  filter(grepl("COLONY", Kinship.Method))%>%
  count(Kinship.Method, sort = T) %>%
  print()

sum(COL_all$n)

## % of all studies 

(sum(COL_all$n)/83)*100

## Same for GERUD 

GER_all <- data %>%
  filter(grepl("GERUD", Kinship.Method))%>%
  count(Kinship.Method, sort = T) %>%
  print()

sum(GER_all$n)

(sum(GER_all$n)/83)*100

## Allele Count 

AC_all <- data %>%
  filter(grepl("Allele counts", Kinship.Method))%>%
  count(Kinship.Method, sort = T) %>%
  print()

sum(AC_all$n)

(sum(AC_all$n)/83)*100

## ----------------------------------------------------------------------------

## Calculate %
estimator <- data.frame(estimator)
estimator <- estimator %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100) 

estimator$Percent <- round(estimator$Percent, 2)

estimator

est_plot <- 
ggplot(estimator, 
       aes(forcats::fct_reorder(Kinship.Method, -n, sum), n, fill = Focus)) +
  geom_col(colour="black") + 
  xlab("Estimator") +
  ylab("No. Papers") +
  scale_fill_brewer(palette = "GnBu") +
  scale_y_continuous(breaks = 1:12) + 
  ggtitle("Summary of estimators used by each research focus") + 
  theme_classic()

est_plot + theme(axis.text.x = element_text(angle = 90), 
                 plot.title = element_text(hjust = 0.5, size = 16))

print(estimator, n = 31)

## -----------------------------------------------------------------------------
## Number of estimators used per study

no_estimators <- data %>% 
  group_by(No..analyses)%>%
  count(No..analyses, sort= T)

print(no_estimators, n = 45)

## Calculate %
no_estimators <- data.frame(no_estimators)

no_estimators <- no_estimators %>% mutate(SUM=sum(n),
                                  Percent=n/SUM*100) 

no_estimators$Percent <- round(no_estimators$Percent, 2)

no_estimators

##------------------------- 

## Need to redo the marker type calculations 

all_fams <- data %>% 
  group_by(Family, Focus)%>% 
  count(Family)


all_fams$Carch <- ifelse(all_fams$Family == "Carcharhinidae", "Carcharhinidae", "Other")

all_fams$Carch

print(all_fams, n = 100)

sp_plot <- ggplot(all_fams, aes(forcats::fct_reorder(Family, -n, sum), n, fill = Family)) + 
  geom_bar(stat = "identity") +
  xlab("Species") +
  ylab("Number of Studies") +
  scale_y_continuous(limits = c(0, 40), n.breaks = 10) + 
  ggtitle("Species studied in relation to research focus") + 
  theme_classic()

sp_plot + theme(axis.text.x = element_text(angle = 90), 
                 plot.title = element_text(hjust = 0.5, size = 16)) 
## Try the Carch vs Other plot 

sp_plot <- ggplot(all_fams, 
                  aes(forcats::fct_reorder(Carch, n, sum), n, fill = Focus)) +
  geom_col() + 
  xlab("Species") +
  ylab("No. Papers") +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_y_continuous(limits = c(0, 48), n.breaks = 14) + 
  ggtitle("Species studied in relation to research focus") + 
  theme_classic()

sp_plot + theme(axis.text = element_text(size = 12), 
                axis.title = element_text(size = 13.5),
                legend.text = element_text(size = 12), 
                legend.title = element_text(size = 12), 
                plot.title = element_text(hjust = 0.5, size = 16)) 


d <- c(16,8,10, 11,7)
sd(d)
mean(d)

## Markers

data$Markers <- recode_factor(data$Markers, "mSats + mtDNA " = "mSats + mtDNA", "SNPS" = "SNPs", "mSats + mtDna_genome" = "mSats + mtDNA", "SNPs + mtDNA_genome" = "SNPs + mtDNA", "SNPs + mtDNA " = "SNPs + mtDNA")

marker <- data.frame(data %>% 
  group_by(Markers, Focus)%>%
  count(Markers))

marker

marker_plot <- ggplot(marker, 
                  aes(forcats::fct_reorder(Markers, -n, sum), n, fill = Focus)) +
  geom_col() + 
  xlab("Marker") +
  ylab("No. Papers") +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 5)) + 
  ggtitle("Marker use grouped by study focus") + 
  theme_classic()

marker_plot + theme(axis.text = element_text(size = 12), 
                    axis.title = element_text(size = 13.5),
                    legend.title = element_text(size = 12.5), 
                    legend.text = element_text(size = 11),
                plot.title = element_text(hjust = 0.5, size = 16), 
                legend.position = "bottom") 

## Now we claculate % to 2 DP for everything...

marker2 <- marker %>% mutate(SUM=sum(n),
                        Percent=n/SUM*100) 

marker2$Percent <- round(marker2$Percent, 2)

marker2

## Number of markers 
no_marker <- data %>% 
  group_by(No..Marker.Types, Focus)%>%
  count(No..Marker.Types, sort= T)

print(no_marker, n = 45)

## Calculate %
no_marker <- data.frame(no_marker)

no_marker <- no_marker %>% mutate(SUM=sum(n),
                                  Percent=n/SUM*100) 

no_marker$Percent <- round(no_marker$Percent, 2)

no_marker$No..Marker.Types <- as.factor(no_marker$No..Marker.Types)

nomarker_plot <- ggplot(no_marker, 
                      aes(forcats::fct_reorder( No..Marker.Types, -n, sum), n, fill = Focus)) +
  geom_col(colour="black") + 
  xlab("No. Marker Types Used") +
  ylab("No. Papers") +
  scale_fill_brewer(palette = "BuGn") +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 65, by = 5)) + 
  ggtitle("Number of marker types used in relation to research focus") + 
  theme_classic()

nomarker_plot + theme(, 
                    plot.title = element_text(hjust = 0.5, size = 16))

## -----------------------------------------------------------------------------

snps <- data %>%
 select(No..SNPs)

mean(snps$No..SNPs, na.rm = T)
summary(snps, na.rm = T)
sd(snps$No..SNPs, na.rm = T)

snps 

levels(data$Markers)
## Need to work out what year the SNP papers were published

snp_year <- data %>%
  group_by(Markers, Year)%>%
  count(Year, sort= T) %>%
  filter(Markers != "mSats", Markers != "mSats + mtDNA") %>%
  na.omit()

snp_year

##-----------------------------------------------------------------------------

## Microsat nitty gritty

microsat <- data %>%
  group_by(No..mSats, Focus) %>%
  count(No..mSats) %>%
  na.omit()

microsat

msat_plot <- ggplot(data = microsat, aes(x = No..mSats, y = Focus, fill = Focus)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0.2, size = 2) +
  scale_fill_brewer(palette = "BuGn")+
  ggtitle("Number of mSat loci used in studies of relatedness") +
  scale_x_continuous(limits = c(0, 26), breaks = seq(0, 26, by = 4)) +
  xlab("Number of mSat loci") +
  theme_classic()


msat_plot <- msat_plot + theme(
                  axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14),
                  plot.title = element_text(hjust = 0.5, size = 16))

print(msat_plot)

## No samples used 

samples <- data %>% 
  group_by(n_samples, Focus) %>% 
  count(n_samples) %>%
  na.omit()

summary(samples)

samples$n_samples <- as.numeric(as.character(samples$n_samples))

n_samples <- ggplot(data = samples, aes(x = n_samples, y = Focus, fill = Focus)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, height = 0.2, size = 2.5, alpha = 0.7) +
  scale_fill_brewer(palette = "BuGn")+
  ggtitle("Sample sizes of studies investigating relatedness") +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 100)) +
  xlab("Number of Individuals") +
  geom_vline(xintercept = 50, linetype = 5, linewidth = 1, col = "red") +
  theme_classic()


n_samples <- n_samples + theme(axis.text.x = element_text(angle = 90),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  plot.title = element_text(hjust = 0.5, size = 18))

print(n_samples)

##-----------------------------------------------------------------------------

library(egg)
library(cowplot)

## Marker power plot 

Marker_PWR <- data %>%
  select(Power_comp, Year, Markers)

Marker_PWR

PWR_comp <- 
  ggplot(data = Marker_PWR, aes(x = Year, y = Power_comp, shape = Markers)) + 
  geom_point(stat = "identity", size = 3, aes(col = Markers, alpha = 0.7)) +
  scale_colour_brewer(palette = "Dark2") +
  ggtitle("Comparison of Marker Power between mSats and SNPs") +
  scale_x_continuous(limits = c(2002, 2024), breaks = seq(2002, 2024, by = 2)) +
  scale_y_continuous(limits = c(0, 9000), breaks = seq(0, 9000, by = 750)) +
  xlab("Year") +
  ylab("Power in terms of no. of SNPs") +
  theme_bw()

PWR_comp <- PWR_comp + theme(
  plot.title = element_text(size = 16, hjust = 0.5), 
  axis.title = element_text(size= 13), 
  axis.text = element_text(size = 12), 
  legend.title = element_text(size = 12.5), 
  legend.text = element_text(size = 11.5)
)

print(PWR_comp)

## Now the inset plot for the mSat pattern 

PWR_comp2 <- 
  ggplot(data = Marker_PWR, aes(x = Year, y = Power_comp, shape = Markers, alpha = 0.7)) + 
  geom_point(stat = "identity", size = 3, aes(col = Markers)) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(2002, 2024), breaks = seq(2002, 2024, by = 4)) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 25)) +
  theme_bw()

PWR_comp2 <- PWR_comp2 + theme(
  plot.title = element_blank(), 
  axis.title = element_blank(), 
  axis.text = element_text(size = 8, angle = 0),
  legend.position = "none"
) 

print(PWR_comp2)

plot.with.inset <-
  ggdraw() +
  draw_plot(PWR_comp) +
  draw_plot(PWR_comp2, x = 0.07, y = 0.52, width = .5, height = .42)

print(plot.with.inset)

## ----------------------------------------------------------------------------

## How many estimators did each study use? 
## Of those with > 1, did they use a diverse range i.e. Ped and Mark? 

colnames(data)

estimator_type <- data %>% 
  group_by(No..analyses, Method_PedvsMark) %>%
  count(No..analyses)

estimator_type

## some lazy calculations 

estimator_type$Two <- ifelse(estimator_type$No..analyses > 1, "T", "F")

estimator_type %>%
  group_by(Two)

tapply(estimator_type$n, estimator_type$Two, FUN=sum)
tapply(estimator_type$n, estimator_type$Method_PedvsMark, FUN=sum)


#Updated at 29/10/2024 ---

## ----------------------------------------------------------------------------

## Reproduction 

repro <- data%>%
  filter(grepl("Reproduction", Focus))%>%
  group_by(Species)%>%
  count(Species)%>%
  print(n = 45)

sum(repro$n)

repro <- data%>%
  filter(grepl("Reproduction", Focus))%>%
  group_by(Order, Family, Species)%>%
  count(Species)

sum(repro$n)


############################## Rep 2.0 ###########################

rep <- data %>%
  filter(grepl("Reproduction", Focus))

rep
## Taxonomy 

rep_tax <- rep %>%
  group_by(Order, Family, Species)%>%
  count(Order)

print(rep_tax, n = 50)

# Carcharhiniformes

carch_rep <- rep %>%
  filter(Order == "Carcharhiniformes") %>%
  select(Order, Family, Species) %>% 
  count(Order, Family, Species)

carch_rep

sum(carch_rep$n)

## Now we need the number of studies by family 

cc_sums <- carch_rep %>%
  filter(Family == "Carcharhinidae")

sum(cc_sums$n)

tri_sums <- carch_rep %>%
  filter(Family == "Triakidae")
sum(tri_sums$n)


## Conservation

IUCN_rep <- rep %>%
  group_by(IUCN.Status)%>%
  count(IUCN.Status)


CE <- (6/sum(IUCN_rep$n))*100 # CE
EN <- (10/sum(IUCN_rep$n))*100 # EN
VU <- (11/sum(IUCN_rep$n))*100 # VU

sum(CE, EN, VU)

(7/sum(IUCN_rep$n))*100 # NT


## Markers 

rep %>%
  group_by(Markers) %>%
  count(Markers)

mean(rep$No..mSats, na.rm = T)
sd(rep$No..mSats, na.rm = T)

## Samples 

## Estimators 

rep <- rep %>%
  select(Kinship.Method)%>%
  count(Kinship.Method)


rep$Kinship.Method <- as.character(rep$Kinship.Method)

rep_expand <- rep %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n)


rep_expand$Kinship.Method <- recode_factor(rep_expand$Kinship.Method, "Allele Counts" = "Allele counts", "Cervus" = "CERVUS")


rep_expand$Kinship.Method <- as.character(rep_expand$Kinship.Method)

rep_expand <- rep_expand %>%
  group_by(Kinship.Method) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

rep_expand

## Sample size rep

rep <- data %>%
  filter(grepl("Reproduction", Focus))

summary(as.numeric(as.character(rep$n_samples)), na.rm = T)
sd(as.numeric(as.character(rep$n_samples)), na.rm = T)

summary(as.numeric(as.character(rep$N_litters)), na.rm = T)
sd(as.numeric(as.character(rep$N_litters)), na.rm = T)

summary(as.numeric(as.character(rep$avg_litter_size)), na.rm = T)
sd(as.numeric(as.character(rep$avg_litter_size)), na.rm = T)

#------------------------------------------------------------------------------

## Popgen 

levels(data$Focus)

## How many of each focus are there?

data %>% 
  filter(grepl("Popgen", Focus)) %>%
  group_by(Focus) %>%
  count(Focus, sort = T)

## COnservation status 

data%>%
  filter(grepl("Popgen", Focus))%>%
  group_by(IUCN.Status)%>%
  count(IUCN.Status)%>%
  print(n = 20)

## Species

popgen_sub <- data%>%
  filter(grepl("Popgen", Focus))%>%
  group_by(Order, Family )%>%
  count(Species, sort = F) ## Count by order, can change to family for simplicities sake

print(popgen_sub)

nlevels(popgen_sub$Family)

sum(popgen_sub$n)

## Marker Use 

data %>%
filter(grepl("Popgen", Focus))%>%
  group_by(Markers) %>%
  count(Markers)

## Marker by year

data%>%
  filter(grepl("Popgen", Focus))%>%
  group_by(Markers, Year)%>%
  count(Markers)

## Focus by year 

Years_gen <- data%>% 
  filter(grepl("Popgen", Focus))%>%
  select(Year)

summary(as.numeric(as.character(Years_gen$Year)))

#mSats

popgen_mSats <- data %>% 
  filter(grepl("Popgen", Focus))%>%
  select(No..mSats)

mean(popgen_mSats$No..mSats, na.rm = T)
sd(popgen_mSats$No..mSats, na.rm = T)

#SNPs

popgen_SNPs <- data %>% 
  filter(grepl("Popgen", Focus))%>%
  select(No..SNPs, Year)

popgen_SNPs

mean(popgen_SNPs$No..SNPs, na.rm = T)
sd(popgen_SNPs$No..SNPs, na.rm = T)

summary(popgen_SNPs$No..SNPs)

## Sample sizes

popgen_N <- data %>% 
  filter(grepl("Popgen", Focus)) %>%
  select(n_samples)

summary(as.numeric(as.character(popgen_N$n_samples)))
sd(as.numeric(as.character(popgen_N$n_samples)))

## Estimators used

Popgen_est <- data %>%
  filter(grepl("Popgen", Focus)) %>%
  group_by(Kinship.Method) %>%
  count(Kinship.Method, sort = T)

Popgen_est$Kinship.Method <- as.character(Popgen_est$Kinship.Method)

popgen_mark_expand <- Popgen_est %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n)

popgen_mark_expand$Kinship.Method <- recode_factor(popgen_mark_expand$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS", "Kingroup" = "KinGroup")

popgen_mark_expand$Kinship.Method <- as.character(popgen_mark_expand$Kinship.Method)

popgen_mark_expand <- popgen_mark_expand %>%
  group_by(Kinship.Method) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

popgen_mark_expand

sum_if_keyword(Popgen_est, "Kinship.Method", "n", "COLONY") 

sum_if_keyword(Popgen_est, "Kinship.Method", "n", "Coancestry")

## -------------------------------------------------------------------------------


## Demography 

data %>% 
  filter(grepl("Demography", Focus)) %>%
  group_by(Order, Family)%>%
  count(Family)

subset_demo <- data %>%
  filter(grepl("Demography", Focus)) 

summary(as.numeric(as.character(subset_demo$n_samples)))
sd(as.numeric(as.character(subset_demo$n_samples)))

# Markers 

demo <- data %>% 
  filter(grepl("Demography", Focus)) %>%
  select(No..SNPs)


mean(demo$No..SNPs)
sd(demo$No..SNPs)

#-------------------------------------------------------------------------------

## Sociality 

soc <- data %>%
  filter(grepl("Social", Focus))
## Taxonomy 

soc %>%
  group_by(Order, Family, Species)%>%
  count(Order)

## Conservation

soc %>%
  group_by(IUCN.Status, Species)%>%
  count(Species)

## Markers 

soc %>%
  group_by(Markers) %>%
  count(Markers)

mean(soc$No..mSats, na.rm = T)
sd(soc$No..mSats, na.rm = T)

## Samples 

## Estimators 

soc <- soc %>%
  select(Kinship.Method)%>%
  count(Kinship.Method)


soc$Kinship.Method <- as.character(soc$Kinship.Method)

soc_expand <- soc %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n)

soc_expand$Kinship.Method <- as.character(soc_expand$Kinship.Method)

soc_expand <- soc_expand %>%
  group_by(Kinship.Method) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")


soc_expand

## Sample size soc

summary(as.numeric(as.character(subset_soc$n_samples)))
sd(as.numeric(as.character(subset_soc$n_samples)))

#-------------------------------------------------------------------------------

## Reproduction: 

## -------------------------------------------------------------------------

## MISC

eda_time <- data%>%
  group_by( Super_Order, Year) %>%
  count(Super_Order)
  

print(eda_time, n = 30)

ggplot(data = eda_time, aes(x = Year, y = n, col = Super_Order)) +
  geom_point()


focus <- data %>%
  filter(grepl("Reproduction", Focus)) %>%
  group_by(Focus, Super_Order, Title) %>%
  
  count(Super_Order)

print(focus, n = 90)


## Quick plot 

ggplot(all_fams, aes(x = Family, y = n, fill = Family)) +
  geom_bar(stat = "identity")

## -- Quick facet grid marker plot

marker_plot_data <- data %>%
  group_by(Markers, Focus, Year) %>%
  count(Markers)

print(marker_plot_data)

ggplot(data = marker_plot_data, aes(x = Markers, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ Focus)

ggplot(data = marker_plot_data, aes(x = Year, y = n, col = Markers)) +
  geom_point(stat = "identity")

data %>%
  group_by(Markers) %>%
  count(Markers)


data %>%
  group_by(Focus) %>%
  count(Focus)


