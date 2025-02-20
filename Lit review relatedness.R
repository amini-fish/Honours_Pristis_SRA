setwd("C:/Users/samue/Desktop")

#install.packages("tidyverse")

library(ggplot2)
library(tidyverse)
library(dplyr)

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)

colnames(data)
View(data)
##------------------Index our results based on incl criteria-------------------##

clean_data <- subset(data, Index == "Keep")

clean_data <- clean_data[, -c(1,2, 4, 20:24)]
head(clean_data)

## Make it a tibble so easier to work with

clean_data <- as_tibble(clean_data); clean_data

summary(clean_data)

clean_data$Species <- recode_factor(clean_data$Species, "Glyphis glyphis " = "Glyphis glyphis",  "Squalus albicaudus\xa0"  =  "Squalus albicaudus")

Taxa_data <- clean_data %>% 
  group_by(Super_Order, Order, Species, IUCN.Status)
  #count(Species, sort = T)


#here we derive the number of papers (i.e. the raw totals) per order

table(Taxa_data$Order)

## The workflow can be done any way you want but I think look at patterns across order

##-------------Lets look at what order has been focussed on most--------------##

#here is where we derived the number of species in a given order 

order_data <- clean_data %>%
  group_by(Order, Species, Focus)%>%
  count(Order, sort = T)

table(order_data$Order)
table(order_data$Species)

##----------------Conservation---------------------------##

IUCN <- clean_data %>% 
  group_by(IUCN.Status) %>%
  count(IUCN.Status, sort = T)

sum(IUCN$n)
## Now we have count data for the number of papers each order is focussed on in included papers..
  
## Plot of research interest in each order...

order_plot <- ggplot(order_data, aes(x = reorder(Order, -n, 
                                             FUN = function(x) sum(x)), y = n, fill = Focus)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Research effort in relatedness within elasmobranchs by Order")+
  labs(y = "No. publications", x = "Order") +
  #coord_flip() +
  scale_y_continuous(breaks = seq(0, 55, by=5)) +
  theme_classic()

order_plot <- order_plot + theme(plot.title = element_text(hjust = 0.5, size = 20), 
                                 axis.text = element_text(size = 13), 
                                 axis.title = element_text(size = 15), 
                                 legend.position = c(0.88, 0.88))
print(order_plot)
## Get the number of species from carcharhiniformes 

###-------------------research trends temporally--------------------------------##

clean_data
Year_paper <- clean_data %>%
  group_by(Year, Super_Order) %>%
  count(Year, sort = F)

Year_paper

year_plot <- ggplot(data = Year_paper, aes(x = Year, y = n, colour = Super_Order)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  ggtitle("Number of Papers Investigating Relatedness in Sharks and Rays Over Time") +
  scale_y_continuous(breaks = seq(0, 10, by=1)) +
  scale_x_continuous(breaks = seq(2002, 2024, by = 2))+
  scale_color_brewer(palette = "Dark2") +
  labs(y = "No. publications") +
  theme_classic()

year_plot <- year_plot +
  theme(legend.position = c(0.075, 0.92), 
                               legend.title = element_text(size = 15),
                               legend.text = element_text(size = 12), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15), 
                               plot.title = element_text(hjust = 0.5, size = 17))

print(year_plot)


##-----------------------Markers over time-----------------------------------##

Year_marker <- clean_data %>%
  group_by(Year, Markers) %>%
  count(Year, sort = F)

Year_marker
colnames(Year_marker)

marker_over_time <- ggplot(data = Year_marker, aes(x = Year, y = cumsum(n), col = Markers)) +
  geom_line(stat = "identity")+
  ggtitle("Marker Use Over Time") +
  scale_y_continuous(breaks = seq(0, 80, by=5)) +
  scale_x_continuous(breaks = seq(2002, 2024, by = 2))+
  scale_color_brewer(palette = "Dark2") +
  labs(y = "No. publications") +
  theme_classic()

marker_over_time <- marker_over_time + theme(legend.position = c(0.075, 0.92), 
                         legend.title = element_text(size = 15),
                         legend.text = element_text(size = 12), 
                         axis.text = element_text(size = 13), 
                         axis.title = element_text(size = 15), 
                         plot.title = element_text(hjust = 0.5, size = 17))

print(marker_over_time)

ggplot() + 
  geom_line(aes(x= Year, y=cumsum(n), col = "green"), Year_marker %>% filter(Markers == 'SNPs')) +
  geom_line(aes(x= Year, y=cumsum(n), col = "red"), Year_marker %>% filter(Markers == 'mSats')) +
  geom_line(aes(x= Year, y=cumsum(n), col = "blue"), Year_marker %>% filter(Markers == "mSats + mtDNA")) +
  geom_line(aes(x= Year, y=cumsum(n), col = "purple"), Year_marker %>% filter(Markers == "SNPs + mtDNA_genome")) +
  theme_bw()
 
###-------------------------No Markers---------------------------##

clean_data$Markers <- recode_factor(clean_data$Markers, "SNPS" = "SNPs", "mSats + mtDNA " = "mSats + mtDNA") #"SNPS + mSats + mtDNA" = "SNPs"

markers_data <- clean_data %>%
  group_by(Markers) %>% # Order
  count(Markers, sort = T)

markers_data

marker_plot <- ggplot(data = markers_data) +
  geom_bar(aes(x = reorder(Markers, -n,
                           FUN = function(x) sum(x)), y = n, fill = Markers), stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = seq(0, 60, by=5)) +
  labs(y = "No. publications", x = "Marker Types Used") +
  ggtitle("Markers Used to Assess Relatedness") +
  theme_classic()

marker_plot  <- marker_plot +  theme(legend.position = c(0.85, 0.85), 
                                     legend.title = element_text(size = 15),
                                     legend.text = element_text(size = 12), 
                                     axis.text = element_text(size = 13), 
                                     axis.title = element_text(size = 15), 
                                     plot.title = element_text(hjust = 0.5, size = 17))
print(marker_plot)

##Colour by order

marker_plot <- ggplot(data = markers_data) +
  geom_bar(aes(x = reorder(Markers, -n,
               FUN = function(x) sum(x)), y = n, fill = Order), stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(breaks = seq(0, 60, by=5)) +
  labs(y = "No. publications", x = "Marker Types Used") +
  ggtitle("Markers Used to Assess Relatedness") +
  theme_classic()

marker_plot  <- marker_plot +  theme(legend.position = c(0.85, 0.85), 
                                     legend.title = element_text(size = 15),
                                     legend.text = element_text(size = 12), 
                                     axis.text = element_text(size = 13), 
                                     axis.title = element_text(size = 15), 
                                     plot.title = element_text(hjust = 0.5, size = 17))
print(marker_plot)

## Look at microsatellites

microsats <- clean_data %>%
  select(No..mSats) %>% 
  filter(No..mSats != "NA")

summary(microsats)
sd(microsats$No..mSats)

snps <- clean_data %>% 
  select(No..SNPs) %>% 
  filter(No..SNPs != "NA")

summary(snps)
sd(snps$No..SNPs)

mtDNA <- clean_data %>% 
  select(No..mtDNA.markers) %>% 
  filter(No..mtDNA.markers != "")

mtDNA

##--- Check out hte data ---## 

View(clean_data)

table(clean_data$Kinship.Method)

## ----------------------------Softwares used-------------------------------##
clean_data$Kinship.Method <- recode_factor(clean_data$Kinship.Method, "GERUD & COLONY " = "GERUD & COLONY", "Allele counts  " = "Allele counts", "IR Values & KINSHIP 1.3" = "Kinship 1.3")


software <- clean_data%>%
  group_by(Kinship.Method, Focus) %>%
  count(Kinship.Method)

software
print(software, n = 34)

software_plot <- ggplot(software, aes(x = reorder(Kinship.Method, -n, 
                                                  FUN = function(x) sum(x)), y = n, fill = Focus)) +
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Dark2") +
  theme_classic()

software_plot <- software_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(software_plot)             


no_estimators <- clean_data %>%
  group_by(No..analyses) %>% # Order
  count(No..analyses, sort = T)


no_estimators                        

## percentages 

(31/(sum(no_estimators$n)))*100 #1 estimator
(30/sum(no_estimators$n))*100
(10/sum(no_estimators$n))*100
(2/sum(no_estimators$n))*100
