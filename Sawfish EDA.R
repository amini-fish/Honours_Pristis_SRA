setwd("C:/Users/samue/Desktop/Honours_Sawfish/analysis")

#install.packages("hierfstat")
#devtools::install_version("ggplot2", "3.4.4")
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)


#save(gl, file = "Sawfish Prelim Data.RData") 

meta <- read.csv("Sawfish_meta2.csv"); meta

View(meta)

?subset

meta <- subset(meta, pop == "Daly", select = c(id, sex, TL_cm))
meta

cols <- c("skyblue", "orange")

ggplot(data = meta, aes(x = sex, y = TL_cm)) +
  geom_boxplot(fill = cols) +
  ggtitle("Size distrbutions (TL) between sexes") +
  xlab("Sex") +
  ylab("Total Length (cm)") +
  theme_minimal()


lapply(meta, FUN = mean)



  ## Lets just explore our data using some simple EDA

## Look at the number inds from each sex by population

# Just save the sex data with populations as a frequency

tab.sex <- as.data.frame(table(meta$sex, meta$pop))
colnames(tab.sex) <- c("sex", "pop", "count"); tab.sex

# I can also simply visualise this as a barplot 
ggplot(data = tab.sex, mapping = aes(x = sex, 
                                     y = count,
                                     fill = sex)) +
  geom_bar(stat = "identity")

# or we can visualise this by population

ggplot(data = tab.sex, mapping = aes(x = sex, 
                                  y = count, 
                                  fill = pop)) +
  geom_bar(stat = "identity", 
           position = "fill") +
  coord_flip() +
  theme_bw()

?geom_bar

## Histogram of the size classes
group_mean<- aggregate(x= meta$TL_cm,
                       # Specify group indicator
                       by = list(meta$sex),      
                       # Specify function (i.e. mean)
                       FUN = median)
print(group_mean)
colnames(group_mean) <- c("sex", "median")

group_mean <- group_mean[-3,]; 

p <- ggplot(data = meta, mapping = aes(TL_cm, fill = sex)) +
  geom_histogram(col = "black", bins = 20) +
  ggtitle("Distribution of sizes between male and female P. pristis") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 13)) +
  labs(y = "Count", x = "Total Length (cm)") +
  geom_vline(data = group_mean, aes(xintercept = median, color= sex),
             linetype="dashed", size = 1.2)

p 

 
p <-ggplot(meta, aes(x= TL_cm, fill =sex)) +
  geom_histogram(col = "black", position="dodge")+
  geom_vline(data= group_mean, aes(xintercept= median, color=sex),
             linetype="dashed", size = 1.2)+
  theme(legend.position="top")
p
p + scale_color_brewer(palette="Dark2")

?subset
meta <- subset(meta, sex != "U")

p <- ggplot(data = meta, aes(x = sex, y = TL_cm, fill = sex)) +
    geom_boxplot(aes(fill = factor(sex))) +
  ggtitle("Distribution of sizes between male and female P. pristis") +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 13)) +
  labs(y = "Total Length (Cm)", x = "Sex")

p

q <- ggplot(data = meta, aes(x = sex, y = TL_cm, fill = sex)) +
  geom_violin() +
  ggtitle("Distribution of sizes between male and female P. pristis") +
  theme(plot.title = element_text(hjust = 0.5, size = 18), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 13)) +
  labs(y = "Total Length (Cm)", x = "Sex")
q

## Plot them next to each other
#library("gridExtra")
?subset
grid.arrange(p,q, ncol = 2)

#Lets calculate some rough descriptive states for relatedness and other predictor variables e.g. TL mean
library(tidyr)
## First we subset the data to only include our Daly samples 
meta
daly <- subset(meta, pop == "Daly", select = c(id, sex, TL_cm, Cohort, pop))
daly # looks good (bare bones metadata)

#lets save it as a data frame..

daly <- as.data.frame(daly); daly

## Replace the NA's with an approximate birth cohort...
daly[is.na(daly)] <- 2008; daly

summary(daly)

group_mean<- aggregate(x= daly$TL_cm,
                       # Specify group indicator
                       by = list(daly$sex),      
                       # Specify function (i.e. mean)
                       FUN = sd)
print(group_mean)

table(daly$sex)

16/29
