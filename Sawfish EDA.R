setwd("C:/Users/samue/Desktop/Honours - Sawfish/pristis_data")

#install.packages("hierfstat")
#devtools::install_version("ggplot2", "3.4.4")
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)


#save(gl, file = "Sawfish Prelim Data.RData") 

meta <- read.csv("Sawfish_meta2.csv"); meta

View(meta)

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

meta

ggplot(data = meta, mapping = aes(TL_cm, 
                                  fill = sex)) +
  geom_histogram() +
  theme_classic()


