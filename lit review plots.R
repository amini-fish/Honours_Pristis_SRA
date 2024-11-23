## Load in relevant packages 

install.packages("ggpubr")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

## Load in the data - review2 is for general taxon + numbers 

setwd("C:/Users/samue/Desktop")

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)

data <- data %>% 
  filter(Index == "Keep")%>%
  droplevels

## Plot 1 - summarise the number of studies of each order, and family of elasmobranch 

# try a dot chart using ggdotchart

# organise the relevant data into a tibble 

tab_1 <- data %>%
  group_by(Order, Family)%>%
  count(Family, sort = T)%>%
  print(n = 50)

# plot it 
tab_1

plot_1 <- ggdotchart(tab_1, x = "Family", y = "n", 
           color = "Order", 
           title = "Number of studies on each family",
           ylab = "No. Studies",
           palette = "Spectral",
           sorting = "descending",
           add = "segments", 
           rotate = F, 
           group = "Order", 
           dot.size = 10, 
           label = round(tab_1$n,1),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),           # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") 

plot_1 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, 
                                  size = 15))


## try this instead
plot_1 <- ggdotchart(tab_1, x = "Family", y = "n",
           color = "Order",                                # Color by groups
           palette = "Spectral", # Custom color palette
           sorting = "descending", 
           title = "Number of studies on each family",
           ylab = "No. Studies",# Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 10,
           label = round(tab_1$n,1), 
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5), # Color y text by groups
           ggtheme = theme_pubr()                        # ggplot2 theme
)+
  theme_cleveland()  

plot_1 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, 
                                    size = 15))


# IDEA - add silhouette of body plans per order to contextualise? 
  # will need to do that in post

## Okay I think that between these I could make somthing work 

# FOCUS AREAS - we need review3 data 

data2 <- read.csv("relatedness_literature_review_working3.csv", stringsAsFactors = T)

data2 <- data2 %>% 
  filter(Index == "Keep")%>%
  droplevels

tab_2 <- data2 %>%
  group_by(Order, Family, Focus)%>%
  count(Family, sort = T)%>%
  print(n = 50)

tab_2

plot_3 <- ggdotchart(tab_2, x = "Family", y = "n",
                     color = "Order",
                     group = "Order",
                     facet.by = "Focus", 
                     palette = "Spectral", # Custom color palette
                     sorting = "descending", 
                     title = "Number of studies on each family",
                     ylab = "No. Studies",# Sort value in descending order
                     rotate = F,
                     add = "segments",
                     dot.size = 10,
                     label = round(tab_2$n,1), 
                     font.label = list(color = "black", size = 9, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_bw()                        # ggplot2 theme
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot_3 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, 
                                         size = 15))
