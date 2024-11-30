## Load in relevant packages 

install.packages("ggpubr")

if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("reshape2")) install.packages("reshape2")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(reshape2)
require(gridExtra)
require(cowplot)

## Load in the data - review2 is for general taxon + numbers 

setwd("C:/Users/samue/Desktop")

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T)

data <- data %>% 
  filter(Index == "Keep")%>%
  droplevels

## Plot 1 - summarise the number of studies of each order, and family of elasmobranch 

# try a dot chart using ggdotchart

# organise the relevant data into a tibble 

tab <- data %>%
  group_by(Order, Family, Species)%>%
  count(Species, sort = F)%>%
  print(n = 50)

write.csv(tab, "Species_Table.csv")


tab_1 <- data %>%
  group_by(Order, Family)%>%
  count(Family, sort = F)%>%
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
           rotate = T, 
           group = "Order", 
           dot.size = 10, 
           label = round(tab_1$n,1),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),           # Adjust label parameters
           ggtheme = theme_bw()                        # ggplot2 theme
)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") 

plot_1 <- plot_1 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, size = 18, margin = margin(10,0,10,0)), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               axis.text = element_text(size = 12), 
               axis.title = element_text(size = 14), 
               axis.title.y = element_text(margin = margin(0,15,0,0)), 
               plot.margin = unit(c(0,1,0,0.5), "cm"))

print(plot_1)

ggsave("plot_1.tiff",
       plot = plot_1,
       width = 30,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)


## try this instead
plot_2 <- ggdotchart(tab_1, x = "Family", y = "n",
           color = "Order",                                # Color by groups
           palette = "", # Custom color palette
           sorting = "descending", 
           title = "Number of studies on each family",
           ylab = "No. of Studies",# Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 10,
           label = round(tab_1$n,1), 
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5), # Color y text by groups
           ggtheme = theme_pubr()                        # ggplot2 theme
)+
  theme_cleveland()  

plot_2 <- plot_2 + theme(legend.position = "bottom", 
                         plot.title = element_text(hjust = 0.5, size = 16, margin = margin(10,0,10,0)), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 12), 
                         axis.title = element_text(size = 14),
                         axis.title.y = element_text(margin = margin(0,10,0,0), angle = 90), 
                         plot.margin = unit(c(0.5,1,1,0.5), "cm"), 
                         panel.background = element_rect(colour = "black", linewidth = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3))

print(plot_2)

ggsave("plot_2.tiff",
       plot = plot_2,
       width = 30,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)


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

  tab_2$Focus <- factor(tab_2$Focus, levels = c("Reproduction", "Popgen", "Demography", "Social"))


plot_3 <- ggdotchart(tab_2, x = "Family", y = "n",
                     color = "Order",
                     group = "Order",
                     facet.by = "Focus", 
                     palette = "Spectral", # Custom color palette
                     sorting = "descending", 
                     title = "Number of studies on each family grouped by research focus",
                     ylab = "No. Studies",# Sort value in descending order
                     rotate = F,
                     add = "segments",
                     dot.size = 8,
                     label = round(tab_2$n,1), 
                     font.label = list(color = "black", size = 10, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_bw()                        # ggplot2 theme
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot_3 <- plot_3 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, size = 18, margin = margin(10,0,10,0)), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               axis.text = element_text(size = 12), 
               axis.title = element_text(size = 14), 
               axis.title.y = element_text(margin = margin(0,10,0,0), angle = 90), 
               axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
               strip.text = element_text(size = 11),
               plot.margin = unit(c(0,1,0,0.5), "cm"))

print(plot_3)

ggsave("plot_3.tiff",
       plot = plot_3,
       width = 30,
       height = 25,
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)

## Now lets look at the markers used (i.e. the distribution between msats and SNPs)

View(data2)



tab_4 <- data2 %>%
  group_by(Markers, Focus) %>%
  count(Focus, sort = T)

tab_4 

tab_4$Focus <- factor(tab_4$Focus, levels = c("Reproduction", "Popgen", "Demography", "Social"))

plot_4 <- ggbarplot(tab_4, x = "Markers", y = "n",
                     fill = "Markers",
                     facet.by = "Focus", 
                    title = "Marker use grouped by research focus",
                    palette = c("grey","orange"), # Custom color palette
                     sorting = "descending", 
                     ylab = "No. Studies",# Sort value in descending order
                     rotate = F,
                     ggtheme = theme_bw()                        
)

plot_4 <- plot_4 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, size = 18, margin = margin(10,0,10,0)), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               axis.text = element_text(size = 12), 
               axis.title = element_text(size = 14), 
               axis.title.y = element_text(margin = margin(0,10,0,0), angle = 90), 
               axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
               strip.text = element_text(size = 11),
               plot.margin = unit(c(0,1,0,0.5), "cm"))


plot_4

ggsave("plot_4.tiff",
       plot = plot_4,
       width = 30,
       height = 25,
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)

## Thingking of alternative ways, you could represent the % of SNP and mSats as stacked in one plot where each bar represents the different focus topics? 

# Give it a whirl

plot_5 <- ggplot(tab_4, aes(x = Focus, y = n, fill = Markers)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Blues") +
  ggtitle("Proportion plot") +
  ylab("Proportion of studies (%)") +
  theme_bw()

plot_5 + theme(legend.position = "bottom", 
                plot.title = element_text(hjust = 0.5, 
                                          size = 16))

## This is okay, but I feel that it can be better - lets try and add a facet grid by year ranges 

# so we need to add year brackets of ~ 5 years? 

data2 <- data2 %>%
  mutate(new_bin = cut(Year, breaks = seq(2000, 2025, by = 5),  # Adjust to your bin size
                       labels = c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025")))


tab_4 <- data2 %>%
  group_by(Markers, Focus, new_bin) %>%
  count(Markers, sort = T)

tab_4 

plot_4 <- ggbarplot(tab_4, x = "Markers", y = "n",
                    fill = "Focus",
                    title ="Changes in Marker Use Over Time",
                    facet.by = "new_bin", 
                    palette = "Spectral", # Custom color palette
                    ylab = "Number of Studies", # Sort value in descending order
                    rotate = F,
                    dot.size = 10,
                    ggtheme = theme_bw()) 

plot_4 <- plot_4 +
  theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, 
                                         size = 16))

print(plot_4)

## I like this, but there could be alterantive ways to present this that look even better? 

# lets try a  dotchart

data2 <- data2 %>%
  mutate(new_bin = cut(Year, breaks = seq(2000, 2025, by = 5),  # Adjust to your bin size
                       labels = c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025")))

tab_4 <- data2 %>%
  group_by(Markers, new_bin) %>%
  count(Markers, sort = T)

tab_4

plot_5 <- ggdotchart(tab_4, x = "Markers", y = "n",
                     color = "Markers",
                     facet.by = "new_bin", 
                     palette = "Spectral", # Custom color palette
                     title = "",
                     ylab = "Number of Studies",# Sort value in descending order
                     rotate = F,
                     add = "segments",
                     dot.size = 8,
                     label = round(tab_4$n,1), 
                     font.label = list(color = "black", size = 8, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_bw()                        # ggplot2 theme
) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

plot_5 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, 
                                         size = 16), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank())

## Sociality

subset_analyses <- data2 %>%
  group_by(No..analyses, Focus) %>%
  count(No..analyses)

subset_analyses

subset_analyses <- subset_analyses %>%
  mutate_at(vars(No..analyses), factor)

plot_6 <- ggbarplot(subset_analyses, x = "No..analyses", y = "n", 
           fill = "No..analyses", 
           facet.by = "Focus", 
           palette = "Spectral", 
           title = "Number of Relatedness Estimators Used By Focus Area", 
           ylab = "Number of Studies",
           xlab = "No. of Analyses",
           ggtheme = theme_bw()) 

plot_6 <- plot_6 + theme(legend.position = "bottom", 
               plot.title = element_text(hjust = 0.5, 
                                         size = 16, 
                                         margin = margin(10,0,10,0)),
               axis.text = element_text(size = 12), 
               axis.title.y = element_text(margin = margin(0,5,0,5)),
               axis.title = element_text(size = 14),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               strip.text.x = element_text(size = 12)
               ) + 
  guides(fill = guide_legend(nrow = 1, title = "No. of Analyses"))

print(plot_6)

ggsave(plot = plot_6, 
       "plot_6.png", 
       dpi = 500, 
       width = 25, 
       height = 25, 
       unit = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots")

## To visualise the estimator use we need to pull in the tibble "estimator" from Lit Review ALL 
## Theres a nice bit of code that fixes all of my shit typos LOL so now we can just load it in
## Bear in mind that you need to use data2 to get accurace faceting based on research focus 

# Load it in

print(estimator, n = 50)

# Visualise a base graphic 

plot(estimator$Kinship.Method, estimator$n, las = 2)


# Load libraries


# Example data

estimator$Kinship.Method <- as.character(estimator$Kinship.Method)

# Extract unique items
unique_items <- unique(unlist(strsplit(paste(estimator$Kinship.Method, collapse = " + "), " \\+ ")))


unique_items

# Initialize matrix
matrix <- matrix(0, nrow = length(unique_items), ncol = length(unique_items),
                 dimnames = list(unique_items, unique_items))


# Populate the matrix
for (i in 1:nrow(estimator)) {
  items <- unlist(strsplit(estimator$Kinship.Method[i], " \\+ "))
  freq <- estimator$n[i]
  
  if (length(items) == 1) {
    # Single factor: Add to diagonal
    matrix[items, items] <- matrix[items, items] + freq
  } else {
    # Combinations: Add to pairwise cells
    for (j in 1:(length(items) - 1)) {
      for (k in (j + 1):length(items)) {
        matrix[items[j], items[k]] <- matrix[items[j], items[k]] + freq
        matrix[items[k], items[j]] <- matrix[items[k], items[j]] + freq # Symmetric
      }
    }
  }
}

# Convert matrix to data frame for ggplot
matrix_df <- melt(matrix)

#matrix_df$value[matrix_df$value == 0] <- NA

# Plot heatmap
plot_7 <- ggplot(matrix_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 1)), color = "black") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  labs(x = "", y = "", fill = "Frequency", 
       title = "Heatmap of Relatedness Estimator Use (Single, Two, and Three Combinations)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, size = 16, margin = margin(10,0,10,0)), 
        axis.text = element_text(size = 12))

print(plot_7)

ggsave("plot_7.tiff",
       plot = plot_7,
       width = 30,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)

## This is a good way to visualise the results but can be a little confusing
## Make sure that there is a detailed explanation of how to interpret it in the figure caption

## Could also just present it as a barplot? 

str(estimator)

unique_items2 <- unique(unlist(strsplit(paste(estimator$Kinship.Method, collapse = " + "), " \\+ ")))

unique_items2

ggdotchart(estimator, 
           x = "Kinship.Method", 
           y = "n",
           color = "Kinship.Method", 
           dot.size = 8, 
           add = "segment", 
           rotate = T, 
           title = "Estimator Use"
           ) +
  theme_bw() +
  theme(legend.position = "bottom") 

## Could try and facet by focus - okay we know this works

data2$Kinship.Method <- recode_factor(data2$Kinship.Method, "Allele counts  " = "Allele counts", "GERUD & COLONY " = "GERUD + COLONY", "COLONY + CERVUS" = "COLONY + CERVUS", "Allele counts & GERUD 1" = "Allele counts + GERUD", "Allele counts & GERUD 2.0" = "Allele counts + GERUD", "GERUD 2.0" = "GERUD", "GERUD 2.0 & COLONY" = "GERUD + COLONY", "Allele counts, GERUD 2.0 & COLONY2" = "Allele counts + GERUD 2.0 + COLONY", "IR Values & KINSHIP 1.3" = "Kinship 1.3", "Kinship 1.3 + Cervus 2.0 " = "Kinship 1.3 + Cervus", "GERUD 1.0, COLONY, STORM & allele counts" = "Allele counts + COLONY + GERUD + STORM", "GERUD 2.0, CERVUS 3.0.7 and COLONY" = "CERVUS + COLONY + GERUD", "Allele counts, GERUD 2.0, COLONY" = "Allele counts + COLONY + GERUD", "GERUD & COLONY" = "COLONY + GERUD", "Allele counts, GERUD 2.0 & COLONY" = "Allele counts + COLONY + GERUD", "Sequoia, COLONY, dartR*" = "COLONY + dartR + Sequoia", "KINFERENCE" = "Kinference", "COLONY2" = "COLONY", "COLONY v.2." = "COLONY", "MLRELATE, COLONY v1.2, KINGROUP 1" = "ML-Relate + COLONY + KinGroup", "KinGroup  " = "KinGroup", "Allele Counts + COLONY + GERUD" = "Allele counts + COLONY + GERUD", "Coancestry + COLONY+ CERVUS " = "Coancestry + COLONY + CERVUS", "COLONY & CERVUS" = "COLONY + CERVUS", "Kinship 1.3 & Cervus" = "Kinship 1.3 + Cervus", "Allele Counts" = "Allele counts", "Allele Counts + ML-Relate" = "Allele counts + ML-Relate")

estimator2 <- data2 %>% 
  group_by(Kinship.Method, Focus)%>%
  count(Kinship.Method, sort= T)


estimator2$Kinship.Method <- recode_factor(estimator2$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS")

estimator2$Kinship.Method <- as.character(estimator2$Kinship.Method)
estimator2$Focus <- as.character(estimator2$Focus)

## Now that we have this data we could try and separate out the combination values as a unique use (rather than a combination...can be used in combination with my heatmap plot)

# Need to use the unique items code - and feed that into a dataframe where each ind occurence = a value 

# ChatGPT gave me some inspo 

expanded_data <- estimator2 %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n, Focus)

expanded_data$Kinship.Method <- recode_factor(expanded_data$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS")

expanded_data$Kinship.Method <- as.character(expanded_data$Kinship.Method)
expanded_data$Focus <- as.character(expanded_data$Focus)

expanded_data <- expanded_data %>%
  group_by(Kinship.Method, Focus) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

expanded_data

## Plot this data 

plot_8 <- ggdotchart(expanded_data, 
           x = "Kinship.Method", 
           y = "Total_Frequency",
           palette = "",
           xlab = "Number of Uses", 
           ylab = "Relatedness Estimator",
           facet.by = "Focus",
           color = "Kinship.Method", 
           dot.size = 7, 
           add = "segment", 
           rotate = T,
           title = "Estimator Use By Focus",  
           label = round(expanded_data$Total_Frequency,1), 
           font.label = list(color = "black", size = 7, 
                             vjust = 0.5), # Color y text by groups
           ggtheme = theme_bw()  
) + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, title = "Relatedness Estimator"))

ggsave("plot_8.tiff",
       plot = plot_8,
       width = 25,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
       )

## All together - no facet 

## Need to reowrk the data to sum by estimator and not focus*estimator

estimator2 <- data2 %>% 
  group_by(Kinship.Method)%>%
  count(Kinship.Method, sort= T)

estimator2$Kinship.Method <- recode_factor(estimator2$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS")

estimator2$Kinship.Method <- as.character(estimator2$Kinship.Method)

expanded_data <- estimator2 %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n)

expanded_data$Kinship.Method <- recode_factor(expanded_data$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS")

expanded_data$Kinship.Method <- as.character(expanded_data$Kinship.Method)

expanded_data <- expanded_data %>%
  group_by(Kinship.Method) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE))

expanded_data

plot_9 <- ggdotchart(expanded_data, 
                     x = "Kinship.Method", 
                     y = "Total_Frequency",
                     palette = "",
                     ylab = "Number of Studies", 
                     xlab = "Relatedness Estimator",
                     color = "Kinship.Method", 
                     dot.size = 10, 
                     add = "segment", 
                     rotate = T,
                     title = "Relatedness Estimator Use",  
                     label = round(expanded_data$Total_Frequency,1), 
                     font.label = list(color = "black", size = 10, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_bw()  
) + 
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16, margin = margin(10,0,10,0)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  guides(colour = guide_legend(nrow = 2, title = "Relatedness Estimator"))

print(plot_9)

## Another really good idea is comparing the marker power for each study where you're standardising the marker powers between SNPs and mSats, one good way to show this is as a boxplot potentially...

## Plot 10 

data$logpwr <- log(data$Power_comp) # note that we log transform to make y axis more easy to interpret

plot_10 <- ggboxplot(data = data, 
                     x = "Markers", 
                     y = "logpwr", 
                     title = "Comparison of Marker Statistical Power",
                     xlab = "Marker Type", 
                     ylab = "log(Marker Power)", 
                     palette = c("grey","orange"),
                     fill = "Markers") +
  theme_bw() +
  geom_jitter(data = data, 
              aes(x = Markers, y = logpwr, size = 7), 
              alpha = 0.5,
              position = position_jitter(height = .2, width = .2)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5, size = 16, margin = margin(10,0,10,0)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(margin = margin(0,10,0,10)),
        axis.title = element_text(size = 15), 
        plot.margin = unit(c(0,1,0,0), "cm")) +
  guides(size = FALSE)

print(plot_10)

ggsave("plot_10.tiff",
       plot = plot_10,
       width = 25,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)


## A cleveland plot could also be great for the overview of estimator use


## Maybe colour by categorical and continuous estimators? 

plot_11 <- ggdotchart(expanded_data, 
                      x = "Kinship.Method", 
                      y = "Total_Frequency",
                     color = "Kinship.Method",                                # Color by groups
                     palette = "", # Custom color palette
                     sorting = "descending", 
                     title = "Relatedness Estimator Use",
                     ylab = "Number of Studies",
                     xlab = "Relatedness Estimator", # Sort value in descending order
                     rotate = TRUE,                                # Rotate vertically
                     dot.size = 10,
                     label = round(expanded_data$Total_Frequency,1), 
                     font.label = list(color = "black", size = 9, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_pubr()                        # ggplot2 theme
)+
  theme_cleveland()  

plot_11 <- plot_11 + theme(legend.position = "bottom", 
                         plot.title = element_text(hjust = 0.5, size = 16, margin = margin(10,0,10,0)), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 12), 
                         axis.title = element_text(size = 14),
                         axis.title.y = element_text(margin = margin(0,10,0,0), angle = 90), 
                         plot.margin = unit(c(0.5,1,1,0.5), "cm"), 
                         panel.background = element_rect(colour = "black", linewidth = 0.5)) +
  guides(colour = FALSE)


print(plot_11)

ggsave("plot_11.tiff",
       plot = plot_11,
       width = 30,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)

## Visualise the markers as a parplot facet wrapped by focus 


## Plot the marker power side by side 

plot_12 <- plot_grid(plot_4, plot_10, 
                     ncol = 2)
  
plot_12

ggsave("plot_12.tiff",
       plot = plot_12,
       width = 40,
       height = 24, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)




plot_13 <- plot_grid(plot_2, plot_11,
                     ncol = 2)

plot_13


## Plot the Conservation status of studies 

# Pull from data

cons_data <- data %>%
  select(IUCN.Status) %>%
  count(IUCN.Status)

cons_data

cons_data$IUCN.Status <- factor(cons_data$IUCN.Status, levels = c("Data deficient", "Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered"))


plot_14 <- ggdotchart(cons_data, x = "IUCN.Status", y = "n",
                      color = "IUCN.Status",
                      palette = c('grey', 'green', 'lightgreen', 'yellow',  'orange', "red"),
                      title = "Number of studies on elasmobranchs by IUCN Status",
                      ylab = "No. of Studies",# Sort value in descending order
                      rotate = F,  
                      sorting = "none",
                      add = "segments",
                      dot.size = 8,
                      label = round(cons_data$n,1), 
                      font.label = list(color = "black", size = 9, 
                                        vjust = 0.5), # Color y text by groups
                      ggtheme = theme_pubr()                        # ggplot2 theme
)

plot_14 <- plot_14 + theme(legend.position = "bottom", 
                         plot.title = element_text(hjust = 0.5, size = 16, margin = margin(10,0,10,0)), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 12), 
                         axis.title = element_text(size = 14),
                         axis.title.y = element_text(margin = margin(0,10,0,0), angle = 90), 
                         plot.margin = unit(c(0.5,1,1,0.5), "cm"), 
                         panel.background = element_rect(colour = "black", linewidth = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size = 3), title = "IUCN"))

print(plot_14)

ggsave("plot_14.tiff",
       plot = plot_14,
       width = 30,
       height = 25, 
       units = "cm", 
       path = "C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots", 
       dpi = 600
)

