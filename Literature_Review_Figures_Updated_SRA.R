## Literature Review Figures 27.03 ##

#### Load packages into R ####
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
library(ggplot2)
library(ggthemes)
library(wesanderson)
library(devEMF)


#### Load data (both data & data_2) ####

data <- read.csv("relatedness_literature_review_working2.csv", stringsAsFactors = T) ## Used when not subsetting by category

data_2 <- read.csv("relatedness_literature_review_working3.csv", stringsAsFactors = T) ## Use when subsetting by category

#### Figure 1 - Taxonomy ####

data %>%
  group_by(Order)%>%
  count(Order, sort = T)%>%
  print(n = 50)

tab_1 <- data %>%
  group_by(Order, Family)%>%
  count(Order, sort = F)%>%
  print(n = 50)

zissou_colors <- wes_palette("Zissou1", type = "continuous", n = length(unique(tab_1$Order))) # set colour palette for taxonomy section

figure_1A <- ggdotchart(tab_1, x = "Family", y = "n",
                     color = "Order", # Custom color palette
                     sorting = "descending", 
                     palette = zissou_colors,
                     ylab = "Number of studies", # Sort value in descending order
                     rotate = TRUE, 
                     dot.size = 7,
                     label = round(tab_1$n,1), 
                     font.label = list(color = "black", size = 8, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_pubr()                        # ggplot2 theme
)+
  theme_cleveland()  

figure_1A <- figure_1A + theme(legend.position = "bottom", 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 10), 
                         axis.title = element_text(size = 10), 
                         axis.title.y = element_text(angle= 90, margin = margin(0,15,0,0)), 
                         plot.margin = unit(c(1,1,1,1), "cm"), 
                         panel.background = element_rect(colour = "black", linewidth = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 2))

print(figure_1A)

### Category Plot ###

tab_2 <- data_2 %>%
  group_by(Order, Family, Focus)%>%
  count(Family, sort = T)%>%
  print(n = 50)

tab_2

tab_2$Focus <- recode_factor(tab_2$Focus, "Popgen" = "Population Genetics", "Social" = "Social Behaviour", "Reproduction" = "Reproductive Behaviour")

tab_2$Focus <- factor(tab_2$Focus, levels = c("Reproductive Behaviour", "Population Genetics", "Demography", "Social Behaviour"))

figure_1B <- ggdotchart(tab_2, x = "Family", y = "n",
                     color = "Order",
                     group = "Order",
                     facet.by = "Focus", 
                     palette = zissou_colors,
                     sorting = "descending",
                     ylab = "Number of studies",
                     rotate = T,
                     dot.size = 5,
                     label = round(tab_2$n,1), 
                     font.label = list(color = "black", size = 8, 
                                       vjust = 0.5), # Color y text by groups
                     ggtheme = theme_bw()                   
) + 
  theme_cleveland()


figure_1B <- figure_1B + theme(legend.position = "bottom", 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.text = element_text(size = 10), 
                         axis.title = element_text(size = 10), 
                         axis.title.y = element_text(margin = margin(0,10,0,0), angle = 90), 
                         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                         strip.text = element_text(size = 11),
                         plot.margin = unit(c(1,1,1,1), "cm")) +
  guides(colour = F)

print(figure_1B)

## Merge plot 2 and 3 

taxa_all <- ggarrange(figure_1A, figure_1B,
                      labels = c("A", "B"),
                      legend = "bottom",
                      common.legend = T,
                      ncol = 2, nrow = 1)

taxa_all

emf("C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/taxa_all_studies.emf", width = 14, height = 11)  # Set the width and height in inches
print(taxa_all)
dev.off()

#### Figure 2 - Markers ####

data2 <- data_2 %>%
  mutate(new_bin = cut(Year, breaks = seq(2000, 2025, by = 5),  # Adjust to your bin size
                       labels = c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025")))

tab_4 <- data2 %>%
  group_by(Markers, new_bin, Focus) %>%
  count(Markers, sort = T)


year_order <- c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025")

tab_4

tab_4$Focus <- recode_factor(tab_4$Focus, "Popgen" = "Population Genetics", "Social" = "Social Behaviour", "Reproduction" = "Reproductive Behaviour")

tab_4$Focus <- factor(tab_4$Focus, levels = c("Reproductive Behaviour", "Population Genetics", "Demography", "Social Behaviour"))

tab_4$new_bin <- factor(tab_4$new_bin, levels = c("2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025"))


tab_4$Markers <- factor(tab_4$Markers, levels = c("RFLP", "AFLP", "mSat", "SNP"))

tab_4

tab_4$new_bin <- factor(tab_4$new_bin, levels = unique(tab_4$new_bin))

figure_2 <- ggplot(tab_4, aes(x = new_bin, y = n, fill = Markers)) +
  geom_bar(color = "black", stat = "identity", position = position_dodge2(preserve = "single"), width = 1) +  # Use stat="identity" to plot pre-summarized data
  facet_wrap(~ Focus) +  # Facet by Focus
  scale_fill_brewer(palette = "Paired") +  
  xlab("Year") + # Set color palette
  labs(y = "Number of studies") +  # Y-axis label
  theme_bw() +
  scale_x_discrete(limits = year_order) + # Base theme
  theme(
    legend.position = "bottom", 
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 11), 
    axis.title.y = element_text(margin = margin(0, 15, 0, 0)), 
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    strip.text = element_text(size = 11)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3), nrow = 1))


# Print the plot
print(figure_2)

emf("C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots/figure_2.emf", width = 12, height = 13)  # Set the width and height in inches
print(figure_2)
dev.off()

### Figure 3 - Estimators ###

data_2 <- read.csv("relatedness_literature_review_working3.csv", stringsAsFactors = T)

data_2 <- data_2 %>% 
  filter(Index == "Keep")%>%
  droplevels

## Fixing estimator names

data_2$Kinship.Method <- recode_factor(data_2$Kinship.Method, "Allele counts  " = "Allele counts", "GERUD & COLONY " = "GERUD + COLONY", "COLONY + CERVUS" = "COLONY + CERVUS", "Allele counts & GERUD 1" = "Allele counts + GERUD", "Allele counts & GERUD 2.0" = "Allele counts + GERUD", "GERUD 2.0" = "GERUD", "GERUD 2.0 & COLONY" = "GERUD + COLONY", "Allele counts, GERUD 2.0 & COLONY2" = "Allele counts + GERUD 2.0 + COLONY", "IR Values & KINSHIP 1.3" = "Kinship 1.3", "Kinship 1.3 + Cervus 2.0 " = "Kinship 1.3 + Cervus", "GERUD 1.0, COLONY, STORM & allele counts" = "Allele counts + COLONY + GERUD + STORM", "GERUD 2.0, CERVUS 3.0.7 and COLONY" = "CERVUS + COLONY + GERUD", "Allele counts, GERUD 2.0, COLONY" = "Allele counts + COLONY + GERUD", "GERUD & COLONY" = "COLONY + GERUD", "Allele counts, GERUD 2.0 & COLONY" = "Allele counts + COLONY + GERUD", "Sequoia, COLONY, dartR*" = "COLONY + dartR + Sequoia", "KINFERENCE" = "Kinference", "COLONY2" = "COLONY", "COLONY v.2." = "COLONY", "MLRELATE, COLONY v1.2, KINGROUP 1" = "ML-Relate + COLONY + KinGroup", "KinGroup  " = "KinGroup", "Allele Counts + COLONY + GERUD" = "Allele counts + COLONY + GERUD", "Coancestry + COLONY+ CERVUS " = "Coancestry + COLONY + CERVUS", "COLONY & CERVUS" = "COLONY + CERVUS", "Kinship 1.3 & Cervus" = "Kinship 1.3 + Cervus", "Allele Counts" = "Allele counts", "Allele Counts + ML-Relate" = "Allele counts + ML-Relate")

estimator2 <- data_2 %>% 
  group_by(Kinship.Method, Focus)%>%
  count(Kinship.Method, sort= T)

estimator2$Kinship.Method <- as.character(estimator2$Kinship.Method)

expanded_data <- estimator2 %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n, Focus)

expanded_data$Kinship.Method <- recode_factor(expanded_data$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS", "VCFtools (Yang)" = "VFCtools", "CKMRsim " = "CKMRsim", "Sequioa " = "Sequoia", "GENAlEX" = "GenAlEx", "KING" = "KING-Robust", "Allele counts  " = "Allele counts", "VFCtools" = "VCFtools", "KinGroup  " = "KinGroup")

expanded_data$Kinship.Method <- as.character(expanded_data$Kinship.Method)
expanded_data$Focus <- as.character(expanded_data$Focus)

expanded_data <- expanded_data %>%
  group_by(Kinship.Method, Focus) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

print(expanded_data, n = 100)

expanded_data$Kinship.Method <- recode_factor(expanded_data$Kinship.Method, "VFCtools" = "VCFtools", "Coancestry" = "COANCESTRY")

expanded_data$Focus <- recode_factor(expanded_data$Focus, "Popgen" = "Population Genetics", "Social" = "Social Behaviour", "Reproduction" = "Reproductive Behaviour")

expanded_data$Focus <- factor(expanded_data$Focus, levels = c("Reproductive Behaviour", "Population Genetics", "Demography", "Social Behaviour"))

## Make Cat_Con_frame: 

Cat_Con_frame <- data.frame(Kinship.Method = c("Allele counts", "CERVUS", "CKMRsim", "COANCESTRY", "COLONY", "demerelate", "Endelman Jannik", "FAMOZ", "GenAlEx", "GERUD", "KINANALYZER", "Kinference", "KING-Robust", "KinGroup", "Kinship 1.3", "ML-Relate", "Sequoia", "SPAGEDI", "VCFtools"), 
                            Cat_Con = c("Cat", "Cat", "Cat", "Con", "Cat", "Con", "Con", "Cat", "Con", "Cat", "Cat", "Cat", "Con", "Cat", "Con", "Both", "Cat", "Con", "Con"))


estimator <- expanded_data %>%
  left_join(Cat_Con_frame, by = "Kinship.Method")


estimator <- estimator %>%
  mutate(Cat_Con = ifelse(is.na(Cat_Con), "Con", Cat_Con))

palette_sa <- c ("orange","skyblue", "grey")

## Use data not data_2 
# For ascending order

print(estimator, n = 40)

estimator$Kinship.Method <- factor(estimator$Kinship.Method, levels=unique(estimator$Kinship.Method))

estimator$Total_Frequency <- as.numeric(as.character(estimator$Total_Frequency))

estimator$Kinship.Method <- factor(estimator$Kinship.Method, 
                                   levels = rev(sort(unique(estimator$Kinship.Method))))

# Now plot
estimatorplot_3B <- ggplot(estimator, aes(x = Kinship.Method, y = Total_Frequency, color = Cat_Con)) +
  geom_point(aes(size = 4), alpha = 0.8) +
  coord_flip() +# Use points as a dot chart
  geom_text(aes(label = Total_Frequency), color = "black", size = 3, vjust = 0.5) +  # Add numeric labels
  facet_wrap(~ Focus, scales = "fixed") +  # Facet by study focus
  scale_color_manual(values = palette_sa) +  # Apply color palette
  scale_size(range = c(3, 8)) +  # Adjust point size range
  labs(x = "Number of Studies", y = "Relatedness Estimator") +
  theme_cleveland() +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +# Clean theme
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    strip.background = element_rect(color = "black", linewidth = 0.5)
  )

# Print the plot
print(estimatorplot_3B)

#### Figure 3 - Estimators ####

## For all focuses 

estimator <- data %>% 
  group_by(Kinship.Method)%>%
  count(Kinship.Method, sort= T)

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

estimator$Kinship.Method <- recode_factor(estimator$Kinship.Method, "CERVUS " = "CERVUS", "Allele counts  " = "Allele counts", "CKMRsim " = "CKMRsim", "KinGroup  " = "KinGroup", "Sequioa " = "Sequoia", "Cervus" = "CERVUS", "VCFtools (Yang)" = "VCFtools", "CERVUS " = "CERVUS")

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

Cat_Con <- c("Cat", "Cat", "Cat", "Con", "Both", "Cat","Cat", "Con", "Cat", "Con", "Con", "Cat", "Cat", "Cat", "Con", "Cat", "Cat", "Con", "Con", "Con")

Cat_Con_frame <- data.frame(Kinship.Method, Cat_Con) ## this is our referemce

glimpse(Cat_Con_frame)

estimator <- estimator %>%
  left_join(Cat_Con_frame, by = "Kinship.Method")

# Print it :)

print(estimator, n = 100)

estimator$Kinship.Method <- factor(estimator$Kinship.Method, 
                                   levels = rev(sort(unique(estimator$Kinship.Method))))

estimatorplot_3A <- ggplot(estimator, aes(x = Kinship.Method, y = Total_Frequency, color = Cat_Con)) +
  geom_point(aes(size = 4), alpha = 0.8) +
  coord_flip() +# Use points as a dot chart
  geom_text(aes(label = Total_Frequency), color = "black", size = 3, vjust = 0.5) +  # Add numeric labels
  scale_color_manual(values = palette_sa) +  # Apply color palette
  scale_size(range = c(3, 8)) +  # Adjust point size range
  labs(x = "Number of Studies", y = "Relatedness Estimator") +
  theme_cleveland() +
  scale_y_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50)) +# Clean theme
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    strip.background = element_rect(color = "black", linewidth = 0.5)
  )

# Print the plot
print(estimatorplot_3A)

### Figure 3B ###

data_2 <- read.csv("relatedness_literature_review_working3.csv", stringsAsFactors = T)

data_2 <- data_2 %>% 
  filter(Index == "Keep")%>%
  droplevels

## Fixing estimator names

data_2$Kinship.Method <- recode_factor(data_2$Kinship.Method, "Allele counts  " = "Allele counts", "GERUD & COLONY " = "GERUD + COLONY", "COLONY + CERVUS" = "COLONY + CERVUS", "Allele counts & GERUD 1" = "Allele counts + GERUD", "Allele counts & GERUD 2.0" = "Allele counts + GERUD", "GERUD 2.0" = "GERUD", "GERUD 2.0 & COLONY" = "GERUD + COLONY", "Allele counts, GERUD 2.0 & COLONY2" = "Allele counts + GERUD 2.0 + COLONY", "IR Values & KINSHIP 1.3" = "Kinship 1.3", "Kinship 1.3 + Cervus 2.0 " = "Kinship 1.3 + Cervus", "GERUD 1.0, COLONY, STORM & allele counts" = "Allele counts + COLONY + GERUD + STORM", "GERUD 2.0, CERVUS 3.0.7 and COLONY" = "CERVUS + COLONY + GERUD", "Allele counts, GERUD 2.0, COLONY" = "Allele counts + COLONY + GERUD", "GERUD & COLONY" = "COLONY + GERUD", "Allele counts, GERUD 2.0 & COLONY" = "Allele counts + COLONY + GERUD", "Sequoia, COLONY, dartR*" = "COLONY + dartR + Sequoia", "KINFERENCE" = "Kinference", "COLONY2" = "COLONY", "COLONY v.2." = "COLONY", "MLRELATE, COLONY v1.2, KINGROUP 1" = "ML-Relate + COLONY + KinGroup", "KinGroup  " = "KinGroup", "Allele Counts + COLONY + GERUD" = "Allele counts + COLONY + GERUD", "Coancestry + COLONY+ CERVUS " = "Coancestry + COLONY + CERVUS", "COLONY & CERVUS" = "COLONY + CERVUS", "Kinship 1.3 & Cervus" = "Kinship 1.3 + Cervus", "Allele Counts" = "Allele counts", "Allele Counts + ML-Relate" = "Allele counts + ML-Relate")

estimator2 <- data_2 %>% 
  group_by(Kinship.Method, Focus)%>%
  count(Kinship.Method, sort= T)

estimator2$Kinship.Method <- as.character(estimator2$Kinship.Method)

expanded_data <- estimator2 %>%
  mutate(Kinship.Method = strsplit(Kinship.Method, " \\+ ")) %>% # Split the combinations
  unnest(Kinship.Method) %>%
  select(Kinship.Method, n, Focus)

expanded_data$Kinship.Method <- recode_factor(expanded_data$Kinship.Method, "Cervus" = "CERVUS", "CERVUS " = "CERVUS", "VCFtools (Yang)" = "VFCtools", "CKMRsim " = "CKMRsim", "Sequioa " = "Sequoia", "GENAlEX" = "GenAlEx", "KING" = "KING-Robust", "Allele counts  " = "Allele counts", "VFCtools" = "VCFtools", "KinGroup  " = "KinGroup")

expanded_data$Kinship.Method <- as.character(expanded_data$Kinship.Method)
expanded_data$Focus <- as.character(expanded_data$Focus)

expanded_data <- expanded_data %>%
  group_by(Kinship.Method, Focus) %>%
  summarise(Total_Frequency = sum(n, na.rm = TRUE), .groups = "drop")

print(expanded_data, n = 100)

expanded_data$Kinship.Method <- recode_factor(expanded_data$Kinship.Method, "VFCtools" = "VCFtools", "Coancestry" = "COANCESTRY")

expanded_data$Focus <- recode_factor(expanded_data$Focus, "Popgen" = "Population Genetics", "Social" = "Social Behaviour", "Reproduction" = "Reproductive Behaviour")

expanded_data$Focus <- factor(expanded_data$Focus, levels = c("Reproductive Behaviour", "Population Genetics", "Demography", "Social Behaviour"))

## Make Cat_Con_frame: 

Cat_Con_frame <- data.frame(Kinship.Method = c("Allele counts", "CERVUS", "CKMRsim", "COANCESTRY", "COLONY", "demerelate", "Endelman Jannik", "FAMOZ", "GenAlEx", "GERUD", "KINANALYZER", "Kinference", "KING-Robust", "KinGroup", "Kinship 1.3", "ML-Relate", "Sequoia", "SPAGEDI", "VCFtools"), 
                            Cat_Con = c("Cat", "Cat", "Cat", "Con", "Cat", "Con", "Con", "Cat", "Con", "Cat", "Cat", "Cat", "Con", "Cat", "Con", "Both", "Cat", "Con", "Con"))


estimator <- expanded_data %>%
  left_join(Cat_Con_frame, by = "Kinship.Method")


estimator <- estimator %>%
  mutate(Cat_Con = ifelse(is.na(Cat_Con), "Con", Cat_Con))

palette_sa <- c ("orange","skyblue", "grey")

## Use data not data_2 
# For ascending order

print(estimator, n = 40)

estimator$Kinship.Method <- factor(estimator$Kinship.Method, levels=unique(estimator$Kinship.Method))

estimator$Total_Frequency <- as.numeric(as.character(estimator$Total_Frequency))

estimator$Kinship.Method <- factor(estimator$Kinship.Method, 
                                   levels = rev(sort(unique(estimator$Kinship.Method))))

# Now plot
estimatorplot_3B <- ggplot(estimator, aes(x = Kinship.Method, y = Total_Frequency, color = Cat_Con)) +
  geom_point(aes(size = 4), alpha = 0.8) +
  coord_flip() +# Use points as a dot chart
  geom_text(aes(label = Total_Frequency), color = "black", size = 3, vjust = 0.5) +  # Add numeric labels
  facet_wrap(~ Focus, scales = "fixed") +  # Facet by study focus
  scale_color_manual(values = palette_sa) +  # Apply color palette
  scale_size(range = c(3, 8)) +  # Adjust point size range
  labs(x = "Number of Studies", y = "Relatedness Estimator") +
  theme_cleveland() +
  scale_y_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +# Clean theme
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), "cm"), 
    strip.background = element_rect(color = "black", linewidth = 0.5)
  )

# Print the plot
print(estimatorplot_3B)


### Figure 3 AB ### 
figure_3AB <- ggarrange(estimatorplot_3A, estimatorplot_3B,
                        labels = c("A", "B"),
                        common.legend = F,
                        ncol = 2, nrow = 1)

emf("C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/New_Plots/estimatorplot_3AB.emf", width = 12, height = 10)  # Set the width and height in inches
print(figure_3AB)
dev.off()

#### Figure 4 #### 

focus_cum <-  data_2 %>%
  group_by(Focus) %>%
  arrange(Year) %>%
  mutate(cumulative_count = row_number())

focus_cum$Focus <- recode_factor(focus_cum$Focus, "Popgen" = "Population Genetics", "Social" = "Sociality")

focus_cum$Focus <- factor(focus_cum$Focus, levels = c("Reproduction", "Population Genetics", "Demography", "Sociality"))

focus_cum

cumplot_1 <- ggplot(data = focus_cum, aes(x = Year, y = cumulative_count, color = Focus, group = Focus)) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Spectral") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Number of studies",
    color = "Focus"
  ) +
  scale_y_continuous(breaks = seq(0, max(focus_cum$cumulative_count), by = 10))+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  geom_vline(xintercept = 2015, linetype = "dashed", size = 1) +
  annotate("text", x = 2016, y = 50, 
           label = "Genomics Era", angle = 90, vjust = -0.5, size = 5)

print(cumplot_1)

emf("C:/Users/samue/Desktop/Honours/Chapter_1_lit_review/Figures_for_pub/figure_4.emf", width = 10, height = 8)  # Set the width and height in inches
print(cumplot_1)
dev.off()


