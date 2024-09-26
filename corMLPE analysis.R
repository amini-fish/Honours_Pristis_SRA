# some codes for corMLPE, MK 23 July 2024
#---------------------------------------------

# https://github.com/nspope/corMLPE 
 devtools::install_github("nspope/corMLPE")

install.packages("nlme")
install.packages("bbmle")
install.packages("ggeffects")
install.packages("usdm")
install.packages("paletteer")
install.packages("performance")
install.packages("emmeans")
install.packages("MuMIn")
install.packages("rcompanion")
install.packages("ggsignif")

## --- START HERE PEOPLES --- ##

setwd("C:/Users/samue/Desktop/Honours/analysis")

library(nlme)
library(corMLPE)
library(bbmle)
library(ggeffects)
library(tidyverse)
library(paletteer)
library(performance)
library(emmeans)
library(lme4)
library(MuMIn)
library(rcompanion)
library(ggsignif)

## --- NOW WE HAVE EVERYTHING LOADED --- ##

data <- read.csv("corMLPE_data.csv", stringsAsFactors=TRUE)
names(data)
summary(data)

data

data2 <- data %>% rename(relatedness = value,
                geogr_dist = dist, 
                ID1 = id_1, 
                ID2 = id)

head(data2)

 nlevels(data2$ID1)
 summary(data2$genet_dist)
 sd(data2$genet_dist)
 
 nrow(data2)
 levels(data2$Billabong_ID1)
 
 levels(data2$billabong_both)
 
 data2$billabong_both <- recode_factor(data2$billabong_both, "Bend:Bank" = "River:Bend", "DalyRiv_Rd:Bend" = "River:Bend", "River:DalyRiv_Rd" = "River:River", "DalyRiv_Rd:Fish" = "River:Fish", "DalyRiv_Rd:Bank" = "River:River", "River:Bank" = "River:River", "Bank:Fish" = "River:Fish")
 
levels(data2$billabong_both)

## Remove the comparisons between individuals with no date of capture...

data2 <- data2[-c(1:55),] ## only use for some models 

ggplot(data = data2, aes(x = geogr_dist, y = relatedness, colour = billabong_diff, shape = year_caught_diff)) +
  geom_jitter(stat = "identity", size = 3, width = 0.3)+
  theme_classic()+
  xlab("Geogrpahic dist (km)")

ggplot(data = data2, aes(x = year_caught_diff, y = relatedness)) +
  geom_boxplot() +
  geom_jitter(aes(colour = geogr_dist), width = 0.25, size = 3)+
  theme_classic()

ggplot(data = data2, aes(x = billabong_diff, y = relatedness)) +
  geom_boxplot() +
  geom_jitter( size = 3)

ggplot(data = data2, aes(x = catch_set_diff, y = relatedness)) +
  geom_boxplot() +
  geom_jitter( size = 3, aes(colour = billabong_diff), width = 0.2, alpha = 0.7)+
  theme_classic()

## doesnt look like there's much going on here...

# do some exploratory plots first ..----
#... i.e. scatter plots e.g. to explore trends btw genetic & geogr.dist - is it linear?

# also check out distribution of genetic distance - approx. normal or strongly skewed?
## if latter, might need to transform first (e.g. log if right skewed)...

hist(data2$relatedness)
hist(log(data2$relatedness))

data <- data2[-c(1:55),]

data <- transform(na.omit(data))

summary(data)

hist(log(data$relatedness))

data2$log_GenD <- log(data2$relatedness)

ggplot(data = data2, aes(x = year_caught_diff, y = log_GenD)) +
  geom_boxplot() +
  geom_jitter(aes(colour = geogr_dist), size = 3)


hist(scale(data2$geogr_dist))
hist(data2$geogr_dist)

# corMLPE models --------------


## log()
m1 <- gls(log(relatedness + 0.001) ~ year_caught_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          na.action = "na.omit", 
          control = list(singular.ok = T))


m2 <- gls(log(relatedness+0.001) ~ scale(geogr_dist), 
                data = data2, 
                correlation = corMLPE(form = ~ID1 + ID2), 
                method = "REML", 
          control = list(singular.ok = T))

m3 <- gls(log(relatedness+0.001) ~ scale(geogr_dist) + year_caught_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ ID1 + ID2), 
          method = "REML", 
          control = list(singular.ok = T))

m4 <- gls(log(relatedness+0.001) ~ billabong_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          na.action = "na.omit", 
          control = list(singular.ok = T))

m5 <- gls(log(relatedness+0.001) ~ catch_set_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          na.action = "na.omit", 
          control = list(singular.ok = T))

m6 <- gls(log(relatedness+0.001) ~ billabong_diff*year_caught_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          control = list(singular.ok = T))

m7 <- gls(log(relatedness+0.001) ~ billabong_both, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          control = list(singular.ok = T))

m8 <- gls(log(relatedness+0.001) ~ year_caught_both, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          control = list(singular.ok = T))

# include year e.g. as random effect if many levels -...
#.... change model to lme() i.e. mixed effect model

tab <- AICctab(m1, m2, m3, m4, m5, m6, m7, m8, weights=TRUE)

tab

write.csv(tab, "AICcTable.csv")

# model summary + plot ---------------

## In order of best fit 

dev.off()

#1
summary(m7) 
anova(m7)
confint(m7)
plot(m7)

plot(m4, resid(., type="normalized") ~ fitted(.), abline = c(0,0))

ggplot(data2, aes(y=resid(m4, type="normalized"), 
                  x= billabong_diff)) + 
  geom_boxplot()

nagelkerke(m4) ## Note that the mdoel was refitted with ML not REML 

## include but don't rely on

# 95% CI's are back-transformed i.e. can ignore warning about SE's being on the log scale


ggpredict(model = m4, terms = c("billabong_diff"), back_transform = TRUE) %>%  
  plot(show_data = TRUE, jitter = 0.07)


## --- Plot for aesthetics --- ## 

plot1 <- ggplot(data = data2, aes(x = catch_set_diff, y = genet_dist, fill = catch_set_diff)) +
  geom_violin(alpha = 0.8) +
  geom_jitter(size = 3, width = 0.3, aes(shape = year_caught_diff, col = geogr_dist)) +
  ggtitle("Comparison of relatedness values between pairs in different vs the same catch set") +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous("Relatedness", breaks = seq(0, 0.3, by = 0.025))+
  scale_colour_viridis_c(option = "magma")+
  theme_classic()

plot1 <- plot1 + theme(plot.title = element_text(hjust = 0.5, size = 15), 
                       axis.title = element_text(size = 12), 
                       axis.text = element_text(size = 11))

plot1

###-------------------------------------------------------------------------

## Clean the palette 
dev.off()

## Plot 2

plot2 <- ggplot(data = data2, aes(x = reorder(billabong_diff, -log(relatedness + 0.001)), y = log(relatedness + 0.001), fill = billabong_both))+
  geom_boxplot()+
  geom_signif(data=data2,comparisons=list( c("River:Fish", "Bend:Fish"),c("River:Fish", "Fish:Fish"), c("River:River","River:Fish")),
              map_signif_level = TRUE, 
              y_position = c(-1.3, -1.0, -0.70, -0.40), 
              annotations = c("*", "***", "***"), 
              textsize = 4.5) +
  ggtitle("Relatedness Between Billabongs") +
  xlab("Pairwise Billabong Comparisons") +
  ylab("Log ( Relatedness + 0.001 )") +
  theme_bw()

plot2 <- plot2 + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="Billabong Comparison"))

print(plot2)

#geom_hline(yintercept = mean(log(data2$relatedness + 0.001)), ) +

## Post hoc comparisons --

emm <- emmeans(m7, ~ billabong_both)
simp <- pairs(emm, simple = "each")
simp

install.packages("multcompView")
library(multcompView)
tukey <- TukeyHSD(aov(mod), "billabong_both")
plot(tukey , las=1 , col="brown", las = 2)

#model effect jtools 

