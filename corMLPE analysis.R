# some codes for corMLPE, MK 23 July 2024
#---------------------------------------------

# https://github.com/nspope/corMLPE 
 devtools::install_github("nspope/corMLPE")

install.packages("nlme")
install.packages("bbmle")
install.packages("ggeffects")
install.packages("usdm")
install.packages("paletteer")

## --- START HERE PEOPLES --- ##

setwd("C:/Users/samue/Desktop/Honours/analysis")

library(nlme)
library(corMLPE)
library(bbmle)
library(ggeffects)
library(tidyverse)
library(paletteer)

## --- NOW WE HAVE EVERYTHING LOADED --- ##

data <- read.csv("corMLPE_data.csv", stringsAsFactors=TRUE)
names(data)
summary(data)

data

data2 <- data %>% rename(genet_dist = value,
                geogr_dist = dist, 
                ID1 = id_1, 
                ID2 = id)

head(data2)
## Remove the comparisons between individuals with no date of capture...

data2 <- data2[-c(1:55),]

ggplot(data = data2, aes(x = geogr_dist, y = genet_dist, colour = billabong_diff, shape = year_caught_diff)) +
  geom_jitter(stat = "identity", size = 3, width = 0.3)+
  theme_classic()+
  xlab("Geogrpahic dist (km)")

ggplot(data = data2, aes(x = year_caught_diff, y = genet_dist)) +
  geom_boxplot() +
  geom_jitter(aes(colour = geogr_dist), width = 0.25, size = 3)+
  theme_classic()

ggplot(data = data2, aes(x = billabong_diff, y = genet_dist)) +
  geom_boxplot() +
  geom_jitter( size = 3)

ggplot(data = data2, aes(x = catch_set_diff, y = genet_dist)) +
  geom_boxplot() +
  geom_jitter( size = 3, aes(colour = billabong_diff), width = 0.2, alpha = 0.7)+
  theme_classic()

## doesnt look like there's much going on here...

# do some exploratory plots first ..----
#... i.e. scatter plots e.g. to explore trends btw genetic & geogr.dist - is it linear?

# also check out distribution of genetic distance - approx. normal or strongly skewed?
## if latter, might need to transform first (e.g. log if right skewed)...

hist(data2$genet_dist)
hist(log(data2$genet_dist))

data <- data2[-c(1:55),]

data <- transform(na.omit(data))

summary(data)

hist(log(data$genet_dist))

data2$log_GenD <- log(data2$genet_dist)

ggplot(data = data2, aes(x = year_caught_diff, y = log_GenD)) +
  geom_boxplot() +
  geom_jitter(aes(colour = geogr_dist), size = 3)


hist(scale(data2$geogr_dist))
hist(data2$geogr_dist)

# corMLPE models --------------


## log()
m1 <- gls(genet_dist ~ year_caught_diff + scale(geogr_dist), 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          na.action = "na.omit", 
          control = list(singular.ok = T))


# use categorical variable instead 
m2 <- gls(genet_dist ~ scale(geogr_dist), 
                data = data2, 
                correlation = corMLPE(form = ~ID1 + ID2), 
                method = "REML", 
          control = list(singular.ok = T))

m3 <- gls(log(genet_dist) ~ year_caught_diff + scale(geogr_dist), 
          data = data2, 
          correlation = corMLPE(form = ~ ID1 + ID2), 
          method = "REML", 
          control = list(singular.ok = T))

m4 <- gls(genet_dist ~ billabong_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          na.action = "na.omit", 
          control = list(singular.ok = T))

m5 <- gls(genet_dist ~ catch_set_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          na.action = "na.omit", 
          control = list(singular.ok = T))

m6 <- gls(genet_dist ~ billabong_diff*year_caught_diff, 
          data = data2, 
          correlation = corMLPE(form = ~ID1 + ID2), 
          method = "REML", 
          control = list(singular.ok = T))

# include year e.g. as random effect if many levels -...
#.... change model to lme() i.e. mixed effect model

#lme1 <- lme(genet_dist ~ catch_set_diff,
          data = data2, 
          random = ~1|#WHAT TO PUT HERE??, 
          corMLPE(form = ~ID1 + ID2|year_caught_diff), 
          method = "REML")

# compare model fit with AICc --------------

tab <- AICctab(m1, m2, m3, m4, m5, m6,  weights=TRUE)

tab

write.csv(tab, "AICcTable.csv")

# model summary + plot ---------------
summary(m5)
anova(m5)

summary(m4)
anova(m4)

# 95% CI's are back-transformed i.e. can ignore warning about SE's being on the log scale

ggpredict(model = m4, terms = c("billabong_diff")) %>%  
  plot(show_data = TRUE, jitter = 0.07)


ggpredict(model = m5, terms = c("catch_set_diff")) %>%  
  plot(show_data = TRUE, jitter = 0.07)

# check model residuals - no patterns across fitted values nor predictors -----

plot(m5, resid(., type="normalized") ~ fitted(.), abline = c(0,0))

ggplot(data2, aes(y=resid(m5, type="normalized"), 
                 x= catch_set_diff)) + 
  geom_boxplot()

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

