#### Install Packages ####

install.packages("ggpubr")
install.packages("caret")
install.packages("interactions")
install.packages("jtools")
install.packages("lsmeans")
install.packages("GGally")
install.packages("DHARMa")
install.packages("pwr")
install.packages("WebPower")
install.packages("simr")
install.packages("splines")
install.packages("brglm2")
install.packages("brms")
install.packages("mclust")
install.packages("devEMF")

#### Load Packages ####

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stats)
library(jtools)
library(bbmle)
library(car)
library(interactions)
library(lsmeans)
library(ggpubr)
library(dplyr)
library(caret)
library(GGally)
library(DHARMa)
library(pwr)
library(WebPower)
library(simr)
library(splines)
library(brglm2)
library(brms)
library(mclust)
library(MASS)     
library(ggeffects)
library(ggplot2)
library(devEMF)
library(MuMIn)
library(jtools)


#### Load data ####

## Set Working directory 

setwd("C:/Users/samue/Desktop/Honours/Daly_ENV")

## Read in data w/ base predictors

pristis <- read.csv("seasonal_predictors_silo.csv", stringsAsFactors = T)

#### Load in additional independent variables - not needed anymore #### 

riverlevel <- read.csv("G8140003_riverlevel.csv")

riverlevel$date <- riverlevel$Timestamp..UTC.09.30.

# Function to estimate the average height of the river in each wet season
average_wet_season_height <- function(data, date_col, level_col) {
  # Convert date column to Date format
  data[[date_col]] <- as.Date(data[[date_col]], format = "%d/%m/%Y")
  
  # List to store results
  avg_height_per_season <- list()
  
  # Loop over the years from 2011 to 2023
  for (year in 2011:2023) {
    # Define start and end dates for the wet season
    start_date <- as.Date(paste0("01/10/", year), format="%d/%m/%Y")
    end_date <- as.Date(paste0("30/04/", year + 1), format="%d/%m/%Y")
    
    # Subset data for the wet season
    season_data <- data %>%
      filter(date >= start_date & date <= end_date)
    
    # Calculate the average river level in the wet season
    avg_height_per_season[[paste0("avg_height_", year)]] <- mean(season_data[[level_col]], na.rm = TRUE)
  }
  
  # Return the list of average heights for each wet season
  return(avg_height_per_season)
}

# Estimate average river height for each wet season
avg_heights <- average_wet_season_height(riverlevel, "date", "Value..m.")

print(avg_heights)

#### Days above X thresholds 

days_above_threshold <- function(data, date_col, level_col, threshold) {
  # Convert date column to Date format
  data[[date_col]] <- as.Date(data[[date_col]], format = "%d/%m/%Y")
  
  # List to store results
  days_above_X_per_season <- list()
  
  # Loop over the years from 2011 to 2023
  for (year in 2011:2023) {
    # Define start and end dates for the wet season
    start_date <- as.Date(paste0("01/10/", year), format="%d/%m/%Y")
    end_date <- as.Date(paste0("30/04/", year + 1), format="%d/%m/%Y")
    
    # Subset data for the wet season
    season_data <- data %>%
      filter(!!sym(date_col) >= start_date & !!sym(date_col) <= end_date) %>%
      filter(!!sym(level_col) > threshold) %>%
      distinct(!!sym(date_col))  # Keep only unique days
    
    # Count the number of unique days
    days_above_X_per_season[[paste0("days_above_", year)]] <- nrow(season_data)
  }
  
  # Return the list of days above threshold for each wet season
  return(days_above_X_per_season)
}

# Calculate the number of days when river level was greater than 12 for each wet season
mod_threshold <- 13.1 
minor_threshold <- 12.6
days_above_threshold <- days_above_threshold(riverlevel, "date", "Value..m.", minor_threshold)
print(days_above_threshold)


#### Exploratory Data Analysis ####

head(df)

eda <- df %>%
  select(ami_standardised, sum_wetseason_rainfall_mm, avg_daily_rainfall_mm, total_wetseason_discharge_ml, average_monthly_discharge_ml, total_days_minorflood, average_river_level)

ggpairs(eda)+
  theme_bw()

ggplot()+
  geom_histogram(data = df, aes(x = sum_wetseason_rainfall_mm),binwidth = 10)

ggplot()+
  geom_histogram(data = df, aes(x = avg_daily_rainfall_mm),binwidth = 0.5)

ggplot()+
  geom_histogram(data = df, aes(x = total_wetseason_discharge_ml),binwidth = 500000) # left skew + outlier

ggplot()+
  geom_histogram(data = df, aes(x = average_monthly_discharge_ml),binwidth = 200000) #left skew + outlier

ggplot()+
  geom_histogram(data = df, aes(x = ami_standardised),binwidth = 0.05)

ggplot()+
  geom_histogram(data = df, aes(x = average_river_level),binwidth = 0.25)

ggplot()+
  geom_histogram(data = df, aes(x = total_days_moderateflood),binwidth = 50) # Right skew

# Parameters
odds_ratio <- 2  # Expected effect size in terms of odds ratio
alpha <- 0.05    # Significance level
power <- 0.70     # Desired power
p <- 0.3        # Probability of event in control group
r2 <- 0.1        # R-squared of other predictors
family <- "Bernoulli"

# Compute sample size

wp.logistic(n = NULL, p0 = p, p1 = p *odds_ratio, alpha = alpha, power = power, family = family)

## we need n = 70 to achieve an odds ratio of 2 


#### Modelling Formulation ####

## Load the data back in - with changes made in excel

df <- read.csv("seasonal_predictors_silo.csv", stringsAsFactors = T)

glimpse(df)

## Total rainfall

m1 <- glm(pristis_presence ~ sum_wetseason_rainfall_mm, 
          data = df, 
          family = "binomial")

m2 <- glm(pristis_presence ~ avg_daily_rainfall_mm,
          data = df, 
          family = "binomial")

m3 <- glm(pristis_presence ~ total_wetseason_discharge_ml, 
          data = df, 
          family = "binomial")

m4 <- glm(pristis_presence ~ average_monthly_discharge_ml, 
          data = df, 
          family = "binomial")

m5 <- glm(pristis_presence ~ sum_wetseason_rainfall_mm, 
          data = df, 
          family = "binomial")

m6 <- glm(pristis_presence ~ ami_standardised, 
          data = df, 
          family = "binomial")

m7 <- glm(pristis_presence ~ average_river_level, 
          data = df, 
          family = "binomial")


m8 <- glm(pristis_presence ~ total_days_minorflood, 
          data = df, 
          family = "binomial")

m9 <- glm(pristis_presence ~ total_days_moderateflood, 
          data = df, 
          family = "binomial")

model_probit <- glm(pristis_presence ~ average_river_level, family = binomial(link = "probit"), data = df)

model_cloglog <- glm(pristis_presence ~ average_river_level, family = binomial(link = "cloglog"), data = df)

# Firth regression good for small n 

model_firth <- glm(pristis_presence ~ average_river_level, family = binomial(link = "cloglog"), data = df, method = "brglmFit") 

#### Model Validation ####

## AIC

AICctab(m1, m2, m3, m4,m5, m6, m7, m8, m9,model_probit, model_cloglog, model_firth, sort=T, base=T, weights=TRUE) # AMI is the best

summary(model_cloglog)

## Check model assumptions m7
cloglogg_resids <- simulateResiduals(model_cloglog)
plot(cloglogg_resids)

m7_darm <- simulateResiduals(m7)

plot(m7_darm)

plot(model_cloglog)

plotResiduals(cloglogg_resids, df$average_river_level)

testDispersion(mod_cont_resids)
testZeroInflation(mod_cont_resids)

ami_resid <- simulateResiduals(m6)
plot(ami_resid)

##

cooksD <- cooks.distance(model_cloglog)
plot(cooksD, type = "h")
abline(h = 4 / length(cooksD), col = "red") 

## Check model assumptions cloglog 
mod_cont_resids <- simulateResiduals(model_cloglog)
plot(mod_cont_resids)

testDispersion(mod_cont_resids)
testZeroInflation(mod_cont_resids)

confint(model_cloglog)

##### Model prediction ####

# Effect size 
predictions <- predict(model_cloglog, newdata = df, type = "link", se.fit = TRUE)

# Calculate the upper and lower bounds for the confidence intervals
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-score for 95% confidence

# Transform the linear predictors (log odds) to probabilities
df$predicted_probabilities <- plogis(predictions$fit)
df$lower_bound <- plogis(predictions$fit - z_value * predictions$se.fit)
df$upper_bound <- plogis(predictions$fit + z_value * predictions$se.fit)

#### ADD VALUES TO ANNOTATE ####

model_summary <- summary(model_cloglog)
coef_table <- model_summary$coefficients  # Extract coefficient table
df_residual <- model_summary$df.residual  # Residual degrees of freedom

# Extract values for annotation
intercept <- coef_table[1, 1]  # Intercept coefficient
slope <- coef_table[2, 1]  # Predictor coefficient
z_value <- coef_table[2, 3]  # z-score for predictor
p_value <- coef_table[2, 4]  # p-value for predictor
odds_ratio <- exp(slope)  # Convert coefficient to odds ratio

annotation_text <- paste0(
  "Intercept = ", round(intercept, 2), "\n",
  "Coefficient = ", round(slope, 2), "\n",
  "Odds Ratio = ", round(odds_ratio, 2), "\n",
  "Z-value = ", round(z_value, 2), "\n",
  "P-value = ", signif(p_value,2), "*\n",
  "Residual df = ", df_residual
)

#### Plot Model ####

plot_cloglog <- ggplot(df, aes(x = average_river_level, 
                               y = pristis_presence)) +
  geom_point(colour = "black", aes(fill = factor(pristis_presence)), 
             size = 3.5, shape = 21) + 
  geom_ribbon(data = df, 
              aes(x = average_river_level, 
                  ymin = lower_bound, 
                  ymax = upper_bound), 
              fill = "grey", 
              alpha = 0.3) + 
  geom_line(data = df, 
            aes(x = average_river_level,
                y = predicted_probabilities), 
            color = "black", size = 1) + 
  annotate("text", 
           x = 2.8, 
           y = 0.9, 
           label = annotation_text, 
           color = "black", 
           size = 3.5, 
           hjust = 0) +
  labs(
    x = "Average Wet Season River Level (m)",
    y = "Probability of Sawfish Presence"
  ) +
  theme_bw() +
  scale_fill_manual(values = c("grey34", "darkolivegreen3"))+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_text(margin = margin(t = 15)),  
        axis.title.y = element_text(margin = margin(r = 15)))


print(plot_cloglog)

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/resence.emf", width = 10, height = 8)  # Set the width and height in inches
print(plot_cloglog)
dev.off()

#### Abundance models ####

m1_abund <- glm(pristis_abund ~ sum_wetseason_rainfall_mm, 
          data = df, 
          family = "gaussian")

m2_abund <- glm(pristis_abund ~ avg_daily_rainfall_mm,
          data = df, 
          family = "gaussian")

m3_abund <- glm(pristis_abund ~ total_wetseason_discharge_ml, 
          data = df, 
          family = "gaussian")

m4_abund <- glm(pristis_abund ~ average_monthly_discharge_ml, 
          data = df, 
          family = "gaussian")

m5_abund <- glm(pristis_abund ~ sum_wetseason_rainfall_mm, 
          data = df, 
          family = "gaussian")

m6_abund <- glm(pristis_abund ~ ami_standardised, 
          data = df, 
          family = "gaussian")

m7_abund <- glm(pristis_abund ~ average_river_level, 
          data = df, 
          family = "gaussian")


m8_abund <- glm(pristis_abund ~ total_days_minorflood, 
          data = df, 
          family = "gaussian")

m9_abund <- glm(pristis_abund ~ total_days_moderateflood, 
          data = df, 
          family = "gaussian")

m10_ab <- lm(pristis_abund ~ ami_standardised, 
             data = df)

m6_pois <- glm(pristis_abund ~ ami_standardised, 
                data = df, 
                family = "poisson")

m6_negbin <- glm.nb(pristis_abund ~ ami_standardised,
                    df)

require(pscl)
require(MASS)
require(boot)

Z_inf <- zeroinfl(pristis_abund ~ ami_standardised,
                  data = df, 
                  dist = "negbin")

summary(m6_negbin)


summ(m6_negbin, exp = T, digits = 3)

m6_resid <- simulateResiduals(m6_negbin)
plot(m6_resid)

plotResiduals(m6_resid, df$ami_standardised)

ggplot(df, aes(x = pristis_abund)) +
  geom_histogram(binwidth = 10)


## Zero inflated works best 

confint(Z_inf)
summary(Z_inf)
### AICc

AICctab(m1_abund, m2_abund, m3_abund, m4_abund, m5_abund, m6_abund, m7_abund, m8_abund, m9_abund,m10_ab,m6_pois, m6_negbin,Z_inf, sort=T, base=T, weights=TRUE)

summary(m10_ab)
confint(m10_ab)

  # For visualization

#### Extract Values from Negative Binomial Model ####

# Create a new data frame for predictions
df$predicted_nb <- predict(m6_negbin, type = "response")  # Get predicted values

# Generate predicted values
preds <- ggpredict(m6_negbin, terms = "ami_standardised")

# Generate confidence intervals
se_fit <- predict(m6_negbin, type = "link", se.fit = TRUE)

# Compute confidence intervals (using log link function)
df$lower <- exp(se_fit$fit - 1.96 * se_fit$se.fit)
df$upper <- exp(se_fit$fit + 1.96 * se_fit$se.fit)

# Extract model summary
model_summary <- summary(m6_negbin)

# Get coefficients (Intercept and Slope)
intercept <- coef(m6_negbin)[1]
slope <- coef(m6_negbin)[2]

# Get p-value for the slope
p_value <- model_summary$coefficients[2, 4]

# Compute pseudo R²
pseudo_r2 <- performance::r2_nagelkerke(m6_negbin)
formatted_p_value <- ifelse(p_value < 0.001, "< 0.001", round(p_value, 3))

#### Plot the model ####

annotation_text <- paste(
  "Model: Negative Binomial Regression\n",
  "y = ", round(intercept, 2), " + ", round(slope, 2), " * x\n",
  "Pseudo R² = ", round(pseudo_r2, 3), "\n",
  "p-value = ", formatted_p_value
)

nb_plot <- ggplot(preds, aes(x = x, y = predicted)) +
  geom_point(data = df, aes(x = ami_standardised, y = pristis_abund), shape = 21, fill = "darkolivegreen3", size = 5, colour = "black") +  # Observed data
  geom_line(color = "black", size = 1) +  # Predicted line
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = "grey") +  # Confidence interval
  labs(
    x = "AMI Standardised",
    y = "Sawfish Abundance") +
  annotate("text", x = -0.2, 
           y = 650, 
           label = annotation_text, 
           hjust = 0, vjust = 1, size = , color = "black") +  
  theme(axis.title = element_text(size = 12)) +
  theme_bw()

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/neg_binom.emf", width = 10, height = 8)  # Set the width and height in inches
print(nb_plot)
dev.off()

# Plot residuals vs. fitted values
ggplot(df, aes(x = predicted, y = residuals(m6_negbin), type = "pearson")) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted values", y = "Pearson Residuals", title = "Residual Plot for NB Model") +
  theme_minimal()

#### Zero Inflated ####

new_data <- data.frame(ami_standardised = seq(min(df$ami_standardised), 
                                              max(df$ami_standardised), length.out = 100))

# Get predicted counts from the count model (non-zero)
new_data$predicted_count <- predict(Z_inf, new_data, type = "response")

ggplot(df, aes(x = ami_standardised, y = pristis_abund)) +
  geom_point(alpha = 0.5) +  # Raw data points
  geom_line(data = new_data, aes(y = predicted_count), color = "blue", size = 1) +
  labs(x = "AMI Standardised", y = "Predicted Pristis Abundance",
       title = "Predicted Pristis Abundance vs. AMI Standardised") +
  theme_bw()

# Get predicted probabilities of being a zero
new_data$predicted_zero <- predict(Z_inf, new_data, type = "zero")

ggplot(new_data, aes(x = ami_standardised, y = predicted_zero)) +
  geom_line(color = "red", size = 1) +
  labs(x = "AMI Standardised", y = "Probability of Excess Zeros",
       title = "Zero-Inflation Probability vs. AMI Standardised") +
  theme_minimal()

## 

new_data <- data.frame(ami_standardised = seq(min(df$ami_standardised), 
                                              max(df$ami_standardised), length.out = 100))

# Get model matrix for the count model
X_count <- model.matrix(~ ami_standardised, data = new_data)

# Extract coefficients for the count model
coef_count <- coef(Z_inf)[1:length(coef(Z_inf, model = "count"))]  # Ensure we only grab the count model coefficients

# Extract variance-covariance matrix for the count model
vcov_Z_inf <- vcov(Z_inf)  # Full variance-covariance matrix
vcov_count <- vcov_Z_inf[names(coef_count), names(coef_count), drop = FALSE]  # Select only count model elements

# Predict on the link scale (log-scale since log link is used)
log_fit <- X_count %*% coef_count
log_se <- sqrt(diag(X_count %*% vcov_count %*% t(X_count)))

# Convert to response scale
new_data$predicted_count <- exp(log_fit)
new_data$lower_CI <- exp(log_fit - 1.96 * log_se)
new_data$upper_CI <- exp(log_fit + 1.96 * log_se)

new_data

# Extract key model results
intercept <- round(coef(Z_inf)["(Intercept)"], 2)
slope <- round(coef(Z_inf)["ami_standardised"], 2)
theta <- round(Z_inf$theta, 2)
p_value_slope <- round(summary(Z_inf)$coefficients$count["ami_standardised", "Pr(>|z|)"], 3)
p_value_zero <- round(summary(Z_inf)$coefficients$zero["ami_standardised", "Pr(>|z|)"], 3)

p_value_slope
p_value_

coef(Z_inf)


z_inflate_plot <- ggplot() +
  geom_point(data = df, aes(x = ami_standardised, y = pristis_abund, fill = pristis_abund),shape = 21,colour = "black", size = 5) +  # Observed points
  geom_line(data = new_data, aes(x = ami_standardised, y = predicted_count), 
            color = "black", size = 1) +  # Model predictions
  geom_ribbon(data = new_data, aes(x = ami_standardised, ymin = lower_CI, ymax = upper_CI), 
              alpha = 0.2, fill = "grey") +  # Confidence interval
  labs(x = "AMI Standardised", y = "Pristis Abundance") +
  scale_fill_gradient(low = "darkgrey", high = "darkolivegreen2") + 
  labs(fill = "Pristis Presence") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 15), size = 15),  
        axis.title.y = element_text(margin = margin(r = 15)), size = 15)

z_inflate_plot <- z_inflate_plot + 
  annotate("text", x = -0.11, y = 69, 
           label = paste0("Intercept: 2.76"), hjust = 1, size = 4) +
  annotate("text", x = -0.11, y = 67, 
           label = paste0("Slope: 2.83 (p = ", p_value_slope, "*)"), hjust = 1, size = 4) +
  annotate("text", -0.11, y = 65, 
           label = paste0("Zero-inflation p: 0.06"), hjust = 1, size = 4) +
  annotate("text", -0.11, y = 63, 
           label = paste0("Theta: ", theta), hjust = 1, size = 4)

print(z_inflate_plot)

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/z_inf_abundance.emf", width = 10, height = 8)  # Set the width and height in inches
print(z_inflate_plot)
dev.off()

#### PLOT the data ####

box_plot <- ggplot(df, aes(x = as.factor(pristis_presence), y = average_river_level, fill = as.factor(pristis_presence))) +
  geom_boxplot() +
  theme_bw() +
  xlab("Sawfish Presence") +
  ylab("Average Wet Season River Level (m)") +
  scale_fill_manual(values = c("grey", "darkolivegreen3"))+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 15), size = 12),  
        axis.title.y = element_text(margin = margin(r = 15)), size = 12)


print(box_plot)

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/presence_boxplot.emf", width = 10, height = 8)  # Set the width and height in inches
print(box_plot)
dev.off()

##### Model prediction ####

# Extracting coefficients (estimates, standard errors, t-values, p-values)
summary_m6 <- summary(m6_abund)

#### Add Citations ####

citation(package = "mvabund")
citation(package = "jtools")
citation(package = "bbmle")
citation(package = "car")
citation(package = "interactions")
citation(package = "lsmeans")
citation(package = "ggpubr")
citation(package = "dplyr")
citation(package = "caret")
citation(package = "CKMRsim")


#### Visualising predictors #### 

glimpse(pristis)

wet_season_labs <- c("Season_2023_2024" = "2023/2024",
                     "Season_2022_2023" = "2023/2024", 
                     "Season_2021_2022" = "2021/2022", 
                     "Season_2020_2021" = "2020/2021", 
                     "Season_2019_2020" = "2019/2020", 
                     "Season_2018_2019" = "2018/2019",
                     "Season_2017_2018" = "2017/2018", 
                     "Season_2016_2017" = "2016/2017",
                     "Season_2015_2016" = "2015/2016",
                     "Season_2014_2015" = "2014/2015", 
                     "Season_2013_2014" = "2013/2014", 
                     "Season_2012_2013" = "2012/2013",
                     "Season_2011_2012" = "2011/2012")

pristis$pristis_presence <- as.factor(pristis$pristis_presence)

ami_bar <- ggplot(pristis, aes(x = Wet_Season, y = ami_standardised
                    , fill = pristis_presence)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("grey", "darkolivegreen3"), 
                    labels = c("0" = "No", "1" = "Yes")) +
  coord_flip()+
  theme_bw() +
  scale_x_discrete(labels = c(wet_season_labs)) +
  ylab("AMI Standardised") +
  xlab("Wet Season") +
  labs(fill = "Pristis Presence") +
  theme(legend.position = c(0.95, 0.2), 
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))

print(ami_bar)

emf("C:/Users/samue/Desktop/Honours/Daly_ENV/ami_barplot.emf", width = 10, height = 8)  # Set the width and height in inches
print(ami_bar)
dev.off()



#### Sawfish lenght & sex EDA ####

library(forcats)
library(dplyr)
library(ggnewscale)


setwd("C:/Users/samue/Desktop/Honours")

sawfish_all <- read.csv("Kyne_Daly River Pristis pristis data_SamAminiHonours.csv")

glimpse(sawfish_all)

sawfish_all <- sawfish_all %>%
  filter(Sawfish_rescue == "Yes")

mean(sawfish_all$TL_mm, na.rm = T)
sd(sawfish_all$TL_mm, na.rm = T)

sawfish_all$Sex <- factor(sawfish_all$Sex, levels = c("Unknown", "M", "F"))

table(sawfish_all$Sex)

Lengths <- ggplot(sawfish_all, aes(x = TL_mm, y = ..density..)) +
  geom_density(alpha = 0.9, fill = "grey", colour = "black") +
  geom_histogram(binwidth = 20, alpha = 0.7, fill = "darkolivegreen3", colour = "black") + 
  geom_vline(xintercept =  1121.229, linetype  = "dashed", linewidth = 1) +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  xlab("Total Length")+
  ylab("Density")

print(Lengths)

emf("C:/Users/samue/Desktop/Honours/Length_Frequency.emf", width = 10, height = 8)  # Set the width and height in inches
print(Lengths)
dev.off()

## Changes in sex proportion over time

Sex_data <- sawfish_all %>%
  group_by(Rescue_year, Sex) %>%
  count(Sex)

sawfish_all %>%
  summarise(
    male_count = sum(Sex == "M", na.rm = TRUE),
    female_count = sum(Sex == "F", na.rm = TRUE),
    sex_ratio = male_count / female_count
  )
  

#position_dodge2(width = 0.9, preserve = "single")

Sex <- ggplot(Sex_data) +
  geom_bar(aes(y = n, x = as.factor(Rescue_year), fill = Sex), stat = "identity", position = "fill", colour = "black") +
  scale_fill_manual(values = c("grey",   "darkolivegreen3", "orange")) +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  xlab("Rescue Year") + 
  ylab("Proportion")

print(Sex)

emf("C:/Users/samue/Desktop/Honours/Sex.emf", width = 10, height = 8)  # Set the width and height in inches
print(Sex)
dev.off()


## Hypothesis 3 - does wet season correlate to the total length of individuals? 

ggplot() +
  geom_boxplot(data = sawfish_all, aes(x = as.factor(Wet_Season), y = TL_mm, fill = Wet_Season)) +
  theme_bw()

ggplot() +
  geom_point(data = sawfish_all, aes(x = ami_pres$ami_standardised, y = TL_mm)) +
  theme_bw()

ggplot(data = sawfish_all, aes(x = as.factor(Rescue_year), y = TL_mm, fill = ami_standardised)) +
  geom_boxplot()
