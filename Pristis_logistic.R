#### Install Packages ####

install.packages("ggpubr")
install.packages("caret")
install.packages("interactions")
install.packages("jtools")
install.packages("lsmeans")
install.packages("GGally")

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

# Function to calculate the number of days when the river level was > X in each wet season
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
      filter(date >= start_date & date <= end_date)
    
    # Count the number of days when the river level is greater than threshold X
    days_above_X_per_season[[paste0("days_above_", year)]] <- sum(season_data[[level_col]] > threshold, na.rm = TRUE)
  }
  
  # Return the list of days above threshold for each wet season
  return(days_above_X_per_season)
}

# Estimate average river height for each wet season
avg_heights <- average_wet_season_height(riverlevel, "date", "Value..m.")

print(avg_heights)

# Calculate the number of days when river level was greater than 12 for each wet season
threshold <- 13.1
days_above_threshold <- days_above_threshold(riverlevel, "date", "Value..m.", threshold)
print(days_above_threshold)



#### Exploratory Data Analysis ####
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

#### Model Validation ####

## AIC

AICctab(m1, m2, m3, m4,m5, m6, m7, m8, m9, sort=T, base=T, weights=TRUE) # AMI is the best

summary(m7)

## Model prediction

# Predicting probabilities
probabilities <- predict(m7, type = "response")

# Convert probabilities to binary predictions (threshold = 0.5)
predictions <- ifelse(probabilities > 0.50, 1, 0)

# View the predicted probabilities and corresponding binary predictions
df$predicted_probabilities <- probabilities
df$predicted_class <- predictions


predictions2 <- predict(m7, newdata = df, type = "link", se.fit = TRUE)

# Calculate the upper and lower bounds for the confidence intervals
alpha <- 0.05  # 95% confidence interval
z_value <- qnorm(1 - alpha / 2)  # Z-score for 95% confidence

# Transform the linear predictors (log odds) to probabilities
df$predicted_probabilities <- plogis(predictions2$fit)
df$lower_bound <- plogis(predictions2$fit - z_value * predictions2$se.fit)
df$upper_bound <- plogis(predictions2$fit + z_value * predictions2$se.fit)

## ADD VALUES TO ANNOTATE ## 

model_summary <- summary(m7)
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
  "P-value = ", signif(p_value,2), "\n",
  "Residual df = ", df_residual
)

#### Quality Control & Checks ####

head(df)

eda <- df %>%
  select(ami_standardised, sum_wetseason_rainfall_mm, avg_daily_rainfall_mm, total_wetseason_discharge_ml, average_monthly_discharge_ml, total_days_minorflood, average_river_level)

ggpairs(eda)+
  theme_bw()


#### Plot Model ####

plot_1 <- ggplot(df, aes(x = average_river_level, y = pristis_presence)) +
  geom_point(aes(color = factor(pristis_presence)), size = 3) + 
  geom_ribbon(data = df, aes(x = average_river_level, ymin = lower_bound, ymax = upper_bound), fill = "grey", alpha = 0.3) + 
  geom_line(data = df, aes(x = average_river_level, y = predicted_probabilities), color = "red", size = 1) + 
  annotate("text", x = 2.8, y = 0.9, label = annotation_text, color = "black", size = 3.5, hjust = 0) +
  labs(
    x = "Standardised Australian Monsoon Index",
    y = "Probability of Sawfish Presence"
  ) +
  theme_bw() +
  scale_color_manual(values = c("grey34", "skyblue"))+
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_text(margin = margin(t = 15)),  
        axis.title.y = element_text(margin = margin(r = 15)))

print(plot_1)


#### General plotting ####


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



