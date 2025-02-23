#### Convert precipitation data ####

setwd("C:/Users/samue/Desktop/Honours/Daly_ENV")

library(dplyr)

#Lat, Long: -13.60 130.55 (DECIMAL DEGREES), 13 36'S 130 33'E

#### Load data ####

precip_data <- read.table("SILO_precipitation.txt", header = T, sep = "", dec = ".")

precip_data <- precip_data %>%
  slice(-1) 

head(precip_data)


#### Function to calculate total wetseason precipition #### 

## Precipitation is denoted as "Rain" in the dataset. To builld this function I need to pull data from 01/10/YearX to 30/04/YearY

subset_by_season <- function(data, date_col) {
  # Convert date column from dd/mm/yyyy format to Date type
  data[[date_col]] <- as.Date(data[[date_col]], format = "%d/%m/%Y")
  
  # Create an empty list to store subsets
  season_subsets <- list()
  
  # Loop from 2011 to 2023 (since last period is Oct 2013 - Apr 2024)
  for (year in 2011:2023) {
    start_date <- as.Date(paste0("01/10/", year), format = "%d/%m/%Y")
    end_date <- as.Date(paste0("30/04/", year + 1), format = "%d/%m/%Y")
    
    # Subset data within the date range
    season_data <- data[data[[date_col]] >= start_date & data[[date_col]] <= end_date, ]
    
    # Store in the list with a descriptive name
    season_subsets[[paste0("Season_", year, "_", year + 1)]] <- season_data
  }
  
  return(season_subsets)
}

# Example usage:
# Assuming your dataset is named `df` and the date column is "Date"
seasonal_data <- subset_by_season(precip_data, "Date2")

# Access a specific season, e.g., 2012-2013
head(seasonal_data[["Season_2012_2013"]])


#### Calculate average daily rainfall for each wet season ####

# Use sapply() to calculate average daily rainfall for each time block

str(seasonal_data[["Season_2011_2012"]])

avg_rainfall_per_season <- sapply(seasonal_data, function(season) {
  if (is.data.frame(season) && nrow(season) > 0) {  # Ensure it's a non-empty dataframe
    if ("Rain" %in% colnames(season)) {  # Check if "Rain" column exists
      season$Rain <- as.numeric(as.character(season$Rain))  # Convert Rain to numeric
      return(mean(season$Rain, na.rm = TRUE))  # Compute mean while ignoring NA
    } else {
      print("Rain column not found in this season!")  # Debug message
    }
  } else {
    print("Season is empty!")  # Debug message
  }
  return(NA)  # Return NA if data is missing
})

# Convert results to a data frame
rainfall_summary <- data.frame(
  Avg_Daily_Rainfall = avg_rainfall_per_season
)

# Print the result
print(rainfall_summary)


#### Calculate the total rainfall for each wet season ####

total_rainfall_per_season <- sapply(seasonal_data, function(season) {
  if (is.data.frame(season) && nrow(season) > 0) {  # Ensure it's a non-empty dataframe
    if ("Rain" %in% colnames(season)) {  # Check if "Rain" column exists
      season$Rain <- as.numeric(as.character(season$Rain))  # Convert Rain to numeric
      return(sum(season$Rain, na.rm = TRUE))  # Compute mean while ignoring NA
    } else {
      print("Rain column not found in this season!")  # Debug message
    }
  } else {
    print("Season is empty!")  # Debug message
  }
  return(NA)  # Return NA if data is missing
})

# Convert results to a data frame
total_rainfall <- data.frame(
  sum_wetseason_rainfall = total_rainfall_per_season
)

# Print the result
print(total_rainfall)


####  Merge two columns and save to a CSV  ####

seasonal_predictors <- bind_cols(total_rainfall, rainfall_summary); seasonal_predictors

write.csv(seasonal_predictors, file = "seasonal_predictors_silo.csv")

#### Now we do the same for discharge from the 8140040 gaguing station ####

oct <- read.csv(file = "G8140040_monthly_total_10.csv", 
                sep = ",", header=TRUE)

nov <- read.csv(file = "G8140040_monthly_total_11.csv", 
                sep = ",", header=TRUE)

dec <- read.csv(file = "G8140040_monthly_total_12.csv", 
                sep = ",", header=TRUE)

jan <- read.csv(file = "G8140040_monthly_total_01.csv", 
                sep = ",", header=TRUE)

feb <- read.csv(file = "G8140040_monthly_total_02.csv", 
                sep = ",", header=TRUE)

mar <- read.csv(file = "G8140040_monthly_total_03.csv")

apr <- read.csv(file = "G8140040_monthly_total_04.csv")

Monthly_streamflow <- list(oct, nov, dec, jan, feb, mar, apr)

df <- Monthly_streamflow %>%
  reduce(inner_join, by = 'Year')

df <- df %>%
  filter(Year >= 2011)

colnames(df) <- c("Year","Oct_flow", "Nov_flow", "Dec_flow", "Jan_flow", "Feb_flow", "Mar_flow", "Apr")

## Calculate the rowsums (unnecessary)

df2 <- df[,2:6]
df2 <- rowSums(df2)
df2 

#### Write to CSV ### 

write.csv(df, "discharge_2011_24.csv")
