## Elasmo conservation status 

library(ggplot2)

## Labels 
Status <- c("DD", "LC", "NT", "VU", "EN", "CR", "EX")
## No species (DD to EX)
count <- c(162, 508, 123, 180, 123, 91, 1)

## Merge into a dataframe

data <- data.frame(cbind(status, count))
 
## Get the proportions 

str(data)

data$fraction <- as.numeric(data$count)/sum(as.numeric(data$count))

View(data)

data$ymax <- cumsum(data$fraction)

data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

colo <- c("grey", "green", 'lightgreen', 'yellow', 'orange', 'red', 'black')


# Make the plot
status <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill = colo)) +
  geom_rect(colour = "black", fill = colo) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=12)) +
  ggtitle("Concervation Status of Elasmobranchs")


status


 scale_color_manual(name='Status',
                     breaks=c('DD','LC','NT', 'VU', 'EN', 'CR', 'EX'),
                     values=c('DD' = 'grey', "LC" = 'green', "NT" = 'lightgreen', "VU"= 'yellow', "EN" = 'orange', "CR" = "red", "EX" = 'black'))

status 

