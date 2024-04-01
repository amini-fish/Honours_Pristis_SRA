setwd("C:/Users/samue/Desktop/Honours_Sawfish/analysis")

install.packages("hierfstat")
devtools::install_version("ggplot2", "3.4.4")
install.packages("related")

library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.sexlinked)
library(viridis)
library(related)


#install.packages("related", repos="http://R-Forge.R-project.org")

## Not run: 
#---Read data into R---#

geno <- gl2related(data.gl, 
                   outfile = "related.txt", , 
                   outpath = "C:/Users/samue/Desktop/Honours_Sawfish/analysis",
                   v = 5)

input <- readgenotypedata("related.txt")

#---Calculate Relatedness---#

output <- coancestry(genotype.data = input$gdata, 
                     trioml = 2, 
                     wang = 2, 
                     quellergt=2, 
                     allow.inbreeding = T,
                     rng.seed = 42)

#---View Point Estimates---#
output$relatedness

#---View 95
output$relatedness.ci95

## End(Not run)	

### Save trioML results into a new dataframe

trioML <- data.frame(output$relatedness[, c(2,3,5)])


mean <- mean(trioML$trioml)
sd <- sd(trioML$trioml)

## Extract the indiviudals with relatedness above x and y -  need to find a solid guideline for this assignment not arbitrary numbers.

pull.kin <- ifelse(trioML$trioml >= 0.15 & trioML$trioml <= 0.25, 
               yes = "hsp", 
               no = ifelse(trioML$trioml >=0.35 & trioML$trioml <= 0.55, 
                           yes = "fsp",
                           "unelated"))

## Combine the two dataframes so we have all putative hsps and fsps identified.

trioML <- data.frame(cbind(trioML, pull.kin))

## from this data we can extract our siblings and put them in a new data frame

trioML.fsp <- trioML[ which(trioML$pull.kin == "fsp"),]
trioML.hsp <- trioML[ which(trioML$pull.kin == "hsp"),]

trioML.sibs<- rbind(trioML.fsp, trioML.hsp)

trioML.sibs

###############
## Visualise ##
###############

rel.hist <- ggplot(data = trioML, aes(x = trioml)) +
  geom_histogram(bins = 100, col = I("black")) +
  labs(title = "Histogram of kinship coefficients (θ) between 2012 & 2013 sawfish") +
  xlab("Relatedness (R)") +
  ylab("Count")
  
print(rel.hist)


############# EXTRA NONESENSE ###################
theme(plot.title = element_text(hjust = 0.5, size = 22), 
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      axis.text = element_text(size = 13)) +
  geom_vline(xintercept = mean, size = 1.1, col = "green")
#geom_vline(xintercept = 0.125, size = 1.1, col = "orange") +
#geom_vline(xintercept = 0.250, size = 1.1, col = "skyblue") +
geom_vline(xintercept = (-sd), col = "black", size = 1, linetype="dotted") +
  geom_vline(xintercept = (sd), col = "black", size = 1,  linetype="dotted") +
  #geom_vline(xintercept = 0.204, col = "black", size = 1, linetype="dotted") +
  #geom_vline(xintercept = 0.296, col = "black", size = 1,  linetype="dotted") +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(n.breaks = 10)









  