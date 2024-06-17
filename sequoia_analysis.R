## install.packages("sequoia")

## Set working directory to wherever you want...
setwd("C:/Users/samue/OneDrive/Desktop/Honours/analysis")

## Load up every package known to man (we really only need a handful of them)

library(sequoia) # analysis
library(SNPRelate)
library(dartRverse) # cleaning/data wrangling 
library(ggplot2) # sexy plots
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.base) 
library(dartR.captive) # becuase dartR is a pain
library(gplots) 
library(graph4lg) # more sexy plots
library(viridis) # why not
library(ggraph) # even more sexy plots

### Now that we've loaded up the packages we need we can inport our CLEANED snp data 

gl <- get(load("C:/Users/samue/OneDrive/Desktop/Honours/analysis/daly_geno_clean.Rdata")) 

gl # always good to have a check and make sure it's what you wanted...

gl.smearplot(gl) ## nice and polymorphic plus no na's YAY!

### If you want to refilter seperately feel free - can also pull datasets filtered at different thresholds (e.g. MAF or callrate etc)

### For now we are going with the filtered dataset because our full dataset will no doubt change the picture by a fair amount

### Using Sequioa ### 

## For the user guide from CRAN see this link... chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/sequoia/vignettes/vignette-main.pdf 


