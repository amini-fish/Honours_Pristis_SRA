setwd("C:/Users/samue/Desktop/Honours_Sawfish/analysis")

#Install some packages

#install_github("green-striped-gecko/dartR.sexlinked@dev")
#install.packages("hierfstat")
#devtools::install_version("ggplot2", "3.4.4")

## Load all required packages
library(dartR.sexlinked)
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)

data <- "Sawfish_SNPs_genotyped.csv"
meta <- "Sawfish_meta2.csv"

## Load all genotype data and metadata

data.gl <- dartR.base::gl.read.dart(filename = data, ind.metafile = meta); data.gl

pop(data.gl) <- data.gl@other$ind.metrics$pop

table(pop(data.gl))

gl.smearplot(data.gl) ## A quick smearplot to check genotype data

## Save unfiltered data 

#save(data.gl, file = "pristis_geno_raw.RData") 

###### FILTERING #####

## Get rid of sex linked loci 

data.gl <- gl.filter.sexlinked(data.gl, system = "xy")

data.gl <- data.gl$autosomal #keeps only our autosomal SNPs

## Get rid of unreliable loci

gl.report.reproducibility(data.gl)

data.gl <- dartR.base::gl.filter.reproducibility(data.gl, threshold = 0.99)

## Callrate 
dartR.base::gl.report.callrate(data.gl)

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "loc", threshold = 0.99)

#Get rid of low and super high read depth loci

dartR.base::gl.report.rdepth(data.gl)

data.gl <- dartR.base::gl.filter.rdepth(data.gl, lower = 10, upper = 75)

## Remove secondary loci

data.gl <- dartR.base::gl.filter.secondaries(data.gl)

#Very low filter – this is only to get rid of your really bad individuals

dartR.base::gl.report.callrate(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "ind", threshold = 0.8)

#Always run this after removing individuals – removes loci that are no longer variable
data.gl <- dartR.base::gl.filter.monomorphs(data.gl)

## remove evidence of DNA contamination ## Important for kin finding 
dartR.base::gl.report.heterozygosity(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.heterozygosity(data.gl,t.lower = 0.2,  t.upper = 0.25)

## Visualise our cleaned data 

gl.smearplot(data.gl)

### Check our filtering steps ###

data.gl@other$history

### SAVE CLEANED SNP DATA ###












