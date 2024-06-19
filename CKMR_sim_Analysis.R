### This may be the solution to simulating expected relatedness vals between dyads whilst accounting for an incredibly low fPos level ###

library(devtools)

#install_github("https://github.com/eriqande/CKMRsim/tree/master") 

library(CKMRsim)
library(dartRverse)
library(tidyverse)

if(system.file("bin", package = "CKMRsim") == "") {
  install_mendel(Dir = system.file(package = "CKMRsim"))
}

#### Good to go ####

setwd("C:/Users/samue/OneDrive/Desktop/Honours/analysis")

gl <- get(load("C:/Users/samue/OneDrive/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl



gl.write.csv(gl, outfile = "outfile.csv", outpath = "C:/Users/samue/OneDrive/Desktop/Honours/analysis", verbose = NULL)

data <- read.csv("outfile.csv")

dim(data)


## -----------------------------------------------------------------------------------------------------


?gl2vcf

position <- gl$other$loc.metrics$ChromPos_WhaleShark_v1_2500len

chrom <- gl$other$loc.metrics$Chrom_WhaleShark_v1_2500len

gl2vcf(gl, 
       plink.bin.path = "C:/Users/samue/OneDrive/Desktop/Honours/analysis/plink", 
       snp.pos = "ChromPos_WhaleShark_v1_2500len", 
       snp.chr = "Chrom_WhaleShark_v1_2500len",
       outfile = "outfile", 
       outpath = getwd())

install.packages("vcfR")
library(vcfR)

geno_dummy <- read.vcfR("outfile.vcf")

genotypes <- geno_dummy@gt
genotypes <- data.frame(genotypes)

genotypes

drop <- c("FORMAT")
genotypes = genotypes[,!(names(genotypes) %in% drop)]
genotypes

references <- data.frame(geno_dummy@fix)

SNPS <- references$ID 


#genotypes$Locus <- SNPS
#genotypes <- genotypes %>% relocate(Locus)

rownames(genotypes) <- SNPS

genotypes


### I want to transpose this genotype data 

genotypes_2 <- t(genotypes)

View(genotypes_2)

genotypes_3 <- data.frame(genotypes_2)

#genotypes_tbl <- tibble(genotypes)
#genotypes_tbl

genotypes_3 <- genotypes_3 %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = "_.")

genotypes_3

genotypes_3$Indiv <- idnames

genotypes_3 <- genotypes_3 %>% relocate(Indiv)

genotypes_3 ## FUCK YES!!!

colnames(genotypes_3)

long_genos <- genotypes_3 %>%
  pivot_longer(
    cols = -Indiv, 
    names_to = c("Locus", "gene_copy"), 
    names_sep = "_.", 
    values_to = "Allele"
  )

long_genos

##---------------------------------------------------------------------------------------------------------------------------------------------------------

### Let's give this a crack shall we...###

gl2plink(gl, plink.bin.path = "C:/Users/samue/OneDrive/Desktop/Honours/analysis/plink", bed.files = T, outfile = "cleaned_geno_plink", outpath = getwd())

### The way that has worked thus far ### 

View(as.matrix(gl))

gl_wide <- as.matrix(gl)

## indnames are a character

### This code has given me a data matrix where the columns are a locus, and rows are individuals...now I want the rownames to be considered a column...

has_rownames(gl_wide) #do we have rownames though?

## NOPE

geno_wide <- as_tibble(gl_wide)

geno_wide$Indiv <- idnames

geno_wide <- geno_wide %>% relocate(Indiv)

geno_wide



