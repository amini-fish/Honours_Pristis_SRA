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
SNPS <- gsub("_", "", SNPS)

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
  separate_wider_delim(everything(), delim = "/", names_sep = "_")

genotypes_3

genotypes_3$Indiv <- idnames

genotypes_3 <- genotypes_3 %>% relocate(Indiv)

genotypes_3 ## FUCK YES!!!

colnames(genotypes_3)

long_genos <- genotypes_3 %>%
  pivot_longer(
    cols = -Indiv, 
    names_to = c("Locus", "gene_copy"), 
    names_sep = "\\_", 
    values_to = "Allele"
  )

long_genos ## Looking MINT!!! We got there! 

locus_names <- unique(long_genos$Locus)
afreqs_ready <- long_genos %>%
  count(Locus, Allele) %>%  
  group_by(Locus) %>%
  mutate(
    Freq = n / sum(n),
    Chrom = "Unk",
    Pos = as.integer(factor(Locus, levels = locus_names))
  ) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele)) %>%
  reindex_markers()

ckmr <- create_ckmr(
  D = afreqs_ready,
  kappa_matrix = kappas[c("PO", "FS", "HS", "FC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

kappas

kappas[c("PO", "FS", "HS", "FC", "U"), ]

##-------------------------------------------------------------------------------------------------------------------------

Qs <- simulate_Qij(
  ckmr, 
  calc_relats = c("PO", "FS", "HS", "FC", "U"),
  sim_relats = c("PO", "FS", "HS", "FC", "U") 
)

PO_U_logls <- extract_logls(
  Qs,
  numer = c(PO = 1),
  denom = c(U = 1)
)

# And we can plot the distribution of those logl ratios for each
# of the different true relationships. 

ggplot(
  PO_U_logls,
  aes(x = logl_ratio, fill = true_relat)
) +
  geom_density(alpha = 0.25)

mc_sample_simple(
  Qs,
  nu = "PO"
)
