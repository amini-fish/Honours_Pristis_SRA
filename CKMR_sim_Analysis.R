### This may be the solution to simulating expected relatedness vals between dyads whilst accounting for an incredibly low fPos level ###

library(devtools)

#install_github("https://github.com/eriqande/CKMRsim/tree/master") 

library(CKMRsim)
library(dartRverse)
library(tidyverse)

if(system.file("bin", package = "CKMRsim") == "") {
  install_mendel(Dir = system.file(package = "CKMRsim"))
}

#install.packages("remotes")
#remotes::install_github("eriqande/CKMRsim")

#### Good to go ####

setwd("C:/Users/samue/OneDrive/Desktop/Honours/analysis")

gl <- get(load("C:/Users/samue/OneDrive/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl



gl.write.csv(gl, outfile = "outfile.csv", outpath = "C:/Users/samue/OneDrive/Desktop/Honours/analysis", verbose = NULL)

data <- read.csv("outfile.csv")

dim(data)


## -----------------------------------------------------------------------------------------------------

## SKIP TO LINE 107 TO START ANALYSIS ## 

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

### I want to transpose this genotype data 

genotypes_2 <- t(genotypes)

genotypes_2

genotypes_3 <- data.frame(genotypes_2)

genotypes_3 <- genotypes_3 %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = "_")

genotypes_3

ind_names <- gl@ind.names

genotypes_3$Indiv <- ind_names

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

## Write the data to a csv file so we can just load it in as a tibble for next time and bypass all the nonesense

#write.csv(long_genos, file = "Pristis_genofor_CKMRsim.csv" )

## START HERE ##

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

##----------------------------------------------------------------------------

### Creating more dataframes for our analysis 

ckmr <- create_ckmr(
  D = afreqs_ready,
  kappa_matrix = kappas[c("FS", "HS", "FC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

## -=-------------------------------------------------------------------------- 

### Select the dyads the want to examine further 
## In ths case it is only FS, HS, and FC

kappas # there are just individual ibd probabilites for each state based on relationship 

kappas[c("FS", "HS", "FC", "U"), ]

##-------------------------------------------------------------------------------------------------------------------------

##  Gives us our true relatedness log likelihoods based on the assumption of no linkage 
##  Good filtering should minimise the effect - can't fully escape it with SNPs


### OUR FPOS RATE 0.00002849002 i.e. 2.849 x 10-5

Qs <- simulate_Qij(
  ckmr, 
  calc_relats = c("FS", "HS", "FC", "U"),
  sim_relats = c("FS", "HS", "FC", "U") 
)

## FS and UP log likelhihood ratio example

FS_U_logls <- extract_logls(
  Qs,
  numer = c(FS = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 
mc_sample_simple(
  Qs,
  nu = "FS", 
  de = "U", 
  method = "IS"
)

## 366 has a low enough fpos rate

## So we want a FPOS rate that is 

## -------------------- Extract this for FSP/UP -------------------------------##

FS_U_gg <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

FS_U_gg

## ------------------------------------------------------------------------------

## Run the pairwise comparisons 

matchers <- find_close_matching_genotypes(
  LG = long_genos,
  CK = ckmr,
  max_mismatch = 500
)

matchers

## Run for all dyads 

pw_4_LRTs <- lapply(
  X = list(
    FSU = c("FS", "U"),
    FSHS = c("FS", "HS"),
    HSU = c("HS", "U"), 
    HSFC = c("HS", "FC"), 
    FCU = c("FC", "U")
  ),
  FUN = function(x) {
    pairwise_kin_logl_ratios(
      D1 = long_genos, 
      D2 = long_genos, 
      CK = ckmr,
      numer = x[1],
      denom = x[2],
      num_cores = 8
    )
  }
) %>%
  bind_rows(
    .id = "lr_type"
  ) %>%
  pivot_wider(names_from = lr_type, values_from = logl_ratio)

pw_4_LRTs

## ------------------------------------------------------------------------------

## First lets look at FULL SIB pairs 

topFS <- pw_4_LRTs %>%
  arrange(desc(FSU)) %>%
  filter(FSU > 0)

topFS

set.seed(42)# for the jittering

FS_U_gg +
  geom_jitter(
    data = topFS,
    mapping = aes(x = FSU, y = -0.002, colour = FSU > 0),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21
  )

## FS vs HS

FS_HS_gg <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(HS = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/HS Logl Ratio")

FS_HS_gg

FS_HS_logls <- extract_logls(
  Qs,
  numer = c(FS = 1),
  denom = c(HS = 1)
)

### Get our false positive and false negative rates to guide selection of T 
mc_sample_simple(
  Qs,
  nu = "FS", 
  de = "HS", 
  method = "IS"
)

#Quick visual inspection
View(pw_4_LRTs)

topFS_HS <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSHS)) %>%
  filter(FSHS > 0)

topFS_HS

set.seed(54) # for the jittering

FS_HS_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSHS, y = -0.002, colour = FSU > 0),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21
  )

## So we have our FSPs but we need to keep going with the rest of the individuals 

remaining_pairs <- pw_4_LRTs %>%
  anti_join(bind_rows(topFS), by = c("D2_indiv", "D1_indiv"))

remaining_pairs # work from this now

##------------------Half siblings and unrelated pairs--------------------------##

HS_UP_gg <- Qs %>%
  extract_logls(numer = c(HS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("HS / UP Logl Ratio")

HS_UP_gg

HS_UP_logls <- extract_logls(
  Qs,
  numer = c(HS = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 
mc_sample_simple(
  Qs,
  nu = "HS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001))

topHS_UP <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSU)) %>%
  filter(HSU > 52)

topHS_UP

seed(54) # for the jittering

HS_UP_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 52),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )

## -----------------Half siblings and first cousins ---------------------------##

HS_FC_gg <- Qs %>%
  extract_logls(numer = c(HS = 1), denom = c(FC = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("HS / FC Logl Ratio")

HS_FC_gg

HS_FC_logls <- extract_logls(
  Qs,
  numer = c(HS = 1),
  denom = c(FC = 1)
)

### Get our false positive and false negative rates to guide selection of T 
mc_sample_simple(
  Qs,
  nu = "HS", 
  de = "FC", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

#Quick visual inspection
topHS_FC <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSFC)) %>%
  filter(HSFC > -0.732)

topHS_FC
topHS_UP

set.seed(54) # for the jittering

HS_FC_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSFC, y = -0.002, colour = HSFC > -0.732),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )


remaining_pairs_2 <- remaining_pairs %>%
  anti_join(bind_rows(topFS, topHS_FC), by = c("D2_indiv", "D1_indiv"))

##---------------------First Cousins vs Unrelated -----------------------------##

all_FCUlogl <- ggplot(remaining_pairs_2, aes(x = FCU)) + 
  geom_histogram(bins = 30)

all_FCUlogl

all_FCU_logsl +
  xlim(-70, NA) +
  ggtitle("Pairs with First COusins > -70")

set.seed(52)

FC_U_gg <- Qs %>%
  extract_logls(numer = c(FC = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FC/UP Logl Ratio")

mc_sample_simple(
  Qs,
  nu = "FC", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)


FC_U_gg + 
  geom_jitter(
    data = remaining_pairs_2, 
    mapping = aes(x = FCU, y = -0.02, colour =  FCU > 2.45),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21, 
    size = 3
  ) +
  coord_cartesian(xlim = c(-50, 125), ylim = c(NA, 0.06))


topFC_UP <- remaining_pairs_2 %>% # remove the PO pairs 
  arrange(desc(FCU)) %>%
  filter(FCU > 2.45)

topFC_UP

## Matches that of EMIBD9 output and related which is cool

## Now remove HSP from the data to leave essentially all unrelated and maybe a cousin

FC_U_gg <- Qs %>%
  extract_logls(numer = c(FC = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FC/U Logl Ratio")

FC_U_gg

remaining

FC_U_gg + 
  geom_jitter(
    data = remaining %>% filter(FCU > -65),
    mapping = aes(x = FCU, y = -0.02),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21
  ) +
  coord_cartesian(xlim = c(-65, 125), ylim = c(NA, 0.06))

topFC <- remaining %>% 
  arrange(desc(FCU)) %>%
  filter(FCU > 0)


topFC
topHS
topFS

## -----------------------------------------------------------------------------

## Check our results against the simulated unrelated pairs

simU_FSHS <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(HS = 1)) %>% 
  filter(true_relat == "U")

ggplot() +
  geom_density(
    data = simU_FSHS,
    mapping = aes(x = logl_ratio),
    fill = "red",
    alpha = 0.3
  ) +
  geom_density(
    data = pw_4_LRTs, 
    mapping = aes(x = FSHS)
  )

## Findings match our expecations and there is a small tail to the right where our FSP are expected to be 

## HS FC

simU_HSFC <- Qs %>%
  extract_logls(numer = c(HS = 1), denom = c(FC = 1)) %>% 
  filter(true_relat == c("U", "FC"))

ggplot() +
  geom_density(
    data = simU_HSFC,
    mapping = aes(x = logl_ratio),
    fill = "red",
    alpha = 0.3
  ) +
  geom_density(
    data = pw_4_LRTs, 
    mapping = aes(x = HSFC)
  )
