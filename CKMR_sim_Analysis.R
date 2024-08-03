### This may be the solution to simulating expected relatedness vals between dyads whilst accounting for an incredibly low fPos level ###

######################################
##                                  ##      
##        GO to line 126            ##
##            dumbass               ##
##                                  ##
##                                  ##
######################################

install.packages("vcfR")

#install.packages("remotes")
#remotes::install_github("eriqande/CKMRsim")

if(system.file("bin", package = "CKMRsim") == "") {
  install_mendel(Dir = system.file(package = "CKMRsim"))
}

library(devtools)

#install_github("https://github.com/eriqande/CKMRsim/tree/master") 

library(CKMRsim)
library(dartRverse)
library(tidyverse)

#### Good to go ####

setwd("C:/Users/samue/Desktop/Honours/analysis")

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

gl.write.csv(gl, outfile = "outfile.csv", outpath = "C:/Users/samue/Desktop/Honours/analysis", verbose = NULL)

data <- read.csv("outfile.csv")

dim(data)

data

library(vcfR)


## -----------------------------------------------------------------------------------------------------


## SKIP TO LINE 107 TO START ANALYSIS ## 

gl2vcf(gl, 
      plink.bin.path = "C:/Users/samue/Desktop/Honours/analysis/plink", 
     outfile = "outfile", 
       outpath = getwd())

geno_dummy <- read.vcfR("outfile.vcf")

geno_dummy

genotypes <- geno_dummy@gt

genotypes <- data.frame(genotypes)

genotypes

drop <- c("FORMAT")
genotypes = genotypes[,!(names(genotypes) %in% drop)]
genotypes

references <- data.frame(geno_dummy@fix)
references
SNPS <- references$ID 
SNPS <- gsub("_", "", SNPS)

genotypes$Locus <- SNPS
genotypes <- genotypes %>% relocate(Locus)

rownames(genotypes) <- SNPS

genotypes

### I want to transpose this genotype data 

genotypes_2 <- t(genotypes)

genotypes_2

genotypes_3 <- data.frame(genotypes_2)

genotypes_3 <- genotypes_3 %>% 
  separate_wider_delim(everything(), delim = "/", names_sep = "_")

length(gl@ind.names)

genotypes_3 <- genotypes_3[-1,]

gen

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

long_genos <- read.csv("Pristis_genofor_CKMRsim.csv"); long_genos

long_genos ## Looking MINT!!! We got there! 

long_genos <- long_genos[,-1]
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
  kappa_matrix = kappas[c("FS", "HS", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

## -=-------------------------------------------------------------------------- 

### Select the dyads the want to examine further 
## In ths case it is only FS, HS, and FC

kappas # there are just individual ibd probabilites for each state based on relationship 

kappas[c("FS", "HS", "HFC", "U"), ]

##-------------------------------------------------------------------------------------------------------------------------

##  Gives us our true relatedness log likelihoods based on the assumption of no linkage 
##  Good filtering should minimise the effect - can't fully escape it with SNPs

### OUR FPOS RATE 0.00002849002 i.e. 2.849 x 10-5

Qs <- simulate_Qij(
  ckmr, 
  calc_relats = c("FS", "HS", "HFC", "U"),
  sim_relats = c("FS", "HS", "HFC", "U") 
)

##------------------------Linkage Model-----------------------------------##

## We will use our linkage model to calibrate our false negative values 

## Lets use Pristis pectinata genome as our model for linkage model

n_chromosomes = 46 #chromosome number

L_genome = 2744*0.001 #Gb genome size

sl_diff = 2723055/149968914	#smallest chromosome / largest chromosome


fake_chromo_lengths <- geometric_chromo_lengths(
  n = n_chromosomes,
  L = L_genome,
  sl = sl_diff
)

fake_chromo_lengths$chrom_length_plot

set.seed(42)

afreqs_link <- sprinkle_markers_into_genome(afreqs_ready, fake_chromo_lengths$chrom_lengths)

ckmr_link <- create_ckmr(
  D = afreqs_link,
  kappa_matrix = kappas[c("FS", "HS", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

## Simulate our linkage simulated true kinship values 

Qs_link_BIG <- simulate_Qij(
  ckmr_link, 
  calc_relats = c("FS", "HS", "HFC", "U"),
  sim_relats = c("FS", "HS", "HFC", "U"),
  unlinked = FALSE, 
  pedigree_list = pedigrees
)

## ------------------------------------------------------------------------------

## Run the pairwise comparisons 

matchers <- find_close_matching_genotypes(
  LG = long_genos,
  CK = ckmr,
  max_mismatch = 500
)

head(long_genos)

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

## --------------------Full sib vs Unrelated ---------------------------------##

FS_U_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FS = 1),
  denom = c(U = 1)
)

set.seed(54) # for the jittering

FS_U_gg + #or FS_U_link_gg
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSU, y = -0.002, colour = FSU > 190),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21
  )

## -------------------- Extract this for FSP/HS------------------------------##

FS_U_gg <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

FS_U_gg

FS_U_link_gg <- Qs_link_BIG %>% 
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FSU/U Linkage Logl Ratio")

FS_U_link_gg

## First lets look at FULL SIB pairs 
### Get our false positive and false negative rates to guide selection of T 
FSU_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "FS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)   

write.csv(FSU_thresholds, "FSU_thresholds_linked.csv")

## Lambda of 178 

topFS_U <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSU)) %>%
  filter(FSU > 166)

topFS_U

topFS_U$rel <- rep("full-sibs")

set.seed(42)# for the jittering

FS_U_link_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSU, y = -0.002, colour = FSU > 190),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 2
  )

##-------------------------Full sibs vs Half Sibs ---------------------------##

#Unlinked 
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
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001)
)


#Linked
FS_HS_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(FS = 1), denom = c(HS = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Linkage FS/HS Logl Ratio")

FS_HS_linked_gg

FS_HS_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FS = 1),
  denom = c(HS = 1)
)

### Get our false positive and false negative rates to guide selection of T 
FSHS_MCMC <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "FS", 
  de = "HS", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001)
)

write.csv(FSHS_MCMC, "FSHS_thresholds_linked.csv")

#Quick visual inspection


topFS_HS <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSHS)) %>%
  filter(FSHS > 27.5) #83.6

topFS_HS

topFS_HS$rel <- rep("full-sibs")


set.seed(54) # for the jittering

FS_HS_linked_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSHS, y = -0.002, colour = FSHS > 27.5),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21
  )

## So we have our FSPs but we need to keep going with the rest of the individuals 

remaining_pairs <- pw_4_LRTs %>%
  anti_join(bind_rows(topFS_HS), by = c("D2_indiv", "D1_indiv"))

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
  filter(HSU > 54)

topHS_UP

set.seed(54) # for the jittering

HS_UP_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 54),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )

## and with the linked model
HS_UP_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Linked HS / UP Logl Ratio")

HS_UP_linked_gg

HS_UP_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HS = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 
HS_U_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001))

write.csv(HS_U_thresholds, "HSU_thresholds_linked.csv")

top_HS_UP <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSU)) %>%
  filter(HSU > 3.48)

topHS_UP

HS_UP_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 3.48),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )

## Super interesting that our linkage model finds the same pairs as EMIBD9 but when unliked it doesn't and makes the dyads all cousins instead

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
  filter(HSFC > 1.10)


set.seed(54) # for the jittering

HS_FC_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSFC, y = -0.002, colour = HSFC > 1.10),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )


### Now we simulate for linkage -------------------------------------------------

HS_FC_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HS = 1), denom = c(FC = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Linked HS / FC Logl Ratio")

HS_FC_linked_gg

HS_FC_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HS = 1),
  denom = c(FC = 1)
)

### Get our false positive and false negative rates to guide selection of T 
HS_FC_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HS", 
  de = "FC", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001))

write.csv(HS_FC_thresholds, "HS_FC_thresholds_linked.csv")


topHS_FC <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSFC)) %>%
  filter(HSFC> -21.1)

topHS_FC

HS_FC_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSFC, y = -0.002, colour = HSFC > -21.1),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )

## Check they match up
topHS_FC

topHS_FC$rel <- rep("half-sib")

remaining_pairs_2 <- remaining_pairs %>%
  anti_join(bind_rows(topFS_HS, topHS_FC), by = c("D2_indiv", "D1_indiv"))

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
  Q_for_fnrs = Qs_link_BIG,
  nu = "FC", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)


FC_U_gg + 
  geom_jitter(
    data = remaining_pairs_2, 
    mapping = aes(x = FCU, y = -0.02, colour =  FCU > 0.06),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21, 
    size = 3
  ) +
  coord_cartesian(xlim = c(-50, 125), ylim = c(NA, 0.06))

## Or with linked model 

FC_UP_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(FC = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("Linked FC / UP Logl Ratio")

FC_UP_linked_gg

FC_UP_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FC = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 
FC_U_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "FC",
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.025, 0.01, 0.001, 0.0001))


write.csv(FC_U_thresholds, "FC_U_thresholds_linked.csv")

FC_UP_linked_gg + 
  geom_jitter(
    data = remaining_pairs_2, 
    mapping = aes(x = FCU, y = -0.02, colour =  FCU > 0.09),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21, 
    size = 3
  ) +
  coord_cartesian(xlim = c(-50, 125), ylim = c(NA, 0.06))


topFC_UP <- remaining_pairs_2 %>% # remove the PO pairs 
  arrange(desc(FCU)) %>%
  filter(FCU > 0.09)

topFC_UP
topFC_UP$rel <- rep("first-cousin") 


##-----------------Stitch all our kin pairs together--------------------------##

sib_groups <- rbind(topFS_HS, topHS_FC, topFC_UP)

sib_groups
