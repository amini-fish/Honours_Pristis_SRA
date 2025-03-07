#### Install packages ####

install.packages("vcfR")

install.packages("remotes")
remotes::install_github("eriqande/CKMRsim")

if(system.file("bin", package = "CKMRsim") == "") {
  install_mendel(Dir = system.file(package = "CKMRsim"))
}


install_github("https://github.com/eriqande/CKMRsim/tree/master") 

#### Load packages ####


library(devtools)
library(vcfR)
library(CKMRsim)
library(dartRverse)
library(tidyverse)
library(devEMF)

#### Load the data ####

setwd("C:/Users/samue/Desktop/Honours/analysis")

gl <- get(load("C:/Users/samue/Desktop/Honours/analysis/daly_geno_clean.Rdata")); gl

gl <- gl.filter.monomorphs(gl)

gl.write.csv(gl, outfile = "outfile.csv", outpath = "C:/Users/samue/Desktop/Honours/analysis", verbose = NULL)

data <- read.csv("outfile.csv")

dim(data)

data

#### Prepare the data ####

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

#### START HERE - Long genos ####

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

afreqs_ready 

#### Create CKMR object ####

### Creating more dataframes for our analysis 

ckmr <- create_ckmr(
  D = afreqs_ready,
  kappa_matrix = kappas[c("FS", "HS", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

ckmr

#### Subset Kappas - IBD probs for desired relats ####

kappas # all possible

kappas[c("FS", "HS", "HFC", "U"), ]

#### Simulate dyads for each to approximate the LLRs and therefore the fpos/fneg ####

##  Gives us our true relatedness log likelihoods based on the assumption of no linkage 
##  Good filtering should minimise the effect - can't fully escape it with SNPs

### Desired fpos rate is 0.00002849002 i.e. 2.849 x 10-5

Qs <- simulate_Qij(
  ckmr, 
  calc_relats = c("FS", "HS", "HFC", "U"),
  sim_relats = c("FS", "HS", "HFC", "U") 
)

#### Linkage Model ####

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

## Our CKMR model that allows for some linkage in our data - i.e. more overlap between LLRs and more stringent kin-assignmnets

ckmr_link <- create_ckmr(
  D = afreqs_link,
  kappa_matrix = kappas[c("FS", "HS", "HFC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

## Simulate our linkage simulated true kinship values 
## This may take some time as uses MENDEL outside of Rstudio...

Qs_link_BIG <- simulate_Qij(
  ckmr_link, 
  calc_relats = c("FS", "HS", "HFC", "U"),
  sim_relats = c("FS", "HS", "HFC", "U"),
  unlinked = FALSE, 
  pedigree_list = pedigrees
)

#### Check the data for matching genotypes ####

# We allow for up to 500 matching loci here...

matchers <- find_close_matching_genotypes(
  LG = long_genos,
  CK = ckmr,
  max_mismatch = 500
)

matchers # have a geeze

#### Run pairwsie estimates of LLRs for our empirical data #### 

pw_4_LRTs <- lapply(
  X = list(
    FSU = c("FS", "U"),
    FSHS = c("FS", "HS"),
    HSHFC = c("HS", "HFC"), 
    HSU = c("HS", "U"),
    HFCU = c("HFC", "U")
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

##### Extract this for FSP/Unrelated ####

## Compute the thresholds we need using MCMC method

FSU_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "FS", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

print(FSU_thresholds)

# What we need - FPR ~ Lamda* = 161

## Visualise the dist of LLRs for non and linked models

FS_U_gg <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

FS_U_gg

FS_U_link_gg <- Qs_link_BIG %>% 
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  #ggtitle("FSU/U Linkage Logl Ratio") +
  scale_fill_brewer(palette = "RdYlGn")

FS_U_link_gg

FS_U_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FS = 1),
  denom = c(U = 1)
)

set.seed(42) # for the jittering

fs_fpos <- 3.49e-137
fs_fneg <-  0.0001

FSU_link_plot <- FS_U_link_gg + 
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x =FSU, y = -0.002, colour = FSU > 161),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 700,
           y = 0.027,  
           label = paste("FPR:", fs_fpos, "\nFNR:", fs_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

print(FSU_link_plot)

emf("C:/Users/samue/Desktop/Honours/FS_UP_LLR_plot.emf", width = 10, height = 8) 
print(FSU_link_plot)
dev.off()

write.csv(FSU_thresholds, "FSU_thresholds_linked.csv")

#### Extract Full Sibling Pairs #### 

topFS_U <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSU)) %>%
  filter(FSU > 161)

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

#### Full sib Half sib ####

#Linked
FS_HS_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(FS = 1), denom = c(HS = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7) +
  scale_fill_brewer(palette = "RdYlGn")

FS_HS_linked_gg

FS_HS_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(FS = 1),
  denom = c(HS = 1)
)

FSHS_MCMC <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "FS", 
  de = "HS", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001)
)

print(FSHS_MCMC)

#Quick visual inspection

topFS_HS <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSHS)) %>%
  filter(FSHS > -7.99)

topFS_HS

topFS_HS$rel <- rep("full-sibs")


set.seed(54) # for the jittering

FSHS_plot <- FS_HS_linked_gg +
  geom_jitter(
    data = pw_4_LRTs,
    mapping = aes(x = FSHS, y = -0.002, colour = FSHS > -7.99),
    width = 0, 
    height = 0.001,
    shape = 21, 
    fill = NA,
    size = 4
  ) +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  xlab("Log Likelihood Ratio") +
  ylab("Density")+
  labs(fill = "True Relationship")

print(FSHS_plot)

emf("C:/Users/samue/Desktop/Honours/FS_HS_LLR_plot.emf", width = 10, height = 8) 
print(FSHS_plot)
dev.off()

## So we have our FSPs but we need to keep going with the rest of the individuals 

remaining_pairs <- pw_4_LRTs %>%
  anti_join(bind_rows(topFS_HS), by = c("D2_indiv", "D1_indiv"))

remaining_pairs # work from this now

#### Half sib Unrelated model ####

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

### What we really care about is Linked Model

HS_UP_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)

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

print(HS_U_thresholds)

top_HS_UP <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSU)) %>%
  filter(HSU > 12.3)

topHS_UP

HSU_fpos <- 1.40e-26
HSU_fneg <- 0.0001 

HSU_plot <- HS_UP_linked_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSU, y = -0.002, colour = HSU > 12.3),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 400,
           y = 0.027,  
           label = paste("FPR:", HSU_fpos, "\nFNR:", HSU_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

plot(HSU_plot)

emf("C:/Users/samue/Desktop/Honours/HSU_LLR_plot.emf", width = 10, height = 8)  # Set the width and height in inches
print(HSU_plot)
dev.off()

remaining_pairs <- remaining_pairs %>%
  anti_join(bind_rows(topHS_UP), by = c("D2_indiv", "D1_indiv"))

### Move on to Half Sib to Half-first cousins ####

## Super interesting that our linkage model finds the same pairs as EMIBD9 but when unliked it doesn't and makes the dyads all cousins instead

HS_HFC_gg <- Qs %>%
  extract_logls(numer = c(HS = 1), denom = c(HFC = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("HS / FC Logl Ratio")

HS_HFC_gg

HS_HFC_logls <- extract_logls(
  Qs,
  numer = c(HS = 1),
  denom = c(HFC = 1)
)

### Get our false positive and false negative rates to guide selection of T 

mc_sample_simple(
  Qs,
  nu = "HS", 
  de = "HFC", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
)

#Quick visual inspection

glimpse(remaining_pairs)

topHS_HFC <- remaining_pairs %>% # remove the PO pairs 
  arrange(desc(HSHFC)) %>%
  filter(HSHFC > 19.9)

topHS_HFC

set.seed(54) # for the jittering

HS_HFC_gg +
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HSHFC, y = -0.002, colour = HSHFC > 1.10),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  )



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

#### HFC and Unrelated ####

HFCU_llr <- ggplot(remaining_pairs, aes(x = HFCU)) + 
  geom_histogram(bins = 30)

HFCU_llr

HFCU_llr +
  xlim(0, NA) +
  ggtitle("Pairs with First COusins > -70")

set.seed(52)

## Get Lamda and thresholds 

mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG,
  nu = "HFC", 
  de = "U", 
  method = "IS", 
  FNRs = c(0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.5))

## Visualise LLRs

HFC_U_gg <- Qs %>%
  extract_logls(numer = c(HFC = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)

HFC_U_gg + 
  geom_jitter(
    data = remaining_pairs, 
    mapping = aes(x = HFCU, y = -0.02, colour =  HFCU > 7.01),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21, 
    size = 3
  ) 
## Or with linked model 

HFC_UP_linked_gg <- Qs_link_BIG %>%
  extract_logls(numer = c(HFC = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.7)

HFC_UP_linked_gg

HFC_UP_linked_logls <- extract_logls(
  Qs_link_BIG,
  numer = c(HFC = 1),
  denom = c(U = 1)
)

### Get our false positive and false negative rates to guide selection of T 
HFC_U_thresholds <- mc_sample_simple(
  Qs,
  Q_for_fnrs = Qs_link_BIG, 
  nu = "HFC",
  de = "U", 
  method = "IS", 
  FNRs = c(0.5, 0.3, 0.2, 0.1, 0.05, 0.025, 0.01, 0.001, 0.0001))

print(HFC_U_thresholds)

hfc_fpos <- 0.0000863
hfc_fneg <-  0.5

HFC_plot <- HFC_UP_linked_gg + 
  geom_jitter(
    data = remaining_pairs,
    mapping = aes(x = HFCU, y = -0.002, colour = HFCU > 7.01),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21, 
    size = 3
  ) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_colour_manual(values = c("red", "black")) +
  theme_bw() +
  labs(fill = "True Relationship")+
  annotate("text", 
           x = 130,
           y = 0.099,  
           label = paste("FPR:", HSU_fpos, "\nFNR:", HSU_fneg), 
           hjust = 0, 
           size = 3, 
           color = "black")+
  xlab("Log Likelihood Ratio") +
  ylab("Density")

plot(HFC_plot)

emf("C:/Users/samue/Desktop/Honours/HSU_LLR_plot.emf", width = 10, height = 8)  # Set the width and height in inches
print(HSU_plot)
dev.off()


topFC_UP <- remaining_pairs_2 %>% # remove the PO pairs 
  arrange(desc(FCU)) %>%
  filter(FCU > 0.09)

topFC_UP
topFC_UP$rel <- rep("first-cousin") 


##-----------------Stitch all our kin pairs together--------------------------##

sib_groups <- rbind(topFS_HS, topHS_FC, topFC_UP)

sib_groups

?simulate_and_calc_Q
