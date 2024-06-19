### This may be the solution to simulating expected relatedness vals between dyads whilst accounting for an incredibly low fPos level ###

library(devtools)

#install_github("https://github.com/eriqande/CKMRsim/tree/master") 

library(CKMRsim)
library(dartRverse)
library(tidyverse)

if(system.file("bin", package = "CKMRsim") == "") {
  install_mendel(Dir = system.file(package = "CKMRsim"))
}

install.packages("remotes")
remotes::install_github("eriqande/CKMRsim")

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

## ------------------------ Fake chromosome stuff ----------------------------##
fake_chromo_lengths <- geometric_chromo_lengths(
  n = 16,
  L = 2.4,
  sl = 0.25
)

fake_chromo_lengths$chrom_length_plot

set.seed(42)

afreqs_link <- sprinkle_markers_into_genome(afreqs_ready, fake_chromo_lengths$chrom_lengths)

ckmr_link <- create_ckmr(
  D = afreqs_link,
  kappa_matrix = kappas[c("PO", "FS", "HS", "FC", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
)

Qs_link_BIG <- simulate_Qij(
  ckmr_link, 
  calc_relats = c("PO", "FS", "HS", "FC", "U"),
  sim_relats = c("PO", "FS", "HS", "FC", "U"),
  unlinked = FALSE, 
  pedigree_list = pedigrees
)

Qs_link_BIG <- simulate_Qij(
  ckmr_link, 
  calc_relats = c("PO", "FS", "HS", "FC", "U"),
  sim_relats = c("PO", "FS", "HS", "FC", "U"),
  unlinked = FALSE, 
  pedigree_list = pedigrees
)

## -------------------- Extract this for POPs/UP -------------------------------##

PO_U_gg <- Qs %>%
  extract_logls(numer = c(PO = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("PO/U Logl Ratio")

PO_U_gg

FS_U_gg <- Qs %>%
  extract_logls(numer = c(FS = 1), denom = c(U = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/U Logl Ratio")

FS_U_gg

## We are running into an error on line 190 - moev on and come back to this 

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
    POU = c("PO", "U"),
    POFS = c("PO", "FS"),
    FSHS = c("FS", "HS"),
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

## PO pairs (there shouldn't be any)

topPO <- pw_4_LRTs %>%
  arrange(desc(POU)) %>%
  filter(POU > 0)

topPO

set.seed(42)# for the jittering

PO_U_gg +
  geom_jitter(
    data = topPO,
    mapping = aes(x = POU, y = -0.002, colour = POFS > 0),
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

topFS <- pw_4_LRTs %>% # remove the PO pairs 
  arrange(desc(FSHS)) %>%
  filter(FSHS > -20)
topFS

set.seed(54) # for the jittering
FS_HS_gg +
  geom_jitter(
    data = topFS,
    mapping = aes(x = FSHS, y = -0.002),
    width = 0, 
    height = 0.001, 
    fill = NA,
    shape = 21
  )

## So we have our FSPs but we need to keep going with the rest of the individuals 

remaining <- pw_4_LRTs %>%
  anti_join(bind_rows(topFS), by = c("D2_indiv", "D1_indiv"))

remaining

all_HSFC_logsl <- ggplot(remaining, aes(x = HSFC)) + 
  geom_histogram(bins = 30)

all_HSFC_logsl

all_HSFC_logsl +
  xlim(-70, NA) +
  ggtitle("Pairs with HSHAN > -70")

set.seed(52)

HS_FC_gg <- Qs %>%
  extract_logls(numer = c(HS = 1), denom = c(FC = 1)) %>%
  ggplot(aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25) +
  ggtitle("FS/HS Logl Ratio")


HS_FC_gg + 
  geom_jitter(
    data = remaining %>% filter(HSFC > -65),
    mapping = aes(x = HSFC, y = -0.02),
    width = 0, 
    height = 0.01, 
    fill = NA,
    shape = 21
  ) +
  coord_cartesian(xlim = c(-65, 125), ylim = c(NA, 0.06))

pw_4_LRTs
topHS <- remaining %>% # remove the PO pairs 
  arrange(desc(HSFC)) %>%
  filter(HSFC > -10)

topHS
topFS

## Matches that of EMIBD9 output and related which is cool

## Now remove HSP from the data to leave essentially all unrelated and maybe a cousin

remaining <- remaining %>%
  anti_join(bind_rows(topFS, topHS), by = c("D2_indiv", "D1_indiv"))

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
