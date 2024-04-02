### DISCLAIMER FORM FLO: These packages are probably out of date atm.

### You need to be able to install R packages from source. GNU/Linux
### OSes can do this out-of-the-box, but Windows users need to add
### Rtools and Mac users need xcode.

## instructions for Rtools installation on Windows:
## https://cran.r-project.org/bin/windows/Rtools/

## instructions for xcode installation on Mac:
## https://clanfear.github.io/CSSS508/docs/compiling.html

### The CRAN packages: Rcpp and remotes
install.packages("Rcpp")
install.packages("remotes")

### The MVB support packages
install.packages(pkgs = c("atease", "mvbutils", "vecless"), repos = "https://markbravington.github.io/Rmvb-repo")

### SNPrelate - only used in part of Session Three
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")

### kinference and gbasics
remotes::install_local("packages/kinference-master.zip", subdir = "kinference")
remotes::install_local("packages/gbasics-master.zip", subdir = "gbasics")
remotes::install_local("packages/genocalldart-master.zip", subdir = "genocalldart")
## Note that you'll have to provide a filepath/to/the/.zip/file that is
## appropriate to your OS.

## Alternative installation for gbasics and kinference
`install.packages(pkgs = c("gbasics", "kinference"), repos = "https://markbravington.github.io/Rmvb-repo")`

## Now, test the installs on this example analysis:
load("data/SHS_cleaned.Rdata") ## you'll need to input the correct filepath for your OS
library(kinference)


pdf(file = "showShaneThisPdf.pdf")
dups <- find_duplicates( SHS, max_diff_loci = 200, showPlot = TRUE)
SHS_b <- SHS[ -c( drop_dups_pairwise_equiv( dups[,2:3])),]
Lfish <- ilglk_geno(SHS_b, list(nclass = 1000, xlim = c(-1500,-1250) ))
SHS_c <-  SHS_b[ ((Lfish > -1450) & (Lfish < -1275)), ]
SHS_d <- hsp_power( SHS_c, k = 0.5)
hetz_rich <- hetzminoo_fancy( SHS_d, "rich")
SHS_e <- SHS_d[ (hetz_rich > 0.215) & (hetz_rich < 0.275), ]
hetz_poor <- hetzminoo_fancy(SHS_e, 'poor')
SHS_f <- SHS_e[ (hetz_poor > 0.175) & (hetz_poor < 0.225), ]
SHS_g <- prepare_PLOD_SPA( SHS_f)
hsps <- find_HSPs( SHS_g, keep_thresh = -20, minbin = -120)
PLOD_loghisto(hsps)
dev.off()


### In the event of errors:
packageVersion("kinference") ## should be 0.0.125
packageVersion("gbasics") ## should be 0.0.93

## If you have had another version of kinference or gbasics installed,
## R _may_ lie to you about versions in your packageVersion() and
## sessionInfo() calls. To be doubly sure that the information is
## correct, restart R after installing new packages. (bug confirmed
## under RStudio on MacOS Mojave, R 3.6.1)