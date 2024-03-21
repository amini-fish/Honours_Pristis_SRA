setwd("C:/Users/samue/Desktop/Honours - Sawfish/pristis_data")

install.packages("hierfstat")
devtools::install_version("ggplot2", "3.4.4")
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)



#save(gl, file = "Sawfish Prelim Data.RData") 

data <- "Sawfish_SNPs_genotyped.csv"
meta <- "Sawfish_meta2.csv"

#Compile them into one gl

raw.gl <- dartR.base::gl.read.dart(filename = data, ind.metafile = meta)
gl.smearplot(raw.gl)
?gl.smearplot


data.gl <- dartR.base::gl.read.dart(filename = data, ind.metafile = meta); data.gl

pop(data.gl) <- data.gl@other$ind.metrics$pop

table(pop(data.gl))

###### FILTERING #####

## Get rid of sex linked loci 
library(devtools)

library(dartRverse)
#install_github("green-striped-gecko/dartR.sexlinked@dev")
library(dartR.sexlinked)

data.gl <- gl.filter.sexlinked(data.gl, system = "xy")

data.gl <- data.gl$autosomal

## Get rid of unreliable loci

gl.report.reproducibility(data.gl)

data.gl <- dartR.base::gl.filter.reproducibility(data.gl, threshold = 0.99)
  
## Callrate 
dartR.base::gl.report.callrate(data.gl)

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "loc", threshold = 0.99)

#Get rid of low and super high read depth loci
#do twice so you can zoom in

dartR.base::gl.report.rdepth(data.gl)

data.gl <- dartR.base::gl.filter.rdepth(data.gl, lower = 10, upper = 75)

data.gl <- dartR.base::gl.filter.secondaries(data.gl)

#Very low filter – this is only to get rid of your really bad individuals
dartR.base::gl.report.callrate(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "ind", threshold = 0.99)

#Always run this after removing individuals – removes loci that are no longer variable
data.gl <- dartR.base::gl.filter.monomorphs(data.gl)

gl.smearplot(data.gl)
### remove evidence of DNA contamination ## Important for kin finding 
dartR.base::gl.report.heterozygosity(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.heterozygosity(data.gl,t.lower = 0.2,  t.upper = 0.25)

###Check our filtering steps ###

data.gl@other$history

# Some popgen EDA 

?gl.dist.pop
dist <- gl.dist.pop(data.gl, method = "euclidean")

gl.plot.heatmap(dist)

## Try a PCoA

pc <- gl.pcoa(data.gl)

pc.plot <- gl.pcoa.plot(glPca = pc, 
                        data.gl)

## No clustering across any of the pops including the Fitzroy :0

## Lets sus out the VDG pops

pc.al <- gl.pcoa(Alligator)
pc.al.plot <- gl.pcoa.plot(glPca = pc.al, 
                          Alligator)

## Nothing at all...

gl.report.heterozygosity(data.gl, 
                         method = "pop")

## Calculate allelic richness 

gi <- gl2gi(data.gl)
hfstat <- genind2hierfstat(gi)

#calculate allelic richness
ar <- allelic.richness(hfstat)
names(ar)
ar$min.all

View(ar)
summary(ar$Ar) 
ar <- as.data.frame(ar$Ar)

View(ar)
mean.ar <- colMeans(ar)

#Change to long format??????? 
par(mar=c(8,3,3,3))

reshape(ar,
        idvar= "Pop", 
        varying = 1:7, 
        timevar = 
        direction = "long")

colo <- viridis::magma(n = 10, begin = 0.2, end = 0.99)

boxplot(ar, ylab="Allelic richness", las=2, col = colo, 
        main = "Allelic Richness of all P. pristis Pops")

#Keep only the daly river individuals 

#Daly <- gl.keep.pop(data.gl, pop.list = "Daly", v = 3)

#save(Daly, file = "Daly Sawfish Clean Prelim Data.RData") 

#Read in the Daly River data 
Daly.gl <- get(load("C:/Users/samue/Desktop/Honours - Sawfish/pristis data for Sam/Daly Sawfish Clean Prelim Data.RData")); Daly.gl

#Now lets run a rough analysis using gl.run.EMIBD

library(devtools)
install_github("green-striped-gecko/dartR.captive@dev_sam")
library(dartR.captive)

daly.rel <- gl.run.EMIBD9(Daly, emibd9.path =  "C:/EMIBD9")

cols <- viridis::magma(n = 100)
dartR.captive::gl.grm.network(daly.rel$rel, 
               x = Daly.gl, 
               method = "gh", 
               link.size = 1.2, 
               title = "Daly River 2012 & 2013 Relatedness")

ibd9Tab <- daly.rel[[2]]; daly.rel

# Kick out self comparisons
ibd9Tab <- ibd9Tab[ibd9Tab$Indiv1 != ibd9Tab$Indiv2,  c(1, 2, 21)]; ibd9Tab

# Add cohorts 
Cohort1 <- Daly@other$ind.metrics$Cohort[as.numeric(ibd9Tab$Indiv1)]
Cohort2 <- Daly@other$ind.metrics$Cohort[as.numeric(ibd9Tab$Indiv2)]
# Flag pairs trapping within G, T and in between (BW)
CC <- ifelse(Cohort1 == Cohort2, 
             yes = "same", 
             no = "different")

# Combine together
ibd9DT <- as.matrix(cbind(ibd9Tab, Cohort1, Cohort2, CC)); ibd9DT

# Compute the mean relatedness
ibd9DT[, mean(as.numeric('r(1,2)')), by = CC]

# Combine together
ibd9DT <- data.frame(cbind(ibd9Tab, Cohort1, Cohort2, CC)); ibd9DT

# Compute the mean relatedness
ibd9DT[, mean(as.numeric(`r(1,2)`)), by=CC]


#Putative siblings
sibs <- ifelse(ibd9DT$r.1.2. >= 0.092 & ibd9DT$r.1.2. <= 0.158, 
               yes = "half-siblings", 
               no = ifelse(ibd9DT$r.1.2. >=0.204 & ibd9DT$r.1.2. <= 0.296, 
                           yes = "fsp or pop",
                           "unelated"))

#Add the sibling assignments to a new data frame 

ibd9DT.2 <- data.frame(cbind(ibd9DT, sibs)); ibd9DT.2

####  Assign kin to sibling network  #####

#hsps - extract

half.sibs <- subset(ibd9DT.2, sibs == "half-siblings"); half.sibs

#fsps - extract 

full.sibs <- subset(ibd9DT.2, sibs == "fsp or pop"); full.sibs

#ALLLLL stitched together now 

sibs.all <- rbind(half.sibs, full.sibs); sibs.all

# Remove duplicated pairs (i.e., Ab & BA)

sibs.all[!duplicated(sibs.all$r.1.2.), ] #this works but I don't trust it for a larger dataset where there may be the same relatedness value for multiple pairs

# keep it for now but work on it with someone

# The next step is making a nice network of all fsp & hsp similar to Patterson et al. 2022
# See supplementary materials (no code but visual aid)
# ggnet2(bip, color = "mode", palette = col, edge.size = "weights") 

gl.report.heterozygosity(Daly.gl, method = "pop")

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



#ALLLLL stitched together now 

sibs.all <- rbind(half.sibs, full.sibs); sibs.all

write.csv(sibs.all, file = "Daly_Sibs.csv")



