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

#ALLLLL stitched together now 

sibs.all <- rbind(half.sibs, full.sibs); sibs.all

write.csv(sibs.all, file = "Daly_Sibs.csv")



