setwd("C:/Users/samue/Desktop/Honours - Sawfish/pristis_data")

install.packages("hierfstat")
devtools::install_version("ggplot2", "3.4.4")
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.sexlinked)
library(viridis)
heat

#save(gl, file = "Sawfish Prelim Data.RData") 

#data <- "Sawfish_SNPs_genotyped.csv"
#meta <- "Sawfish_meta2.csv"

## Compile them into one gl

#data.gl <- dartR.base::gl.read.dart(filename = data, ind.metafile = meta); data.gl

#pop(data.gl) <- data.gl@other$ind.metrics$pop

#table(pop(data.gl))

## Keep the Daly Inds

Daly.gl <- gl.keep.pop(data.gl, pop.list = "Daly")

## Make a smearplot to show before and after effects of filtering on data quality 

gl.smearplot(Daly.gl, ind.labels = T, plot.theme = theme_classic(), plot.file = "unfiltered_smear_daly", plot.dir = "/pristis_data")


################################################################################

## Filtering steps

Daly.gl <- dartR.sexlinked::gl.filter.sexlinked(Daly.gl, system = "xy")

data.gl <- Daly.gl$autosomal

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

data.gl <- dartR.base::gl.filter.maf(data.gl, t = 0.1)

#Very low filter – this is only to get rid of your really bad individuals

dartR.base::gl.report.callrate(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.callrate(data.gl, method = "ind", threshold = 0.9)

#Always run this after removing individuals – removes loci that are no longer variable

data.gl <- dartR.base::gl.filter.monomorphs(data.gl)

### remove evidence of DNA contamination ## Important for kin finding 

dartR.base::gl.report.heterozygosity(data.gl, method = "ind")

data.gl <- dartR.base::gl.filter.heterozygosity(data.gl,t.lower = 0.2,  t.upper = 0.3)


###Check our filtering steps ###

data.gl@other$history

gl.smearplot(data.gl, ind.labels = T)

unfiltered_smear_daly <- readRDS("C:/Users/samue/AppData/Local/Temp/RtmpgTwQhk/unfiltered_smear_daly.RDS"); unfiltered_smear_daly

data.gl

################################################################################

## Run some preliminary kin finding analyses

install_github("green-striped-gecko/dartR.captive@dev_sam")
library(dartR.captive)
library(gplots)

daly.rel <- gl.run.EMIBD9(data.gl, emibd9.path =  "C:/EMIBD9")

ibd9Tab <- daly.rel[[2]]; daly.rel

daly.rel$rel
colnames(daly.rel$rel)

# Kick out self comparisons
ibd9Tab <- ibd9Tab[ibd9Tab$Indiv1 != ibd9Tab$Indiv2,  c(1, 2, 21)]; ibd9Tab

# Add cohorts 
Cohort1 <- Daly@other$ind.metrics$cohort[as.numeric(ibd9Tab$Indiv1)]
Cohort2 <- Daly@other$ind.metrics$cohort[as.numeric(ibd9Tab$Indiv2)]
# Flag pairs trapping within G, T and in between (BW)
CC <- ifelse(Cohort1 == Cohort2, 
             yes = "same", 
             no = "different")

# Combine together
ibd9DT <- as.matrix(cbind(ibd9Tab, Cohort1, Cohort2, CC)); ibd9DT

# Combine together
ibd9DT <- data.frame(cbind(ibd9Tab, Cohort1, Cohort2, CC)); ibd9DT

mean(as.numeric(ibd9DT$r.1.2.)) ## Mean is 0.0129 
sd(as.numeric(ibd9DT$r.1.2.)) ## SD is +- 0.0323


## Test for within and between cohort means, note this isn't perfect
group_mean<- aggregate(x= as.numeric(ibd9DT$r.1.2.),
                       # Specify group indicator
                       by = list(ibd9DT$CC),      
                       # Specify function (i.e. mean)
                       FUN = mean)

group_mean

group_sd <- aggregate(x= as.numeric(ibd9DT$r.1.2.),
                       # Specify group indicator
                       by = list(ibd9DT$CC),      
                       # Specify function (i.e. mean)
                       FUN = sd)

group_sd
#A quick look at the spread of rel values...

hist(as.numeric(ibd9DT$r.1.2.), breaks = 250)

## Simulate 10 HSP & FSP to get rough relatedness estimates - use > 10 LOL

fsp.sim <- dartR.captive::gl.sim.relatedness(data.gl, rel = "full.sib", nboots = 10, emibd9.path = "C:/EMIBD9")


hsp.sim <- dartR.captive::gl.sim.relatedness(data.gl, rel = "half.sib", nboots = 10, emibd9.path = "C:/EMIBD9")


## Putative siblings - we use the 95% CI from the simulations to groundtruth our sibling relationships as mean rel is 0.01
sibs <- ifelse(ibd9DT$r.1.2. >= 0.092 & ibd9DT$r.1.2. <= 0.158, 
               yes = "hsp", 
               no = ifelse(ibd9DT$r.1.2. >=0.204 & ibd9DT$r.1.2. <= 0.296, 
                           yes = "fsp",
                           "unelated"))

#Add the sibling assignments to a new data frame 

ibd9DT.2 <- data.frame(cbind(ibd9DT, sibs)); ibd9DT.2

####  Assign kin to sibling network  #####

#hsps - extract

half.sibs <- subset(ibd9DT.2, sibs == "hsp"); half.sibs

#fsps - extract 

full.sibs <- subset(ibd9DT.2, sibs == "fsp"); full.sibs

#ALLLLL stitched together now 

sibs.all <- rbind(half.sibs, full.sibs); sibs.all

# Remove duplicated pairs (i.e., Ab & BA)

sibs.all[!duplicated(sibs.all$r.1.2.), ] 
sibs.all

# Kick out self comparisons
ibd9Tab <- ibd9Tab[ibd9Tab$Indiv1 != ibd9Tab$Indiv2,  c(1, 2, 21)]; ibd9Tab

str(ibd9DT)

mean.rel <- mean(as.numeric(ibd9DT$r.1.2.)) ## Mean is 0.0129 
sd.rel <- sd(as.numeric(ibd9DT$r.1.2.)) ## SD is +- 0.0323
med <- median(as.numeric(ibd9DT$r.1.2.)) ## 0.00575

rel.hist <- ggplot(data = ibd9DT, aes(x = as.numeric(r.1.2.))) +
  geom_histogram(bins = 200, col = I("black")) +
  labs(title = "Histogram of kinship coefficients (θ) between 2012 & 2013 sawfish") +
  xlab("Kinship (θ)") +
  ylab("Count") +
    theme(plot.title = element_text(hjust = 0.5, size = 22), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        axis.text = element_text(size = 13)) +
  geom_vline(xintercept = med, size = 1.1, col = "green") +
  geom_vline(xintercept = 0.125, size = 1.1, col = "orange") +
  geom_vline(xintercept = 0.250, size = 1.1, col = "skyblue") +
  geom_vline(xintercept = 0.092, col = "black", size = 1, linetype="dotted") +
  geom_vline(xintercept = 0.158, col = "black", size = 1,  linetype="dotted") +
  geom_vline(xintercept = 0.204, col = "black", size = 1, linetype="dotted") +
  geom_vline(xintercept = 0.296, col = "black", size = 1,  linetype="dotted") +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(n.breaks = 10) + 
  scale_color_manual(name = "kinships", values = c("mean" = "green", "half-sib" = "orange", "full-sib" = "skyblue")) 


rel.hist

# #geom_density(alpha=.2, fill="#FF6666") Keep for a rainy day.
colo <- viridisLite::viridis(n = 20)
heatmap.2(daly.rel$rel, scale = "none", col = colo, 
          trace = "none", density.info = "none", 
          main = "Heatmap of relatedness (θ) in sawfish from the Daly River", 
          dendrogram = c("none"))


### Lets reload our related individuals updated with names and make the heatmap...

sibs  <- read.csv("Daly_Sibs.csv")
meta <- read.csv("Daly_meta.csv")
sibs

meta <- meta[meta$id %in% c(sibs$id1, sibs$id2),]
data <- sibs
kinNWdata <- data %>% 
  dplyr::mutate(Cohort_gap = Cohort1 - Cohort2) %>%
  dplyr::select(id1, id2, sibs, Cohort_gap)

#Calculate the gap between cohorts, and add it into a new dataframe

kinNWdata$Cohort_gap[kinNWdata$Cohort_gap<0] <- 1 # this replaces all of our "negative" cohort gaps so that we have a tidy plot

kinNWdata
#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE) 

df <- data.frame(id = igraph::V(network)$name)

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, sex, cohort)
vertices <- as_data_frame(vertices, what = "vertices")


network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')
attributes(layout)

kin_network1 <- ggraph::ggraph(network, layout = layout) + 
  ggraph::geom_edge_link( 
    aes(width = Cohort_gap,
        edge_colour = factor(sibs)),
    #arrow = arrow(length = unit(3, 'mm')), 
    end_cap = ggraph::circle(2, 'mm'),
    edge_alpha = 1) +
  ggraph::scale_edge_width(range = c(1,2), breaks = c(0,1), name = "Cohort Gap") +
  ggraph::scale_edge_colour_manual(values = c("orange", "skyblue"),
                                   name = "Kin-Type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  ggraph::geom_node_point(aes(shape = sex),
                          size = 4) +
  ggraph::geom_node_text( aes(label = df$id), repel = TRUE, 
                          size = 5, color = "black") +
  ggplot2::scale_color_manual(values = adegenet::funky(9), na.value = "grey50") +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "right",
    plot.margin = unit(rep(1,4), "cm")
  ) 


print(kin_network1)

?scale_edge_width
