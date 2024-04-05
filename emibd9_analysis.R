setwd("C:/Users/samue/Desktop/Honours_Sawfish/analysis")

### LOAD REQUIRED PACKAGES ###

#install.packages("hierfstat")
#install.packages("graph4lg")
devtools::install_version("ggplot2", "3.4.4")
library(dartRverse)
library(ggplot2)
library(hierfstat)
library(dplyr)
library(devtools)
library(dartR.sexlinked)
library(dartR.captive)
library(gplots)
library(graph4lg)

################################################################################

### LOAD IN CLEAN GENOTYPE DATA ###

gl <- get(load("C:/Users/samue/Desktop/Honours_Sawfish/analysis/daly_geno_clean.Rdata")); gl

## Run our analysis - using EMIBD9 implementation 

daly.rel <- gl.run.EMIBD9(gl, 
                          Inbreed = 1, 
                          emibd9.path =  "C:/EMIBD9")

## Lets extract the relatedness data from our output file

diag(daly.rel$rel) <- 0 # removes self comparisons hooray

emibd.rel <- pw_mat_to_df(daly.rel$rel) # convert the square matrix to an edge based dataframe

emibd.rel <- emibd.rel[,c(1,2,4)] # remove your third column (not needed and is confusing)

## Check your results: 

View(emibd.rel)

# A neat bit of code to kick out self comparisons

#ibd9Tab <- ibd9Tab[ibd9Tab$Indiv1 != ibd9Tab$Indiv2,  c(1, 2, 21)]; ibd9Tab

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
emibd.sibs <- ifelse(emibd.rel$value >= 0.092 & emibd.rel$value <= 0.158, 
               yes = "hsp", 
               no = ifelse(emibd.rel$value >=0.204 & emibd.rel$value <= 0.296, 
                           yes = "fsp",
                           "unelated"))

#Add the sibling assignments to a new data frame 

emibd.results <- data.frame(cbind(emibd.rel, emibd.sibs)); emibd.results

####  Assign kin to sibling network  #####

#hsps - extract

half.sibs <- subset(emibd.results, sibs == "hsp"); half.sibs

#fsps - extract 

full.sibs <- subset(emibd.results, sibs == "fsp"); full.sibs

# ALLLLL stitched together now 

emibd.siblings <- rbind(half.sibs, full.sibs); sibs.all

# Remove duplicated pairs (i.e., Ab & BA)

sibs.all[!duplicated(sibs.all$r.1.2.), ] 
sibs.all

## Calculate our summary statistics for relatedness...

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

## write our sibling results as csv 

#write.csv(emibd.siblings, "emibd_sibs_daly.csv")


sibs  <- read.csv("emibd_sibs_daly.csv")
meta <- read.csv("Daly_meta.csv")

meta <- meta[meta$id %in% c(sibs$id_1, sibs$id_2),]
data <- sibs
kinNWdata <- data %>% 
  dplyr::mutate(Cohort_gap = cohort1 - cohort2) %>%
  dplyr::select(id_1, id_2, relatedness, Cohort_gap)

meta

#Calculate the gap between cohorts, and add it into a new dataframe

kinNWdata$Cohort_gap[kinNWdata$Cohort_gap<0] <- 1 # this replaces all of our "negative" cohort gaps so that we have a tidy plot

kinNWdata

#This makes our data frame from which the pariwise network plot between select individuals will be drawn 

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE) 

df <- data.frame(id = igraph::V(network)$name)
df

vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id, sex, Cohort)
vertices <- as_data_frame(vertices, what = "vertices")

vertices

network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')
attributes(layout)

kin_network1 <- ggraph::ggraph(network, layout = layout) + 
  ggraph::geom_edge_link( 
    aes(width = Cohort_gap,
        edge_colour = factor(relatedness)),
    #arrow = arrow(length = unit(3, 'mm')), 
    #end_cap = ggraph::circle(2, 'mm'),
    edge_alpha = 1) +
  ggraph::scale_edge_width(range = c(1,2), breaks = c(0,1), name = "Cohort Gap") +
  ggraph::scale_edge_colour_manual(values = c("skyblue", "orange"),
                                   name = "Kin Type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  ggraph::geom_node_point(aes(shape = sex),
                          size = 4) +
  ggraph::geom_node_text( aes(label = df$id), repel = TRUE, 
                          size = 4, color = "black") +
  ggplot2::scale_color_manual(values = adegenet::funky(9), na.value = "grey50") +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "right",
    plot.margin = unit(rep(1,4), "cm")) 


print(kin_network1)

