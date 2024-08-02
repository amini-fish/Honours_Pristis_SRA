setwd("C:/Users/samue/Desktop/Honours/analysis/sims")


sim_data <- readLines("sims.ibd9")

strt <- which(grepl("^IBD", sim_data)) + 2
stp <- which(grepl("Indiv genotypes", sim_data)) - 4
linez_headings <- sim_data[strt]
linez_data <- sim_data[(strt + 1):stp]
tmp_headings <- unlist(stringr::str_split(linez_headings, " "))
tmp_data <- stringr::str_split(linez_data, " ")

#Raw data 
tmp_data_raw_1 <- lapply(tmp_data, "[", c(2:22))
tmp_data_raw_2 <- do.call("rbind", tmp_data_raw_1)
tmp_data_raw_3 <- as.data.frame(tmp_data_raw_2)

tmp_data_raw_3

tmp_data_raw_3$V3 <- lapply(tmp_data_raw_3$V3, as.numeric)

tmp_data_raw_3


colnames(tmp_data_raw_3) <- tmp_headings[2:22]

tmp_data_raw_3

df <- data.frame(ind1=tmp_data_raw_3$Indiv1, ind2=tmp_data_raw_3$Indiv2,rel= tmp_data_raw_3$`r(1,2)`)

df$rel <- as.numeric(df$rel); df

#Relatedness

res <- matrix(df$rel, nrow = 380, ncol = 380)

for (i in 1:nrow(df)) {
  res[df[i, 1], df[i, 2]] <- df[i, 3]
}


View(res)

ggplot(data = df, aes(x = rel))+
  geom_histogram(bins = 500)+
  xlim(0, 0.25)## this is inclusive of self comparisons


gl.plot.heatmap(res)

colnames(res) <- tmp_data_raw_3$Indiv1
rownames(res) <- tmp_data_raw_3$Indiv1

df
summary(df)
