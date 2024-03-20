setwd("C:/Users/samue/Desktop/Honours - Sawfish//pristis data for Sam")


meta2 <- readr::read_csv("Daly_Sibs.csv")
cohorts <- readr::read_csv("Glyphis garricki_Metadata_with_cohort.csv")
meta <- meta[meta$id %in% c(pwdata2$iname, pwdata2$jname),]
data <- pwdata2
kinNWdata <- data %>% 
  dplyr::mutate(Cohort_gap = Cohort_j - Cohort_i) %>%
  dplyr::select(iname, jname, relation, Cohort_gap)
network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE) 


df <- data.frame(id = igraph::V(network)$name)
vertices <- dplyr::left_join(df, meta, by = "id") %>%
  dplyr::select(id,sample_location, sex, HSP, yearcollect, Age, 
                Cohort2, HT)
vertices$HT[vertices$HT == ""] <- NA
vertices <- as_data_frame(vertices, what = "vertices")


network <- igraph::graph_from_data_frame(d = kinNWdata, directed = TRUE,
                                         vertices = vertices ) 

layout <- ggraph::create_layout(network, layout = 'igraph', 
                                circular = FALSE, algorithm = 'fr')
attributes(layout)

kin_network1 <- ggraph::ggraph(network, layout = layout) + 
  ggraph::geom_edge_link( 
    aes(width = Cohort_gap,
        edge_colour = factor(relation)), 
    # arrow = arrow(length = unit(3, 'mm')), 
    #                end_cap = ggraph::circle(2, 'mm'),
    edge_alpha = 1) +
  ggraph::scale_edge_width(range = c(0.5, 3), name = "Cohort_gap") +
  ggraph::scale_edge_colour_manual(values = c("blue", "black", "red"),
                                   name = "kin-type",
                                   aesthetics = "edge_colour",
                                   na.value = "grey50") +
  ggraph::geom_node_point(aes(color = sample_location),
                          size = 4) +
  ggraph::geom_node_text( aes(label = HT), repel = TRUE, 
                          size = 5, color = "black") +
  ggplot2::scale_color_manual(values = adegenet::funky(9), na.value = "grey50") +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "right",
    plot.margin = unit(rep(1,4), "cm")
  ) 
print(kin_network1)

ggplot2::ggsave(filename = "0.Kinship_network_all_relationships2.png",
                bg = "white", plot = kin_network, device = "png", 
                width = 40, height = 30, units = "cm")