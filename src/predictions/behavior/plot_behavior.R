

plot_behavior_sum <- function(behavior_sum, rows = 3, cols = 2) {
  # Plots a grid of mean paths of 3 models for nrow * ncol conditions
  # Takes list of lists of numeric vectors, integer scalars and saves an image
  
  behavior_sum %>% mutate(graph = cut(condition, 
                                      breaks = 3,
                                      labels = FALSE)) %>%
    group_by(graph) %>% 
    mutate(subgraph = cut(condition, breaks = 6, labels = FALSE)) %>% 
    do({        
        qp_graph <- ggplot(data = .,
                           mapping = aes(x = prob, y = quant, 
                                         group = interaction(quant_rank,
                                                             subgraph,
                                                             model),
                                         color = model)) +
          facet_wrap(facets = ~subgraph, nrow = 3, ncol = 2, 
                     scales = c("free"), as.table = TRUE) +
          geom_line(alpha = 1, size = .2) +
          geom_point(size = 1) +
          theme_solarized_2() + 
          ylab("Reaction Time (sec)") +
          xlab("Choice Probability")
                    
        ggsave(qp_graph, 
               filename = paste0("results/behavior/qp-plot-corr",
                                 unique(.$graph), "-", Sys.Date(), ".pdf"))
  })
}


