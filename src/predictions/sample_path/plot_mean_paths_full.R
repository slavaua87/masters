

plot_mean_paths <- function(paths, cores = 1, nrow = 6, 
                            ncol = 4, nplots = 22, speeds = 2) {
  # Plots a grid of mean paths of 3 models for nrow * ncol conditions
  # Takes list of lists of numeric vectors, integer scalars and saves an image
  library(package = "magrittr")
  library(package = "reshape2")
  library(package = "dplyr")
  library(package = "ggplot2")
  source("src/predictions/sample_path/calculate_mean_path_full.R")
  source("src/predictions/sample_path/wiener_parameters.R")
  source("src/predictions/sample_path/combine_parameters.R")
  
  ind_param <- combine_param(nu = nu, wiener = wiener,
                             rho = rho, omega = omega)
   
  thresholds <- data.frame(upper = rep(unique(ind_param$alpha),
                                       each = nplots / speeds),
                          lower = 0, 
                          condition = seq_len(9 * nplots),
                          graphs = rep(1:9, each = nplots))
  
  calc_model_paths(paths, cores = cores) %>% melt %>% 
    rename(state = value, condition = L2, model = L1) %>%
    mutate(graph = cut(x = condition, labels = FALSE, 
                       breaks = max(condition) / nplots)) %>%
    group_by(condition, model) %>% 
    mutate(time = seq(from = 1, by = 1, length.out = n())) %>% 
    group_by(graph) %>% 
    do({graph_index <- select(., graph) %>% unique %>% as.data.frame(.)
        graph_index <- as.numeric(graph_index)
        alphas <- filter(thresholds, graphs == graph_index)
        
        ts.graph <- ggplot(data = .) +
                    aes(x = time, y = state, 
                        group = interaction(condition, model),
                        color = model) + 
                    facet_wrap(facets = ~condition, nrow = nrow, ncol = ncol, 
                               scales = c("free")) + 
                    geom_line() +
                    geom_hline(mapping = aes(yintercept = c(lower, upper)),
                               data = alphas) +
                    theme_bw() + ylab("Evidence") + xlab("Time (ms)")
        ggsave(ts.graph, 
               filename = paste0("results/sample_path/path-plot-corr",
                                 unique(.$graph), "-", Sys.Date(), ".pdf"))
  })
}


