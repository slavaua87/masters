

plot_mean_paths <- function(paths, cores = 1, rows = 6, 
                            cols = 2) {
  # Plots a grid of mean paths of 3 models for nrow * ncol conditions
  # Takes list of lists of numeric vectors, integer scalars and saves an image
  library(package = "magrittr")
  library(package = "reshape2")
  library(package = "dplyr")
  library(package = "ggplot2")
  library(package = "ggthemes")
  library(package = "doParallel")
  source("src/predictions/sample_path/calculate_path_stats.R")
  source("src/predictions/sample_path/wiener_parameters.R")
  source("src/predictions/sample_path/combine_parameters.R")
  
  ind_param <- combine_param(nu = nu, wiener = wiener,
                             rho = rho, omega = omega) %>%
               slice(1:72) %>% mutate(condition = 1:72)
  acc_index <- filter(ind_param, alpha == 0.221) %>% select(condition)
  
  bounds <- c(0, rev(unique(ind_param$alpha)))
  thresholds <- data.frame(upper = rep(unique(ind_param$alpha),
                                       each = 6),
                           lower = 0, 
                           plots = seq_len(72),
                           graphs = rep(1:6, each = 12))
  
  model_paths <- calc_paths_stats(paths, alpha = bounds, cores = cores) %>% 
    melt %>% 
    rename(value = value, sum_stat = L4, response = L3,
           condition = L2, model = L1) %>%
    mutate(graph = cut(x = condition, labels = FALSE, 
                       breaks = 6)) %>%
    group_by(sum_stat, response, condition, model) %>% 
    mutate(time = seq(from = 1, by = 1, length.out = n())) %>%
    group_by(condition) %>%
    mutate(instruction = rep(unique(condition) %in% acc_index$condition + 1,
                                    length(condition))) %>% ungroup %>%
    group_by(instruction) %>% 
    do({cutoff <- select(., instruction) %>% distinct %>% 
          unlist %>% '['(c(150, 1000), .)
        filter(., time < cutoff)
    })
  
  model_paths %>% group_by(graph) %>% 
    do({graph_index <- select(., graph) %>% distinct %>% as.numeric
        #y_limits <- filter(thresholds, graphs == graph_index) 
        plot_order <- group_by(., instruction) %>% select(condition) %>% unique
        plot_order <- ungroup(plot_order) %>% select(condition) %>% unlist %>%
          matrix(nrow = 2, ncol = 6, byrow = TRUE) %>% c
        reordered <- mutate(., plots = factor(x = .$condition, 
                            levels = as.character(plot_order)))
                
        conf_int <- group_by(reordered, response, plots, model) %>%  
                    do({means <- filter(., sum_stat == "mean_path") %>% 
                                 select(value) %>% unlist
                        ci <- filter(., sum_stat == "path_sd") %>% 
                              select(value) %>% 
                              unlist %>% 
                              multiply_by(2)
                        ci <- ci %>% divide_by(sqrt(length(ci)))
                        data.frame(means = means, ci = ci) %>%
                        mutate(time = seq(from = 1, by = 1, length.out = n()))
                    }) %>%           
                    ungroup

        #y_limits <- y_limits[match(plot_order, unlist(y_limits$plots)), ]
        
        ts.graph <- ggplot() +
                    facet_wrap(facets = ~plots, nrow = 6, ncol = 2, 
                               scales = c("free"), as.table = TRUE) +
                    geom_line(aes(x = time, y = value, 
                                  group = interaction(response, plots, model),
                                  color = model),
                              data = filter(reordered, sum_stat == "mean_path"),
                              alpha = 1, size = .2) +                    
          theme_fivethirtyeight() + ylab("Evidence") + xlab("Time (ms)") +
                    geom_ribbon(aes(x = time, ymin = means - ci, 
                                    ymax = means + ci, 
                                    fill = model, 
                                    group = interaction(response, plots, model)),
                                alpha = 0.5,
                                data = conf_int)
                    
        ggsave(ts.graph, 
               filename = paste0("results/sample_path/path-plot-corr",
                                 unique(.$graph), "-", Sys.Date(), ".pdf"))
  })
}


