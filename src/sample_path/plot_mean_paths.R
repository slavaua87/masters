



process_paths <- function(paths) {
  # Purpose: reshape simulation output, save results and then return results
  # Inputs: list of lists of double vectors paths
  # Outputs: data_frame model_paths
  
  ind_param <- combine_param(nu = nu, wiener = wiener,
                             rho = rho, omega = omega) %>%
    slice(1:72) %>% 
    mutate(condition = 1:72)
  acc_index <- filter(ind_param, alpha == 0.221) %>% dplyr::select(condition)
  
  alpha <- c(0, rev(unique(ind_param$alpha)))
  thresholds <- data_frame(upper = rep(unique(ind_param$alpha),
                                       each = 6),
                           lower = 0, 
                           plots = seq_len(72),
                           graphs = rep(1:6, each = 12))
  model_paths <- paths %>% 
    melt %>%
    rename(value = value, sum_stat = L4, response = L3,
           condition = L2, model = L1) %>%
    mutate(graph = cut(x = condition, labels = FALSE, 
                       breaks = 6)) %>%
    group_by(sum_stat, response, condition, model) %>% 
    mutate(time = seq(from = 1, by = 1, length.out = n())) %>%
    group_by(condition) %>%
    mutate(instruction = rep(unique(condition) %in% acc_index$condition + 1,
                             length(condition))) %>% ungroup
  write.table(model_paths, "results/sample_path/plot_values.txt",
              row.names = FALSE)
}


plot_mean_paths <- function(paths, cores = 1, rows = 6, cols = 2) {
  # Purpose: plots a grid of mean paths of 3 models for nrow * ncol conditions
  # Inputs: list of lists of double vectors paths, 
  #         integer scalars cores, rows, cols
  # Output: pdf images
  
  library("reshape2")
  library("ggplot2")
  library("ggthemes")
  source("src/sample_path/load_dependencies.R")
  
  process_paths(paths) %>% 
    group_by(instruction) %>% 
    do({cutoff <- dplyr::select(., instruction) %>% distinct %>% 
          unlist %>% '['(c(250, 3000), .)
        filter(., time < cutoff)
    }) %>% 
    group_by(graph) %>% 
    do({graph_index <- dplyr::select(., graph) %>% 
          distinct %>%
          as.numeric
        plot_order <- group_by(., instruction) %>%
          dplyr::select(condition) %>%
          unique
        plot_order <- ungroup(plot_order) %>%
          dplyr::select(condition) %>% 
          unlist %>%
          matrix(nrow = 2, ncol = 6, byrow = TRUE) %>%
          c
        reordered <- mutate(., plots = factor(x = .$condition, 
                            levels = as.character(plot_order)))
        
        ts.graph <- ggplot() +
                    facet_wrap(facets = ~plots, nrow = 6, ncol = 2, 
                               scales = c("free"), as.table = TRUE) +
                    geom_line(aes(x = time, y = value, 
                                  group = interaction(response, plots,
                                                      model, sum_stat),
                                  color = model),
                              data = filter(reordered, 
                                            sum_stat %in% c("fast_mean",
                                                            "middle_mean",
                                                            "slow_mean")),
                              alpha = 1, size = .2) +                    
                   theme_solarized_2() + 
                   ylab("Relative Evidence") + 
                   xlab("Time (ms)") 
        ggsave(ts.graph, 
               filename = paste0("results/sample_path/path-plot-corr",
                                 unique(.$graph), "-", Sys.Date(), ".pdf"))
  })
}


