
combine_nu_wiener <- function(nu, wiener) {
  # Purpose: combines bounded wiener process parameters into a single matrix
  # Inputs: double vector drifts, double matrix wiener
  # Outputs: double matrix nu_wiener
  
  n_combos <- dim(nu)[2] * 2
  wiener_index <- rep(c("acc", "spd"), each = n_combos / 2)
  nu_wiener <- matrix(0, nrow = n_combos, ncol = 7)
  nu <- c(as.matrix(nu["acc", ]), as.matrix(nu["spd", ]))
  for (item in seq_len(n_combos)) {
    nu_wiener[item, ] <- c(nu[item], as.matrix(wiener[wiener_index[item], ]))
  }
  nu_wiener
}
  
combine_wiener_copula <- function(nu_wiener, rho, omega) {
  # Purpose: combines all wiener process parameters with copula parameters
  # Inputs: double matrix nu_wiener, double vector rho, double scalar omega
  # Outputs: double matrix all_params
  
  n_corrs <- dim(rho)[1]
  n_wiener <- dim(nu_wiener)[1]
  n_combos <- n_wiener * n_corrs
  all_params <- cbind(matrix(rep(t(nu_wiener), times = n_corrs), 
                             ncol = 7, nrow = n_combos, byrow = TRUE),
                      as.matrix(rho[rep(seq_len(n_corrs), each = n_wiener), ]),
                      rep(omega, times = n_combos))
  all_params
}

combine_param <- function(nu, wiener, rho, omega) {
  # Purpose: wraps both combine functions and cleans up the output
  # Inputs: double vector nu, double matrix wiener,
  #        double vector rho, double scalr omega
  # Outputs: double matrix all_params
  
  nu_wiener <- combine_nu_wiener(nu = nu, wiener = wiener)
  all_params <- as.data.frame(combine_wiener_copula(nu_wiener = nu_wiener,
                                                    rho = rho, omega = omega))
  colnames(all_params) <- c("nu", colnames(wiener), colnames(rho), "omega")
  rownames(all_params) <- seq_len(dim(all_params)[1])
  all_params <- all_params[, c(2, 1, seq(3, 11))]
  all_params
}






