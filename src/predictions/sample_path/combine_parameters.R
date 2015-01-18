
### Combines bounded wiener process parameters into a single data.frame
combine_nu_wiener <- function(nu, wiener) {
  # Combines drifts to other wiener process parameters
  # Takes 2 data.frames and returns a matrix
  n_combos <- dim(nu)[2] * 2
  wiener_index <- rep(c("acc", "spd"), each = n_combos / 2)
  nu_wiener <- matrix(0, nrow = n_combos, ncol = 7)
  nu <- c(as.matrix(nu["acc", ]), as.matrix(nu["spd", ]))
  for (item in seq_len(n_combos)) {
    nu_wiener[item, ] <- c(nu[item], as.matrix(wiener[wiener_index[item], ]))
  }
  return(nu_wiener)
}
  
combine_wiener_copula <- function(nu_wiener, rho, omega) {
  # Combines wiener process parameters with copula parameters
  # Takes a matrix, a data.frame and a scalar to return a matrix
  n_corrs <- dim(rho)[1]
  n_wiener <- dim(nu_wiener)[1]
  n_combos <- n_wiener * n_corrs
  all_params <- cbind(matrix(rep(t(nu_wiener), n_corrs), 
                             ncol = 7, nrow = n_combos, byrow = TRUE),
                      as.matrix(rho[rep(seq_len(n_corrs), each = n_wiener), ]),
                      rep(omega, times = n_combos))
  return(all_params)
}

simul_conditions <- function(nu, wiener, rho, omega) {
  # Wraps both combine functions and cleans up the output
  # Takes and outputs data.frames
  nu_wiener <- combine_nu_wiener(nu = nu, wiener = wiener)
  all_params <- as.data.frame(combine_wiener_copula(nu_wiener = nu_wiener,
                                                    rho = rho, omega = omega))
  colnames(all_params) <- c("nu", colnames(wiener), colnames(rho), "omega")
  rownames(all_params) <- seq_len(dim(all_params)[1])
  all_params <- all_params[, c(2, 1, seq(3, 11))]
  return(all_params)
}






