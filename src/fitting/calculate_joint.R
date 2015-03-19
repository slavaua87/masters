setwd("~/Dropbox/Slava/Masters/")

source("src/fitting/integrate_pdf.R")
source("src/fitting/calculate_weibull.R")

dtmvt <- tmvtnorm::dtmvt

lkj_logdensity <- function(rho, omega) {
  logkernel <- (1 - sum(rho ^ 2) + 2 * prod(rho)) ^ (omega - 1)
  return(logkernel)
}

joint_logdensity <- function(behav_data, theta, model) {
  if (model == "independent") {
    group_logdensity <- log(c(dunif(x = theta[1:2], mix = 0.001, max = 0.786),
                              dunif(x = theta[3:6], mix = 0, max = 1.172),
                              dunif(x = theta[7:8], mix = 0, max = 50),
                              dunif(x = theta[9:10], mix = 0, max = 5),
                              dunif(x = theta[11:12], mix = 0, max = 1),
                              dunif(x = theta[13:14], mix = 0, max = 1.884),
                              dtmvt(x = theta[15:19], sigma = as.matrix(25),
                                    df = 1, lower = 0, upper = Inf), 
                              dunif(x = theta[20], mix = 0, max = .5),
                              dtmvt(x = theta[21], sigma = as.matrix(1), df = 1,
                                    lower = 0, upper = Inf)))
    person_logdensity <- log(c(
      dgamma(x = theta[22:24], 
             shape = (theta[1] / theta[15]) ^ 2, 
             scale = theta[15] ^ 2 / theta[1]),  
      dgamma(x = theta[25:27],            
             shape = (theta[2] / theta[15]) ^ 2,             
             scale = theta[15] ^ 2 / theta[2]),         
      dgamma(x = theta[28:30],        
             shape = (theta[3] / theta[16]) ^ 2, 
             scale = theta[16] ^ 2 / theta[3]),
      dgamma(x = theta[31:33],               
             shape = (theta[4] / theta[16]) ^ 2,   
             scale = theta[16] ^ 2 / theta[4]),         
      dgamma(x = theta[34:36],               
             shape = (theta[4] / theta[16]) ^ 2,                
             scale = theta[16] ^ 2 / theta[4]),         
      dgamma(x = theta[37:39],             
             shape = (theta[6] / theta[17]) ^ 2, 
             scale = theta[17] ^ 2 / theta[6]),       
      dgamma(x = theta[40:42],               
             shape = (theta[7] /theta[18]) ^ 2,   
             scale = theta[18] ^ 2 / theta[7]),         
      dgamma(x = theta[43:45],               
             shape = (theta[8] / theta[18]) ^ 2, 
             scale = theta[18] ^ 2 / theta[8]),         
      dgamma(x = theta[46:48],
             shape = (theta[9] / theta[19]) ^ 2,             
             scale = theta[19] ^ 2 / theta[9]),         
      dgamma(x = theta[49:51],               
             shape = (theta[10] / theta[19]) ^ 2,   
             scale = theta[19] ^ 2 / theta[10]),              
      dbeta(x = theta[52:54],               
            shape1 = shape1(theta[11], theta[20]),
            shape2 = shape2(theta[11], theta[20])),
      dbeta(x = theta[55:57],               
            shape1 = shape1(theta[12], theta[20]),
            shape2 = shape2(theta[12], theta[20])),
      dgamma(x = theta[58:60],               
             shape = (theta[13] / theta[21]) ^ 2, 
             scale = theta[21] ^ 2 / theta[13]),       
      dgamma(x = theta[61:63],               
             shape = (theta[14] / theta[21]) ^ 2, 
             scale = theta[21] ^ 2 / theta[14]),             
      dunif(x = theta[64:66], mix = 0, max = 0.658),
      dunif(x = theta[67:69], mix = 0, max = 0.860),
      dunif(x = theta[70:72], mix = 0, max = 1.260)))
  }
  if (model == "normal") {
    group_logdensity <- log(c(dunif(x = theta[1:2], mix = 0.001, max = 0.786),
                              dunif(x = theta[3:6], mix = 0, max = 1.172),
                              dunif(x = theta[7:8], mix = 0, max = 50),
                              dunif(x = theta[9:10], mix = 0, max = 5),
                              dunif(x = theta[11:12], mix = 0, max = 1),
                              dunif(x = theta[13:14], mix = 0, max = 1.884),
                              dtmvt(x = theta[15:19], sigma = as.matrix(25),
                                    df = 1, lower = 0, upper = Inf), 
                              dunif(x = theta[20], mix = 0, max = .5),
                              dtmvt(x = theta[21], sigma = as.matrix(25), 
                                    df = 1, lower = 0, upper = Inf),
                              rtmvt(x = theta[22], sigma = as.matrix(625),
                                    df = 1, lower = 0, upper = Inf)))
    person_logdensity <- log(c(
      dgamma(x = theta[23:25], 
             shape = (theta[1] / theta[15]) ^ 2, 
             scale = theta[15] ^ 2 / theta[1]),  
      dgamma(x = theta[26:28],            
             shape = (theta[2] / theta[15]) ^ 2,             
             scale = theta[15] ^ 2 / theta[2]),         
      dgamma(x = theta[29:31],        
             shape = (theta[3] / theta[16]) ^ 2, 
             scale = theta[16] ^ 2 / theta[3]),
      dgamma(x = theta[32:34],               
             shape = (theta[4] / theta[16]) ^ 2,   
             scale = theta[16] ^ 2 / theta[4]),         
      dgamma(x = theta[35:37],               
             shape = (theta[4] / theta[16]) ^ 2,                
             scale = theta[16] ^ 2 / theta[4]),         
      dgamma(x = theta[38:40],             
             shape = (theta[6] / theta[17]) ^ 2, 
             scale = theta[17] ^ 2 / theta[6]),       
      dgamma(x = theta[41:43],               
             shape = (theta[7] /theta[18]) ^ 2,   
             scale = theta[18] ^ 2 / theta[7]),         
      dgamma(x = theta[44:46],               
             shape = (theta[8] / theta[18]) ^ 2, 
             scale = theta[18] ^ 2 / theta[8]),         
      dgamma(x = theta[47:49],
             shape = (theta[9] / theta[19]) ^ 2,             
             scale = theta[19] ^ 2 / theta[9]),         
      dgamma(x = theta[50:52],               
             shape = (theta[10] / theta[19]) ^ 2,   
             scale = theta[19] ^ 2 / theta[10]),              
      dbeta(x = theta[53:55],               
            shape1 = shape1(theta[11], theta[20]),
            shape2 = shape2(theta[11], theta[20])),
      dbeta(x = theta[56:58],               
            shape1 = shape1(theta[12], theta[20]),
            shape2 = shape2(theta[12], theta[20])),
      dgamma(x = theta[59:61],               
             shape = (theta[13] / theta[21]) ^ 2, 
             scale = theta[21] ^ 2 / theta[13]),       
      dgamma(x = theta[62:64],               
             shape = (theta[14] / theta[21]) ^ 2, 
             scale = theta[21] ^ 2 / theta[14]),             
      dunif(x = theta[65:67], mix = 0, max = 0.658),
      dunif(x = theta[68:70], mix = 0, max = 0.860),
      dunif(x = theta[71:73], mix = 0, max = 1.260),
      lkj_logkernel(rho = theta[74:76], omega = theta[22]),
      lkj_logkernel(rho = theta[77:79], omega = theta[22]),
      lkj_logkernel(rho = theta[80:82], omega = theta[22])))
    
    theta <- inits[c(seq(23, 64, 3)[seq(1, 14, 2)], seq(65, 73, 3), 74:76)]
    theta <- c(theta[1], weibull(17, theta[2], theta[3], theta[4], theta[5]),
               theta[8], theta[6], theta[9], theta[7], theta[10], theta[11:13])
    
    # make sure length is flexible to account for subsampling
    dat_mat <- group_by(obs, subj) %>% 
        transmute(rt = rt, choice = resp, alpha = , nu = , eta = ,
                  lambda = , gamma = , chi = , phi = , rho = )
    })
    
    
    
    c(seq(24, 73, 3), 77:79)
    c(seq(25, 73, 3), 80:82)
    
    group_by(rr_data, )
    behav_logdensity <- log(c(
      integrate_density_vec(rt = , choice = , sigma = .1, 
                            params = params, model = model, tol = 1e-2)))
                              
    # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
    # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt
    
  }
}
  