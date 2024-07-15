generateContiVars <- function(distr = "rnorm", center = NA, var = NA, 
                              minimum = NA, maximum = NA,
                              size){ 
  # generate continuous variables using given parameters.
  # One can add some constraints to the max. and min. of the data.
  
  variable_vec <- switch(distr,
     "rnorm" = rnorm(n = size, mean = center, sd = sqrt(var)), # only center and var
     "rtruncnorm" = rtruncnorm(n = size, a = minimum, b = maximum, mean = center, sd = sqrt(var)),
                         # center, var, max, min
     "runif" = runif(n = size, min = minimum, max = maximum), # only max and min
     "rbeta" = {
       rbetaGen <- function(){ # max, min, center, var. will only be centered approximately.
         total_range <- maximum - minimum 
         scaled_mean <- (center - minimum) / total_range
         scaled_var <- var / (total_range ^ 2)
         
         alpha <- ((1 - scaled_mean) / scaled_var - 1 / scaled_mean) * scaled_mean ^ 2
         beta <- alpha * (1 / scaled_mean - 1)
         beta_distribution <- rbeta(size, alpha, beta)
         return(minimum + beta_distribution * total_range)
       }
       rbetaGen()},
     "rlnorm" = { # only center and var
                           rlnormGen <- function(){
                             sdlog <- sqrt(log(var/center^2 + 1))
                             meanlog <- log(center) - 0.5 * sdlog^2
                             lognormal_distribution <- rlnorm(size, meanlog, sdlog)
                             return(lognormal_distribution)
                           }
                           rlnormGen()
                         },
     "rgamma" = rgamma(size, shape = (center^2) / var, scale = var / center), # only center and var
     "rchisq" = rchisq(size,df = center), # only center
     stop("???")
  )
  return(variable_vec)
}