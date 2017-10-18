
update_mu <- function(A, B=0, mu=0, mu_previous) {
  return(A%*%mu_previous)
}

update_sigma <- function(A, R, sigma_previous) {
  return(A%*%sigma_previous%*%t(A) + R)
}

calculate_kalman_gain <- function(sigma, C, Q) {
  sigma%*%t(C) %*% solve(C%*%sigma%*%t(C) + Q)
}

scale_mu_with_kalman_gain <- function(mu_bar, K, z, C) {
  mu_bar + K%*%(z-C%*%mu_bar)
}

scale_sigma_with_kalman_gain <- function(K, C, sigma_bar) {
  dimensions <- dim(K)
  I <- diag(1, nrow=dim(K)[1], ncol=dim(C)[2])
  
  return((I-K%*%C)%*%sigma_bar)
}

kalman_filter <- function(A, B, C, R, Q, mu_0, sigma_0, z) {
  n <- length(z)
  mu_list <- list()
  sigma_list <- list()
  
  mu_list[[1]] <- mu_0
  sigma_list[[1]] <- sigma_0
  
  for (t in 2:(n+1)) {
    mu_bar <- update_mu(A=A, mu_previous=mu_list[[t-1]])
    sigma_bar <- update_sigma(A=A, R=R, sigma_list[[t-1]])
    K <- calculate_kalman_gain(sigma=sigma_bar, C=C, Q=Q)
    mu_list[[t]] <- scale_mu_with_kalman_gain(mu_bar, K, z[t-1], C)
    sigma_list[[t]] <- scale_sigma_with_kalman_gain(K, C, sigma_bar)
  }
  return(structure(list(
    mu <- mu_list,
    sigma <- sigma_list
  )))
}


# Test
A <- matrix(c(1, 1, 0, 1), byrow=TRUE, nrow=2)
B <- matrix(c(0, 0, 0, 0), byrow=TRUE, nrow=2)
C <- matrix(c(1, 0), byrow=TRUE, nrow=1)
R <- matrix(c(0.035, 0.035+3.06*10^(-12), 3.06*10^(-12), 3.06*10^(-12)), byrow=TRUE, nrow=2)
Q <- 0.035
mu_0 <- matrix(c(10, 0), byrow=TRUE, nrow=2)
sigma_0 <- matrix(c(10^2, 0, 0, 10^2), byrow=TRUE, nrow=2)

load("../data/Radiation_data.Rda")
z <- Radiation_data$dose

debugonce(kalman_filter)
kalman_values <- kalman_filter(A, B, C, R, Q, mu_0, sigma_0, z)
kalman_mean <- unlist(lapply(kalman_values[[1]], function(x) x[[1]]))

plot(Radiation_data$dose, pch=16)
lines(kalman_mean, col="red")
