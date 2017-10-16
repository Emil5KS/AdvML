# Define the transition, emission and initialization probabilities --------

emission_probs <- matrix(c(.2, .2, .2, 0, 0, 0, 0, 0, .2, .2,
                           .2, .2, .2, .2, 0, 0, 0, 0, 0, .2,
                           .2, .2, .2, .2, .2, 0, 0, 0, 0, 0,
                           0, .2, .2, .2, .2, .2, 0, 0, 0, 0,
                           0, 0, .2, .2, .2, .2, .2, 0, 0, 0,
                           0, 0, 0, .2, .2, .2, .2, .2, 0, 0,
                           0, 0, 0, 0, .2, .2, .2, .2, .2, 0,
                           0, 0, 0, 0, 0, .2, .2, .2, .2, .2,
                           .2, 0, 0, 0, 0, 0, .2, .2, .2, .2,
                           .2, .2, 0, 0, 0, 0, 0, .2, .2, .2), byrow=TRUE, nrow=10)

transition_probs <- matrix(c(.5, .5, 0, 0, 0, 0, 0, 0, 0, 0,
                             0, .5, .5, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, .5, .5, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, .5, .5, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, .5, .5, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, .5, .5, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, .5, .5, 0, 0,
                             0, 0, 0, 0, 0, 0, 0, .5, .5, 0,
                             0, 0, 0, 0, 0, 0, 0, 0, .5, .5,
                             .5, 0, 0, 0, 0, 0, 0, 0, 0, .5), byrow=TRUE, nrow=10)

tProbDensity <- function(zt, zt_1) {
  return(transition_probs[zt_1, zt])
}

eProbDensity <- function(xt, zt) {
  return(emission_probs[zt, xt])
}

initProbDensity <- function(z0) {
  return(dunif(z0, min=1, max=10))
}


# Simulate data -----------------------------------------------------------

library(HMM)

robotHmm <- HMM::initHMM(
  States = 1:10,
  Symbols = 1:10,
  transProbs = transition_probs,
  emissionProbs = emission_probs
)

simHMM <- function(hmm, length) {
  simulation <- HMM::simHMM(hmm, length)
  return(structure(simulation, class="HmmSimulation"))
}

nSim <- 100
robotSimultation <- simHMM(hmm=robotHmm, length=nSim)

X <- robotSimultation$observation
Z <- robotSimultation$states

# Implement Viterbi -------------------------------------------------------

possibleStates <- 1:10
get_omega <- function(Z, Omega, Z_next, x_next) {
  sapply(Z_next, function(z_next) {
    term1 <- log(eProbDensity(x_next, z_next))
    
    term2 <- sapply(Z, function(z) {
      log(tProbDensity(z_next, z))
    }) + Omega
    
    return(term1+ max(term2))
  })
}

get_phi <- function(Z, Z_next, Omega) {
  sapply(Z_next, function(z_next) {
    term <- sapply(Z, function(z) {
      log(tProbDensity(z_next, z))
    }) + Omega
    return(Z[which.max(term)])
  }) 
}

viterbi <- function(observations, possibleStates) {
  cardinality <- length(possibleStates)
  t_total <- length(observations)
  
  omega_0 <- vector("numeric", length = cardinality)
  for (i in 1:cardinality) {
    omega_0[i] <- log(initProbDensity(possibleStates[i])) + log(eProbDensity(observations[1], possibleStates[i]))
  }
  
  
  omega <- matrix(NA, nrow=t_total, ncol=cardinality)
  phi <- matrix(NA, nrow=t_total, ncol=cardinality)
  omega[1, ] <- omega_0
  
  for (i in 1:(t_total-1)) {
    omega[i+1, ] <- get_omega(possibleStates, omega[i, ], possibleStates, observations[i+1])
    phi[i+1, ] <- get_phi(possibleStates, possibleStates, omega[i, ])
  }
  
  mpp <- rep(NA, t_total)
  mpp[t_total] <- possibleStates[which.max(omega[t_total, ])]
  for (t in (t_total - 1):1) {
    mpp[t] <- phi[t + 1, possibleStates[mpp[t + 1]] == possibleStates]
  }
  
  return(list(path = mpp, omega = omega, phi = phi))
  
}

results <- viterbi(X, possibleStates)
results$path

results_HMM <- HMM::viterbi(robotHmm, X)

cbind(results$path, results_HMM)
