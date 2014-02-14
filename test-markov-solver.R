##########################################################################################
## Copyright (c) 2013 Llorenç Cerdà-Alabern,
## http://personals.ac.upc.edu/llorenc This file is free software: you
## can redistribute it and/or modify it under the terms of the GNU
## Affero Public License as published by the Free Software Foundation,
## either version 3 of the License, or (at your option) any later
## version.  test-markov-solver.R is distributed in the hope that it
## will be useful, but WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU Affero Public License for more details.  You should
## have received a copy of the GNU Affero Public License along with
## test-markov-solver.R.  If not, see <http://www.gnu.org/licenses/>.
##########################################################################################
## Implementation of the algorithm described in the paper:
## Closed Form Transient Solution of Continuous Time Markov Chains Through Uniformization
## Presented at Valuetools'13
## The method is also described in the public technical report:
## Transient Solution of Markov Chains Using the Uniformized Vandermonde Method
## UPC-DAC-RR-XCSD-2010-2
## https://www.ac.upc.edu/app/research-reports/html/research_center_index-XCSD-2010,en.html
## http://www.ac.upc.edu/RR/2010/53.pdf
##########################################################################################
## Markov Transient solution examples.

library(Matrix)

source("oo-markov-undetermined-coefficients.R")

build.ctmc.m.m.1.N.matrix.sparse <- function(N, lambda, mu) {
  return(bandSparse(n=N+1,
                    k=c(-1, 0, 1),
                    diag=list(rep(mu, N),
                      c(-lambda, rep(-lambda-mu, N-1), -mu),
                      rep(lambda, N))))
}

##
## M/M/M/1/N 
##
## Example where the chain converges very fast to stationary regime. The method is very accurate.
lambda <- 1e-3
mu <- 1
rho <- lambda/mu
N <- 500
i <- 1 # initial state.
j <- 1:min(N+1, 9)  # target states
ei.max.mult <- 5    # maximum multiplicity of the eigenvalues (0 for unlimited)

Q <- build.ctmc.m.m.1.N.matrix.sparse(N, lambda, mu)
uc.mat <- UCmatrix(Q)
uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult)

t <- time.log.scale(dec=c(-6,7), num=20)
uc.s$plot(t, logy=TRUE, logx=TRUE)

## Comparison with the R exponential matrix function, expm
tt <- 1 # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## Comparison with the R exponential matrix function, expm
tt <- 10 * N * uc.s$qn # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

# Some UCs information
uc.s
# Some eigenvalues information
uc.s$ei
# The computed UCs
uc.s$coef
# The computed UCs
uc.s$coef
# Single eigenvalues
uc.s$ei$one
# Confluent eigenvalues
uc.s$ei$confl
# Multiplicity of confluent eigenvalues
uc.s$ei$mult

##
## M/M/M/1/N 
##
## Example with initial state equal to a probability vector
lambda <- 1e-3
mu <- 1
rho <- lambda/mu
N <- 500
i <- 1 
i <- c(rep(1/5, 5), rep(0, N+1-5)) # initial state as a probability vector
j <- 1:min(N+1, 9)  # target states
ei.max.mult <- 5    # maximum multiplicity of the eigenvalues (0 for unlimited)

Q <- build.ctmc.m.m.1.N.matrix.sparse(N, lambda, mu)
uc.mat <- UCmatrix(Q)
uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult)

t <- time.log.scale(dec=c(-6,7), num=20)
uc.s$plot(t, logy=TRUE, logx=TRUE)

## Comparison with the R exponential matrix function, expm
tt <- 1 # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- i%*%expm(Q*tt)[,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-as.numeric(p.expm))^2))

## Comparison with the R exponential matrix function, expm
tt <- 10 * N * uc.s$qn # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- i%*%expm(Q*tt)[,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-as.numeric(p.expm))^2))

##
## M/M/M/1/N 
##
## Example where the chain does not converge very fast to stationary regime.
lambda <- 0.5
mu <- 1
rho <- lambda/mu
N <- 500
i <- 1 # initial state. Can be also a probability vector, e.g. i <- rep(1/(N+1), N+1)
j <- 1:min(N+1, 9)  # target states
ei.max.mult <- 5    # maximum multiplicity of the eigenvalues (0 for unlimited)

Q <- build.ctmc.m.m.1.N.matrix.sparse(N, lambda, mu)
uc.mat <- UCmatrix(Q)
uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult)

t <- time.log.scale(dec=c(-6,7), num=20)
uc.s$plot(t, logy=TRUE, logx=TRUE)

## Comparison with the R exponential matrix function, expm
tt <- 1 # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## Comparison with the R exponential matrix function, expm
tt <- 10 * N * uc.s$qn # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## The error in stationary regime can be reduced adding the option
## use.lim=TRUE. In this case, the script solves the stationary
## solution, and it is added to the system that computes the UCs.
uc.mat <- UCmatrix(Q)
uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult, use.lim=TRUE)

t <- time.log.scale(dec=c(-6,7), num=20)
uc.s$plot(t, logy=TRUE, logx=TRUE)

## Comparison with the R exponential matrix function, expm
tt <- 1 # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## Comparison with the R exponential matrix function, expm
tt <- 10 * N * uc.s$qn # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## For lambda=1 the method diverges for N > around 40
## M/M/M/1/N where the method diverge (for large values of t)
##
lambda <- 1
mu <- 1  # For mu=1 the method diverges for N > around 40
rho <- lambda/mu
N <- 100
i <- 1 # initial state
j <- 1:min(N+1, 9)  # target states
ei.max.mult <- 5    # maximum multiplicity of the eigenvalues (0 for unlimited)

Q <- build.ctmc.m.m.1.N.matrix.sparse(N, lambda, mu)
uc.mat <- UCmatrix(Q)

uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult)

t <- time.log.scale(dec=c(-6,7), num=20)
uc.s$plot(t, logy=TRUE, logx=TRUE)

## Comparison with the R exponential matrix function, expm
tt <- 1 # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## Comparison with the R exponential matrix function, expm
tt <- 10 * N * uc.s$qn # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## Option choose.samples=TRUE may help the method to converge. See section 7 in the report
## http://www.ac.upc.edu/RR/2010/53.pdf
uc.mat <- UCmatrix(Q)
uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult, choose.samples=TRUE)

t <- time.log.scale(dec=c(-6,7), num=20)
uc.s$plot(t, logy=TRUE, logx=TRUE)

## Comparison with the R exponential matrix function, expm
tt <- 1 # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))

## Comparison with the R exponential matrix function, expm
tt <- 10 * N * uc.s$qn # sample t value
p.uvand <- uc.s$prob(tt) # probabilities obtained with the closed form computed by the method.
p.expm <- expm(Q*tt)[i,j] # probabilities obtained with the R exponential matrix function, expm.
p.uvand 
p.expm
message("error: ", sum((p.uvand-p.expm)^2))
