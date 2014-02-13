##########################################################################################
## Copyright (c) 2013 Llorenç Cerdà-Alabern, http://personals.ac.upc.edu/llorenc
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## oo-uc.R is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero Public License for more details.
## You should have received a copy of the GNU Affero Public License
## along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
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
lambda <- 0.5
mu <- 1
rho <- lambda/mu
N <- 500
i <- 1 # initial state. Can be also a probability vector, e.g. i <- rep(1/(N+1), N+1)
j <- 1:min(N+1, 9)  # target states
ei.max.mult <- 5    # maximum multiplicity of the eigenvalues (0 for unlimited)

uc.mat <- UCmatrix(Q=build.ctmc.m.m.1.N.matrix.sparse(N, lambda, mu))

uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
#                        use.lim=T,
#                        choose.samples=T,
                        ei.max.mult=ei.max.mult)

t <- time.log.scale(dec=c(-6,7), num=50)
uc.s$plot(t, logy=TRUE, logx=TRUE)

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

## For lambda=1 the method diverges for N > around 40
## M/M/M/1/N where the method diverge
##
lambda <- 1
mu <- 1  # For mu=1 the method diverges for N > around 40
rho <- lambda/mu
N <- 100
i <- 1 # initial state
j <- 1:min(N+1, 9)  # target states
ei.max.mult <- 5    # maximum multiplicity of the eigenvalues (0 for unlimited)

uc.mat <- UCmatrix(Q=build.ctmc.m.m.1.N.matrix.sparse(N, lambda, mu))

uc.s <- uc.mat$solve.uc(method='vand', unif=TRUE, j=j, i=i,
                        ei.max.mult=ei.max.mult)

t <- time.log.scale(dec=c(-6,7), num=50)
uc.s$plot(t, logy=TRUE, logx=TRUE)

