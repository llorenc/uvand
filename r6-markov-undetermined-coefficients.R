## Copyright (c) 2013 Llorenç Cerdà-Alabern,
## http://personals.ac.upc.edu/llorenc This file is free software: you
## can redistribute it and/or modify it under the terms of the GNU
## Affero Public License as published by the Free Software Foundation,
## either version 3 of the License, or (at your option) any later
## version.  oo-markov-undetermined-coefficients.R is distributed in
## the hope that it will be useful, but WITHOUT ANY WARRANTY; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.  See the GNU Affero Public License for more
## details.  You should have received a copy of the GNU Affero Public
## License along with oo-markov-undetermined-coefficients.R.  If not,
## see <http://www.gnu.org/licenses/>.
##
## Compute the transient solution using the undetermined coeficient
##
## save an object: save(obj, fname)
## load an object: obj <- Object$load(fname)
library(R6)
library(Matrix)
source("r6-eigen.R")
source("r6-uc.R")

##########################################################################################
debug <- TRUE
UCmatrix.Q <- Matrix()
Vmatrix.Q <- NULL
Bmatrix.Q <- NULL
UCmatrix.states <- integer()
UCmatrix.type <- character()
UCmatrix.qii.max <- numeric()

solve.uc <- function(Q=NULL, qii.max=NULL, j=NULL, i=1, method='vand',
                     unif=FALSE, qn=NULL, choose.samples=FALSE, alpha=10.0,
                     lim.dist=NULL, use.lim=FALSE,
                     states.names=NULL, ei.max=0, ei.max.mult=0,
                     use.eigen.from.file=NULL, ei=NULL) {
  if(!is.null(Q)) {
    UCmatrix.Q <<- Q
  }
  stopifnot(!is.null(UCmatrix.Q)) ;
  UCmatrix.qii.max <<- qii.max
  UCmatrix.states <<- nrow(UCmatrix.Q)
  UCmatrix.type <<- guess.Q.type(UCmatrix.Q)
  msg(class="", "solve.uc", "solving for ", UCmatrix.type)
  #
  j <- choose.j(j, nrow(UCmatrix.Q))
  if(!is.null(lim.dist)) use.lim <- TRUE
  if(unif) {
    if(UCmatrix.type == 'P') {
      if(is.null(UCmatrix.qii.max))
        stop('The matrix seems uniformized and qii.max is unknown')
    } else {
      if(is.null(UCmatrix.qii.max)) UCmatrix.qii.max <<- abs(min(diag(UCmatrix.Q)))
    }
  }
  uc <- UC$new(type=UCmatrix.type,method=method, states=UCmatrix.states, j=j, i=i,
               states.names=states.names, unif=unif, qii.max=UCmatrix.qii.max, qn=qn, 
               choose.samples=choose.samples, alpha=alpha, lim=use.lim,
               ei.max=ei.max, ei.max.mult=ei.max.mult,
               use.eigen.from.file=use.eigen.from.file, ei=ei)
  t.tot <-
    switch(method,
           vand =  system.time(solve.uc.vand(uc, lim=lim.dist)),
           evec = system.time(solve.uc.evec(uc)),
           stop('Unknown method (must be "vand" or "evec")')
           )
  ##
  uc$add.time(list(t.tot=t.tot))
  if(!empty(uc$coef) && isTRUE(all(Re(uc$coef) == 0))) uc$reset.coef() # the system failed
  ##
  if(debug) {
    msg(class="", "solve.uc", 'done')
    Q.print()
    uc$print()
  }
  return(uc)
}

Q.print <- function() {
  if(debug) cat('---------------\n')
  msg("", "Matrix", "summary")
  msg("", " type", UCmatrix.type)
  msg("", " Number of states", UCmatrix.states)
}

##
## Return the uniformization matrix.
##
uniformize <- function(qn=NULL) {
  ## Compute Q = I + 1/qn Q
  msg("", "uniformize", "uniformizing Q")
  stopifnot(UCmatrix.type == 'Q')
  if(is.null(qn)) {
    stopifnot(!is.null(UCmatrix.qii.max))
    UCmatrix.Q <<- UCmatrix.Q/UCmatrix.qii.max
  } else UCmatrix.Q <<- UCmatrix.Q/qn
  diag(UCmatrix.Q) <<- 1 + diag(UCmatrix.Q)
  UCmatrix.type <<- 'P'
}

##
## Rescale the uniformization matrix.
##
rescale <- function(qn) {
  ## Compute Q = I + 1/qn Q
  msg(class="", "rescale", "rescaling Q")
  stopifnot(UCmatrix.type == 'P')
  stopifnot(qn > UCmatrix.qii.max)
  scale <- UCmatrix.qii.max/qn
  UCmatrix.Q <<- UCmatrix.Q * scale
  diag(UCmatrix.Q) <<- (1-scale) + diag(UCmatrix.Q)
}

##
## Solve a CTMC using uniformization
##
limit.distribution <- function(i, absorbing) {
  msg(class="", "limit.distribution", "Compute the limit dist.")
  if(absorbing > 0) absorbing.prob(i, absorbing)
  else stationary.prob()
}

########################################
## UCmatrix solving methods
########################################
##
## Using QR decomposition of a Vanderemonde matrix
##
confl.vand.mat <- function(uc, nr, powers) {
    msg(class="", "confl.vand.mat", 'Build the Vandermonde matrix')
    nc <- length(uc$ei$one)
    if(length(uc$ei$confl) > 0) {
        nc <- nc + sum(uc$ei$mult)
    }
    if(is.null(powers)) vand.matrix(uc$ei$one, nr, nc)
    else vand.matrix.powers(uc$ei$one, powers, nc)
    ##
    if(length(uc$ei$confl) > 0) {
        multi <- length(uc$ei$one)
        if(is.null(powers)) {
            for(i in 1:length(uc$ei$confl)) {
                Vmatrix.Q[,(multi+1):(multi+uc$ei$mult[i])] <<-
                    switch(uc$ei$type,
                           Q = vand.ctmc.confl.block(uc$ei$confl[i], uc$ei$mult[i], nr),
                           P = vand.dtmc.confl.block(uc$ei$confl[i], uc$ei$mult[i], nr))
                multi <- multi+uc$ei$mult[i]
            }
        } else {
            for(i in 1:length(uc$ei$confl)) {
                Vmatrix.Q[,(multi+1):(multi+uc$ei$mult[i])] <<-
                    switch(uc$ei$type,
                           Q = vand.ctmc.confl.block.powers(uc$ei$confl[i], powers, uc$ei$mult[i]),
                           P = vand.dtmc.confl.block.powers(uc$ei$confl[i], powers, uc$ei$mult[i]))
                multi <- multi+uc$ei$mult[i]
            }
        }
    }
}

##
## Coefficients are computed solving a Vandermonde system
##
compute.qn <- function(uc) {
  msg(class="", "compute.qn",'Computing the uniformization parameter')
  if(is.null(uc$qn)) {
    ## most negative eigenvalue
    abs.min.ei <- abs(min(Re(uc$ei$confl), Re(uc$ei$one)))
    ## second largest eigenvalue in modulus
    ei2 <- abs(uc$ei$second.largest.re()$value)
    N <- uc$ei$number()
    if(uc$ei$type == 'P') {
      min.qn <- (uc$qii.max - uc$qii.max * ei2)/(1-.Machine$double.eps^(1/N))
    } else {
      min.qn <- ei2/(1-.Machine$double.eps^(0.8/as.numeric(N)))
      ## min.qn <- ei2/(1-.Machine$double.eps^(0.5/as.numeric(N)))
    }
    uc$set.qn(max(uc$qii.max, abs.min.ei, min.qn))
    if(uc$qn > ei2*N) {
      uc$set.qn(max(uc$qii.max, ei2*N))
    }
  } else msg(class="", "compute.qn", 'qn already initialized')
}

##
## Coefficients are computed using the eigenvectors of Q.
##
solve.uc.eigenvectors <- function(uc) {
    msg(class="", "solve.uc.eigenvectors",'Solving the eigenvectors')
    if(length(uc$i) == 1) {
        x0 <- rep(0,len=nrow(uc$ei$evec)) ; x0[uc$i] <- 1
    } else {
        x0 <- uc$i
    }
    s <- try.solve.kappa(uc$ei$evec, x0)
    uc$set.kappa(s$kappa)
    uc$set.Vmatrix.characteristics(dim=dim(uc$ei$evec),
                                   zeros=length(which(uc$ei$evec == 0)),
                                   bytes=object.size(uc$ei$evec))
    if(!empty(s$x)) {
        uc$set.coef(sapply(uc$j, function(n) uc$ei$evec[n,] * s$x),
                    colnames(uc$coef) <- uc$get.states.names())
    }
    msg(class="", 'solve.uc.eigenvectors', 'Removing eigenvectors')
    uc$ei$remove.evec() # eigenvectors are not needed anymore
    ## correct the eigenvalues if computed with an uniformized matrix
    if(uc$unif && (uc$ei$type == 'P')) uc$ei$toggle.Q.P.map(uc$qii.max)
}

##
## Coefficients are computed solving a Vandermonde system
##
solve.uc.vandermonde <- function(uc, lim=NULL, powers=NULL) {
  confl.vand.mat(uc, nrow(Bmatrix.Q), powers)
  if(!is.null(lim)) {
      msg(class="", "solve.uc.vandermonde",'add the limit distribution')
      Vmatrix.Q <<- rbind(Vmatrix.Q, c(1, rep(0, ncol(Vmatrix.Q)-1)))
      Bmatrix.Q <<- rbind(Bmatrix.Q, lim)
  }
  msg(class="", "solve.uc.vandermonde",'Solve the Vandermonde system')
  s <- try.solve.kappa(Vmatrix.Q, Bmatrix.Q)
  if(!empty(s$x)) {
      uc$set.coef(s$x, uc$get.states.names())
  }
  uc$set.kappa(s$kappa)
  uc$set.Vmatrix.characteristics(dim(Vmatrix.Q),
                                 length(which(Vmatrix.Q == 0)),
                                 object.size(Vmatrix.Q))
}

solve.uc.vand <- function(uc, lim=NULL) {
  msg(class="", "solve.uc.vand", "Solve using QR decomposition of Vandermonde system")
  ## compute the eigenvalues
  if(!is.null(uc$use.eigen.from.file) && !is.null(uc$ei$type)) {
      msg(class="", "solve.uc.vand", "Using eigenvalues from file", uc$use.eigen.from.file)
      t.e <- 0.0
  } else {
      t.e <- system.time(uc$compute.eigenvalues(type=UCmatrix.type,
                                                ei.max=uc$ei.max,
                                                ei.max.mult=uc$ei.max.mult))
  }
  powers <- NULL
  ##
  if(uc$unif) {
    compute.qn(uc)
    if(uc$choose.samples) powers <- uc$compute.samples()
    if(UCmatrix.type == 'P') {
      ## the matrix is yet uniformized
      if(UCmatrix.qii.max != uc$qn) {
        ## rescale the uniformized matrix
        rescale(uc$qn)
        ## rescale the eigenvalues
        uc$ei$rescale(UCmatrix.qii.max, uc$qn)
      }
      ei.tmp <- NULL
    } else {
      ## uniformize the matrix
      uniformize(uc$qn)
      ## save a copy of the current eigenvalues
      ei.tmp <- Eigen$new()
      ei.tmp$copy(uc$ei)
      ## compute the eigenvalues of the uniformized matrix
      uc$ei$toggle.Q.P.map(uc$qn)
    }
    t.c <- uc$time.constant()
    if(t.c > nrow(UCmatrix.Q)) {
      if(!uc$choose.samples) {
        add.w <- "  Try the option choose.samples=TRUE\n" ;
      } else {
        add.w <- ""
      }
      warning(paste(sep='',
                    sprintf("Transient time in slots (%6.1e) larger than nuber of states (%6.1e)\n" ,
                            t.c, nrow(UCmatrix.Q)),
                    sprintf("  The solution can be innacurate for t > %6.2e\n" ,
                            nrow(UCmatrix.Q)/uc$qn),
                    add.w
                    ))
    }
  }
  if(is.null(powers)) t.b <- system.time(compute.b(uc$i, uc$j, nr=uc$ei$number()))
  else  t.b <- system.time(compute.b.powers(uc$i, uc$j, powers))
  t.s <-
    system.time(if(uc$lim) {
      msg(class="", "solve.uc.vand", "Solve using the limit distribution")
      if(is.null(lim)) {
        lim <- limit.distribution(uc$i, uc$ei$absorbing)[uc$j]
      }
    })
  ## Q is not needed anymore
  # msg(class="", "solve.uc.vand", "removing Q")
  remove(UCmatrix.Q, pos = ".GlobalEnv") # <<- NULL
  call.gc("", "solve.uc.vand")
  t.c <- system.time(solve.uc.vandermonde(uc, lim, powers))
  remove(Bmatrix.Q, pos = ".GlobalEnv") # <<- NULL
  remove(Vmatrix.Q, pos = ".GlobalEnv") # <<- NULL
  call.gc("", "solve.uc.vand")
  if(uc$unif) {
    ## compute the UC of the ctmc
    uc$correct.confl.coef()
    ## set the eigenvalues of the ctmc
    if(!is.null(ei.tmp)) uc$ei$copy(ei.tmp)
    else uc$ei$toggle.Q.P.map(uc$qn)
  }
  uc$init.time(list(ei=t.e, lim=t.s, b=t.b, coef=t.c))
}

##
## Using eigenvectors
##
solve.uc.evec <- function(uc) {
  msg(class="", "solve.uc.evec", "Solve using eigenvectors")
  t.e <- system.time(uc$compute.eigenvalues(type=UCmatrix.type, evec=TRUE))
  ## Q is not needed anymore
  msg(class="", "solve.uc.evec", "removing Q")
  remove(UCmatrix.Q, pos = ".GlobalEnv") # <<- NULL
  call.gc("", 'solve.uc.evec')
  t.c <- system.time(uc$solve.uc.eigenvectors())
  uc$init.time(list(ei=t.e, coef=t.c))
}

########################################
## UCmatrix b matrix
########################################
compute.b <- function(i, j, nr) {
  msg(class="", "compute.b", "Compute the matrix b")
  n.j <- 1:length(j)
  if(length(i) == 1) {
    r <- matrix(UCmatrix.Q[i,], nr=1, nc=UCmatrix.states)
    Bmatrix.Q <<- matrix(nc=length(n.j), nr=nr)
    for(n in n.j) { # first row
      if(j[n] == i) { Bmatrix.Q[1,n] <<- 1 }
      else          { Bmatrix.Q[1,n] <<- 0 }
    }
  } else {
    r <- i %*% UCmatrix.Q
    Bmatrix.Q <<- matrix(nc=length(n.j), nr=nr)
    Bmatrix.Q[1,n.j] <<- i[j]
  }
  Bmatrix.Q[2,n.j] <<- r[1,j]
  if(nr > 2) {
    for(n in 3:nr) {
      r <- r %*% UCmatrix.Q
      Bmatrix.Q[n,n.j] <<- r[1,j]
    }
  }
}

########################################
## UCmatrix b matrix for the given powers
########################################
compute.b.powers <- function(i, j, powers) {
  msg(class="", "compute.b.powers", "Compute the matrix b", ", max power=", max(powers))
  n.j <- 1:length(j)
  if(length(i) == 1) {
    r <- matrix(UCmatrix.Q[i,], nr=1, nc=UCmatrix.states)
    Bmatrix.Q <<- matrix(nc=length(n.j), nr=length(powers)+1)
    for(n in n.j) { # first row
      if(j[n] == i) { Bmatrix.Q[1,n] <<- 1 }
      else          { Bmatrix.Q[1,n] <<- 0 }
    }
  } else {
    r <- i %*% UCmatrix.Q
    Bmatrix.Q <<- matrix(nc=length(n.j), nr=length(powers)+1)
    Bmatrix.Q[1,n.j] <<- i[j]
  }
  curr.power <- 1
  n <- 2
  for(pow in powers) {
    while(curr.power < pow) {
      r <- r %*% UCmatrix.Q
      curr.power <- curr.power+1
    }
    Bmatrix.Q[n,n.j] <<- r[1,j]
    n <- n+1
  }
}

#######################################
## UCmatrix limit distribution methods
#######################################
#
# Return the probabilities to be absorbed starting from state i
#
absorbing.prob <- function(i, absorbing) {
  msg(class="", "absorbing.prob", "Compute the absorbing probs.")
  stopifnot(absorbing > 0)
  ## Assume absorbing states and Q in canonical form (absorbing states have
  ## the highest indexes)
  UCmatrix.states <<- nrow(UCmatrix.Q)
  trans <- UCmatrix.states -  absorbing
  if(i > trans) stop('State is not transient:', i, '\n')
  if(absorbing == 1) {
    return(c(rep(0, trans), 1))
  }
  ## Remove last equation and normalize
  R <- as.matrix(UCmatrix.Q[1:trans,(trans+1):UCmatrix.states])
  abs.p <-
    switch(UCmatrix.type,
           Q = try.solve(UCmatrix.Q[1:trans,1:trans], -R),
           P = try.solve(diag(nr=trans)-UCmatrix.Q[1:trans,1:trans], R))
  if(!is.null(abs.p)) return(c(rep(0, trans), abs.p[i,]))
  return(NULL)
}

stationary.prob <- function() {
  ## Remove the first equation and normalize
  msg(class="", "stationary.prob", "Compute the stationary dist.")
  n <- nrow(UCmatrix.Q)
  stat <-
    switch(UCmatrix.type,
           Q = try.solve(t(UCmatrix.Q[2:n,2:n]), -as.matrix(UCmatrix.Q[1,2:n])),
           P = try.solve(t(UCmatrix.Q[2:n,2:n]-Diagonal(n-1, 1)), -as.matrix(UCmatrix.Q[1,2:n])))
  if(!is.null(stat)) {
    stat[stat < 0] <- 0
    stat <- rbind(1, stat)/(sum(stat)+1)
  }
  return(stat)
}

#############################################################################
## Supporting functions
#############################################################################
msg <- function(class, fname, msg=NULL, ...) {
  if(debug) {
    if(!empty(class)) cat(sep='', class, '$')
    cat(sep='', fname)
    if(!is.null(msg)) cat(sep='', ': ', msg, ...)
    cat("\n")
  }
}

empty <- function(obj) {
  is.null(obj) || (length(obj) == 0 || obj == '')
}

choose.j <- function(j, dim.Q) {
  if(is.null(j)) {
    j <- 1:min(10, dim.Q) # take all j up to 10
  } else {
    if((min(j) < 1) || (max(j) > dim.Q)) stop("j out of bounds")
  }
  return(j)
}

#############################################################################
## Matrix functions
#############################################################################
try.solve <- function(A, b) {
  sol <- try(solve(qr(A, LAPACK=TRUE), b))
  if(inherits(sol, "try-error")) {
    msg('', 'try.solve', 'failure solving the system')
    return(NULL)
  }
  msg('', 'try.solve', 'success')
  return(sol)
}

call.gc <- function(o.class=NULL, o.func='call.gc') {
  msg(o.class, o.func, 'garbage collection')
  if(debug) print(gc(reset=TRUE))
  else gc(reset=TRUE)
}

## try.solve.kappa <- function(A, b) {
##   sol <- try(solve(A, b))
##   if(inherits(sol, "try-error")) {
##     msg('', 'try.solve.kappa', 'failure solving the system')
##     return(list(x=NULL, kappa=NULL))
##   }
##   kap <- try(kappa(A, exact=TRUE))
##   if(inherits(kap, "try-error")) {
##     msg('', 'try.solve.kappa', 'failure solving kappa')
##     kap <- NULL
##   }
##   msg('', 'try.solve.kappa', 'success')
##   return(list(x=as.matrix(sol), kappa=kap))
## }

try.solve.kappa <- function(A, b) {
  qr.A <- try(qr(A, LAPACK=TRUE))
  if(inherits(qr.A, "try-error")) {
    msg('', 'try.solve.kappa', 'failure computing QR')
    return(list(x=NULL, kappa=NULL))
  }
  sol <- try(solve.qr(qr.A, b))
  if(inherits(sol, "try-error")) {
    msg('', 'try.solve.kappa', 'failure solving the system')
    return(list(x=NULL, kappa=NULL))
  }
  kap <- try(kappa(qr.A, exact=TRUE))
  if(inherits(kap, "try-error")) {
    msg('', 'try.solve.kappa', 'failure solving kappa')
    kap <- NULL
  }
  msg('', 'try.solve.kappa', 'success')
  return(list(x=as.matrix(sol), kappa=kap))
}

vand.matrix <- function(lambda, nr, nc) {
  ## nc <- length(lambda)
  ## V <- Matrix(data=rep(1, nc), byrow=FALSE, nr=nr, nc=nc, forceCheck=TRUE)
  Vmatrix.Q <<- matrix(nr=nr, nc=nc)
  if(nc > 0) {
    Vmatrix.Q[1,] <<- rep(1, nc)
    for(i in 2:nr) {
      Vmatrix.Q[i,1:length(lambda)] <<- Vmatrix.Q[i-1,1:length(lambda)] * lambda
    }
  }
}

vand.matrix.powers <- function(lambda, powers, nc) {
  ## powers is a vector with the powers>0 to evaluate
  stopifnot(length(powers) > 0)
  Vmatrix.Q <<- matrix(nr=length(powers)+1, nc=nc)
  Vmatrix.Q[1,1:length(lambda)] <<- rep(1, length(lambda))
  ## V <- Matrix(rep(1, length(lambda)), byrow=TRUE, nr=length(powers)+1, nc=length(lambda))
  for(n in 1:length(powers))
    Vmatrix.Q[n+1,1:length(lambda)] <<- lambda^powers[n]
}

##
## Build a dtmc confluent Vandermonde block.
##
vand.dtmc.confl.block.powers <- function(lambda, powers, m) {
    ## powers is a vector with the powers>0 to evaluate
    nr <- length(powers)+1
    stopifnot(nr > 1)
    if(lambda == 0) return(diag(1, nr=nr, nc=m))
    ## A <- Matrix(rep(lambda, nr-1)^c(1:(nr-1)), byrow=FALSE, nr=nr-1, nc=m)
    A <- matrix(nr=nr, nc=m)
    ## First row
    A[1,] <- c(1, rep(0, m-1))
    for(n in 1:(nr-1)) {
        A[n+1,] <- rep(lambda, len=m)^powers[n]*powers[n]^(0:(m-1))
    }
    return(A)
}

vand.ctmc.confl.block.powers <- function(lambda, powers, m) {
    stop("vand.dtmc.confl.block.powers: not yet implemented")
}
##
## Build a ctmc confluent Vandermonde block.
##
vand.ctmc.confl.block <- function(lambda, m, nr) {
  if(lambda == 0) return(diag(1, nr=nr, nc=m))
  ## A <- Matrix(c(1, (rep(lambda, (nr-1))^c(1:(nr-1)))), byrow=FALSE, nr=nr, nc=m)
  A <- matrix(nr=nr, nc=m)
  ## First col
  A[,1] <- c(1, (rep(lambda, (nr-1))^c(1:(nr-1))))
  for(j in 1:(m-1)) {
    A[,(j+1)] <- 
      c(rep(0, j), sapply(j:(nr-1), function(n) prod(n:(n-j+1))) *
        rep(lambda, (nr-j))^c(0:(nr-j-1)))
  }
  return(A)
}

##
## Build a dtmc confluent Vandermonde block.
##
vand.dtmc.confl.block <- function(lambda, m, nr) {
  if(lambda == 0) return(diag(1, nr=nr, nc=m))
  ## A <- Matrix(rep(lambda, nr-1)^c(1:(nr-1)), byrow=FALSE, nr=nr-1, nc=m)
  A <- matrix(nr=nr-1, nc=m)
  ## First col
  A[,1] <- rep(lambda, nr-1)^c(1:(nr-1))
  for(j in 2:m) A[,j] <-  A[,j-1] * c(1:(nr-1))
  return(rbind(c(1, rep(0, m-1)), A))
}

##
## Björck and Pereyra's Algorithm
## 
solve.vandermonde.bp <- function(alpha, b) {
  n <- length(alpha)
  x <- b
  for(k in 1:(n-1)) {
    for(j in n:(k+1)) {
      x[j] <- x[j] - alpha[k] * x[j-1]
    }
  }
  for(k in (n-1):1) {
    for(j in (k+1):n) {
      x[j] <- x[j]/(alpha[j] - alpha[j-k])
    }
    for(j in k:(n-1)) {
      x[j] <- x[j] - x[j+1]
    }
  }
  return(x)
}

guess.Q.type <- function(Q) {
  stopifnot(!is.null(Q))
  if(isTRUE(all.equal(apply(Q,1,sum), rep(0,nrow(Q))))) {
    ## Is an infinitesimal generator
    type <- 'Q'
  } else if(isTRUE(all.equal(apply(Q,1,sum), rep(1,nrow(Q))))) {
    ## Is a stochastic matrix
    type <- 'P'
  } else {
    stop("Nor stochastic or infinitesimal generator matrix?")
  }
  return(type)
}

time.log.scale <- function(dec, num) {
  c(sapply(c(dec[1]:dec[2]), function(t) 10^(t + seq(0,0.9, len=num))), 10^dec[2])
}

## return n integer samples evenly spaced in log scale
int.log.scale <- function(max, N) {
  n <- N/ceiling(log10(max))
  c(1:(n-1), round(10^(seq(log10(n),log10(max), len=N-n+1))))
}
