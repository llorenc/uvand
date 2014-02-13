## Copyright (c) 2013 Llorenç Cerdà-Alabern, http://personals.ac.upc.edu/llorenc
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU Affero Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## oo-markov-undetermined-coefficients.R is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero Public License for more details.
## You should have received a copy of the GNU Affero Public License
## along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
##
## Compute the transient solution using the undetermined coeficient
##
## save an object: save(obj, fname)
## load an object: obj <- Object$load(fname)
library(R.oo)
library(Matrix)
library(R.matlab)
source("oo-eigen.R")
source("oo-uc.R")

debug <- TRUE
##########################################################################################
## UC Matrix class
##########################################################################################
setConstructorS3("UCmatrix", function(Q=NULL, qii.max=NULL) {
  if(!is.null(Q)) {
    call.gc("UCmatrix", "constructor")
    stopifnot(inherits(Q, "Matrix"))
    states <- dim(Q)[1]
    type <- guess.Q.type(Q)
  } else {
    states <- NULL
    type <- NULL
  }
  ##
  extend(Object(), "UCmatrix",
         Q=Q,             # the matrix
         type=type,       # matrix type 'Q', 'P'
         qii.max=qii.max, # max_i |Qii| if the matrix is uniformized
         states=states    # states of the chain
         );
})

matlab <- NULL
#######################################
## UC Matrix public methods
#######################################
##
## Solve a CTMC/DTMC
##
setMethodS3("solve.uc", "UCmatrix", appendVarArgs=FALSE, function(this,
                                      j=NULL, i=1, method='vand',
                                      unif=FALSE, qn=NULL, choose.samples=FALSE,
                                      lim.dist=NULL, use.lim=FALSE, rm.lim.eq=FALSE,
                                      states.names=NULL, ei.max=0, ei.max.mult=0,
                                      use.matlab=FALSE) {
  msg(class=class(this)[1], "solve.uc", "solving for ", this$type)
  if(use.matlab) {
    msg(class=class(this)[1], "solve.uc", "using matlab")
    ## Start the matlab server on the same machine
    Matlab$startServer()
    ## Create a Matlab client
    matlab <<- Matlab(host="localhost")
    ## Connect to the Matlab server
    if (!open(matlab))
      throw("Matlab server is not running: waited 30 seconds.")
  }
  j <- choose.j(j, dim(this$Q)[1])
  if(!is.null(lim.dist)) use.lim <- TRUE
  if(unif) {
    if(this$type == 'P') {
      if(is.null(this$qii.max))
        stop('The matrix seems uniformized and qii.max is unknown')
    } else {
      if(is.null(this$qii.max)) this$qii.max <- abs(min(diag(this$Q)))
    }
  }
  uc <- UC(type=this$type, method=method, states=this$states, j=j, i=i,
           states.names=states.names, unif=unif, qii.max=this$qii.max, qn=qn, 
           choose.samples=choose.samples, lim=use.lim, rm.lim.eq=rm.lim.eq,
           ei.max=ei.max, ei.max.mult=ei.max.mult, use.matlab=use.matlab)
  t.tot <-
    switch(method,
           vand =  system.time(this$solve.uc.vand(uc, lim=lim.dist)),
           evec = system.time(this$solve.uc.evec(uc)),
           stop('Unknown method (must be "vand" or "evec")')
           )
  ##
  uc$time <- c(uc$time, list(t.tot=t.tot))
  if(!empty(uc$coef) && isTRUE(all(Re(uc$coef) == 0))) uc$coef <- NULL # the system failed
  ##
  if(debug) {
    msg(class=class(this)[1], "solve.uc", 'done')
    this$print()
    uc$print()
  }
  if(use.matlab) {
    ## close the Matlab client and shutdown the server
    close(matlab)
  }
  return(uc)
})

setMethodS3("print", "UCmatrix", function(this) {
  if(debug) cat('---------------\n')
  msg("", "UCmatrix", "summary")
  msg("", " Matrix type", this$type)
  msg("", " Number of states", this$states)
})

##
## Return the uniformization matrix.
##
setMethodS3("uniformize", "UCmatrix", appendVarArgs=FALSE, function(this, qn=NULL) {
  ## Compute Q = I + 1/qn Q
  msg(class=class(this)[1], "uniformize", "uniformizing Q")
  stopifnot(this$type == 'Q')
  if(is.null(qn)) {
    stopifnot(!is.null(this$qii.max))
    this$Q <- this$Q/this$qii.max
  } else this$Q <- this$Q/qn
  diag(this$Q) <- 1 + diag(this$Q)
  this$type <- 'P'
})

##
## Rescale the uniformization matrix.
##
setMethodS3("rescale", "UCmatrix", appendVarArgs=FALSE, function(this, qn) {
  ## Compute Q = I + 1/qn Q
  msg(class=class(this)[1], "rescale", "rescaling Q")
  stopifnot(this$type == 'P')
  stopifnot(qn > this$qii.max)
  scale <- this$qii.max/qn
  this$Q <- this$Q * scale
  diag(this$Q) <- (1-scale) + diag(this$Q)
})

##
## Solve a CTMC using uniformization
##
setMethodS3("limit.distribution", "UCmatrix", appendVarArgs=FALSE, function(this,
                                                              i, absorbing) {
  msg(class=class(this)[1], "limit.distribution", "Compute the limit dist.")
  if(absorbing > 0) this$absorbing.prob(i, absorbing)
  else this$stationary.prob()
})

#############################################################################
## UCmatrix private methods
#############################################################################
########################################
## UCmatrix solving methods
########################################
##
## Using QR decomposition of a Vanderemonde matrix
##
setMethodS3("solve.uc.vand", "UCmatrix", appendVarArgs=FALSE, priv=TRUE, function(this, uc, lim=NULL) {
  msg(class=class(this)[1], "solve.uc.vand", "Solve using QR decomposition of Vandermonde system")
  ## compute the eigenvalues
  t.e <- system.time(uc$ei$init(Q=this$Q, type=this$type, ei.max=uc$ei.max,
                                ei.max.mult=uc$ei.max.mult, use.matlab=uc$use.matlab))
  powers <- NULL
#browser()
  if(uc$unif) {
    uc$compute.qn()
    if(uc$choose.samples) powers <- uc$compute.samples()
    if(this$type == 'P') {
      ## the matrix is yet uniformized
      if(this$qii.max != uc$qn) {
        ## rescale the uniformized matrix
        this$rescale(uc$qn)
        ## rescale the eigenvalues
        uc$ei$rescale(this$qii.max, uc$qn)
      }
      ei.tmp <- NULL
    } else {
      ## uniformize the matrix
      this$uniformize(uc$qn)
      ## save a copy of the current eigenvalues
      ei.tmp <- clone(uc$ei)
      ## compute the eigenvalues of the uniformized matrix
      uc$ei$toggle.Q.P.map(uc$qn)
    }
    t.c <- uc$time.constant(alpha=0.5)
    if(t.c > dim(this$Q)[1]) {
      if(!uc$choose.samples) {
        add.w <- "  Try the option choose.samples=TRUE\n" ;
      } else {
        add.w <- ""
      }
      warning(paste(sep='',
                    sprintf("Transient time in slots (%6.1e) larger than nuber of states (%6.1e)\n" ,
                            t.c, dim(this$Q)[1]),
                    sprintf("  The solution can be innacurate for t > %6.2e\n" ,
                            dim(this$Q)[1]/uc$qn),
                    add.w
                    ))
    }
  }
  ## if(uc$rm.lim.eq) nr <- this$states-uc$ei$lim.mult
  ## else nr <- this$states-uc$ei$lim.mult+1
  if(is.null(powers)) t.b <- system.time(b <- this$compute.b(uc$i, uc$j, nr=uc$ei$number()))
  else  t.b <- system.time(b <- this$compute.b.powers(uc$i, uc$j, powers))
  t.s <-
    system.time(if(uc$lim) {
      msg(class=class(this)[1], "solve.uc.vand", "Solve using the limit distribution")
      if(is.null(lim)) {
        lim <- this$limit.distribution(uc$i, uc$ei$absorbing)[uc$j]
      }
      if(uc$rm.lim.eq) {
        if(uc$ei$type == 'Q') {
          b[1,] <- b[1,] - lim
        } else {
          b <- b - matrix(rep(lim, dim(b)[1]), nr=dim(b)[1], byr=TRUE)
        }
      }
    })
  ## Q is not needed anymore
  msg(class=class(this)[1], "solve.uc.vand", "removing Q")
  this$Q <- NULL
  call.gc(class(this)[1], "solve.uc.vand")
  t.c <- system.time(uc$solve.uc.vandermonde(b, lim, powers))
  if(uc$unif) {
    ## compute the UC of the ctmc
    uc$correct.confl.coef()
    ## set the eigenvalues of the ctmc
    if(!is.null(ei.tmp)) uc$ei <- ei.tmp
    else uc$ei$toggle.Q.P.map(uc$qn)
  }
  uc$time <- list(ei=t.e, lim=t.s, b=t.b, coef=t.c)
})

##
## Using eigenvectors
##
setMethodS3("solve.uc.evec", "UCmatrix", appendVarArgs=FALSE, priv=TRUE, function(this, uc,
                                                            use.matlab=FALSE) {
  msg(class=class(this)[1], "solve.uc.evec", "Solve using eigenvectors")
  t.e <- system.time(uc$ei$init(this$Q, this$type, evec=TRUE, use.matlab=use.matlab))
  ## Q is not needed anymore
  msg(class=class(this)[1], "solve.uc.evec", "removing Q")
  this$Q <- NULL
  call.gc(class(this)[1], 'solve.uc.evec')
  t.c <- system.time(uc$solve.uc.eigenvectors())
  uc$time <- list(ei=t.e, coef=t.c)
})

########################################
## UCmatrix b matrix
########################################
setMethodS3("compute.b", "UCmatrix", appendVarArgs=FALSE, priv=TRUE, function(this, i, j, nr) {
  msg(class=class(this)[1], "compute.b", "Compute the matrix b")
  n.j <- 1:length(j)
  if(length(i) == 1) {
    r <- matrix(this$Q[i,], nr=1, nc=this$states)
    b <- matrix(nc=length(n.j), nr=nr)
    for(n in n.j) { # first row
      if(j[n] == i) { b[1,n] <- 1 }
      else          { b[1,n] <- 0 }
    }
  } else {
    r <- i %*% this$Q
    b <- matrix(nc=length(n.j), nr=nr)
    b[1,n.j] <- i[j]
  }
  b[2,n.j] <- r[1,j]
  if(nr > 2) {
    for(n in 3:nr) {
      r <- r %*% this$Q
      b[n,n.j] <- r[1,j]
    }
  }
  return(b)
})

########################################
## UCmatrix b matrix for the given powers
########################################
setMethodS3("compute.b.powers", "UCmatrix", appendVarArgs=FALSE, priv=TRUE, function(this, i, j, powers) {
  msg(class=class(this)[1], "compute.b.powers", "Compute the matrix b", ", max power=", max(powers))
  n.j <- 1:length(j)
  if(length(i) == 1) {
    r <- matrix(this$Q[i,], nr=1, nc=this$states)
    b <- matrix(nc=length(n.j), nr=length(powers)+1)
    for(n in n.j) { # first row
      if(j[n] == i) { b[1,n] <- 1 }
      else          { b[1,n] <- 0 }
    }
  } else {
    r <- i %*% this$Q
    b <- matrix(nc=length(n.j), nr=length(powers)+1)
    b[1,n.j] <- i[j]
  }
  curr.power <- 1
  n <- 2
  for(pow in powers) {
    while(curr.power < pow) {
      r <- r %*% this$Q
      curr.power <- curr.power+1
    }
    b[n,n.j] <- r[1,j]
    n <- n+1
  }
  return(b)
})

#######################################
## UCmatrix limit distribution methods
#######################################
#
# Return the probabilities to be absorbed starting from state i
#
setMethodS3("absorbing.prob", "UCmatrix", appendVarArgs=FALSE, priv=TRUE, function(this,
                                                             i, absorbing) {
  msg(class=class(this)[1], "absorbing.prob", "Compute the absorbing probs.")
  stopifnot(absorbing > 0)
  ## Assume absorbing states and Q in canonical form (absorbing states have
  ## the highest indexes)
  states <- dim(this$Q)[1]
  trans <- states -  absorbing
  if(i > trans) stop('State is not transient:', i, '\n')
  if(absorbing == 1) {
    return(c(rep(0, trans), 1))
  }
  ## Remove last equation and normalize
  R <- as.matrix(this$Q[1:trans,(trans+1):states])
  abs.p <-
    switch(this$type,
           Q = try.solve(this$Q[1:trans,1:trans], -R),
           P = try.solve(diag(nr=trans)-this$Q[1:trans,1:trans], R))
  if(!is.null(abs.p)) return(c(rep(0, trans), abs.p[i,]))
  return(NULL)
})

setMethodS3("stationary.prob", "UCmatrix", appendVarArgs=FALSE, priv=TRUE, function(this) {
  ## Remove the first equation and normalize
  msg(class=class(this)[1], "stationary.prob", "Compute the stationary dist.")
  n <- dim(this$Q)[1]
  stat <-
    switch(this$type,
           Q = try.solve(t(this$Q[2:n,2:n]), -as.matrix(this$Q[1,2:n])),
           P = try.solve(t(this$Q[2:n,2:n]-Diagonal(n-1, 1)), -as.matrix(this$Q[1,2:n])))
  if(!is.null(stat)) {
    ##:ess-bp-start::browser:##
    ##:ess-bp-end:##
    stat <- as.vector(stat)
    stat[stat < 0] <- 0
    stat <- c(1, stat)/(sum(stat)+1)
  }
  return(stat)
})

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
  V <- matrix(nr=nr, nc=nc)
  if(nc > 0) {
    V[1,] <- rep(1, nc)
    for(i in 2:nr) {
      V[i,1:length(lambda)] <- V[i-1,1:length(lambda)] * lambda
    }
  }
  return(V)
}

vand.matrix.powers <- function(lambda, powers, nc) {
  ## powers is a vector with the powers>0 to evaluate
  stopifnot(length(powers) > 1)
  V <- matrix(nr=length(powers)+1, nc=nc)
  V[1,1:length(lambda)] <- rep(1, length(lambda))
  ## V <- Matrix(rep(1, length(lambda)), byrow=TRUE, nr=length(powers)+1, nc=length(lambda))
  for(n in 1:length(powers))
    V[n+1,1:length(lambda)] <- lambda^powers[n]
  return(V)
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
  if(isTRUE(all.equal(apply(Q,1,sum), rep(0,dim(Q)[1])))) {
    ## Is an infinitesimal generator
    type <- 'Q'
  } else if(isTRUE(all.equal(apply(Q,1,sum), rep(1,dim(Q)[1])))) {
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
