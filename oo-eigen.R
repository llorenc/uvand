## Copyright (c) 2013 Llorenç Cerdà-Alabern,
## http://personals.ac.upc.edu/llorenc This file is free software: you
## can redistribute it and/or modify it under the terms of the GNU
## Affero Public License as published by the Free Software Foundation,
## either version 3 of the License, or (at your option) any later
## version.  oo-eigen.R is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero Public License for more details.  You should have
## received a copy of the GNU Affero Public License along with
## oo-eigen.R.  If not, see <http://www.gnu.org/licenses/>.
##
## Compute confluent eigenvalues
##
setConstructorS3("Eigen", function() {
  extend(Object(), "Eigen",
         type=NULL,               # eigenvalues type: 'Q', 'P'
         one=vector('numeric'),   # single eigenvalues
         confl=vector('numeric'), # confluent eigenvalues
         mult=vector('integer'),  # multiplicity of confluent eigenvalues
         cmult=vector('integer'), # multplicity before applying ei.max.mult
         evec=NULL,               # right eigenvectors
         lim.mult=integer(),      # multiplicity of the limit eigenvalue
         absorbing=integer(),     # number of absorbing states
         tol=.Machine$double.eps^0.5 # tolerance in some comparisons
         );
})

##########################################################################################
## Eigen methods
##########################################################################################
setMethodS3("print", "Eigen", appendVarArgs=FALSE, function(this, ...) {
  if(debug) cat('---------------\n')
  msg("", "Eigen", "summary")
  msg("", " Matrix type", this$type)
  msg("", " number of single eigenvalues", length(this$one))
  msg("", " number of confluent eigenvalues", length(this$confl))
  if(length(this$confl)>0) {
      msg("", " computed multiplicities", paste(sep='', this$cmult, collapse=','))
      msg("", " applied multiplicities", paste(sep='', this$mult, collapse=','))      
  }
  msg("", " right eigenvectors", !empty(this$evec))
  msg("", " multiplicity of the limit eigenvalue", this$lim.mult)
  msg("", " number of absorbing states", this$absorbing)
  msg("", " tolerance", this$tol)
})

##
## Compute the eigenvalues and its multiplicity
##
setMethodS3("init", "Eigen", appendVarArgs=FALSE, function(this, Q, type, evec=FALSE,
                               ei.max=0, ei.max.mult=0, use.matlab=FALSE) {
  msg(class=class(this)[1], "init", "Compute the eigenvalues")
  N <- nrow(Q)
  stopifnot(N>0)
  stopifnot((ei.max.mult==0)||(ei.max.mult>1))
  this$tol <- (((1+.Machine$double.eps)^N)-1)^0.5
  this$type <- type
  if(use.matlab) {
    if(evec) {
      setVariable(matlab, Q=t(Q))
      res <- evaluate(matlab, paste(sep='', '[V,ei]=eigs(Q,', N,');'))
      this$one <- getVariable(matlab, c("ei"))$ei
      this$evec <- getVariable(matlab, c("V"))$V
    } else {
      setVariable(matlab, Q=Q)
      res <- evaluate(matlab, paste(sep='', 'ei=eigs(Q,', N,');'))
      this$one <- getVariable(matlab, c("ei"))$ei
    }
  } else {
    if(evec) {
      ei <- eigen(t(Q), only.values=FALSE, sym=FALSE)
      this$evec <- ei$vectors
      this$one <- ei$values
    } else {
      this$one <- eigen(Q, only.values=TRUE, sym=FALSE)$values
    }
  }
  this$lim.mult <- 1
  if(this$type == 'Q') {
    ## First eigenvalue must be 0
    this$one <- this$one[length(this$one):1] 
    if(!is.null(this$evec)) this$evec <- this$evec[,ncol(this$evec):1]
    stopifnot(isTRUE(all.equal(abs(this$one[1]), 0, tol=this$tol)))
    this$one[1] <- 0
    ## Check for absorbing states
    abs.s <- which(sapply(1:N, function(i) all(Q[i,]==0)))
  } else if(this$type == 'P') {
    ## First eigenvalue must be 1
    stopifnot(isTRUE(all.equal(abs(this$one[1]), 1, tol=this$tol)))
    this$one[1] <- 1
    ## Check for absorbing states
    abs.s <- which(sapply(1:N, function(i) Q[i,i]==1))
  } else {
    stop("Nor stochastic or infinitesimal generator matrix?")
  }
  this$absorbing <- length(abs.s)
  if((this$absorbing > 0) && (length(intersect(c(1:(abs.s[1]-1)), abs.s)) > 0))
    warning('There are absorbing states and Q is not in canonical form.\n')
  if(evec == FALSE) this$compute.mult(ei.max, ei.max.mult)
})

##
## Split ei.one in single and multiple eigenvalues
##
setMethodS3("compute.mult", "Eigen", appendVarArgs=FALSE, private=TRUE, function(this,
                                                        ei.max, ei.max.mult) {
  msg(class=class(this)[1], "compute.mult", "check for confluent eigenvalues")
  # browser()
  ##
  ## Supporting functions
  ##
  ## swap an eigenvalue pair
  ei.swap <- function(i, j) {
    # browser()
    if(i != j) { 
      tmp <- this$one[i]
      this$one[i] <- this$one[j]
      this$one[j] <- tmp
    }
  } 
  ## Make real a conjugate pair of eigenvalue i
  set.real <- function(i) {
    # browser()
    if(i < length(this$one)) {
      ei.conj <- Conj(this$one[i])
      for(j in (i+1):length(this$one)) {
        if(ei.conj == this$one[j]) {
          msg(class=class(this)[1], 'set.real', 'set Im part to zero: ',
              this$one[i], ', ', this$one[j])
          this$one[i] <- Re(this$one[i])
          this$one[j] <- Re(this$one[j])
          return()
        }
      }
    }
    warning("Eigen$compute.mult: Conjugate pair not found: ", this$one[i])
  }
  ## Swap the conjugate of eigenvalue i from position >= j to position j
  swap.conjugate <- function(i, j) {
    # browser()
    if(i < length(this$one)) {
      ei.conj <- Conj(this$one[i])
      for(k in j:length(this$one)) {
        if(ei.conj == this$one[k]) ei.swap(k, j)
        return()
      }
    }
    warning("Eigen$compute.mult: Conjugate pair not found: ", this$one[i])
  }
  ## Find the multiplicity of a real eigenvalue.
  ## Eigenvalues are assumed sorted by its real part.
  find.real.ei.multiplicity <- function(i) {
    # browser()
    m <- 1    # multiplicity
    if(i < length(this$one)) {
      for(j in (i+1):length(this$one)) {
        if(abs(Re(this$one[i])-Re(this$one[j])) < this$tol) {
          if(Im(this$one[j]) == 0) {
            ## found multiplicity
            ei.swap(i+m, j) # put confluent eigenvalues consecutive
            m <- m+1
          }
        } else {
          if(m>1)
            msg(class=class(this)[1], 'find.real.ei.multiplicity',
                'found confl. eigenvalue: ', this$one[i])
          break
        }
      }
    }
    return(m)
  }
  ## Find the multiplicity of complex eigenvalues
  ## Eigenvalues are assumed sorted by its real part.
  find.complex.ei.multiplicity <- function(i) {
    # browser()
    m <- 1    # multiplicity
    if(i < length(this$one)) {
      for(j in (i+1):length(this$one)) {
        if(abs(Re(this$one[i])-Re(this$one[j])) < this$tol) {
          if(abs(Im(this$one[i])-Im(this$one[j])) < this$tol) {
            ## found multiplicity
            ei.swap(i+m, j)
            m <- m+1
          }
        } else {
          if(m>1)
            msg(class=class(this)[1], 'find.complex.ei.multiplicity',
                'found confl. eigenvalue: ', this$one[i])
          break
        }
      }
    }
    return(m)
  }
  set.ei.mult <- function(i, m) {
    # browser()
    this$confl[mult.i] <- this$one[i]
    ## remove the confl. eig. from  this$one
    if(i+m > length(this$one)) length(this$one) <- i-1
    else this$one <- c(this$one[1:(i-1)], this$one[(i+m):length(this$one)])
    ## adjust the multiplicity
    if(ei.max.mult > 0) {
      this$cmult[mult.i] <- m
      if(m > ei.max.mult) {
        msg(class=class(this)[1], 'set.ei.mult', 'reducing multiplicity of eigenvalue ',
            this$confl[mult.i], ': ', m, ' to ', ei.max.mult)
        m <- ei.max.mult
      }
    }
    msg(class=class(this)[1], 'set.ei.mult', 'setting confluent eigenvalue ',
        this$confl[mult.i], ': ', m)
    this$mult[mult.i] <- m
    mult.i <<- mult.i + 1
  }
  ## find the conjugate of eigenvalue i 
  find.conjugate <- function(ei, i) {
    # browser()
    ei.conj <- Conj(ei[i])
    id <- which(ei.conj == ei)
    if(length(id) == 0) stop('find.conjugate: ', ei.conj, ' not found')
    return(id)
  }
  ##
  ## body method
  ##
  if(this$absorbing > 0) {
    ## remove the absorbing eigenvalues
    this$one <- this$one[this$absorbing:length(this$one)]
  }
  ## Remove Im part of almost real eigenvalues
  for(i in 2:(length(this$one)-1)) {
    if((Im(this$one[i]) != 0) && (abs(Im(this$one[i])) < this$tol)) set.real(i)
  }
  ## sort by real part to facilitate finding multiplicities
  this$one <- this$one[order(Re(this$one), decreasing=TRUE)]
  ## Check multiplicity
  curr.i <- 2 # current ei (do not check for the lim. eigenvalue)
  mult.i <- 1 # index to next free element in this$confl
    # browser()
  repeat {
    if(Im(this$one[curr.i]) == 0) {
      ## look for real comfluent eigenvalue
      m <- find.real.ei.multiplicity(curr.i)
      if(m > 1) set.ei.mult(curr.i, m) # process the confluent eigenvalue
      else curr.i <- curr.i+1
    } else {
      ## look for complex comfluent eigenvalue
      m <- find.complex.ei.multiplicity(curr.i)
      if(m > 1) { # process the confluent eigenvalue
        ## Move the conjugate pairs to the next positions, since they must
        ## have the same multiplicity.
        for(j in curr.i:(curr.i+m-1)) swap.conjugate(j, j+m)
        set.ei.mult(curr.i, m)
        set.ei.mult(curr.i, m)
      } else curr.i <- curr.i+1
    }
    if(curr.i >= length(this$one)) break
  }
    # browser()
  ## choose at most ei.max eigenvalues
  if((ei.max > 0) && (this$number() > ei.max)) {
    msg(class=class(this)[1], 'compute.mult', 'reducing mult. to ', ei.max)
    single <- length(this$one)
    mult <- length(this$confl)
    choose.ei.one <- c(1)
    choose.ei.confl <- c()
    sndl <- this$second.largest.re()
    if(sndl$confl) rand.samples <- single+sndl$id
    else rand.samples <- sndl$id
    set.seed(1)
    rand.samples <- c(rand.samples, sample(2:(single+mult)))
    num.of.ei <- 1
    for(curr.sample in rand.samples) {
      if(curr.sample <= single) {
        if(!any(choose.ei.one == curr.sample)) {
          if(Im(this$one[curr.sample]) == 0) {
            num.of.ei <- num.of.ei + 1
            choose.ei.one <- c(choose.ei.one, curr.sample)
          } else {
            num.of.ei <- num.of.ei + 2
            choose.ei.one <- c(choose.ei.one, curr.sample,
                               find.conjugate(this$one, curr.sample))
          }
        }
      } else {
        curr.sample <- curr.sample-single
        if(!any(choose.ei.confl == curr.sample)) {
          if(Im(this$confl[curr.sample]) == 0) {
            num.of.ei <- num.of.ei + this$mult[curr.sample]
            choose.ei.confl <- c(choose.ei.confl, curr.sample)
          } else {
            num.of.ei <- num.of.ei + 2 * this$mult[curr.sample]
            choose.ei.confl <- c(choose.ei.confl, curr.sample,
                               find.conjugate(this$confl, curr.sample))
          }
        }
      }
      if(num.of.ei >= ei.max) break
    }
    if(length(this$one) > length(choose.ei.one)) {
      msg(class=class(this)[1], 'compute.mult', 'reducing single eigenvalues from ',
          length(this$one), ' to ', length(choose.ei.one))
      this$one <- this$one[choose.ei.one]
    }
    if(length(this$confl) > length(choose.ei.confl)) {
      msg(class=class(this)[1], 'compute.mult', 'reducing confl eigenvalues from ',
          length(this$confl), ' to ', length(choose.ei.confl))
      this$confl <- this$confl[choose.ei.cofl]
      this$mult <- this$mult[choose.ei.confl]
    }
    msg(class=class(this)[1], 'compute.mult', ' final number of eigenvalues ', num.of.ei)
  }
    # browser()
})

##
## Return the second smallest eigenvalue in modulus (id, value and confl indicator)
##
setMethodS3("second.smallest.abs", "Eigen", appendVarArgs=FALSE, function(this) {
  find.min.abs <- function(x) {
    id <- which.min(abs(x))
    list(id=id, value=x[id])
  }
  ##
  if(length(this$one) > 2) {
    min.one <- find.min.abs(this$one[2:length(this$one)])
    min.one$id <- min.one$id+1
  } else min.one <- NULL
  ##
  if(length(this$confl) > 0) {
    min.confl <- find.min.abs(this$confl)
    if(is.null(min.one) || (abs(min.confl$value) < abs(min.one$value)))
      return(c(min.confl, confl=TRUE))
  }
  return(c(min.one, confl=FALSE))
})

##
## Return the second largest real part (id, value and confl indicator)
##
setMethodS3("second.largest.re", "Eigen", appendVarArgs=F, function(this) {
  find.max.re <- function(x) {
    id <- which.max(Re(x))
    list(id=id, value=x[id])
  }
  ##
  if(length(this$one) > 2) {
    max.one <- find.max.re(this$one[2:length(this$one)])
    max.one$id <- max.one$id+1
  } else max.one <- NULL
  ##
  if(length(this$confl) > 0) {
    max.confl <- find.max.re(this$confl)
    if(is.null(max.one) || (Re(max.confl$value) > Re(max.one$value)))
      return(c(max.confl, confl=TRUE))
  }
  return(c(max.one, confl=FALSE))
})

##
## Number of eigenvalues
##
setMethodS3("number", "Eigen", appendVarArgs=F, function(this) {
  if(length(this$confl) == 0) {
    num <- length(this$one)
  } else {
    num <- length(this$one) + sum(this$mult)
  }
  return(num)
})

##
## Number of eigenvalues
##
setMethodS3("plot", "Eigen", appendVarArgs=F, function(this, lim=FALSE) {
  ## If lim=TRUE plot the limit eigenvalue
  ei <- vector('complex')
  if(length(this$one) > 0) {
    ei <- this$one
    if(!lim) ei <- ei[-1]
  }
  if(length(this$confl) > 0) {
    ei <- c(ei, this$confl)
  }
  plot(Re(ei), Im(ei))
})


##
## Toggle P/Q-eigenvalues to Q/P-eigenvalues
##
setMethodS3("toggle.Q.P.map", "Eigen", appendVarArgs=F, function(this, qn) {
  msg(class=class(this)[1], "toggle.Q.P.map")
  stopifnot(!is.null(qn))
  single <- length(this$one)   # number of single eigenvalues
  mult <- length(this$confl) # number of non single eigenvalues
  switch(this$type,
         Q = {
           msg(class=class(this)[1], "toggle.Q.P.map","toggle Q-eigenvalues to P-eigenvalues")
           this$type <- 'P'
           if(single) this$one <- 1 + (this$one/qn)
           if(mult) this$confl <- 1 + (this$confl/qn)     
         },
         P = {
           msg(class=class(this)[1], "toggle.Q.P.map","toggle P-eigenvalues to Q-eigenvalues")
           this$type <- 'Q'
           if(single) this$one <- qn * (this$one-1)
           if(mult) this$confl <- qn * (this$confl-1)
         },
           stop('Unknown type: ', this$type)
         )
})

##
## rescale uniformized eigenvalues
##
setMethodS3("rescale", "Eigen", appendVarArgs=F, function(this, q.old, q.new) {
  msg(class=class(this)[1], "rescale", "rescale uniformized eigenvalues")
  stopifnot(this$type == 'P')
  rescale <- function(ei, q1, q2) (q2 - q1 + q1 * ei)/q2
  single <- length(this$one)   # number of single eigenvalues
  mult <- length(this$confl) # number of non single eigenvalues
  if(single > 1) this$one[2:single] <- rescale(this$one[2:single], q.old, q.new)
  if(mult > 0)  this$confl <- rescale(this$confl, q.old, q.new)
})
