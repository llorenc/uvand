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

library(R6)
##
## Compute confluent eigenvalues
##
Eigen <-
    R6Class("Eigen",
            portable = TRUE,
            cloneable = FALSE,
            private = list(
                ##
                ## Split ei.one in single and multiple eigenvalues
                ##
                compute.mult = function(ei.max, ei.max.mult) {
                    msg(class="Eigen", "compute.mult", "check for confluent eigenvalues")
                                        # browser()
                    ##
                    ## Supporting functions
                    ##
                    ## swap an eigenvalue pair
                    ei.swap <- function(i, j) {
                                        # browser()
                        if(i != j) { 
                            tmp <- self$one[i]
                            self$one[i] <- self$one[j]
                            self$one[j] <- tmp
                        }
                    } 
                    ## Make real a conjugate pair of eigenvalue i
                    set.real <- function(i) {
                                        # browser()
                        if(i < length(self$one)) {
                            ei.conj <- Conj(self$one[i])
                            for(j in (i+1):length(self$one)) {
                                if(ei.conj == self$one[j]) {
                                    msg(class="Eigen", 'set.real', 'set Im part to zero: ',
                                        self$one[i], ', ', self$one[j])
                                    self$one[i] <- Re(self$one[i])
                                    self$one[j] <- Re(self$one[j])
                                    return()
                                }
                            }
                        }
                        warning("Eigen$compute.mult: Conjugate pair not found: ", self$one[i])
                    }
                    ## Swap the conjugate of eigenvalue i from position >= j to position j
                    swap.conjugate <- function(i, j) {
                                        # browser()
                        if(i < length(self$one)) {
                            ei.conj <- Conj(self$one[i])
                            for(k in j:length(self$one)) {
                                if(ei.conj == self$one[k]) ei.swap(k, j)
                                return()
                            }
                        }
                        warning("Eigen$compute.mult: Conjugate pair not found: ", self$one[i])
                    }
                    ## Find the multiplicity of a real eigenvalue.
                    ## Eigenvalues are assumed sorted by its real part.
                    find.real.ei.multiplicity <- function(i) {
                                        # browser()
                        m <- 1    # multiplicity
                        if(i < length(self$one)) {
                            for(j in (i+1):length(self$one)) {
                                if(abs(Re(self$one[i])-Re(self$one[j])) < self$tol) {
                                    if(Im(self$one[j]) == 0) {
                                        ## found multiplicity
                                        ei.swap(i+m, j) # put confluent eigenvalues consecutive
                                        m <- m+1
                                    }
                                } else {
                                    if(m>1)
                                        msg(class="Eigen", 'find.real.ei.multiplicity',
                                            'found confl. eigenvalue: ', self$one[i])
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
                        if(i < length(self$one)) {
                            for(j in (i+1):length(self$one)) {
                                if(abs(Re(self$one[i])-Re(self$one[j])) < self$tol) {
                                    if(abs(Im(self$one[i])-Im(self$one[j])) < self$tol) {
                                        ## found multiplicity
                                        ei.swap(i+m, j)
                                        m <- m+1
                                    }
                                } else {
                                    if(m>1)
                                        msg(class="Eigen", 'find.complex.ei.multiplicity',
                                            'found confl. eigenvalue: ', self$one[i])
                                    break
                                }
                            }
                        }
                        return(m)
                    }
                    set.ei.mult <- function(i, m) {
                                        # browser()
                        self$confl[mult.i] <- self$one[i]
                        ## remove the confl. eig. from  self$one
                        if(i+m > length(self$one)) length(self$one) <- i-1
                        else self$one <- c(self$one[1:(i-1)], self$one[(i+m):length(self$one)])
                        ## adjust the multiplicity
                        if(ei.max.mult > 0) {
                            self$cmult[mult.i] <- m
                            if(m > ei.max.mult) {
                                msg(class="Eigen", 'set.ei.mult', 'reducing multiplicity of eigenvalue ',
                                    self$confl[mult.i], ': ', m, ' to ', ei.max.mult)
                                m <- ei.max.mult
                            }
                        }
                        msg(class="Eigen", 'set.ei.mult', 'setting confluent eigenvalue ',
                            self$confl[mult.i], ': ', m)
                        self$mult[mult.i] <- m
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
                    if(self$absorbing > 0) {
                        ## remove the absorbing eigenvalues
                        self$one <- self$one[self$absorbing:length(self$one)]
                    }
                    ## Remove Im part of almost real eigenvalues
                    for(i in 2:(length(self$one)-1)) {
                        ## if((Re(self$one[i]) != 0) && (Im(self$one[i]) != 0)
                        ##    && (abs(Im(self$one[i])) < self$tol)) set.real(i)
                        if((Re(self$one[i]) < self$tol) && (Im(self$one[i]) < self$tol)
                           && (abs(Im(self$one[i])) < self$tol)) set.real(i)
                    }
                    ## sort by real part to facilitate finding multiplicities
                    self$one <- self$one[order(Re(self$one), decreasing=TRUE)]
                    ## Check multiplicity
                    curr.i <- 2 # current ei (do not check for the lim. eigenvalue)
                    mult.i <- 1 # index to next free element in self$confl
                                        # browser()
                    repeat {
                        if(Im(self$one[curr.i]) == 0) {
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
                        if(curr.i >= length(self$one)) break
                    }
                                        # browser()
                    ## choose at most ei.max eigenvalues
                    if((ei.max > 0) && (self$number() > ei.max)) {
                        msg(class="Eigen", 'compute.mult', 'reducing mult. to ', ei.max)
                        single <- length(self$one)
                        mult <- length(self$confl)
                        choose.ei.one <- c(1)
                        choose.ei.confl <- c()
                        sndl <- self$second.largest.re()
                        if(sndl$confl) rand.samples <- single+sndl$id
                        else rand.samples <- sndl$id
                        set.seed(1)
                        rand.samples <- c(rand.samples, sample(2:(single+mult)))
                        num.of.ei <- 1
                        for(curr.sample in rand.samples) {
                            if(curr.sample <= single) {
                                if(!any(choose.ei.one == curr.sample)) {
                                    if(Im(self$one[curr.sample]) == 0) {
                                        num.of.ei <- num.of.ei + 1
                                        choose.ei.one <- c(choose.ei.one, curr.sample)
                                    } else {
                                        num.of.ei <- num.of.ei + 2
                                        choose.ei.one <- c(choose.ei.one, curr.sample,
                                                           find.conjugate(self$one, curr.sample))
                                    }
                                }
                            } else {
                                curr.sample <- curr.sample-single
                                if(!any(choose.ei.confl == curr.sample)) {
                                    if(Im(self$confl[curr.sample]) == 0) {
                                        num.of.ei <- num.of.ei + self$mult[curr.sample]
                                        choose.ei.confl <- c(choose.ei.confl, curr.sample)
                                    } else {
                                        num.of.ei <- num.of.ei + 2 * self$mult[curr.sample]
                                        choose.ei.confl <- c(choose.ei.confl, curr.sample,
                                                             find.conjugate(self$confl, curr.sample))
                                    }
                                }
                            }
                            if(num.of.ei >= ei.max) break
                        }
                        if(length(self$one) > length(choose.ei.one)) {
                            msg(class="Eigen", 'compute.mult', 'reducing single eigenvalues from ',
                                length(self$one), ' to ', length(choose.ei.one))
                            self$one <- self$one[choose.ei.one]
                        }
                        if(length(self$confl) > length(choose.ei.confl)) {
                            msg(class="Eigen", 'compute.mult', 'reducing confl eigenvalues from ',
                                length(self$confl), ' to ', length(choose.ei.confl))
                            self$confl <- self$confl[choose.ei.cofl]
                            self$mult <- self$mult[choose.ei.confl]
                        }
                        msg(class="Eigen", 'compute.mult', ' final number of eigenvalues ', num.of.ei)
                    }
                    if(self$type == 'Q') {
                        if((length(self$confl) > 0) &&
                           (abs(Re(self$confl[1])) < self$tol) && (abs(Im(self$confl[1])) < self$tol)) {
                            msg(class="Eigen", "compute.mult", "Removing eigenvalue=0 with mult. ", self$mult[1])
                            self$confl <- self$confl[-1]
                            self$zero.mult <- self$mult[1]
                            self$mult <- self$mult[-1]
                        }
                    }
                }
            ),
            public = list(
                type=NULL,   # eigenvalues type: 'Q', 'P'
                one=vector('numeric'),     # single eigenvalues
                confl=vector('numeric'), # confluent eigenvalues
                mult=vector('numeric'),   # multiplicity of confluent eigenvalues
                cmult=vector('numeric'), # multplicity before applying ei.max.mult
                evec=NULL,           # right eigenvectors
                lim.mult=integer(),   # multiplicity of the limit eigenvalue
                absorbing=integer(), # number of absorbing states
                zero.mult=0,
                tol=.Machine$double.eps^0.5,      # tolerance in some comparisons
                ##
                ## Compute the eigenvalues and its multiplicity
                ##
                initialize = function() {
                    msg(class="Eigen", "initialize", "new Object")
                },
                compute.eigenvalues = function(type, evec=FALSE, ei.max=0, ei.max.mult=0) {
                    msg(class="Eigen", "init", "Compute the eigenvalues")
                    N <- nrow(UCmatrix.Q)
                    stopifnot(N>0)
                    stopifnot((ei.max.mult==0)||(ei.max.mult>1))
                    self$tol <- (((1+.Machine$double.eps)^N)-1)^0.5
                    self$type <- type
                    if(evec) {
                        ei <- eigen(t(UCmatrix.Q), only.values=FALSE, sym=FALSE)
                        self$evec <- ei$vectors
                        self$one <- ei$values
                    } else {
                        self$one <- eigen(UCmatrix.Q, only.values=TRUE, sym=FALSE)$values
                    }
                    self$lim.mult <- 1
                    if(self$type == 'Q') {
                        if(any(Re(self$one) > 0)) {
                            gtz.idx <- which(Re(self$one) > 0)
                            msg(class="Eigen", "init", "Found ",
                                length(gtz.idx), " eigenvalues with real part > 0: ",
                                paste(sep='', Re(self$one[gtz.idx]), collapse=','))
                            msg(class="Eigen", "init", "Setting eigenvalues with real part > 0 to 0 ")
                            self$one[gtz.idx] <- 0
                        }
                        ## First eigenvalue must be 0
                        self$one <- self$one[length(self$one):1] 
                        if(!is.null(self$evec)) self$evec <- self$evec[,ncol(self$evec):1]
                        if(!isTRUE(all.equal(abs(self$one[1]), 0, tol=self$tol))) {
                            warning('Dominant eigenvalue is not 0:', self$one[1])
                        }
                        self$one[1] <- 0
                        ## Check for absorbing states
                        abs.s <- which(sapply(1:N, function(i) all(UCmatrix.Q[i,]==0)))
                    } else if(self$type == 'P') {
                        ## First eigenvalue must be 1
                        if(!isTRUE(all.equal(abs(self$one[1]), 1, tol=self$tol))) {
                            warning('Dominant eigenvalue is not 1:', self$one[1])
                        }
                        self$one[1] <- 1
                        ## Check for absorbing states
                        abs.s <- which(sapply(1:N, function(i) UCmatrix.Q[i,i]==1))
                    } else {
                        stop("Nor stochastic or infinitesimal generator matrix?")
                    }
                    self$absorbing <- length(abs.s)
                    if((self$absorbing > 0) && (length(intersect(c(1:(abs.s[1]-1)), abs.s)) > 0))
                        warning('There are absorbing states and Q is not in canonical form.\n')
                    if(evec == FALSE) private$compute.mult(ei.max, ei.max.mult)
                },
                print = function(...) {
                    if(debug) cat('---------------\n')
                    msg("", "Eigen", "summary")
                    msg("", " Matrix type", self$type)
                    msg("", " number of single eigenvalues", length(self$one))
                    msg("", " number of confluent eigenvalues", length(self$confl))
                    if(length(self$confl)>0) {
                        msg("", " computed multiplicities", paste(sep='', self$cmult, collapse=','))
                        msg("", " applied multiplicities", paste(sep='', self$mult, collapse=','))      
                    }
                    msg("", " right eigenvectors", !empty(self$evec))
                    msg("", " multiplicity of the limit eigenvalue", self$lim.mult)
                    msg("", " number of absorbing states", self$absorbing)
                    msg("", " tolerance", self$tol)
                },
                copy = function(ei) {
                    self$type <- ei$type
                    self$one <- ei$one
                    self$confl <- ei$confl
                    self$mult <- ei$mult
                    self$cmult <- ei$cmult
                    self$evec <- ei$evec
                    self$lim.mult <- ei$lim.mult
                    self$absorbing <- ei$absorbing
                    self$zero.mult <- ei$zero.mult
                    self$tol <- ei$tol
                },
                ##
                ## Return the second smallest eigenvalue in modulus (id, value and confl indicator)
                ##
                second.smallest.abs = function() {
                    find.min.abs <- function(x) {
                        id <- which.min(abs(x))
                        list(id=id, value=x[id])
                    }
                    ##
                    if(length(self$one) == 2) {
                        min.one <- list(value=abs(self$one[2]), id=2)
                    } else if(length(self$one) > 2) {
                        min.one <- find.min.abs(self$one[2:length(self$one)])
                        min.one$id <- min.one$id+1
                    } else min.one <- NULL
                    ##
                    if(length(self$confl) > 0) {
                        min.confl <- find.min.abs(self$confl)
                        if(is.null(min.one) || (abs(min.confl$value) < abs(min.one$value)))
                            return(c(min.confl, confl=TRUE))
                    }
                    return(c(min.one, confl=FALSE))
                },
                ##
                ## Return the second largest real part (id, value and confl indicator)
                ##
                second.largest.re = function() {
                    find.max.re <- function(x) {
                        id <- which.max(Re(x))
                        list(id=id, value=x[id])
                    }
                    ##
                    if(length(self$one) == 2) {
                        max.one = list(value=self$one[2], id=2)
                    } else if(length(self$one) > 2) {
                        max.one <- find.max.re(self$one[2:length(self$one)])
                        max.one$id <- max.one$id+1
                    } else max.one <- NULL
                    ##
                    if(length(self$confl) > 0) {
                        max.confl <- find.max.re(self$confl)
                        if(is.null(max.one) || (Re(max.confl$value) > Re(max.one$value)))
                            return(c(max.confl, confl=TRUE))
                    }
                    return(c(max.one, confl=FALSE))
                },
                ##
                ## Number of eigenvalues
                ##
                number = function() {
                    if(length(self$confl) == 0) {
                        num <- length(self$one)
                    } else {
                        num <- length(self$one) + sum(self$mult)
                    }
                    return(num)
                },
                ##
                ## Number of eigenvalues
                ##
                plot = function(lim=FALSE) {
                    ## If lim=TRUE plot the limit eigenvalue
                    ei <- vector('complex')
                    if(length(self$one) > 0) {
                        ei <- self$one
                        if(!lim) ei <- ei[-1]
                    }
                    if(length(self$confl) > 0) {
                        ei <- c(ei, self$confl)
                    }
                    plot(Re(ei), Im(ei))
                },
                ##
                ## Toggle P/Q-eigenvalues to Q/P-eigenvalues
                ##
                toggle.Q.P.map = function(qn) {
                    msg(class="Eigen", "toggle.Q.P.map")
                    stopifnot(!is.null(qn))
                    single <- length(self$one)   # number of single eigenvalues
                    mult <- length(self$confl) # number of non single eigenvalues
                    switch(self$type,
                           Q = {
                               msg(class="Eigen", "toggle.Q.P.map","toggle Q-eigenvalues to P-eigenvalues")
                               self$type <- 'P'
                               if(single) self$one <- 1 + (self$one/qn)
                               if(mult) self$confl <- 1 + (self$confl/qn)     
                           },
                           P = {
                               msg(class="Eigen", "toggle.Q.P.map","toggle P-eigenvalues to Q-eigenvalues")
                               self$type <- 'Q'
                               if(single) self$one <- qn * (self$one-1)
                               if(mult) self$confl <- qn * (self$confl-1)
                           },
                           stop('Unknown type: ', self$type)
                           )
                },
                ##
                ## rescale uniformized eigenvalues
                ##
                rescale = function(q.old, q.new) {
                    msg(class="Eigen", "rescale", "rescale uniformized eigenvalues")
                    stopifnot(self$type == 'P')
                    rescale <- function(ei, q1, q2) (q2 - q1 + q1 * ei)/q2
                    single <- length(self$one)   # number of single eigenvalues
                    mult <- length(self$confl) # number of non single eigenvalues
                    if(single > 1) self$one[2:single] <- rescale(self$one[2:single], q.old, q.new)
                    if(mult > 0)  self$confl <- rescale(self$confl, q.old, q.new)
                }
            ))
## end of UC class

