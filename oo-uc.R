## Copyright (c) 2013 Llorenç Cerdà-Alabern,
## http://personals.ac.upc.edu/llorenc This file is free software: you
## can redistribute it and/or modify it under the terms of the GNU
## Affero Public License as published by the Free Software Foundation,
## either version 3 of the License, or (at your option) any later
## version.  oo-uc.R is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero Public License for more details.  You should have
## received a copy of the GNU Affero Public License along with
## oo-uc.R.  If not, see <http://www.gnu.org/licenses/>.
##########################################################################################
## UC class
##########################################################################################
library(reshape) # melt
library(ggplot2)

setConstructorS3("UC", function(type=NULL,  method=NULL, states=NULL, j=NULL, i=NULL,
                                states.names=NULL,
                                unif=FALSE, qii.max=NULL, qn=NULL,
                                choose.samples=FALSE, lim=FALSE, rm.lim.eq=FALSE,
                                ei.max=0, ei.max.mult=0, use.matlab=FALSE) {
  if(!is.null(j) && !is.null(states.names) && (length(states.names) != length(j))) {
    warning("UC: Incorrect number of states.names")
    states.names <- NULL
  } else {
    states.names <- as.character(states.names)
  }
  if(!is.null(type)) {
    mc.type = switch(type,
      P = 'dtmc',
      Q = 'ctmc',
      stop('Unknown type')
      )
  } else {
    mc.type = NULL ;
  }
  if(rm.lim.eq && lim==FALSE) lim <- TRUE
  ##
  extend(Object(), "UC",
         method=method,          # UC solution method
         states=states,          # states of the chain
         mc.type=mc.type,        # Markov chain type: 'ctmc', 'dtmc'
         matrix.type=type,       # Matrix type: 'Q', 'P'
         states.names=states.names, # Names of the Markov chain states
         unif=unif,              # true if the UC are solved using uniformization
         lim=lim,                # true if the UC are solved using the lim. dist.
         rm.lim.eq=rm.lim.eq,    # true if the lim. dist. eq. is removed.
         i=i,                    # initial state
         j=j,                    # target states
         qii.max=qii.max,        # max_i |Qii|
         qn=qn,                  # uniformization parameter
         choose.samples=choose.samples, # if TRUE, samples are choosen evely spaced in log scale.
         max.power=NULL,         # maximum power used in the Vand. system when choose.samples=TRUE
         ei=Eigen(),             # eigenvalues
         ei.max=ei.max,          # Take the ei.max largest Re eigenvalues (0 for unlimited)
         ei.max.mult=ei.max.mult,# maximum multiplicity of the eigenvalues (0 for unlimited)
         use.matlab=use.matlab,  # TRUE it the eigenvalues are computed with matlab
         coef=NULL,              # undetermined coefficients (UC)
         kappa=NA,               # condition number of the matix used to solve the coefs.
         matrix.dim=NULL,        # matrix dimensions
         matrix.zeros=NULL,      # Number of matrix elements = 0 in the equations.
         matrix.bytes=NULL,      # Matrix storage in Mb.
         time=list()             # computation times
         );
})


##########################################################################################
## UC public methods
##########################################################################################
##
## read matlab data
##
setMethodS3("read.uc.from.matlab.file", "UC", appendVarArgs=FALSE, function(this, fname) {
  msg(class=class(this)[1], "read.uc.from.matlab.file", "reading file ", fname)
  mat <- readMat(fname)
  for(v in rownames(mat$uc)) {
    if((v != 'ei') && any(names(this) == v)) {
      this[[v]] <- mat$uc[[v,1,1]] ;
    }
  }
  for(v in rownames(mat$uc[['ei',1,1]])) {
    if(any(names(this$ei) == v)) {
      this$ei[[v]] <- mat$uc[['ei',1,1]][[v,1,1]] ;
    }
  }
  ## for(v in names(mat)) {
  ##   switch(gsub("^(..).*$", "\\1", v),
  ##          uc = {
  ##            v.f <- gsub("^uc.", "", v)
  ##            this[[v.f]] <- as.matrix(mat[[v]])   
  ##          },
  ##          ei = {
  ##            v.f <- gsub("^ei.", "", v)
  ##            this$ei[[v.f]] <- as.matrix(mat[[v]])
  ##          },
  ##          stop('Unknown prefix')
  ##          )
  ## }
  ##
  if(debug) {
    msg(class=class(this)[1], "solve.uc", 'done')
    this$print()
  }
})

##
## UC probability methods
##
## Evaluate pn(t)
setMethodS3("prob", "UC", appendVarArgs=FALSE, function(this, t) {
  nc <- ncol(this$coef)
  if(is.null(this$coef)) {
    prob <- data.frame()
    warning("empty coefficients, ", this$title())
    prob <- data.frame()
  } else {
    prob <-
      switch(this$mc.type,
             ctmc = data.frame(matrix(this$p.ctmc(t), nr=length(t), nc=nc)),
             dtmc = data.frame(matrix(this$p.dtmc(t), nr=length(t), nc=nc)))
    colnames(prob) <- this$get.states.names()
  }
  return(prob)
})

setMethodS3("print", "UC", appendVarArgs=FALSE, function(this, t) {
  if(debug) cat('---------------\n')
  msg("", "UC", "summary")
  msg("", " UC solution method", this$method)
  if(empty(this$coef)) msg("", " UC coef", "failure")
  else {
    ei.max <- max(Re(this$coef))
    ei.min <- min(Re(this$coef))
    msg("", " UC coef", "max(Re):", ei.max, ", min(Re):", ei.min)
  }
  msg("", " kappa", this$kappa)
  msg("", " states of the chain", this$states)
  msg("", " Markov chain type", this$mc.type)
  msg("", " Matrix type", this$matrix.type)
  if(!is.null(this$states.names))
    msg("", " states names", paste(sep='', this$states.names, collapse=','))
  msg("", " matrix dimensions", paste(sep='', this$matrix.dim, collapse=','))
  msg("", " matrix zeros", this$matrix.zeros)
  msg("", " matrix storage", this$matrix.bytes, ' bytes')
  msg("", " solved using uniformization", this$unif)
  if(this$unif) {
    msg("", " uniformization parameter", this$qn)
    msg("", " qii.max", this$qii.max)
  }
  msg("", " choose samples", this$choose.samples)
  if(!is.null(this$max.power)) msg("", " max power", this$max.power)
  msg("", " time constant", this$time.constant())
  msg("", " solved using the lim. dist.", this$lim)
  msg("", " remove lim. eq", this$rm.lim.eq)
  if(length(this$i) > 1) {
    if(length(this$i) <= 5) info <- paste(sep='', sprintf('%.2e',this$i), collapse=',')
    else info <- paste(sep='', c(sprintf('%.2e', this$i[1:5]), '...'), collapse=',')
    msg("", " initial state", info)
  } else  msg("", " initial state", this$i)
  msg("", " target states", paste(sep='', this$j, collapse=','))
  msg("", " running times",)
  t.names <- c('sys.self', 'elapsed')
  if(!is.null(this$time)) {
    for(t in names(this$time)) {
      cat(sep='', '  ', t, ': ',
          paste(t.names, format(this$time[[t]][t.names]), collapse=', '), '\n')
      ##    for(n in c('sys.self', 'elapsed')) cat(n, this$time[[t]][[n]]))
      ## cat('\n')
      ## cat(paste(sep='', format(this$time[[t]]), collapse=', '), '\n')
    }
  }
  msg("", " eigenvalues maximum multiplicity", this$ei.max.mult)
  msg("", " eigenvalues maximum number", this$ei.max)
  msg("", " eigenvalues computed with matlab", this$use.matlab)
  this$ei$print()
})

##
## Coefficients are computed using the eigenvectors of Q.
##
setMethodS3("solve.uc.eigenvectors", "UC", appendVarArgs=FALSE, priv=TRUE, function(this) {
  msg(class=class(this)[1], "solve.uc.eigenvectors",'Solving the eigenvectors')
  if(length(this$i) == 1) {
    x0 <- rep(0,len=nrow(this$ei$evec)) ; x0[this$i] <- 1
  } else {
    x0 <- this$i
  }
  s <- try.solve.kappa(this$ei$evec, x0)
  this$kappa <- s$kappa
  this$matrix.dim <- dim(this$ei$evec)
  this$matrix.zeros <- length(which(this$ei$evec == 0))
  this$matrix.bytes <- object.size(this$ei$evec)
  if(!empty(s$x)) {
    this$coef <- sapply(this$j, function(n) this$ei$evec[n,] * s$x)
    colnames(this$coef) <- this$get.states.names()
  }
  msg(class=class(this)[1], 'solve.uc.eigenvectors', 'Removing eigenvectors')
  this$ei$evec <- NULL # eigenvectors are not needed anymore
  ## correct the eigenvalues if computed with an uniformized matrix
  if(this$unif && (uc$ei$type == 'P')) uc$ei$toggle.Q.P.map(uc$qii.max)
})

##
## Coefficients are computed solving a Vandermonde system
##
setMethodS3("compute.qn", "UC", appendVarArgs=FALSE, priv=TRUE, function(this) {
  msg(class=class(this)[1], "compute.qn",'Computing the uniformization parameter')
  if(is.null(this$qn)) {
#browser()
    ## most negative eigenvalue
    abs.min.ei <- abs(min(Re(this$ei$confl), Re(this$ei$one)))
    ## second largest eigenvalue in modulus
    ei2 <- abs(this$ei$second.largest.re()$value)
    N <- this$ei$number()
    if(this$ei$type == 'P')
      min.qn <- (this$qii.max - this$qii.max * ei2)/(1-.Machine$double.eps^(1/N))
    else
      min.qn <- ei2/(1-.Machine$double.eps^(1/as.numeric(N)))
    this$qn <- max(this$qii.max, abs.min.ei, min.qn)
    if(this$qn > ei2*N) {
      this$qn <- max(this$qii.max, ei2*N)
    }
  } else msg(class=class(this)[1], "compute.qn", 'qn already initialized')
})

##
## Return the eigenvalue with the second lagest real part and its coefficient
## (id, value and confl indicator, coefficient)
##
setMethodS3("second.largest.re", "UC", appendVarArgs=FALSE, function(this) {
  max.re <- this$ei$second.largest.re()
  if(max.re$confl) {
    mult <- this$ei$mult[max.re$id]
    id <- length(this$ei$one) + sum(this$ei$mult[1:max.re$id])
    coef <- this$coef[(id-mult+1):id]
    return(c(max.re, mult=mult, coef=coef))
  }
  return(c(max.re, coef=this$ei$one[max.re$id]))
})

##
## Return the eigenvalue with the second lagest real part and its coefficient
## (id, value and confl indicator, coefficient)
##
setMethodS3("time.constant", "UC", appendVarArgs=FALSE, function(this, alpha=0.5) {
  switch(this$ei$type,
         P = alpha/(1-abs(this$ei$second.largest.re()$value)),
         Q = (alpha*this$qn)/abs(this$ei$second.smallest.abs()$value),
         stop())
})
  
##
## Return the eigenvalue with the second lagest real part and its coefficient
## (id, value and confl indicator, coefficient)
##
setMethodS3("compute.samples", "UC", appendVarArgs=FALSE, function(this) {
  msg(class=class(this)[1], "compute.samples",'Compute the samples')
  int.log.scale <- function(max, N) {
    n <- ceiling(N/log10(max))
    c(1:(n-1), round(10^(seq(log10(n),log10(max), len=N-n+1))))
  }
  if(length(this$ei$confl) > 0) {
    warning(paste(class(this)[1],
                  ': option choose.samples not yet implemented with confluent eigenvalues'))
    return(NULL)
  }
  if(is.null(this$max.power)) 
    this$max.power <- max(ceiling(this$time.constant()), this$ei$number()-1)
  if(this$max.power < this$ei$number()) return(NULL)
  else return(int.log.scale(this$max.power, this$ei$number()-1))
})

##
## Coefficients are computed solving a Vandermonde system
##
setMethodS3("solve.uc.vandermonde", "UC", appendVarArgs=FALSE, priv=TRUE, function(this,
                                                             b, lim=NULL, powers=NULL) {
  ##:ess-bp-start::browser:##
  #browser()##:ess-bp-end:##
  msg(class=class(this)[1], "solve.uc.vandermonde",'Solving the vandermonde system')
  if(this$rm.lim.eq) {
    msg(class=class(this)[1], "solve.uc.vandermonde",'remove the limit eigenv equation')
    tmp <- this$ei$one[1]
    this$ei$one <- this$ei$one[-1]
    V <- this$confl.vand.mat(nrow(b), powers)
    this$ei$one <- c(tmp, this$ei$one)
  } else {
    V <- this$confl.vand.mat(nrow(b), powers)
    if(!is.null(lim)) {
      msg(class=class(this)[1], "solve.uc.vandermonde",'add the limit distribution')
      V <- rbind(V, c(1, rep(0, ncol(V)-1)))
      b <- rbind(b, lim)
    }
  }
  msg(class=class(this)[1], "solve.uc.vandermonde",'Solve the Vandermonde system')
  if(this$use.matlab) {
    msg(class=class(this)[1], "solve.uc.vandermonde",'Solving with matlab')
    # Start the matlab server on the same machine
    ## nz <- which(V!=0, arr.ind=T)
    ## setVariable(matlab, i=nz[,1])
    ## setVariable(matlab, j=nz[,2])
    ## setVariable(matlab, s=V[nz])
    ## setVariable(matlab, m=nrow(V))
    ## setVariable(matlab, n=ncol(V))
    ## setVariable(matlab, nzmax=length(nz))
    setVariable(matlab, b=b)
    setVariable(matlab, V=V)
    ## res <- evaluate(matlab, 'VS=sparse(i,j,s,m,n,nzmax);')
    res <- evaluate(matlab, 'V=sparse(V);')
    res <- evaluate(matlab, paste(sep='', 's=b\\V;'))
    this$coef <- getVariable(matlab, c("s"))$s
    ## colnames(this$coef) <- this$get.states.names()
  } else {
    s <- try.solve.kappa(V, b)
    if(!empty(s$x)) {
      this$coef <- s$x
      colnames(this$coef) <- this$get.states.names()
    }
    this$kappa <- s$kappa
  }
  this$matrix.dim <- dim(V)
  this$matrix.zeros <- length(which(V == 0))
  this$matrix.bytes <- object.size(V)
  if(!empty(this$coef) && this$rm.lim.eq) {
    this$coef <- rbind(lim, this$coef)
  }
})

##
## Compute the confluent Vandermonde matrix.
##
## setMethodS3("confl.vand.mat", "UC", appendVarArgs=FALSE, priv=TRUE, function(this, nr, powers) {
##   msg(class=class(this)[1], "confl.vand.mat", 'Build the Vandermonde matrix')
##   if(is.null(powers)) V <- vand.matrix(this$ei$one, nr)
##   else V <- vand.matrix.powers(this$ei$one, powers)
##   ##
##   if(length(this$ei$confl) > 0) {
##     stopifnot(is.null(powers)) # not yet implemented
##     for(i in 1:length(this$ei$confl)) {
##       V <-
##         cbind(V, switch(this$ei$type,
##                         Q = vand.ctmc.confl.block(this$ei$confl[i], this$ei$mult[i], nr),
##                         P = vand.dtmc.confl.block(this$ei$confl[i], this$ei$mult[i], nr)))
##     }
##   }
##   return(V)
## })
setMethodS3("confl.vand.mat", "UC", appendVarArgs=FALSE, priv=TRUE, function(this, nr, powers) {
  msg(class=class(this)[1], "confl.vand.mat", 'Build the Vandermonde matrix')
  nc <- length(this$ei$one)
  if(length(this$ei$confl) > 0) {
    nc <- nc + sum(this$ei$mult)
  }
  if(is.null(powers)) V <- vand.matrix(this$ei$one, nr, nc)
  else V <- vand.matrix.powers(this$ei$one, powers, nc)
  ##
  ##browser()
  if(length(this$ei$confl) > 0) {
    stopifnot(is.null(powers)) # not yet implemented
    multi <- length(this$ei$one)
    for(i in 1:length(this$ei$confl)) {
      V[,(multi+1):(multi+this$ei$mult[i])] <-
        switch(this$ei$type,
               Q = vand.ctmc.confl.block(this$ei$confl[i], this$ei$mult[i], nr),
               P = vand.dtmc.confl.block(this$ei$confl[i], this$ei$mult[i], nr))
      multi <- multi+this$ei$mult[i]
    }
  }
  return(V)
})

## Correct the coefficients of the confluent eigenvalues
setMethodS3("correct.confl.coef", "UC", appendVarArgs=FALSE, function(this) {
  if(empty(this$ei$confl) || empty(this$coef)) return()
  ##
  msg(class=class(this)[1], "correct.confl.coef",
      "Correcting the coefficients of the confluent eigenvalues")
  single <- length(this$ei$one) # number of single eigenvalues
  j <- single+1
  for(n.ei in 1:length(this$ei$confl)) {
    mult <- this$ei$mult[n.ei]
    if(this$ei$confl[n.ei] == 0) {
      this$coef[(j+1):(j+mult-1),] <-
        (this$coef[(j+1):(j+mult-1),] * this$qn^(1:(mult-1))) / factorial(1:(mult-1))
    } else {
      for(m in 1:(mult-1)) {
        if(m < (mult-1)) {
          for(k in (m+1):(mult-1)) {
            this$coef[j+m,] <- this$coef[j+m,] + this$coef[j+k,] * q.n.coefs(k, m)
          }
        }
        this$coef[j+m,] <- this$coef[j+m,] * (this$qn * this$ei$confl[n.ei])^m
      }
    }
    j <- j + mult
  }
})

##
## return the coefficients q_i^(m) for a solution with confluent eigenvalues.
##
q.n.coefs.matrix <- list()
q.n.coefs <- function(m, i) {
  if((length(q.n.coefs.matrix) < m) ||
     is.null(q.n.coefs.matrix[[m]])) {
    ## cat('create', n, '\n')
    q.n.coefs.matrix[[m]] <<- vector(len=m)
  }
  if(q.n.coefs.matrix[[m]][i] == FALSE) {
    if(i == 1) {
      q.n.coefs.matrix[[m]][i] <<- 1
    } else {
      sum <- 0
      for(k in (i-1):(m-1))  {
        sum <- sum + choose(m-1, k) * q.n.coefs(k , i-1)
      }
      q.n.coefs.matrix[[m]][i] <<- sum
    }
  }
  return(q.n.coefs.matrix[[m]][i])
}

setMethodS3("plot", "UC", appendVarArgs=FALSE, function(this,
                            t, comp=NULL, x.comp=NULL, logy=TRUE, logx=FALSE,
                            title=NULL, ylim=NULL) {
  if(is.null(title)) title <- this$title()
  plot.facet.grid.row.col(data.main=list(list(prob=this$prob(t))), x=t,
                          data.comp=comp, x.comp=x.comp,
                          field.x='t', field.y='prob', field.legend='n',
                          logy=logy, logx=logx, title=title, #type='single',
                          ylim=ylim)
})

##
## UC info methods
##
setMethodS3("title", "UC", appendVarArgs=FALSE, function(this, show.states=TRUE) {
  if(this$unif) type <- 'CTMC, Unif.'
  else type <- switch(this$mc.type, ctmc = 'CTMC', dtmc = 'DTMC')
  method <- switch(this$method,
                   vand = 'VAND',
                   evec = 'Eigenvectors')
  title <- paste(type, method)
  if(this$lim) title <- paste(sep='', title,  ', limit dist')
  if(this$rm.lim.eq) title <- paste(sep='', title,  ', rm limit eq')
  if(show.states) title <- paste(sep='', title, ', N=', this$states)
  return(title)
})

setMethodS3("get.states.names", "UC", appendVarArgs=FALSE, function(this) {
  if(empty(this$states.names)) as.character(this$j)
  else this$states.names
})

#############################################################################
## Private methods for UC
#############################################################################
########################################
## UC probability methods
########################################
## Evaluate p1n(t) of a CTMC Markov chain. ei are the eigenvalues of Q as
## returned by eigen.mult. t and j can be vectors.
setMethodS3("p.ctmc", "UC", appendVarArgs=FALSE, private=TRUE, function(this, t) {
  single <- length(this$ei$one) # number of single eigenvalues
  mult <- length(this$ei$confl) # number of non single eigenvalues
  p <- sapply(1:ncol(this$coef), function(j) {
    sapply(t, function(tt) sum(this$coef[1:single,j] * exp(tt * this$ei$one)))
  })
  if(mult > 0) { # evaluate non single eigenvalues
    p <- p + sapply(1:ncol(this$coef), function(j) {
      sapply(t, function(tt) {
        i <- single+1
        sum(sapply(1:mult, function(n.ei) {
          s <- sum(this$coef[i:(i+this$ei$mult[n.ei]-1),j] *
              tt^c(0:(this$ei$mult[n.ei]-1))) *
                exp(tt * this$ei$confl[n.ei])
          i <<- i + this$ei$mult[n.ei]
          return(s)
        }))
      })
    })
  }
  return(Re(p))
})

## Evaluate p1n(t) of a DTMC Markov chain. ei are the eigenvalues of P as
## returned by eigen.mult.
setMethodS3("p.dtmc", "UC", appendVarArgs=FALSE, private=TRUE, function(this, t) {
  single <- length(this$ei$one) # number of single eigenvalues
  mult <- length(this$ei$confl) # number of non single eigenvalues
  p <- sapply(1:ncol(this$coef), function(j) {
    sapply(t, function(t) sum(this$coef[1:single,j] * this$ei$one^t))
  })
  if(mult > 0) { # evaluate non single eigenvalues
    p <- p + sapply(1:ncol(this$coef), function(j) {
      sapply(t, function(t) {
        i <- single+1
        sum(sapply(1:mult, function(n.ei) {
          s <- if(abs(this$ei$confl[n.ei]) == 0) {
            if(t < this$ei$mult[n.ei]) this$coef[i+t,j]
            else 0
          } else {
            sum(this$coef[i:(i+this$ei$mult[n.ei]-1),j] *
                t^c(0:(this$ei$mult[n.ei]-1))) *
                  this$ei$confl[n.ei]^t
          }
          i <<- i + this$ei$mult[n.ei]
          return(s)
        }))
      })
    })
  }
  return(Re(p))
})

##
## helping functions
##
##
## Build a data.frame suitable for plotting.
##
build.df <- function(data.l, field.y, x=NULL, field.x=NULL, field.legend='n', type='single',
                     ylim=NULL, logy=TRUE, logx=FALSE, mark.na=NULL, depth=0) {
  build.df.single <- function(data, x, field.x, field.y, field.legend) {
    if(length(data[[field.y]]) == 0) return(NULL)
    ## msg("", "build.df.single", 'building df for: ', names(data))
    df <- data.frame(cbind(t=x, p=melt(data[[field.y]])))
    df <- rename(df, c(t=field.x, p.variable=field.legend, p.value=field.y))
  }
  build.df.row.col <- function(data, x, field.x, field.y, field.legend) {
    if(is.null(data$facet.row) || is.null(data$facet.col)) {
      stop('is.null(data$facet.row) || is.null(data$facet.col):',
           'names(data)=', paste(names(data), col='-', sep=''),
           ', field.x=', field.x, ', field.y=', field.y,
           ', facet.row=', data$facet.row,', facet.col=', data$facet.col, '\n')
    }
    if(length(data[[field.y]]) == 0) {
      warning('Field ', field.y, 'empty',
              ', facet.row=', data$facet.row,', facet.col=', data$facet.col)
      return(NULL)
    }
    msg("", "build.df", 'Found data for: ',
         'facet.row=', data$facet.row,', facet.col=', data$facet.col)
    df <- data.frame(cbind(t=x,
                           p=melt(data[[field.y]]),
                           facet.row=data$facet.row,
                           facet.col=data$facet.col))
    ##browser()
    df <- rename(df, c(t=field.x, p.variable=field.legend, p.value=field.y))
  }
  stopifnot(depth<100)
  ## msg("", "build.df", 'building df for:\n',
  ##     ' names(data)=', paste(names(data), col='-', sep=''),
  ##     ', field.x=', field.x, ', field.y=', field.y,
  ##     ', facet.row=', data$facet.row,', facet.col=', data$facet.col, ', type=', type)
  res <- list(df=data.frame(), title=NULL)
  if(is.null(field.legend)) field.legend <- 'n'
  repeat {
    ##    if(!is.null(data.l[[as.character(field.y)]])) {
    if(any(names(data.l) == field.y) && (class(data.l[[field.y]]) == 'data.frame')) {
      if(is.null(x)) x <- data.l[[field.x]]
      ## if(logy) data.l[[field.y]] <- log10(data.l[[field.y]])
      ## if(logx) x <- log10(x)
      ## if(!is.null(ylim)) {
      ##   data.l[[field.y]][data.l[[field.y]] < ylim[1]] <- ylim[1]
      ##   data.l[[field.y]][data.l[[field.y]] > ylim[2]] <- ylim[2]
      ## }
      df <- switch(type,
                   single = build.df.single(data.l, x, field.x, field.y, field.legend),
                   facet  = build.df.row.col(data.l, x, field.x, field.y, field.legend),
                   stop("Unknown type"))
      if(is.null(mark.na)) df <- na.omit(df)
      else {
        id.na <- which(is.na(df[,field.y]))
        id.ok <- which(!is.na(df[,field.y]))
        df[,field.y][id.ok] <- NA
        df[,field.y][id.na] <- mark.na
      }
      return(list(df=na.omit(df), title=data.l$title))
    } else {
      ##browser()
      r.try <- NULL
      if((class(data.l) == 'list') && (length(data.l) > 0)) {
        ## cat(depth, '-', length(data.l), '-', class(data.l), '-', names(data.l), '\n')
        r.try <- build.df(data.l=data.l[[1]], field.y=field.y, x=x, field.x=field.x,
                          field.legend=field.legend, type=type, ylim=ylim, logy=logy,
                          logx=logx, mark.na=mark.na, depth=depth+1)
        if(!is.null(r.try)) {
          if(is.null(res$title)) res$title <- r.try$title
          else if(!is.null(r.try$title))  res$title <- paste(sep='', res$title, '\n', r.try$title)
          res <- list(df=rbind(res$df, r.try$df), title=res$title)
        }
        data.l <- data.l[-1]
      } else return(res)
    }
  }
}

build.aes <- function(logx, field.x, logy, field.y, field.legend=NULL, linetype=F)
{
  set.log <- function(log, var) {
    if(log) paste('log10(', var, ')', sep='')
    else var
    ##var
  }
  ##
  if(is.null(field.legend))
    eval(parse(text=paste('aes(x=', set.log(logx, field.x),
                 ', y=', set.log(logy, field.y), ')', sep='')))
  else {
    if(linetype) aestype <- 'linetype'
    else aestype <- 'colour'
    eval(parse(text=paste('aes(x=', set.log(logx, field.x),
                 ', y=', set.log(logy, field.y),
                 ', ', aestype, '=', field.legend, ')', sep='')))
  }
}

plot.facet.grid.row.col <-
  function(data.main, x=NULL, data.comp=NULL, x.comp=NULL, field.x, field.y, field.legend=NULL,
           grid=TRUE, logy=TRUE, logx=FALSE, title=NULL, xlab='time (s)', ylab='Probability',
           xlim=NULL, ylim=NULL, mark.na=NULL, file=NULL, font.size=12, line.size=NULL,
           kpos='right', kjus=c(1, 1), w=12, h=10, geom='line', scales="free_y", gray=F,
           linetype=F, addplot=NULL)
{
  msg("", "plot.facet.grid.row.col", 'building main df')
  stopifnot(class(data.main) == 'list')
  stopifnot(length(data.main) > 0)
  if(length(data.main) == 1) type='single'
  else type='facet'
  msg("", "plot.facet.grid.row.col", 'plot type=', type)
  df <- build.df(data.l=data.main, field.y=field.y, x=x, field.x=field.x,
                 field.legend=field.legend, type=type, ylim=ylim, logy=logy, logx=logx)
  if(!is.null(mark.na))
    df.na <- build.df(data.l=data.main, field.y=field.y, x=x, field.x=field.x,
                      field.legend=field.legend, type=type, ylim=ylim,
                      logy=logy, logx=logx, mark.na=1)
  ##browser()
  if(logx) xlab <- paste(xlab, ', log10 scale', sep='')
  if(logy) ylab <- paste(ylab[1], ', log10 scale', ylab[-1], sep='')
  ##
  g <- ggplot(df$df) +
    build.aes(logx, field.x, logy, field.y, field.legend, linetype) + xlab(xlab) + ylab(ylab)
  if(geom == 'line') {
    if(is.null(line.size)) g <- g + geom_line()
    else g <- g + geom_line(size=line.size)
  } else {
    if(is.null(line.size)) g <- g + geom_point()
    else g <- g + geom_point(size=line.size)
  }
  if(type == 'facet') {
    g <- g + facet_grid(facet.row~facet.col, scales=scales)
  }
  if(!is.null(mark.na) && (length(df.na[[1]][[1]]) > 0))
    g <- g + geom_point(data=df.na$df)#, colour='red') 
  ##
  if(!is.null(data.comp)) {
    if(is.null(x.comp)) x.comp <- x
    msg("", "plot.facet.grid.row.col", 'building comp df')
    df.comp <- build.df(data.l=data.comp, field.y=field.y, x=x.comp, field.x=field.x,
                        field.legend=field.legend, type=type,
                        ylim=ylim, logy=logy, logx=logx)$df
    msg("", "plot.facet.grid.row.col", 'done')
    if(linetype) {
      g <- g + geom_line(data=df.comp, size=1)
    } else {
      g <- g + geom_line(data=df.comp, linetype="dashed", size=1)
    }
  }
  if(!is.null(xlim)) g <- g + scale_x_continuous(lim=xlim)
  if(!is.null(ylim)) g <- g + scale_y_continuous(lim=ylim)
  if(is.null(title)) title <- df$title
  if(!is.null(title))
    g <- g + opts(title=title,
                  plot.title=theme_text(size=8,face="bold"))
  msg("", "plot.facet.grid.row.col", 'print ggplot')
  ## browser()
  if(gray && !is.null(field.legend)) {
    colnames <- levels(df$df[[field.legend]])
    lcolors <- grey.colors(length(colnames), start=.3, end=0, gamma=0.8)
    names(lcolors) <- colnames
    lshapes <- c(0:(length(lcolors)-1))
    ltypes <- c(1:2, 4:(length(lcolors)-1))
    names(lshapes) <- names(lcolors)
    g <- g + scale_colour_manual(values=lcolors) + scale_shape_manual(values=lshapes) +
      scale_linetype_manual(values=ltypes)
  }
  if(!is.null(addplot)) {
    g <- g + eval(parse(text=addplot))
  }
#  if(!is.null(file))
  ## g <- g + theme_bwll(base_size=font.size, kpos=kpos, kjus=kjus)
  if(!grid)
    g <- g + opts(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank())
  print(g)
  if(!is.null(file))
    ggsave(file, w=inches(w), h=inches(h), paper="Letter", textspecial=TRUE)
  # return(df)
}
