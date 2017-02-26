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
library(R6)
library(reshape) # melt
library(ggplot2)

UC <-
    R6Class("UC",
            portable = TRUE,
            cloneable = FALSE,
            private = list(
                ##
                ## Private methods for UC
                ##
                ## UC probability methods
                ## Evaluate p1n(t) of a CTMC Markov chain. ei are the eigenvalues of Q as
                ## returned by eigen.mult. t and j can be vectors.
                p.ctmc = function(t) {
                    single <- length(self$ei$one) # number of single eigenvalues
                    mult <- length(self$ei$confl) # number of non single eigenvalues
                    p <- sapply(1:ncol(self$coef), function(j) {
                        sapply(t, function(tt) sum(self$coef[1:single,j] * exp(tt * self$ei$one)))
                    })
                    if(mult > 0) { # evaluate non single eigenvalues
                        p <- p + sapply(1:ncol(self$coef), function(j) {
                            sapply(t, function(tt) {
                                i <- single+1
                                sum(sapply(1:mult, function(n.ei) {
                                    s <- sum(self$coef[i:(i+self$ei$mult[n.ei]-1),j] *
                                             tt^c(0:(self$ei$mult[n.ei]-1))) *
                                        exp(tt * self$ei$confl[n.ei])
                                    i <<- i + self$ei$mult[n.ei]
                                    return(s)
                                }))
                            })
                        })
                    }
                    return(Re(p))
                },
                ## Evaluate p1n(t) of a DTMC Markov chain. ei are the eigenvalues of P as
                ## returned by eigen.mult.
                p.dtmc = function(t) {
                    single <- length(self$ei$one) # number of single eigenvalues
                    mult <- length(self$ei$confl) # number of non single eigenvalues
                    p <- sapply(1:ncol(self$coef), function(j) {
                        sapply(t, function(t) sum(self$coef[1:single,j] * self$ei$one^t))
                    })
                    if(mult > 0) { # evaluate non single eigenvalues
                        p <- p + sapply(1:ncol(self$coef), function(j) {
                            sapply(t, function(t) {
                                i <- single+1
                                sum(sapply(1:mult, function(n.ei) {
                                    s <- if(abs(self$ei$confl[n.ei]) == 0) {
                                             if(t < self$ei$mult[n.ei]) self$coef[i+t,j]
                                             else 0
                                         } else {
                                             sum(self$coef[i:(i+self$ei$mult[n.ei]-1),j] *
                                                 t^c(0:(self$ei$mult[n.ei]-1))) *
                                                 self$ei$confl[n.ei]^t
                                         }
                                    i <<- i + self$ei$mult[n.ei]
                                    return(s)
                                }))
                            })
                        })
                    }
                    return(Re(p))
                }
            ),
            ##
            public = list(
                method=NULL,            # UC solution method
                states=NULL,            # states of the chain
                mc.type=NULL,           # Markov chain type: 'ctmc', 'dtmc'
                matrix.type=NULL,       # Matrix type: 'Q', 'P'
                states.names=NULL,      # Names of the Markov chain states
                lim=NULL,               # true if the UC are solved using the lim. dist.
                qn=NULL,                # uniformization parameter
                choose.samples=NULL,    # if TRUE, samples are choosen evely spaced in log scale.
                alpha=NULL,
                max.power=NULL,         # maximum power used in the Vand. system when choose.samples=TRUE
                use.eigen.from.file=NULL,
                ei.max=NULL,            # Take the ei.max largest Re eigenvalues (0 for unlimited)
                ei.max.mult=NULL,       # maximum multiplicity of the eigenvalues (0 for unlimited)
                coef=NULL,              # undetermined coefficients (UC)
                matrix.dim=NULL,        # matrix dimensions
                matrix.zeros=NULL,      # Number of matrix elements = 0 in the equations.
                matrix.bytes=NULL,      # Matrix storage in Mb.
                time=NULL,              # computation times
                unif=NULL,              # true if the UC are solved using uniformization
                qii.max=NULL,           # max_i |Qii|
                i=NULL,                 # initial state
                j=NULL,                 # target states
                ei=NULL,                # eigenvalues
                kappa=NULL,             # condition number of the matix used to solve the coefs.
                ##
                initialize = function(type=NULL, method=NULL, states=NULL, j=NULL, i=NULL,
                                      states.names=NULL,
                                      unif=FALSE, qii.max=NULL, qn=NULL,
                                      choose.samples=FALSE, alpha=10.0, lim=FALSE,
                                      ei.max=0, ei.max.mult=0,
                                      use.eigen.from.file=NULL, uc=NULL, ei=NULL) {                    
                    self$method <- method          # UC solution method
                    self$states <- states          # states of the chain
                    self$matrix.type <- type       # Matrix type: 'Q', 'P'
                    self$states.names <- states.names # Names of the Markov chain states
                    self$unif <- unif              # true if the UC are solved using uniformization
                    self$lim <- lim                # true if the UC are solved using the lim. dist.
                    self$i <- i                    # initial state
                    self$j <- j                    # target states
                    self$qii.max <- qii.max        # max_i |Qii|
                    self$qn <- qn                  # uniformization parameter
                    self$choose.samples <- choose.samples # if TRUE, samples are choosen evely spaced in log scale.
                    self$alpha <- alpha
                    self$ei <- ei                  # eigenvalues
                    self$use.eigen.from.file <- use.eigen.from.file
                    self$ei.max <- ei.max          # Take the ei.max largest Re eigenvalues (0 for unlimited)
                    self$ei.max.mult <- ei.max.mult# maximum multiplicity of the eigenvalues (0 for unlimited)
                    self$coef <- coef              # undetermined coefficients (UC)
                    self$time <- time               # computation times
                    if(!is.null(type)) {
                        self$mc.type = switch(type,
                                                 P = 'dtmc',
                                                 Q = 'ctmc',
                                                 stop('Unknown type')
                                                 )
                    } else {
                        self$mc.type = NULL ;
                    }
                    if(!is.null(j) && !is.null(states.names) && (length(states.names) != length(j))) {
                        warning("UC: Incorrect number of states.names")
                        self$states.names <- NULL
                    } else {
                        self$states.names <- as.character(states.names)
                    }
                    ##
                    self$kappa <- kappa            # condition number of the matix used to solve the coefs.
                },
                compute.eigenvalues = function(type, evec=FALSE, ei.max, ei.max.mult) {
                    self$ei <- Eigen$new()
                    self$ei$compute.eigenvalues(type=type, evec=evec, ei.max=ei.max, ei.max.mult=ei.max.mult)
                },
                ##
                ## Evaluate pn(t)
                prob = function(t) {
                    nc <- ncol(self$coef)
                    if(is.null(self$coef)) {
                        prob <- data.frame()
                        warning("empty coefficients, ", self$title())
                        prob <- data.frame()
                    } else {
                        prob <-
                            switch(self$mc.type,
                                   ctmc = data.frame(matrix(private$p.ctmc(t), nr=length(t), nc=nc)),
                                   dtmc = data.frame(matrix(private$p.dtmc(t), nr=length(t), nc=nc)))
                        colnames(prob) <- self$get.states.names()
                    }
                    return(prob)
                },
                ##
                print = function(t) {
                    if(debug) cat('---------------\n')
                    msg("", "UC", "summary")
                    msg("", " UC solution method", self$method)
                    if(!is.null(self$use.eigen.from.file)) {
                        msg("", "Use eigenvalues from file", self$use.eigen.from.file)
                    }
                    if(empty(self$coef)) msg("", " UC coef", "failure")
                    else {
                        ei.max <- max(Re(self$coef))
                        ei.min <- min(Re(self$coef))
                        msg("", " UC coef", "max(Re):", ei.max, ", min(Re):", ei.min)
                    }
                    msg("", " kappa", self$kappa)
                    msg("", " states of the chain", self$states)
                    msg("", " Markov chain type", self$mc.type)
                    msg("", " Matrix type", self$matrix.type)
                    if(!is.null(self$states.names))
                        msg("", " states names", paste(sep='', self$states.names, collapse=','))
                    msg("", " matrix dimensions", paste(sep='', self$matrix.dim, collapse=','))
                    msg("", " matrix zeros", self$matrix.zeros)
                    msg("", " matrix storage", self$matrix.bytes, ' bytes')
                    msg("", " solved using uniformization", self$unif)
                    if(self$unif) {
                        msg("", " uniformization parameter", self$qn)
                        msg("", " qii.max", self$qii.max)
                    }
                    msg("", " choose samples", self$choose.samples)
                    if(!is.null(self$max.power)) msg("", " max power", self$max.power)
                    msg("", " time constant", self$time.constant())
                    msg("", " time constant, alpha", self$alpha)
                    msg("", " solved using the lim. dist.", self$lim)
                    if(length(self$i) > 1) {
                        if(length(self$i) <= 5) info <- paste(sep='', sprintf('%.2e',self$i), collapse=',')
                        else info <- paste(sep='', c(sprintf('%.2e', self$i[1:5]), '...'), collapse=',')
                        msg("", " initial state", info)
                    } else  msg("", " initial state", self$i)
                    msg("", " target states", paste(sep='', self$j, collapse=','))
                    msg("", " running times",)
                    t.names <- c('sys.self', 'elapsed')
                    if(!is.null(self$time)) {
                        for(t in names(self$time)) {
                            cat(sep='', '  ', t, ': ',
                                paste(t.names, format(self$time[[t]][t.names]), collapse=', '), '\n')
                            ##    for(n in c('sys.self', 'elapsed')) cat(n, self$time[[t]][[n]]))
                            ## cat('\n')
                            ## cat(paste(sep='', format(self$time[[t]]), collapse=', '), '\n')
                        }
                    }
                    msg("", " eigenvalues maximum multiplicity", self$ei.max.mult)
                    msg("", " eigenvalues maximum number", self$ei.max)
                    msg("", " eigenvalues computed with matlab", self$use.matlab)
                    self$ei$print()
                },
                ## Return the eigenvalue with the second lagest real part and its coefficient
                ## (id, value and confl indicator, coefficient)
                ##
                second.largest.re = function() {
                    max.re <- self$ei$second.largest.re()
                    if(max.re$confl) {
                        mult <- self$ei$mult[max.re$id]
                        id <- length(self$ei$one) + sum(self$ei$mult[1:max.re$id])
                        coef <- self$coef[(id-mult+1):id]
                        return(c(max.re, mult=mult, coef=coef))
                    }
                    return(c(max.re, coef=self$ei$one[max.re$id]))
                },
                ##
                ## Return the eigenvalue with the second lagest real part and its coefficient
                ## (id, value and confl indicator, coefficient)
                ##
                time.constant = function() {
                    switch(self$ei$type,
                           P = self$alpha/(1-abs(self$ei$second.largest.re()$value)),
                           Q = (self$alpha*self$qn)/abs(self$ei$second.smallest.abs()$value),
                           stop())
                },
                ##
                ##
                compute.samples = function() {
                    msg(class="UC", "compute.samples",'Compute the samples')
                    int.log.scale <- function(maxp, N) {
                        maxp <- max(1, maxp)
                        n <- max(1, ceiling(N/log10(maxp)))
                        if((N > n) && (maxp > n)) {
                            return(c(1:(n-1), round(10^(seq(log10(n), log10(maxp), len=N-n+1)))))
                        } else {
                            return(1:N)
                        }
                    }
                    if(is.null(self$max.power)) 
                        self$max.power <- max(ceiling(self$time.constant()), self$ei$number()-1)
                    if(self$max.power < self$ei$number()) return(NULL)
                    else return(int.log.scale(self$max.power, self$ei$number()-1))
                },
                ## Correct the coefficients of the confluent eigenvalues
                correct.confl.coef = function() {
                    if(empty(self$ei$confl) || empty(self$coef)) return()
                    ##
                    msg(class="UC", "correct.confl.coef",
                        "Correcting the coefficients of the confluent eigenvalues")
                    single <- length(self$ei$one) # number of single eigenvalues
                    j <- single+1
                    for(n.ei in 1:length(self$ei$confl)) {
                        mult <- self$ei$mult[n.ei]
                        if(self$ei$confl[n.ei] == 0) {
                            self$coef[(j+1):(j+mult-1),] <-
                                (self$coef[(j+1):(j+mult-1),] * self$qn^(1:(mult-1))) / factorial(1:(mult-1))
                        } else {
                            for(m in 1:(mult-1)) {
                                if(m < (mult-1)) {
                                    for(k in (m+1):(mult-1)) {
                                        self$coef[j+m,] <- self$coef[j+m,] + self$coef[j+k,] * q.n.coefs(k, m)
                                    }
                                }
                                self$coef[j+m,] <- self$coef[j+m,] * (self$qn * self$ei$confl[n.ei])^m
                            }
                        }
                        j <- j + mult
                    }
                },
                plot = function(t, comp=NULL, x.comp=NULL, logy=TRUE, logx=FALSE,
                            title=NULL, ylim=NULL) {
                    if(is.null(title)) title <- self$title()
                    plot.facet.grid.row.col(data.main=list(list(prob=self$prob(t))), x=t,
                                            data.comp=comp, x.comp=x.comp,
                                            field.x='t', field.y='prob', field.legend='n',
                                            logy=logy, logx=logx, title=title, #type='single',
                                            ylim=ylim)
                },
                title = function(show.states=TRUE) {
                    if(self$unif) type <- 'CTMC, Unif.'
                    else type <- switch(self$mc.type, ctmc = 'CTMC', dtmc = 'DTMC')
                    method <- switch(self$method,
                                     vand = 'VAND',
                                     evec = 'Eigenvectors')
                    title <- paste(type, method)
                    if(self$lim) title <- paste(sep='', title,  ', limit dist')
                    if(show.states) title <- paste(sep='', title, ', N=', self$states)
                    return(title)
                },
                get.states.names = function() {
                    if(empty(self$states.names)) as.character(self$j)
                    else self$states.names
                },
                ##
                set.Vmatrix.characteristics = function(dim, zeros, bytes) {
                    self$matrix.dim <- dim
                    self$matrix.zeros <- zeros
                    self$matrix.bytes <- bytes
                },
                ##
                set.coef = function(coef, cnames=NULL) {
                    self$coef <- coef
                    if(!is.null(cnames)) colnames(self$coef) <- cnames
                },
                ##
                init.time = function(t.list) {
                    self$time <- t.list
                },
                ##
                add.time = function(t.list) {
                    self$time <- c(self$time, t.list)
                },
                reset.coef = function() {
                    self$coef <- NULL
                },
                ##
                set.ei.one = function(one) {
                    self$ei$one <- one
                },
                ##
                set.kappa = function(kappa) {
                    self$kappa <- kappa
                },
                ##
                set.qn = function(qn) {
                    self$qn <- qn
                }
            ))
## end of UC class


##########################################################################################
## UC public methods
##########################################################################################
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
    g <- g + ggtitle(title) +
        theme(plot.title=element_text(size=8,face="bold"))
  msg("", "plot.facet.grid.row.col", 'print ggplot')
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
    g <- g + theme(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank())
  print(g)
  if(!is.null(file))
    ggsave(file, w=inches(w), h=inches(h), paper="Letter", textspecial=TRUE)
  # return(df)
}
