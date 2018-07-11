## Copyright 2017 Jason Vander Woude and Nathan Ryder
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.


#source("sim-functions.R")
library(ggplot2)

## Generate an RIL population based on the recombination probabilities and the
## generation. For each pair of sequential sites, plot the recombination rate and the
plot.gen.trans.from.recomb <- function(n.gen, n.mem, recomb.probs, k) {

    ## The theoretical equations are based on the probability of an odd number
    ## of physical recombinations
    if (n.gen == 2) {
        f1 <- function(x){x <- 2*x; 1/8*x^2 - 1/2*x + 1/2}
        f2 <- function(x){x <- 2*x; 1/4*x^2 - 1/2*x + 1/2}
        f3 <- function(x){x <- 2*x; -1/2*x^2 + x}
        f4 <- function(x){x <- 2*x; 1/8*x^2}
    } else if (n.gen == 3) {
        f1 <- function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 3/8*x^2 - 3/4*x + 3/4}
        f2 <- function(x){x <- 2*x; 1/16*x^4 - 1/4*x^3 + 1/2*x^2 - 1/2*x + 1/4}
        f3 <- function(x){x <- 2*x; -1/8*x^4 + 1/2*x^3 - x^2 + x}
        f4 <- function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 1/8*x^2 + 1/4*x}
    } else if (n.gen == 4) {
        f1 <- function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 5/16*x^3 + 19/32*x^2 - 7/8*x + 7/8}
        f2 <- function(x){x <- 2*x; 1/64*x^6 - 3/32*x^5 + 9/32*x^4 - 1/2*x^3 + 9/16*x^2 - 3/8*x + 1/8}
        f3 <- function(x){x <- 2*x; -1/32*x^6 + 3/16*x^5 - 9/16*x^4 + x^3 - 9/8*x^2 + 3/4*x}
        f4 <- function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 3/16*x^3 - 1/32*x^2 + 1/2*x}
    } else if (n.gen == 5) {
        f1 <- function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 19/64*x^4 - 1/2*x^3 + 3/4*x^2 - 15/16*x + 15/16}
        f2 <- function(x){x <- 2*x; 1/256*x^8 - 1/32*x^7 + 1/8*x^6 - 5/16*x^5 + 17/32*x^4 - 5/8*x^3 + 1/2*x^2 - 1/4*x + 1/16}
        f3 <- function(x){x <- 2*x; -1/128*x^8 + 1/16*x^7 - 1/4*x^6 + 5/8*x^5 - 17/16*x^4 + 5/4*x^3 - x^2 + 1/2*x}
        f4 <- function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 15/64*x^4 - 1/8*x^3 - 1/4*x^2 + 11/16*x}
    } else {
        stop("bad generation")
    }

    types <- c("AA", "HH", "AH", "AB")
    mat <- matrix(c("AA", "AH", "AB",
                    "AH", "HH", "AH",
                    "AB", "AH", "AA"), byrow = TRUE, nrow=3)

    pop <- pop.to.geno.matrix(create.ril.pop(n.gen, n.mem, recomb.probs))[-(1:2), ]
    temp <- mapply(function(f, s){mat[f+1,s+1]}, pop[,1], pop[,ncol(pop)])
    print(sapply(types, function(t) {perc(temp==t)}))

    inter.df <- function(i, k) {
        first = pop[, i]
        second = pop[, i+1+k]
        each.type = mapply(function(f, s){mat[f+1,s+1]}, first, second)

        percs <- sapply(types, function(t) {perc(each.type == t)})
        if (k >= 0) {
            r1 <- recomb.probs[i]
            recomb <- r1
        }
        if (k >= 1) {
            r2 <- recomb.probs[i+1]
            recomb <- -4*r1*r2 + r1 + r2
        }
        if (k >= 2) {
            r3 <- recomb.probs[i+2]
            recomb <- 16*r1*r2*r3 - 4*r1*r2 - 4*r1*r3 - 4*r2*r3 + r1 + r2 + r3
        }
        if (k >= 3) {
            r4 <- recomb.probs[i+3]
            recomb <- -64*r1*r2*r3*r4 + 16*r1*r2*r3 + 16*r1*r2*r4 + 16*r1*r3*r4 + 16*r2*r3*r4 - 4*r1*r2 - 4*r1*r3 - 4*r2*r3 - 4*r1*r4 - 4*r2*r4 - 4*r3*r4 + r1 + r2 + r3 + r4
        }

        data.frame(recomb = recomb,
                   percs = percs,
                   type = types)
    }


    df <- NULL
    if (k >= 0) {
        df <- rbind(df, do.call(rbind, lapply(1:(length(recomb.probs)), inter.df, 0)))
    }
    if (k >= 1) {
        df <- rbind(df, do.call(rbind, lapply(1:(length(recomb.probs)-1), inter.df, 1)))
    }
    if (k >= 2) {
        df <- rbind(df, do.call(rbind, lapply(1:(length(recomb.probs)-2), inter.df, 2)))
    }
    if (k >= 3) {
        df <- rbind(df, do.call(rbind, lapply(1:(length(recomb.probs)-3), inter.df, 3)))
    }


    xlims <- c(0, 0.5)
    ggplot(df, aes(x=recomb, y=percs, color=type)) + geom_point() +
        stat_function(aes(color="AA"), fun=f1, xlim=xlims) +
        stat_function(aes(color="HH"), fun=f2, xlim=xlims) +
        stat_function(aes(color="AH"), fun=f3, xlim=xlims) +
        stat_function(aes(color="AB"), fun=f4, xlim=xlims) +
        coord_cartesian(xlim = xlims) +
        scale_y_continuous(breaks=seq(0, 1, 0.1))


}

## 3-site: (((2*r1)*(1-2*r2) + (1-2*r1)*(2*r2))/2).expand()
## 4-site: (((2*r1)*(2*r2)*(2*r3) + (2*r1)*(1-2*r2)*(1-2*r3) + (1-2*r1)*(2*r2)*(1-2*r3) + (1-2*r1)*(1-2*r2)*(2*r3))/2).expand()
## 5-site: (((2*r1)*(2*r2)*(2*r3)*(1-2*r4) + (2*r1)*(1-2*r2)*(1-2*r3)*(1-2*r4) + (1-2*r1)*(2*r2)*(1-2*r3)*(1-2*r4) + (1-2*r1)*(1-2*r2)*(2*r3)*(1-2*r4) + (1-2*r1)*(2*r2)*(2*r3)*(2*r4) + (1-2*r1)*(1-2*r2)*(1-2*r3)*(2*r4) + (2*r1)*(2*r2)*(1-2*r3)*(2*r4) + (2*r1)*(1-2*r2)*(2*r3)*(2*r4))/2).expand()


## rec <- rep(0.4, 4) + runif(4, -0.05, 0.05); print(rec); k <- 3; plot.gen.trans.from.recomb(2, 100, rec, k)


props.from.vcf <- function(vcf, snp1, snp2, n.gen) {
    GT <- getGT(vcf)
    a <- sapply(GT[snp1,], function(gt) {sum(gt.to.num(gt))})
    b <- sapply(GT[snp2,], function(gt) {sum(gt.to.num(gt))})

    valid <- !is.na(a) & !is.na(b)

    a <- a[valid]
    b <- b[valid]
    print(length(a))

    types <- c("AA", "HH", "AH", "AB")
    mat <- matrix(c("AA", "AH", "AB",
                    "AH", "HH", "AH",
                    "AB", "AH", "AA"), byrow = TRUE, nrow=3)

    each.type = mapply(function(f, s){mat[f+1,s+1]}, a, b)
    props <- sapply(types, function(t) {perc(each.type == t)})

    if (n.gen == 2) {
        f1 <- function(x){x <- 2*x; 1/8*x^2 - 1/2*x + 1/2}
        f2 <- function(x){x <- 2*x; 1/4*x^2 - 1/2*x + 1/2}
        f3 <- function(x){x <- 2*x; -1/2*x^2 + x}
        f4 <- function(x){x <- 2*x; 1/8*x^2}
    } else if (n.gen == 3) {
        f1 <- function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 3/8*x^2 - 3/4*x + 3/4}
        f2 <- function(x){x <- 2*x; 1/16*x^4 - 1/4*x^3 + 1/2*x^2 - 1/2*x + 1/4}
        f3 <- function(x){x <- 2*x; -1/8*x^4 + 1/2*x^3 - x^2 + x}
        f4 <- function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 1/8*x^2 + 1/4*x}
    } else if (n.gen == 4) {
        f1 <- function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 5/16*x^3 + 19/32*x^2 - 7/8*x + 7/8}
        f2 <- function(x){x <- 2*x; 1/64*x^6 - 3/32*x^5 + 9/32*x^4 - 1/2*x^3 + 9/16*x^2 - 3/8*x + 1/8}
        f3 <- function(x){x <- 2*x; -1/32*x^6 + 3/16*x^5 - 9/16*x^4 + x^3 - 9/8*x^2 + 3/4*x}
        f4 <- function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 3/16*x^3 - 1/32*x^2 + 1/2*x}
    } else if (n.gen == 5) {
        f1 <- function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 19/64*x^4 - 1/2*x^3 + 3/4*x^2 - 15/16*x + 15/16}
        f2 <- function(x){x <- 2*x; 1/256*x^8 - 1/32*x^7 + 1/8*x^6 - 5/16*x^5 + 17/32*x^4 - 5/8*x^3 + 1/2*x^2 - 1/4*x + 1/16}
        f3 <- function(x){x <- 2*x; -1/128*x^8 + 1/16*x^7 - 1/4*x^6 + 5/8*x^5 - 17/16*x^4 + 5/4*x^3 - x^2 + 1/2*x}
        f4 <- function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 15/64*x^4 - 1/8*x^3 - 1/4*x^2 + 11/16*x}
    } else {
        stop("bad generation")
    }

    df <- do.call(rbind, lapply(seq(from=0, to=0.5, length.out=30), function(recomb) {
        data.frame(recomb = recomb,
                   props = props,
                   type = types)
    }))

    ggplot(df, aes(x=recomb, y=props, color=type)) + geom_point() +
        stat_function(aes(color="AA"), fun=f1) +
        stat_function(aes(color="HH"), fun=f2) +
        stat_function(aes(color="AH"), fun=f3) +
        stat_function(aes(color="AB"), fun=f4)
}



## plot.real.est <- function(n.gen, n.mem, recomb.probs) {

##     ## The theoretical equations are based on the probability of an odd number
##     ## of physical recombinations
##     if (n.gen == 2) {
##         f1 <- function(x){x <- 2*x; 1/8*x^2 - 1/2*x + 1/2}
##         f2 <- function(x){x <- 2*x; 1/4*x^2 - 1/2*x + 1/2}
##         f3 <- function(x){x <- 2*x; -1/2*x^2 + x}
##         f4 <- function(x){x <- 2*x; 1/8*x^2}
##     } else if (n.gen == 3) {
##         f1 <- function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 3/8*x^2 - 3/4*x + 3/4}
##         f2 <- function(x){x <- 2*x; 1/16*x^4 - 1/4*x^3 + 1/2*x^2 - 1/2*x + 1/4}
##         f3 <- function(x){x <- 2*x; -1/8*x^4 + 1/2*x^3 - x^2 + x}
##         f4 <- function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 1/8*x^2 + 1/4*x}
##     } else if (n.gen == 4) {
##         f1 <- function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 5/16*x^3 + 19/32*x^2 - 7/8*x + 7/8}
##         f2 <- function(x){x <- 2*x; 1/64*x^6 - 3/32*x^5 + 9/32*x^4 - 1/2*x^3 + 9/16*x^2 - 3/8*x + 1/8}
##         f3 <- function(x){x <- 2*x; -1/32*x^6 + 3/16*x^5 - 9/16*x^4 + x^3 - 9/8*x^2 + 3/4*x}
##         f4 <- function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 3/16*x^3 - 1/32*x^2 + 1/2*x}
##     } else if (n.gen == 5) {
##         f1 <- function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 19/64*x^4 - 1/2*x^3 + 3/4*x^2 - 15/16*x + 15/16}
##         f2 <- function(x){x <- 2*x; 1/256*x^8 - 1/32*x^7 + 1/8*x^6 - 5/16*x^5 + 17/32*x^4 - 5/8*x^3 + 1/2*x^2 - 1/4*x + 1/16}
##         f3 <- function(x){x <- 2*x; -1/128*x^8 + 1/16*x^7 - 1/4*x^6 + 5/8*x^5 - 17/16*x^4 + 5/4*x^3 - x^2 + 1/2*x}
##         f4 <- function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 15/64*x^4 - 1/8*x^3 - 1/4*x^2 + 11/16*x}
##     } else {
##         stop("bad generation")
##     }

##     types <- c("AA", "HH", "AH", "AB")
##     mat <- matrix(c("AA", "AH", "AB",
##                     "AH", "HH", "AH",
##                     "AB", "AH", "AA"), byrow = TRUE, nrow=3)

##     pop <- pop.to.geno.matrix(create.ril.pop(n.gen, n.mem, recomb.probs))[-(1:2), ]
##     temp <- mapply(function(f, s){mat[f+1,s+1]}, pop[,1], pop[,ncol(pop)])
##     print(sapply(types, function(t) {perc(temp==t)}))

##     inter.df <- function(i, k) {
##         first = pop[, i]
##         second = pop[, i+1+k]
##         each.type = mapply(function(f, s){mat[f+1,s+1]}, first, second)

##         counts <- sapply(types, function(t) {sum(each.type == t)})
##         percs <- sapply(types, function(t) {perc(each.type == t)})

##         r1 <- recomb.probs[i]
##         recomb <- r1

##         est1 <- optim(0.25, function(r){(f1(r) - percs[1])^2 + (f2(r) - percs[2])^2 + (f3(r) - percs[3])^2 + (f4(r) - percs[4])^2}, lower=0, upper=0.5)$par

##         est2 <- optim(0.25, function(r){multinomial.pmf(counts, c(f1(r), f2(r), f3(r), f4(r)))}, lower=0, upper=0.5)$par

##         data.frame(recomb = c(rep(recomb, 4), rep(est1, 4), rep(est2, 4)),
##                    percs = rep(percs, 3),
##                    type = rep(types, 3),
##                    opt = as.factor(c(rep(19, 4), rep(25, 4), rep(22, 4))))
##     }

##     df <- do.call(rbind, lapply(1:(length(recomb.probs)), inter.df, 0))

##     xlims <- c(0, 0.5)
##     ggplot(df, aes(x=recomb, y=percs, color=type, shape=opt)) + geom_point(size=3) +
##         stat_function(aes(color="AA"), fun=f1, xlim=xlims) +
##         stat_function(aes(color="HH"), fun=f2, xlim=xlims) +
##         stat_function(aes(color="AH"), fun=f3, xlim=xlims) +
##         stat_function(aes(color="AB"), fun=f4, xlim=xlims) +
##         coord_cartesian(xlim = xlims) +
##         scale_y_continuous(breaks=seq(0, 1, 0.1))


## }


plot.real.est <- function(n.gen, n.mem, recomb.probs, n.recombs) {

    pop <- pop.to.geno.matrix(create.ril.pop(n.gen, n.mem, recomb.probs))[-(1:2), ]

    inter.df <- function(i) {

        important.markers <- pop[,i:(i+n.recombs), drop=FALSE]  # matrix; cols are snps/markers

        fun <- objective.fun(important.markers, n.gen)
        ## vals <- sapply(seq(0,0.5,0.01), fun)
        ## print(paste0("min: ", min(vals)))
        ## print(paste0("max: ", max(vals)))

        result <- optim(rep(0, n.recombs),
                        fun,
                        method="L-BFGS-B",
                        lower=0,
                        upper=0.5)

        ## browser()
        ## print(recomb.probs[i])
        ## print(result$par)
        ## x <- seq(0,0.5,0.01)
        ## qplot(x, sapply(x, fun))


        ## recombs <- recomb.probs[1:n.recombs]

        ## data.frame(recomb = c(rep(recomb, 4), rep(est1, 4), rep(est2, 4)),
        ##            percs = rep(percs, 3),
        ##            type = rep(types, 3),
        ##            opt = as.factor(c(rep(19, 4), rep(25, 4), rep(22, 4))))
        c(rep(NA, i-1), result$par, rep(NA, (length(recomb.probs) - n.recombs - i + 1)))
    }

    df <- do.call(rbind, lapply(1:(length(recomb.probs) - n.recombs + 1), inter.df))

    ## xlims <- c(0, 0.5)
    ## ggplot(df, aes(x=recomb, y=percs, color=type, shape=opt)) + geom_point(size=3) +
    ##     stat_function(aes(color="AA"), fun=f1, xlim=xlims) +
    ##     stat_function(aes(color="HH"), fun=f2, xlim=xlims) +
    ##     stat_function(aes(color="AH"), fun=f3, xlim=xlims) +
    ##     stat_function(aes(color="AB"), fun=f4, xlim=xlims) +
    ##     coord_cartesian(xlim = xlims) +
    ##     scale_y_continuous(breaks=seq(0, 1, 0.1))

    df
}


log.multinomial.pmf <- function(xs, ps) {
    ps[ps==0] <- 1e-10
    sum(xs * log(ps))
}


## m is a matrix with columns as SNPs. NAs are allowed
objective.fun <- function(m, gen) {
    props <- prop.funs(gen)

    function(indiv.recombs) {
        results <- sapply(1:ncol(m), function(col) {
            sapply(1:ncol(m), function(index) {
                if(col <= index) {
                    0
                } else {
                    r <- calc.recomb(col, index, indiv.recombs)
                    marker.pair.log.liklihood(m[,col], m[,index], r, props)
                }
            })
        })
        -1 * sum(results)
    }
}


prop.funs <- function(n.gen) {
    if (n.gen == 2) {
        c(function(x){x <- 2*x; 1/8*x^2 - 1/2*x + 1/2},
          function(x){x <- 2*x; 1/4*x^2 - 1/2*x + 1/2},
          function(x){x <- 2*x; -1/2*x^2 + x},
          function(x){x <- 2*x; 1/8*x^2})
    } else if (n.gen == 3) {
        c(function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 3/8*x^2 - 3/4*x + 3/4},
          function(x){x <- 2*x; 1/16*x^4 - 1/4*x^3 + 1/2*x^2 - 1/2*x + 1/4},
          function(x){x <- 2*x; -1/8*x^4 + 1/2*x^3 - x^2 + x},
          function(x){x <- 2*x; 1/32*x^4 - 1/8*x^3 + 1/8*x^2 + 1/4*x})
    } else if (n.gen == 4) {
        c(function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 5/16*x^3 + 19/32*x^2 - 7/8*x + 7/8},
          function(x){x <- 2*x; 1/64*x^6 - 3/32*x^5 + 9/32*x^4 - 1/2*x^3 + 9/16*x^2 - 3/8*x + 1/8},
          function(x){x <- 2*x; -1/32*x^6 + 3/16*x^5 - 9/16*x^4 + x^3 - 9/8*x^2 + 3/4*x},
          function(x){x <- 2*x; 1/128*x^6 - 3/64*x^5 + 9/64*x^4 - 3/16*x^3 - 1/32*x^2 + 1/2*x})
    } else if (n.gen == 5) {
        c(function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 19/64*x^4 - 1/2*x^3 + 3/4*x^2 - 15/16*x + 15/16},
          function(x){x <- 2*x; 1/256*x^8 - 1/32*x^7 + 1/8*x^6 - 5/16*x^5 + 17/32*x^4 - 5/8*x^3 + 1/2*x^2 - 1/4*x + 1/16},
          function(x){x <- 2*x; -1/128*x^8 + 1/16*x^7 - 1/4*x^6 + 5/8*x^5 - 17/16*x^4 + 5/4*x^3 - x^2 + 1/2*x},
          function(x){x <- 2*x; 1/512*x^8 - 1/64*x^7 + 1/16*x^6 - 5/32*x^5 + 15/64*x^4 - 1/8*x^3 - 1/4*x^2 + 11/16*x})
    } else {
        stop("bad generation")
    }
}


## compute the liklihood of the given recombination rate between two markers
marker.pair.log.liklihood <- function(a, b, recomb, prop.funs) {
    types <- c("AA", "HH", "AH", "AB")
    mat <- matrix(c("AA", "AH", "AB",
                    "AH", "HH", "AH",
                    "AB", "AH", "AA"), byrow = TRUE, nrow=3)

    valid <- !is.na(a) & !is.na(b)

    each.type = mapply(function(x, y){mat[x+1,y+1]}, a[valid], b[valid])
    counts <- sapply(types, function(t) {sum(each.type == t)})
    theory.props <- sapply(prop.funs, function(f) {f(recomb)})

    log.multinomial.pmf(counts, theory.props)
}


calc.recomb <- function(i, j, recomb.probs) {
    a <- min(i,j)
    b <- max(i,j)
    k <- b - a  # number of gaps; n.sites - 1
    if (k >= 1) {
        r1 <- recomb.probs[a]
        recomb <- r1
        if (k >= 2) {
            r2 <- recomb.probs[a+1]
            recomb <- -4*r1*r2 + r1 + r2
            if (k >= 3) {
                r3 <- recomb.probs[a+2]
                recomb <- 16*r1*r2*r3 - 4*r1*r2 - 4*r1*r3 - 4*r2*r3 + r1 + r2 + r3
                if (k >= 4) {
                    r4 <- recomb.probs[a+3]
                    recomb <- -64*r1*r2*r3*r4 + 16*r1*r2*r3 + 16*r1*r2*r4 + 16*r1*r3*r4 + 16*r2*r3*r4 - 4*r1*r2 - 4*r1*r3 - 4*r2*r3 - 4*r1*r4 - 4*r2*r4 - 4*r3*r4 + r1 + r2 + r3 + r4
                    if (k >=5) {
                        stop("invalid number of recombination probabilities")
                    }
                }
            }
        }
    } else {
        stop("must contain at least one recomb.prob")
    }

    recomb
}
