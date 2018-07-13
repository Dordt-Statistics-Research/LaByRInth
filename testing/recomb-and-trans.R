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


est.recombs <- function(pop.ad, n.gen, n.recombs, rerr) {

    init <- rep(0.1, ncol(pop.ad) - 1)

    inter.df <- function(i) {
        important.markers <- pop.ad[,i:(i+n.recombs), , drop=FALSE]  # matrix; cols are snps/markers

        fun <- objective.fun(important.markers, n.gen, rerr)

        result <- optim(init[i:(i+n.recombs-1)],
                        fun,
                        method="L-BFGS-B",
                        lower=0,
                        upper=0.5)

        c(rep(NA, i-1), result$par, rep(NA, (ncol(pop.ad) - n.recombs - i)))
    }

    df <- do.call(rbind, lapply(1:(ncol(pop.ad) - n.recombs), inter.df))
    return(apply(df, 2, mean, na.rm=T))
}


plot.pop <- function(pop, gen, recombs, n.recombs, rerr) {
    inter.df <- function(i) {
        important.markers <- pop[,i:(i+n.recombs), , drop=FALSE]  # matrix; cols are snps/markers
        counts <- get.weighted.counts(pop[,i,], pop[,i+1,], gen, rerr)
        if (counts["AA"] < counts["AB"]) {
            message(colnames(pop)[i:(i+1)])
            message(i-1, ",", i)
            print(counts)
            message("")
            ## browser()
        }
        s <- sum(counts)

        data.frame(
            type = c("AA", "HH", "AH", "AB"),
            recomb = recombs[i],
            props = counts / s
        )
    }

    df <- do.call(rbind, lapply(1:(ncol(pop) - n.recombs), inter.df))

    funs <- prop.funs(gen)

    xlims <- c(0, 0.5)
    ggplot(df, aes(x=recomb, y=props, color=type)) + geom_point() +
        stat_function(aes(color="AA"), fun=funs[[1]], xlim=xlims) +
        stat_function(aes(color="HH"), fun=funs[[2]], xlim=xlims) +
        stat_function(aes(color="AH"), fun=funs[[3]], xlim=xlims) +
        stat_function(aes(color="AB"), fun=funs[[4]], xlim=xlims) +
        coord_cartesian(xlim = xlims) +
        scale_y_continuous(breaks=seq(0, 1, 0.1))
}


log.multinomial.pmf <- function(xs, ps) {
    ps[ps==0] <- 1e-10
    sum(xs * log(ps))
}


## m is a matrix with columns as SNPs. NAs are allowed
objective.fun <- function(m, gen, rerr) {
    props <- prop.funs(gen)

    function(indiv.recombs) {
        results <- sapply(1:ncol(m), function(col) {
            sapply(1:ncol(m), function(index) {
                if(col <= index) {
                    0
                } else {
                    r <- calc.recomb(col, index, indiv.recombs)
                    marker.pair.log.liklihood(m[,col,], m[,index,], r, props, gen, rerr)
                }
            })
        })
        -1 * sum(results)
    }
}


prop.funs <- function(n.gen) {
    ## if (n.gen == 5) {
    ##     c(
    ##         function(r) {1/2*r^8 - 2*r^7 + 4*r^6 - 5*r^5 + 19/4*r^4 - 4*r^3 + 3*r^2 - 15/8*r + 15/16},
    ##         function(r) {r^8 - 4*r^7 + 8*r^6 - 10*r^5 + 17/2*r^4 - 5*r^3 + 2*r^2 - 1/2*r + 1/16},
    ##         function(r) {-2*r^8 + 8*r^7 - 16*r^6 + 20*r^5 - 17*r^4 + 10*r^3 - 4*r^2 + r},
    ##         function(r) {1/2*r^8 - 2*r^7 + 4*r^6 - 5*r^5 + 15/4*r^4 - r^3 - r^2 + 11/8*r}
    ##     )
    ## } else {
    ##     stop("temporary bad generation")
    ## }

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
marker.pair.log.liklihood <- function(a, b, recomb, prop.funs, gen, rerr) {
    counts <- get.weighted.counts(a, b, gen, rerr)
    theory.props <- sapply(prop.funs, function(f) {f(recomb)})

    log.multinomial.pmf(counts, theory.props)
}


get.counts <- function(a, b) {
    types <- c("AA", "HH", "AH", "AB")
    mat <- matrix(c("AA", "AH", "AB",
                    "AH", "HH", "AH",
                    "AB", "AH", "AA"), byrow = TRUE, nrow=3)

    valid <- !is.na(a) & !is.na(b)

    each.type = mapply(function(x, y){mat[x+1,y+1]}, a[valid], b[valid])
    counts <- sapply(types, function(t) {sum(each.type == t)})
}


## a and b should be matrices with an arbitrary but equal number of rows and two
## columns which represen the number of parent 1 allele reads and the number of
## parent 2 allele reads at each site
get.weighted.counts <- function(a, b, gen, rerr) {

    ## probability of each of 3 genetic states given the reads
    perc.het <- 0.5 ^ (gen - 1)
    perc.hom <- 1 - perc.het
    state.probs.3 <- function(ad) {
        p1.reads <- ad[1]
        p2.reads <- ad[2]

        p.p1  <- perc.hom/2 * rerr^p2.reads * (1-rerr)^p1.reads
        p.p2  <- perc.hom/2 * rerr^p1.reads * (1-rerr)^p2.reads
        p.het <- perc.het   * 0.5^(p1.reads + p2.reads)

        probs <- c(p.p1, p.het, p.p2)
        probs / sum(probs)
    }

    ## types <- c("AA", "HH", "AH", "AB")
    ## mat <- matrix(c("AA", "AH", "AB",
    ##                 "AH", "HH", "AH",
    ##                 "AB", "AH", "AA"), byrow = TRUE, nrow=3)

    valid <- mapply(function(quant.a, quant.b) {
        quant.a != 0 && quant.b != 0
    }, apply(a, 1, sum), apply(b, 1, sum))




    ## valid <- rep(T, nrow(a))



    mat.probs <- lapply(which(valid), function(i) {
        col.vec <- matrix(state.probs.3(a[i, ]), ncol=1)
        row.vec <- matrix(state.probs.3(b[i, ]), nrow=1)
        col.vec %*% row.vec
    })

    bound <- do.call(abind3, mat.probs)
    probs <- apply(bound, 1:2, sum)
    counts <- c(
        "AA" = probs[1,1] + probs[3,3],
        "HH" = probs[2,2],
        "AH" = probs[1,2] + probs[2,1] + probs[2,3] + probs[3,2],
        "AB" = probs[1,3] + probs[3,1]
    )
    counts
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
            recomb <- -2*r1*r2 + r1 + r2
            if (k >= 3) {
                r3 <- recomb.probs[a+2]
                recomb <- 4*r1*r2*r3 - 2*r1*r2 - 2*r1*r3 - 2*r2*r3 + r1 + r2 + r3
                if (k >= 4) {
                    r4 <- recomb.probs[a+3]
                    recomb <- -8*r1*r2*r3*r4 + 4*r1*r2*r3 + 4*r1*r2*r4 + 4*r1*r3*r4 + 4*r2*r3*r4 - 2*r1*r2 - 2*r1*r3 - 2*r2*r3 - 2*r1*r4 - 2*r2*r4 - 2*r3*r4 + r1 + r2 + r3 + r4
                    if (k >= 5) {
                        r5 <- recomb.probs[a+4]
                        recomb <- 16*r1*r2*r3*r4*r5 - 8*r1*r2*r3*r4 - 8*r1*r2*r3*r5 - 8*r1*r2*r4*r5 - 8*r1*r3*r4*r5 - 8*r2*r3*r4*r5 + 4*r1*r2*r3 + 4*r1*r2*r4 + 4*r1*r3*r4 + 4*r2*r3*r4 + 4*r1*r2*r5 + 4*r1*r3*r5 + 4*r2*r3*r5 + 4*r1*r4*r5 + 4*r2*r4*r5 + 4*r3*r4*r5 - 2*r1*r2 - 2*r1*r3 - 2*r2*r3 - 2*r1*r4 - 2*r2*r4 - 2*r3*r4 - 2*r1*r5 - 2*r2*r5 - 2*r3*r5 - 2*r4*r5 + r1 + r2 + r3 + r4 + r5
                        if (k >= 6) {
                            r6 <- recomb.probs[a+5]
                            recomb <- -32*r1*r2*r3*r4*r5*r6 + 16*r1*r2*r3*r4*r5 + 16*r1*r2*r3*r4*r6 + 16*r1*r2*r3*r5*r6 + 16*r1*r2*r4*r5*r6 + 16*r1*r3*r4*r5*r6 + 16*r2*r3*r4*r5*r6 - 8*r1*r2*r3*r4 - 8*r1*r2*r3*r5 - 8*r1*r2*r4*r5 - 8*r1*r3*r4*r5 - 8*r2*r3*r4*r5 - 8*r1*r2*r3*r6 - 8*r1*r2*r4*r6 - 8*r1*r3*r4*r6 - 8*r2*r3*r4*r6 - 8*r1*r2*r5*r6 - 8*r1*r3*r5*r6 - 8*r2*r3*r5*r6 - 8*r1*r4*r5*r6 - 8*r2*r4*r5*r6 - 8*r3*r4*r5*r6 + 4*r1*r2*r3 + 4*r1*r2*r4 + 4*r1*r3*r4 + 4*r2*r3*r4 + 4*r1*r2*r5 + 4*r1*r3*r5 + 4*r2*r3*r5 + 4*r1*r4*r5 + 4*r2*r4*r5 + 4*r3*r4*r5 + 4*r1*r2*r6 + 4*r1*r3*r6 + 4*r2*r3*r6 + 4*r1*r4*r6 + 4*r2*r4*r6 + 4*r3*r4*r6 + 4*r1*r5*r6 + 4*r2*r5*r6 + 4*r3*r5*r6 + 4*r4*r5*r6 - 2*r1*r2 - 2*r1*r3 - 2*r2*r3 - 2*r1*r4 - 2*r2*r4 - 2*r3*r4 - 2*r1*r5 - 2*r2*r5 - 2*r3*r5 - 2*r4*r5 - 2*r1*r6 - 2*r2*r6 - 2*r3*r6 - 2*r4*r6 - 2*r5*r6 + r1 + r2 + r3 + r4 + r5 + r6
                            if (k >= 7) {
                                r7 <- recomb.probs[a+6]
                                recomb <- 64*r1*r2*r3*r4*r5*r6*r7 - 32*r1*r2*r3*r4*r5*r6 - 32*r1*r2*r3*r4*r5*r7 - 32*r1*r2*r3*r4*r6*r7 - 32*r1*r2*r3*r5*r6*r7 - 32*r1*r2*r4*r5*r6*r7 - 32*r1*r3*r4*r5*r6*r7 - 32*r2*r3*r4*r5*r6*r7 + 16*r1*r2*r3*r4*r5 + 16*r1*r2*r3*r4*r6 + 16*r1*r2*r3*r5*r6 + 16*r1*r2*r4*r5*r6 + 16*r1*r3*r4*r5*r6 + 16*r2*r3*r4*r5*r6 + 16*r1*r2*r3*r4*r7 + 16*r1*r2*r3*r5*r7 + 16*r1*r2*r4*r5*r7 + 16*r1*r3*r4*r5*r7 + 16*r2*r3*r4*r5*r7 + 16*r1*r2*r3*r6*r7 + 16*r1*r2*r4*r6*r7 + 16*r1*r3*r4*r6*r7 + 16*r2*r3*r4*r6*r7 + 16*r1*r2*r5*r6*r7 + 16*r1*r3*r5*r6*r7 + 16*r2*r3*r5*r6*r7 + 16*r1*r4*r5*r6*r7 + 16*r2*r4*r5*r6*r7 + 16*r3*r4*r5*r6*r7 - 8*r1*r2*r3*r4 - 8*r1*r2*r3*r5 - 8*r1*r2*r4*r5 - 8*r1*r3*r4*r5 - 8*r2*r3*r4*r5 - 8*r1*r2*r3*r6 - 8*r1*r2*r4*r6 - 8*r1*r3*r4*r6 - 8*r2*r3*r4*r6 - 8*r1*r2*r5*r6 - 8*r1*r3*r5*r6 - 8*r2*r3*r5*r6 - 8*r1*r4*r5*r6 - 8*r2*r4*r5*r6 - 8*r3*r4*r5*r6 - 8*r1*r2*r3*r7 - 8*r1*r2*r4*r7 - 8*r1*r3*r4*r7 - 8*r2*r3*r4*r7 - 8*r1*r2*r5*r7 - 8*r1*r3*r5*r7 - 8*r2*r3*r5*r7 - 8*r1*r4*r5*r7 - 8*r2*r4*r5*r7 - 8*r3*r4*r5*r7 - 8*r1*r2*r6*r7 - 8*r1*r3*r6*r7 - 8*r2*r3*r6*r7 - 8*r1*r4*r6*r7 - 8*r2*r4*r6*r7 - 8*r3*r4*r6*r7 - 8*r1*r5*r6*r7 - 8*r2*r5*r6*r7 - 8*r3*r5*r6*r7 - 8*r4*r5*r6*r7 + 4*r1*r2*r3 + 4*r1*r2*r4 + 4*r1*r3*r4 + 4*r2*r3*r4 + 4*r1*r2*r5 + 4*r1*r3*r5 + 4*r2*r3*r5 + 4*r1*r4*r5 + 4*r2*r4*r5 + 4*r3*r4*r5 + 4*r1*r2*r6 + 4*r1*r3*r6 + 4*r2*r3*r6 + 4*r1*r4*r6 + 4*r2*r4*r6 + 4*r3*r4*r6 + 4*r1*r5*r6 + 4*r2*r5*r6 + 4*r3*r5*r6 + 4*r4*r5*r6 + 4*r1*r2*r7 + 4*r1*r3*r7 + 4*r2*r3*r7 + 4*r1*r4*r7 + 4*r2*r4*r7 + 4*r3*r4*r7 + 4*r1*r5*r7 + 4*r2*r5*r7 + 4*r3*r5*r7 + 4*r4*r5*r7 + 4*r1*r6*r7 + 4*r2*r6*r7 + 4*r3*r6*r7 + 4*r4*r6*r7 + 4*r5*r6*r7 - 2*r1*r2 - 2*r1*r3 - 2*r2*r3 - 2*r1*r4 - 2*r2*r4 - 2*r3*r4 - 2*r1*r5 - 2*r2*r5 - 2*r3*r5 - 2*r4*r5 - 2*r1*r6 - 2*r2*r6 - 2*r3*r6 - 2*r4*r6 - 2*r5*r6 - 2*r1*r7 - 2*r2*r7 - 2*r3*r7 - 2*r4*r7 - 2*r5*r7 - 2*r6*r7 + r1 + r2 + r3 + r4 + r5 + r6 + r7
                                if (k >= 8) {
                                    stop("invalid number of recombination probabilities")
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        stop("must contain at least one recomb.prob")
    }

    recomb
}


test.stuff <- function() {

    test.rec.1 <- runif(7, 0, 0.3); test.pop.1 <- pop.to.geno.matrix(create.ril.pop(5, 100, test.rec.1))[-(1:2), ]
    sum((apply(plot.real.est(5, 0, test.rec.1, 4, test.pop.1), 2, mean, na.rm=T) - test.rec.1)^2)

    test.rec.2 <- runif(7, 0, 0.3); test.pop.2 <- pop.to.geno.matrix(create.ril.pop(5, 1000, test.rec.2))[-(1:2), ]

    test.rec.3 <- runif(7, 0, 0.3);
    test.pop.3 <- pop.to.geno.matrix(create.ril.pop(5, 10000, test.rec.3))[-(1:2), ]

    set.seed(0)
    test.rec.4 <- runif(50, 0, 0.3);
    test.pop.4 <- pop.to.geno.matrix(create.ril.pop(5, 10000, test.rec.3))[-(1:2), ]

    vcf <- read.vcfR("~/Desktop/LaByRInth-master/out2.vcf.gz")
    lf.pop <- vcf.to.p1.gt(vcf, c("LAKIN", "FULLER"))
    rerr <- 0
    pop.ad <- get.ad.array.w.p1.ref(vcf, c("LAKIN", "FULLER"))
    recomb.ests <- est.recombs(pop.ad, 5, 1, 0)
    plot.pop(lf.pop, 5, recomb.ests, 1, 0)

    sim.pop <- pop.to.geno.matrix(create.ril.pop(n.gen, n.mem, recomb.probs))[-(1:2), ]

}


sse <- function(x, y) {
    sum((x-y)^2)
}


vcf.to.gt.mat <- function(vcf) {
    mat <- t(getGT(vcf))
    mat <- apply(mat, 2, gt.to.num)
    mat
}


plot.theory <- function(gen) {
    funs <- prop.funs(gen)

    xlims <- c(0, 0.5)
    ggplot(data.frame(x=xlims), aes(x)) +
        stat_function(aes(color="AA"), fun=funs[[1]], xlim=xlims) +
        stat_function(aes(color="HH"), fun=funs[[2]], xlim=xlims) +
        stat_function(aes(color="AH"), fun=funs[[3]], xlim=xlims) +
        stat_function(aes(color="AB"), fun=funs[[4]], xlim=xlims) +
        coord_cartesian(xlim = xlims) +
        scale_y_continuous(breaks=seq(0, 1, 0.1))
}


depths.to.gt <- function(depths) {
    if (depths[1] == 0 && depths[2] == 0) {
        NA
    } else if (depths[1] == 0 && depths[2] != 0) {
        2
    } else if (depths[1] != 0 && depths[2] == 0) {
        0
    } else if (depths[1] != 0 && depths[2] != 0) {
        1
    } else {
        stop("logic error")
    }
}


vcf.to.p1.gt <- function(vcf, parents) {
    ad.arr <- get.ad.w.p1.ref(vcf, parents)
    t(apply(ad.arr, 1:2, depths.to.gt))
}


get.ad.array <- function(vcf) {
    ads <- apply(getAD(vcf), 1:2, ad.to.num)
    abind(ads[1, , ], ads[2, , ], along=3)
}


get.ad.array.w.p1.ref <- function(vcf, parents) {
    str.ad <- getAD(vcf)[ , parents[1]]
    p1.is.ref <- sapply(str.ad, function(str) {
        ## split the string and check if reference read is nonzero
        ad.to.num(str)[1] != 0
    })

    str.ad <- getAD(vcf)[ , parents[2]]
    p2.is.ref <- sapply(str.ad, function(str) {
        ## split the string and check if reference read is nonzero
        ad.to.num(str)[1] != 0
    })

    if (any(p1.is.ref == p2.is.ref))
        stop("problem in finding reference parent")

    ad.arr <- get.ad.array(vcf)

    for (i in seq_along(p2.is.ref)) {
        if (p2.is.ref[i]) {
            ad.arr[i, , ] <- cbind(ad.arr[i, , 2], ad.arr[i, , 1])
        }
    }

    abind3(t(ad.arr[,,1]), t(ad.arr[,,2]))
}


abind3 <- function(...) {
    abind(..., along=3)

}
