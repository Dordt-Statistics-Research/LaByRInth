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

        ## result <- optim(init[i:(i+n.recombs-1)],
        ##                 fun,
        ##                 method="L-BFGS-B",
        ##                 lower=0,
        ##                 upper=0.5)
        ## result <- result$par

        x <- seq(0, 0.5, 0.01)
        result <- x[which.max(sapply(x, fun))]

        c(rep(NA, i-1), result, rep(NA, (ncol(pop.ad) - n.recombs - i)))
    }

    df <- do.call(rbind, lapply(1:(ncol(pop.ad) - n.recombs), inter.df))
    ## print(df)
    return(apply(df, 2, mean, na.rm=T))
}


plot.pop <- function(pop, gen, recombs, n.recombs, rerr, chroms) {
    if (length(chroms) != ncol(pop))
        stop("bad chrom vector")


    inter.df <- function(i) {
        important.markers <- pop[,i:(i+n.recombs), , drop=FALSE]  # matrix; cols are snps/markers
        ## counts <- get.weighted.counts(pop[,i,], pop[,i+1,], gen, rerr)
        counts <- get.counts(apply(pop[,i,], 1, depths.to.gt),
                             apply(pop[,i+1,], 1, depths.to.gt))
        s <- sum(counts)
        props <- counts / s

        ## if (counts["AA"] < counts["AB"]) {
        if (props["AB"] > 0.8) {
            message(colnames(pop)[i:(i+1)])
            message(i-1, ",", i)
            print(rbind(counts=counts, props=round(props, 3)))
            message("")
            browser()
        }

        data.frame(
            id = i,
            type = c("AA", "HH", "AH", "AB"),
            recomb = recombs[i],
            props = props
        )
    }


    dfs <- lapply(unique(chroms), function(chrom) {
        indices <- which(chroms==chrom)
        begin <- min(indices)
        end <- max(indices) - 1
        do.call(rbind, lapply(begin:end, inter.df))
    })

    df <- do.call(rbind, dfs)

    funs <- prop.funs(gen)

    color1 <- "#37474F"
    color2 <- '#CACACA'
    xlims <- c(0, 0.5)
    ylims <- c(0, 1)
    gentext <- ifelse(gen==2, "Generation", "Generations")
    thickness <- 2
    ggplot(df, aes(x=recomb, y=props, color=type)) + geom_point(size=1.5) +
        stat_function(aes(color="AA"), fun=funs[[1]], xlim=xlims, size=thickness) +
        stat_function(aes(color="HH"), fun=funs[[2]], xlim=xlims, size=thickness) +
        stat_function(aes(color="AH"), fun=funs[[3]], xlim=xlims, size=thickness) +
        stat_function(aes(color="AB"), fun=funs[[4]], xlim=xlims, size=thickness) +
        coord_cartesian(xlim = xlims) +
        coord_cartesian(ylim = ylims) +
        scale_y_continuous(breaks=seq(0, 1, 0.2)) +
        ggtitle(paste0("Transitional Characteristics for\n", gen-1, " ", gentext, " of Inbreeding (F", gen, ")")) +
        xlab("Probability of Recombinations (Odd Number)") +
        ylab("Proportion of Population") +
        scale_color_discrete(name="Type") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(fill = '#37474F', colour = '#37474F'),
              panel.background = element_rect(fill = color2, color = color2),
              text = element_text(size=20),
              legend.key = element_rect(fill = color1, color=color1),
              legend.background = element_rect(fill = color1, color=color1),
              axis.text.x=element_text(color = color2),
              axis.text.y=element_text(color = color2),
              axis.title.x=element_text(color = color2, vjust = -3),
              axis.title.y=element_text(color = color2, vjust = 3),
              legend.title=element_text(color = color2),
              legend.text=element_text(color = color2),
              title=element_text(color = color2),
              plot.margin = unit(c(1,1,1,1), "cm"),
              legend.key.size = unit(3, 'lines'))

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
                    get.obj.fun(m[,col,], m[,index,], rerr)(r)
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
    } else if (n.gen == 6) {
        c(
            function(r){1/2*r^10 - 5/2*r^9 + 25/4*r^8 - 10*r^7 + 45/4*r^6 - 39/4*r^5 + 59/8*r^4 - 21/4*r^3 + 109/32*r^2 - 31/16*r + 31/32},
            function(r){r^10 - 5*r^9 + 25/2*r^8 - 20*r^7 + 45/2*r^6 - 37/2*r^5 + 45/4*r^4 - 5*r^3 + 25/16*r^2 - 5/16*r + 1/32},
            function(r){-2*r^10 + 10*r^9 - 25*r^8 + 40*r^7 - 45*r^6 + 37*r^5 - 45/2*r^4 + 10*r^3 - 25/8*r^2 + 5/8*r},
            function(r){1/2*r^10 - 5/2*r^9 + 25/4*r^8 - 10*r^7 + 45/4*r^6 - 35/4*r^5 + 31/8*r^4 + 1/4*r^3 - 59/32*r^2 + 13/8*r}
        )
    } else if (n.gen == 7) {
        c(
            function(r){1/2*r^12 - 3*r^11 + 9*r^10 - 35/2*r^9 + 195/8*r^8 - 51/2*r^7 + 21*r^6 - 59/4*r^5 + 311/32*r^4 - 99/16*r^3 + 117/32*r^2 - 63/32*r + 63/64},
            function(r){r^12 - 6*r^11 + 18*r^10 - 35*r^9 + 195/4*r^8 - 51*r^7 + 41*r^6 - 51/2*r^5 + 195/16*r^4 - 35/8*r^3 + 9/8*r^2 - 3/16*r + 1/64},
            function(r){-2*r^12 + 12*r^11 - 36*r^10 + 70*r^9 - 195/2*r^8 + 102*r^7 - 82*r^6 + 51*r^5 - 195/8*r^4 + 35/4*r^3 - 9/4*r^2 + 3/8*r},
            function(r){1/2*r^12 - 3*r^11 + 9*r^10 - 35/2*r^9 + 195/8*r^8 - 51/2*r^7 + 20*r^6 - 43/4*r^5 + 79/32*r^4 + 29/16*r^3 - 81/32*r^2 + 57/32*r}
        )
    } else if (n.gen == 8) {
        c(
            function(r){1/2*r^14 - 7/2*r^13 + 49/4*r^12 - 28*r^11 + 371/8*r^10 - 469/8*r^9 + 931/16*r^8 - 93/2*r^7 + 1003/32*r^6 - 617/32*r^5 + 743/64*r^4 - 219/32*r^3 + 487/128*r^2 - 127/64*r + 127/128},
            function(r){r^14 - 7*r^13 + 49/2*r^12 - 56*r^11 + 371/4*r^10 - 469/4*r^9 + 931/8*r^8 - 92*r^7 + 931/16*r^6 - 469/16*r^5 + 371/32*r^4 - 7/2*r^3 + 49/64*r^2 - 7/64*r + 1/128},
            function(r){-2*r^14 + 14*r^13 - 49*r^12 + 112*r^11 - 371/2*r^10 + 469/2*r^9 - 931/4*r^8 + 184*r^7 - 931/8*r^6 + 469/8*r^5 - 371/16*r^4 + 7*r^3 - 49/32*r^2 + 7/32*r},
            function(r){1/2*r^14 - 7/2*r^13 + 49/4*r^12 - 28*r^11 + 371/8*r^10 - 469/8*r^9 + 931/16*r^8 - 91/2*r^7 + 859/32*r^6 - 321/32*r^5 - 1/64*r^4 + 107/32*r^3 - 389/128*r^2 + 15/8*r}
        )
    } else if (n.gen == 9) {
        c(
            function(r){1/2*r^16 - 4*r^15 + 16*r^14 - 42*r^13 + 161/2*r^12 - 119*r^11 + 140*r^10 - 267/2*r^9 + 1675/16*r^8 - 277/4*r^7 + 163/4*r^6 - 23*r^5 + 417/32*r^4 - 233/32*r^3 + 249/64*r^2 - 255/128*r + 255/256},
            function(r){r^16 - 8*r^15 + 32*r^14 - 84*r^13 + 161*r^12 - 238*r^11 + 280*r^10 - 267*r^9 + 1667/8*r^8 - 267/2*r^7 + 70*r^6 - 119/4*r^5 + 161/16*r^4 - 21/8*r^3 + 1/2*r^2 - 1/16*r + 1/256},
            function(r){-2*r^16 + 16*r^15 - 64*r^14 + 168*r^13 - 322*r^12 + 476*r^11 - 560*r^10 + 534*r^9 - 1667/4*r^8 + 267*r^7 - 140*r^6 + 119/2*r^5 - 161/8*r^4 + 21/4*r^3 - r^2 + 1/8*r},
            function(r){1/2*r^16 - 4*r^15 + 16*r^14 - 42*r^13 + 161/2*r^12 - 119*r^11 + 140*r^10 - 267/2*r^9 + 1659/16*r^8 - 257/4*r^7 + 117/4*r^6 - 27/4*r^5 - 95/32*r^4 + 149/32*r^3 - 217/64*r^2 + 247/128*r}
        )
    } else if (n.gen == 10) {
        c(
            function(r){1/2*r^18 - 9/2*r^17 + 81/4*r^16 - 60*r^15 + 261/2*r^14 - 441/2*r^13 + 1197/4*r^12 - 333*r^11 + 4923/16*r^10 - 3811/16*r^9 + 5011/32*r^8 - 361/4*r^7 + 1549/32*r^6 - 827/32*r^5 + 899/64*r^4 - 121/16*r^3 + 2017/512*r^2 - 511/256*r + 511/512},
            function(r){r^18 - 9*r^17 + 81/2*r^16 - 120*r^15 + 261*r^14 - 441*r^13 + 1197/2*r^12 - 666*r^11 + 4923/8*r^10 - 3803/8*r^9 + 4923/16*r^8 - 333/2*r^7 + 1197/16*r^6 - 441/16*r^5 + 261/32*r^4 - 15/8*r^3 + 81/256*r^2 - 9/256*r + 1/512},
            function(r){-2*r^18 + 18*r^17 - 81*r^16 + 240*r^15 - 522*r^14 + 882*r^13 - 1197*r^12 + 1332*r^11 - 4923/4*r^10 + 3803/4*r^9 - 4923/8*r^8 + 333*r^7 - 1197/8*r^6 + 441/8*r^5 - 261/16*r^4 + 15/4*r^3 - 81/128*r^2 + 9/128*r},
            function(r){1/2*r^18 - 9/2*r^17 + 81/4*r^16 - 60*r^15 + 261/2*r^14 - 441/2*r^13 + 1197/4*r^12 - 333*r^11 + 4923/16*r^10 - 3795/16*r^9 + 4835/32*r^8 - 305/4*r^7 + 845/32*r^6 - 55/32*r^5 - 377/64*r^4 + 91/16*r^3 - 1855/512*r^2 + 251/128*r}
        )
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
    test.rec.4 <- runif(200, 0, 0.5);
    test.pop.4 <- create.ril.pop(5, 300, test.rec.4)
    sample.pop.4 <- sample.population(test.pop.4, 0, 4)[-(1:2), , ]

    set.seed(0)
    test.rec.5 <- runif(200, 0, 0.5);
    test.pop.5 <- create.ril.pop(5, 300, test.rec.5)
    sample.pop.5 <- sample.population(test.pop.5, 0, 2)[-(1:2), , ]

    set.seed(0)
    test.rec.6 <- runif(200, 0, 0.5);
    test.pop.6 <- create.ril.pop(5, 300, test.rec.6)
    sample.pop.6 <- sample.population(test.pop.6, 0, 1)[-(1:2), , ]

    vcf <- read.vcfR("~/Desktop/LaByRInth-master/out2.vcf.gz")
    lf.pop <- vcf.to.p1.gt(vcf, c("LAKIN", "FULLER"))
    rerr <- 0
    pop.ad <- get.ad.array.w.p1.ref(vcf, c("LAKIN", "FULLER"))
    recomb.ests <- est.recombs(pop.ad, 5, 1, 0)
    plot.pop(lf.pop, 5, recomb.ests, 1, 0)

    sim.pop <- pop.to.geno.matrix(create.ril.pop(n.gen, n.mem, recomb.probs))[-(1:2), ]

    ## LaByRInth recomb estimates vs physical distance
    vcf <- read.vcfR("~/Desktop/LaByRInth-master/out2.vcf.gz")
    rerr <- 0
    pop.ad <- get.ad.array.w.p1.ref(vcf, c("LAKIN", "FULLER"))
    lf.indices <- which(rownames(lb.pop.ad) %in% c("LAKIN","FULLER"))
    all.chroms <- getCHROM(vcf)
    chroms <- unique(all.chroms)
    recomb.ests <- est.recombs(pop.ad, 5, 1, 0)
    dists <- diff(getPOS(vcf))
    dists[dists < 0] <- 0


    ## LB-testing
    vcf <- read.vcfR("~/Desktop/LaByRInth-master/out2.vcf.gz")
    lb.pop.ad <- get.ad.array.w.p1.ref(vcf, c("LAKIN", "FULLER"))
    lf.indices <- which(rownames(lb.pop.ad) %in% c("LAKIN","FULLER"))
    all.chroms <- getCHROM(vcf)
    chroms <- unique(all.chroms)
    lb.pop.ad <- lb.pop.ad[-lf.indices, all.chroms=="1A", ]
    dist.to.r <- function(d) 0.5*(1 - exp(-d / 1e6))
    lb.recs <- dist.to.r(diff(getPOS(vcf)[all.chroms=="1A"]))
    recs.list <- lapply(chroms, function(chrom) {c(dist.to.r(diff(getPOS(vcf)[chrom==all.chroms])), 0.5)})
    lb.recs <- do.call(c, recs.list)
    lb.recs <- lb.recs[1:(length(lb.recs)-1)]

    lb.recs.bak <- lb.recs
    lb.pop.ad.bak <- lb.pop.ad

    inter.df <- function(i) {
        important.markers <- lb.pop.ad[,i:(i+1), , drop=FALSE]  # matrix; cols are snps/markers
        counts <- get.counts(apply(lb.pop.ad[,i,], 1, depths.to.gt),
                             apply(lb.pop.ad[,i+1,], 1, depths.to.gt))
        s <- sum(counts)

        data.frame(
            type = c("AA", "HH", "AH", "AB"),
            recomb = lb.recs[i],
            props = counts / s
        )
    }
#BEGIN
    df <- do.call(rbind, lapply(1:(ncol(lb.pop.ad) - 1), inter.df))

    funs <- list(function(x) 0.94*(1-2*x),  # AA
                 function(x) 0.06*(1-2*x),  # HH
                 function(x) 1.06*x,        # AH
                 function(x) 0.94*x)        # AB

    color1 <- "#37474F"
    color2 <- '#CACACA'
    xlims <- c(0, 0.5)
    ylims <- c(0, 1)
    thickness <- 2
    ggplot(df, aes(x=recomb, y=props, color=type)) + geom_point(size=1) +
        stat_function(aes(color="AA"), fun=funs[[1]], xlim=xlims, size=thickness) +
        stat_function(aes(color="HH"), fun=funs[[2]], xlim=xlims, size=thickness) +
        stat_function(aes(color="AH"), fun=funs[[3]], xlim=xlims, size=thickness) +
        stat_function(aes(color="AB"), fun=funs[[4]], xlim=xlims, size=thickness) +
        coord_cartesian(xlim = xlims) +
        coord_cartesian(ylim = ylims) +
        scale_y_continuous(breaks=seq(0, 1, 0.2)) +
        ggtitle(paste0("Transitional Characteristics for F?")) +
        xlab("Probability of Recombinations (Odd Number)") +
        ylab("Proportion of Population") +
        scale_color_discrete(name="Type") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(fill = '#37474F', colour = '#37474F'),
              panel.background = element_rect(fill = color2, color = color2),
              text = element_text(size=30),
              legend.key = element_rect(fill = color1, color=color1),
              legend.background = element_rect(fill = color1, color=color1),
              axis.text.x=element_text(color = color2),
              axis.text.y=element_text(color = color2),
              axis.title.x=element_text(color = color2, vjust = -3),
              axis.title.y=element_text(color = color2, vjust = 3),
              legend.title=element_text(color = color2),
              legend.text=element_text(color = color2),
              title=element_text(color = color2),
              plot.margin = unit(c(1,1,1,1), "cm"),
              legend.key.size = unit(3, 'lines'))






    
#END
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
    ylims <- c(0, 1)
    gentext <- ifelse(gen==2, "Generation", "Generations")
    thickness <- 2
    ggplot(data.frame(x=xlims), aes(x)) +
        stat_function(aes(color="AA"), fun=funs[[1]], xlim=xlims, size=thickness) +
        stat_function(aes(color="HH"), fun=funs[[2]], xlim=xlims, size=thickness) +
        stat_function(aes(color="AH"), fun=funs[[3]], xlim=xlims, size=thickness) +
        stat_function(aes(color="AB"), fun=funs[[4]], xlim=xlims, size=thickness) +
        coord_cartesian(xlim = xlims) +
        coord_cartesian(ylim = ylims) +
        scale_y_continuous(breaks=seq(0, 1, 0.2)) +
        ggtitle(paste0("Transitional Characteristics for\n", gen-1, " ", gentext, " of Inbreeding (F", gen, ")")) +
        xlab("Probability of Recombinations (Odd Number)") +
        ylab("Proportion of Population") +
        scale_color_discrete(name="Type") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(fill = '#37474F', colour = '#37474F'),
              panel.background = element_rect(fill = color2, color = color2),
              text = element_text(size=20),
              legend.key = element_rect(fill = color1, color=color1),
              legend.background = element_rect(fill = color1, color=color1),
              axis.text.x=element_text(color = color2),
              axis.text.y=element_text(color = color2),
              axis.title.x=element_text(color = color2, vjust = -3),
              axis.title.y=element_text(color = color2, vjust = 3),
              legend.title=element_text(color = color2),
              legend.text=element_text(color = color2),
              title=element_text(color = color2),
              plot.margin = unit(c(1,1,1,1), "cm"),
              legend.key.size = unit(3, 'lines'))

        ## guides(fill=guide_legend(title="Transition Type"))
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


get.ad.array.w.p1.ref <- function(vcf, parents) {

    get.ad.array <- function(vcf) {
        ads <- apply(getAD(vcf), 1:2, ad.to.num)
        abind(ads[1, , ], ads[2, , ], along=3)
    }

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






## TODO(Jason): currently analyzing each SNP twice
get.obj.fun <- function(snp.a.reads, snp.b.reads, rerr) {
    ## return length 3 vector indicating probability of reads for each of the 3
    ## states (P1, H, P2)
    read.probs.given.states <- function(reads) {
        n1 <- reads[1]
        n2 <- reads[2]
        n  <- n1 + n2

        k <- choose(n, n1)  # = choose(n, n2)

        p1.read.emm  <- k * (1-rerr)^n1 * (rerr)^n2
	p2.read.emm  <- k * (1-rerr)^n2 * (rerr)^n1
	het.read.emm <- k * (0.5)^n

        c("P1" = p1.read.emm,
          "H1" = het.read.emm,
          "H2" = het.read.emm,
          "P2" = p2.read.emm
          )
    }

    ## list of 3x3 matrices representing probability of reads given the genetic
    ## state of each site (P1, H, P2)
    trans.prob.mats <- lapply(seq(nrow(snp.a.reads)), function(i) {
        ## construct column and row vectors for matrix multiplication
        a.read.probs <- read.probs.given.states(snp.a.reads[i,])
        b.read.probs <- read.probs.given.states(snp.b.reads[i,])

        col.vec <- matrix(a.read.probs, ncol=1)
        row.vec <- matrix(b.read.probs, nrow=1)

        ## probably not necessary but helpful for debugging
        names(col.vec) <- names(row.vec) <- names(a.read.probs)

        ## matrix multiplication; element i,j of the result gives the
        ## probability of the reads given that snp.a is c(P1, H, P2)[i] and
        ## snp.b is c(P1, H, P2)[j]
        col.vec %*% row.vec
    })

    function(r) {
        ## get.site.pair.trans is sourced from a generation specific file at runtime
        genetic.probs <- get.site.pair.trans(r)
        site.probs <- sapply(trans.prob.mats, function(trans.probs) {
            ## sum entries of element by element multiplication of the matrices
            ## to get the probabilty of the reads at this site given the genetic
            ## transition probabilities
            sum(trans.probs * genetic.probs)
        }
)
        ## return the sum of the log of the site.probs instead of the product of
        ## all the site probs for numerical stability
        sum(log(site.probs))
    }

}


quick.test <- function(pop.ad, i,j,k) {
    g <- function(a,b) {
        f <- get.obj.fun(pop.ad[,a,], pop.ad[,b,], 0); x <- seq(0, 0.5, 0.001); y <- sapply(x, f);
        x[which.max(y)]
   }

    gij <- g(i,j)
    gjk <- g(j,k)
    print(gij)
    print(gjk)
    print(g(i,k))
    print(-2*gij*gjk + gij + gjk)
}


h <- function(a,b) {
    f <- get.obj.fun(pop.ad[,a,], pop.ad[,b,], 0)

    result <- optim(,
                    fun,
                    method="L-BFGS-B",
                    lower=0,
                    upper=0.5)

}


g <- function(pop.ad,a,b) {
    f <- get.obj.fun(pop.ad[,a,], pop.ad[,b,], 0); x <- seq(0, 0.5, 0.01); y <- sapply(x, f);
    x[which.max(y)]
}


g.plot <- function(pop.ad,a,b) {
    f <- get.obj.fun(pop.ad[,a,], pop.ad[,b,], 0)
    x <- seq(0, 0.5, 0.01); y <- exp(sapply(x, f))
    message("max:   ", max(y),
            "\n",
            "min:   ", min(y),
            "\n",
            "ratio: ", max(y) / min(y))
    qplot(x, y) + geom_line()
}

estimate <- function(pop.ad, a,b,c) {
    f1 <- get.obj.fun(pop.ad[,a,], pop.ad[,b,], 0)
    f2 <- get.obj.fun(pop.ad[,b,], pop.ad[,c,], 0)
    f3 <- get.obj.fun(pop.ad[,a,], pop.ad[,c,], 0)
    obj.fun <- function(r1, r2) {
        f1(r1) + f2(r2) + f3(-2*r1*r2 + r1 + r2)
    }

    rs <- seq(0, 0.5, length.out=100)

    mat <- matrix(rep(0, 100*100), nrow=100)

    for (i in 1:100) {
        for (j in 1:100) {
            mat[i,j] <- obj.fun(rs[i], rs[j])
        }
    }

    w <- as.vector(which(mat==max(mat), arr.ind=TRUE))

    c(rs[w[1]], rs[w[2]])
}


masked.deep.calls.vcf <- function(vcf, parents, depth) {
    ad <- get.ad.array.w.p1.ref(vcf, parents)
    depths <- apply(ad, 1:2, sum)
    locations <- which(depths >= depth, arr.ind=TRUE)
    parent.rows <- which(getSAMPLES(vcf) %in% parents)
    non.parent.loc <- locations[! locations[,"row"] %in% parent.rows, ]
    # update vcf@fix based on the non.parent.loc to have not gt or ad or dp
    browser()
    print("DONE")
}
