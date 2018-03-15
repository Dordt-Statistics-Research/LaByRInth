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

GENOME <- "genome"

get.sites <- function(file) {
    read.csv(file, header=FALSE)[,1]
}

get.recomb.prob.fun <- function(peak, trough, flatness, start, end) {
    ## function(x) {
    ##     (peak - trough) * (2*x - 1)^(2*flatness) + trough
    ## }
    integral <- function(x) {
        (peak - trough) * (2*x - 1)^(2*flatness + 1) / (4*flatness + 2) + trough*x
    }
    function(a, b) {
        abs( (1e-6) * (integral(a) - integral(b)) )
    }
}

get.recomb.prob.fun <- function(peak, trough, d, sites) { # dist in megabase
    start <- sites[1]
    end <- sites[length(sites)]
    peak <- peak-trough
    y <- 0.0001  # allowed error at peak of logistic function on x axis
    x_shift <- d / 2  # logistic function shift param
    k <- log(1/y - 1) / (-x_shift)  # logistic function steepness param
    l <- end - start  # length

    logistic <- function(x) {
        -peak / (1 + exp(-k * (-x + x_shift)))
    }

    f(pos) <- function(pos) {
        if (pos < l / 2) {
            logistic(pos)
        } else {
            logistig(end - pos)
        }
    }

    ## implicitly return this function
    function(a, b) {
        if (a > b) {
            temp <- a
            a <- b
            b <- temp
        }
        a <- a - start  # normalize
        b <- b - start  # normalize

        megabases <- (b-a)/(1e-6)  # base distance to megabase distance
        cM <- (f(b) + f(a)) / 2 * (b - a)
    }
}

sapply.pairs <- function(vec, fun) {
    if (length(vec) < 2) {
        stop("first argument must have length >= 2")
    }
    sapply(1:(length(vec)-1), function(i){
        fun(vec[i], vec[i+1])
    })
}

gamete <- function(parent, fun) {
    if (! "genome" %in% class(parent)) {
        stop("parent must be of type genome")
    }
    gamete <- sample(1:4, size=1)  # uniform randomly choose 1 which gamete
    if (gamete==1) {
        parent$cA
    } else if (gamete==2) {
        parent$cB
    } else {
        if (gamete==3) {
            base.gamete <- parent$cA
            recomb.gamete <-parent$cB
        } else {  # gamete==4
            base.gamete <- parent$cB
            recomb.gamete <-parent$cA
        }
        ## prob of recomb (in map units/centimorgans) is returned by
        ## sapply.pairs, but must be multiplied by 2 to account for the fact
        ## that we have already determined that the gamete will be one of the
        ## two chromatids that potentially cross over. The 'raw' map units
        ## cannot be used because if we randomly started with one of the two
        ## homologous chromasomes and swapped alleles based on the centimorgan
        ## probability, it would be very unlikely to return one of the two
        ## original chromosomes even though one of those should be returned 50%
        ## of the time. This is why we randomize equally over the four possible
        ## chromatids before worrying about the possibility of recombination.
        recomb.probs <- sapply.pairs(parent$sites, fun) * 2

        ## vector indicating if a physical recombination will occur between each
        ## pair of adjacent sites
        recomb.bools <- runif(length(recomb.probs)) < recomb.probs

        gamete.result <- base.gamete  # initialization
        odd.recombs <- FALSE  # indicate even or odd num of recombs up to point

        for (i in seq_along(recomb.bools)) {
            if (recomb.bools[i]) {
                # a recombination has occurred between index i and i+1 of
                # gamete.result
                odd.recombs <- !odd.recombs
            }
            if (odd.recombs) {
                ## the +1 is necessary because recomb.bools is one element
                ## smaller than base.gamete, recomb.gamete, and gamete.result
                ## and because index 1 of gamete.result is implied to be the
                ## same as index 1 of base.gamete by design
                gamete.result[i+1] <- recomb.gamete[i+1]
            }
        }
        gamete.result  # implicit return
    }


}

cross <- function(p1, p2, fun) {
    if (! "genome" %in% class(p1)) {
        stop("p1 must be of type genome")
    }
    if (! "genome" %in% class(p2)) {
        stop("p2 must be of type genome")
    }
    if (length(p1) != length(p2)) {
        stop("p1 and p2 have different lengths")
    }
    if (! identical(p1$sites, p2$sites)) {
        stop("site vectors differ for p1 and p2")
    }
    if (FALSE) {
        ## TODO(Jason): warn if probability of recombination is greater than 50% across
    }

    offspring <- list()
    class(offspring) <- c("genome", class(offspring))

    offspring$cA <- gamete(p1, fun)
    offspring$cB <- gamete(p2, fun)

    offspring  # implicit return
}


create.ril.pop <- function(n.generations, sites) {
    n.sites <- length(sites)

    ## Use booleans to keep track of reference and alternate allele
    ## Each parent has two homologous chromosomes (cA and cB)
    p1 <- p2 <- list()  # empty object initialization
    p1$cA <- p1$cB <- rep(TRUE, times=n.sites)  # parent 1
    p2$cA <- p2$cB <- rep(FALSE, times=n.sites)  # parent 2

    recomb.fun <- get.recomb.prob.fun(5, 0.0001, 10000, sites)

    f1 <- cross(p1, p2, recomb.fun)

    f2s <- replicate(n, cross(f1, f1)) # n.sites x n matrix
}
