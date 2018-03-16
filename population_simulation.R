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

get.sites <- function(file) {
    read.csv(file, header=FALSE)[,1]
}

## THIS IS DEPRICATED
## get.recomb.prob.fun <- function(peak, trough, flatness, start, end) {
##     ## function(x) {
##     ##     (peak - trough) * (2*x - 1)^(2*flatness) + trough
##     ## }
##     integral <- function(x) {
##         (peak - trough) * (2*x - 1)^(2*flatness + 1) / (4*flatness + 2) + trough*x
##     }
##     function(a, b) {
##         abs( (1e-6) * (integral(a) - integral(b)) )
##     }
## }


##' Return a scaled, mirrored, shifted logistic / sigmoid function
##'
##' The returned function will be as follows. The asymptotic maximum and minimum
##' values of function are specified by peak and trough. The function will be
##' shifted upward as necessary so that the minimum asymptotic value is the value
##' of trough rather than 0, and it will be shifted left so that f(0) is very
##' close to the peak value, though not exact because the function will
##' asymptotically approach the peak value as x approaches negative
##' infinity. The function will then be mirrored accross the y-axis so that it
##' is monotonically decreasing.
##' @title
##' @param peak the maximum value of the function
##' @param trough the minimum value of the function
##' @param d the distance on the x-axis from peak to trough in number of bases
##' @param sites numeric vector of site positions in base pairs
##' @param error allowed error at intersection of logistic function with y-axis
##' @return the semi-logistic function
##' @author Jason
get.recomb.profile.fun <- function(peak, trough, d, sites, error=0.0001) {
    start <- sites[1]
    end <- sites[length(sites)]
    x_shift <- d / 2  # logistic function shift param
    k <- log(1/error - 1) / (-x_shift)  # logistic function steepness param
    l <- end - start  # length

    ## implicitly return this function
    logistic <- function(x) {
        -(peak-trough) / (1 + exp(-k * (-x + x_shift))) + peak
    }

    function (positions) {
        positions <- positions - start  # offset
        ifelse(positions < l/2, logistic(positions), logistic(l-positions))
    }
}


get.recomb.probs <- function(peak, trough, d, sites) {
    f <- get.recomb.profile.fun(peak, trough, d, sites)

    sapply.pairs(sites, function(s1, s2) {
        megabases <- abs((s2-s1)/(1e6))  # base distance to megabase distance
        mid <- (s1 + s2) / 2
        cM.per.Mb <- f(mid)  # "average" centimorgan per megabase value

        ## return the centimorgan value which is probability of observing a
        ## recombination in the progeny (it is half of the probability of a
        ## physical recombination occurring between two of the chromotids
        ## between these sites and is thus theoretically limited to 50%)
        pmin(cM.per.Mb * megabases, 0.5)  # implicit return
    })
}

sapply.pairs <- function(vec, fun) {
    if (length(vec) < 2) {
        stop("first argument must have length >= 2")
    }
    v1 <- vec[-length(vec)]  # remove last entry
    v2 <- vec[-1]  # remove first entry

    mapply(fun, v1, v2)
}

gamete <- function(parent, recomb.probs) {
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
        recomb.probs <- 2 * recomb.probs  # convert from cM to prob of physical

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

cross <- function(p1, p2, recomb.probs) {
    if (! "genome" %in% class(p1)) {
        stop("p1 must be of type genome")
    }
    if (! "genome" %in% class(p2)) {
        stop("p2 must be of type genome")
    }
    if (length(p1) != length(p2)) {
        stop("p1 and p2 have different lengths")
    }
    if (FALSE) {
        ## TODO(Jason): check that cA and cB for p1 and p2 are all same length
        ## and length of recomb.probs is one less
    }
    if (FALSE) {
        ## TODO(Jason): warn if probability of recombination is greater than 50% across
    }

    offspring <- list()
    class(offspring) <- c("genome", class(offspring))

    offspring$cA <- gamete(p1, recomb.probs)
    offspring$cB <- gamete(p2, recomb.probs)

    offspring  # implicit return
}


create.ril.pop <- function(n.gen, n.mem, sites, recomb.probs) {
    if (n.generations < 2) {
        stop("n.generations must be at least 2")
    }
    n.sites <- length(sites)

    ## Use booleans to keep track of reference and alternate allele because it
    ## is more efficient than using 0 and 1
    ## Each parent has two homologous chromosomes (denoted here as cA and cB)
    ## and the parents are homozygous within and polymorphic between.
    p1 <- p2 <- list()  # empty object initialization
    p1$cA <- p1$cB <- rep(TRUE, times=n.sites)  # parent 1
    p2$cA <- p2$cB <- rep(FALSE, times=n.sites)  # parent 2
    class(p1) <- class(p2) <- c("genome", class(p1))

    f1 <- cross(p1, p2, recomb.fun)
    f2s <- replicate(n.members, cross(f1, f1, recomb.probs)) # n.sites x n.members matrix

    if (n.generations == 2) {
        return k(f2s)
    }

    fis <- f2s
    for (i in 3:n.generations) {
        fis <- apply(fis, 2, function(column){cross(column, column, recomb.fun)})
    }
    fis  # implicit return
}
