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



## This file contains functions related to generating simulated
## populations. Currently, only single chromosome pair (homologous chromosomes)
## populations can be generated, but this is not a hinderance because these
## functions are used only for testing and because every homologous chromosome
## pair is imputed seperately anyway, there is no need to be concerned with
## populations of members with more than 2 chromosomes.


get.sites <- function(file) {
    read.csv(file, header=FALSE)[,1]
}


perc <- function(vec) {
    sum(vec) / length(vec)
}


get.member.homozygosity <- function(genome) {
    perc(genome$cA == genome$cB)
}


get.pop.homozygosity <- function(population) {
    mean(sapply(population, get.member.homozygosity))
}


list.of.cols <- function(matrix) {
    lapply(1:ncol(matrix), function(col) {matrix[, col]})
}



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
##' @author Jason Vander Woude
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


##' Breed a population from two full homozygous parents
##'
##' TRUE and FALSE are used as the two alleles. Both homologous chromosomes for
##' parent 1 are set to TRUE at all sites, and both homologous chromosomes for
##' parent 2 are set to FALSE at all sites. These parents are crossed to produce
##' an F1 genome which will necessarily have one homologous chromosome be TRUE
##' at all sites and the other will be FALSE at all sites. The F1 genome will be
##' selfed n.mem times to produce n.mem members of the F2 population. If n.gen
##' is 2, this population will be returned in a list, otherwise each member of
##' the F2 population will be selfed to produce the F3 generation, then each of
##' the F3 members will be selfed, etc. until an F{n.gen} population has been
##' generated which will then be returned in a list.
##' @title
##' @param n.gen the generation number of the population
##' @param n.mem the number of members in the population
##' @param recomb.probs vector of genetic map units between sites i.e. a vector
##' of probabilities of observing a recombination between each pair sites in the
##' population. In theory each probability should be less than 50%, but the code
##' will still work as expected if this is not the case
##' @return a list of elements of class genome with the genome length being 1
##' more than the length of recomb.probs
##' @author Jason Vander Woude
create.ril.pop <- function(n.gen, n.mem, recomb.probs) {
    if (n.gen < 2) {
        stop("n.gen must be at least 2")
    }
    if (n.mem < 1) {
        stop("n.mem must be at least 1")
    }
    n.sites <- length(recomb.probs) + 1

    ## Use booleans to keep track of reference and alternate allele because it
    ## is more efficient than using 0 and 1
    ## Each parent has two homologous chromosomes (denoted here as cA and cB)
    ## and the parents are homozygous within and polymorphic between.
    p1 <- p2 <- list()  # empty object initialization
    p1$cA <- p1$cB <- rep(TRUE, times=n.sites)  # parent 1
    p2$cA <- p2$cB <- rep(FALSE, times=n.sites)  # parent 2
    class(p1) <- class(p2) <- c("genome", class(p1))

    f1 <- cross(p1, p2, recomb.fun)

    f2s <- lapply(seq(n.mem), function(i){cross(f1, f1, recomb.probs)})

    if (n.gen == 2) {
        return (f2s)
    }

    fis <- f2s
    for (i in 3:n.gen) {
        fis <- lapply(fis, function(member){cross(member, member, recomb.probs)})
    }
    fis <- c(list(p1), list(p2), fis)  # prepend parents to list of progeny
    class(fis) <- c("population", class(fis))
    fis  # implicit return
}


## Choose the alleles that will be sampled
## This assumes that the parents are the first two entries in the population and
## will not be sampled
## TODO(Jason): allow for false reads
SamplePopulation <- function(pop, coverage=1, type="uniform", parents=TRUE) {
    if (! "population" %in% class(pop)) {
        stop("pop must be of class population")
    }
    if (! type %in% c("uniform")) {
        stop("currently only type supported is uniform")
    }

    n.variants <- length(pop)  # how many members
    n.sites <- length(pop[[1]]$cA)  # check length of one chrom of first variant
    ## 2 parents are not sampled and 2 chromosomes per member
    ## TODO(Jason): verify that this is how coverage is computed
    total <- (n.variants - 2) * n.sites * 2
    n.coverage.sites <- total * coverage
    variants <- sample(3:n.variants, size=n.coverage.sites, replace=TRUE)
    sites    <- sample(1:n.sites,    size=n.coverage.sites, replace=TRUE)
    chroms   <- sample(1:2,          size=n.coverage.sites, replace=TRUE)

    ## also sample both chromosomes of both parents at every site
    p.variants <- rep(1:2, each=2*n.sites)  # 1,1,...,1,2,...,2,2
    p.chroms <- rep(rep(1:2, each=n.sites), times=2)
    p.sites <- rep(1:n.sites, times=2*2)  # 1,2,3,...,n,1,2,3,...

    variants <- c(p.variants, variants)
    sites <- c(p.sites, sites)
    chroms <- c(p.chroms, chroms)

    result <- list()
    class(result) <- c("sample", class(result))
    result$population <- pop
    result$sample <- list.of.cols(rbind(variants, sites, chroms))

    result  # implicit return
}



EmptyVCF <- function(pop, sites) {
    if (! "population" %in% class(pop)) {
        stop("pop must be of class population")
    }

    ## TODO(Jason): either skip this or make it more efficient
    n.sites <- length(sites)
    for (member in pop) {
        if (length(member$cA) != n.sites ||
            length(member$cB) != n.sites) {
            stop("not all members of population have same length chromosomes as the sites vector")
        }
    }

    n.variants <- length(pop)  # how many members
    chr.name <- "UN"

    vcf <- list()
    class(vcf) <- "vcf"

    ## first two variants are parents, the rest are progeny
    vcf$variant.names <- c("P1", "P2", paste0("V", 1:(n.variants - 2)))
    vcf$chrom.names <- chr.name

    ## initialize the three main components of vcf
    vcf$GT <- array(NA_integer_, dim=c(n.sites, n.variants, 2))
    vcf$AD <- array(0, dim=c(n.sites, n.variants, 2))
    vcf$DP <- matrix(0, nrow=n.sites, ncol=n.variants)

    ## shorthand notation to allow naming all rows and columns together
    `cn<-` <- `colnames<-`
    `rn<-` <- `rownames<-`

    ## name all rows and columns of each matrix/array
    cn(vcf$AD) <- cn(vcf$DP) <- cn(vcf$GT) <- vcf$variant.names
    rn(vcf$AD) <- rn(vcf$DP) <- rn(vcf$GT) <- paste0(chr.name, ":", sites)

    ## create vcf header
    fields <- c("CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT")
    ## single string of all fields and variant names seperated by tabs
    tab.sep.line <- paste0("#", paste0(c(fields, vcf$variant.names), collapse="\t"))
    vcf$header <- c("##fileformat=VCFv4.0", tab.sep.line)

    ## set field data
    vcf$variants <- cbind(rep(chr.name, n.sites),        # CHROM
                          as.character(sites),           # POS
                          paste0(chr.name, "_", sites),  # ID
                          rep("A", n.sites),             # REF
                          rep("C", n.sites),             # ALT
                          rep("0", n.sites),             # QUALITY
                          rep("PASS", n.sites),          # FILTER
                          rep("", n.sites),              # INFO
                          rep("GT:AD:DP", n.sites))      # FORMAT
    cn(vcf$variants) <- fields

    vcf  # implicit return
}


FullPopulationVCF <- function(pop, sites){
    vcf <- EmptyVCF(pop, sites)

    ## properly set the DP field to 2 at every location
    vcf$DP <- apply(vcf$DP, 1:2, function(i){2})

    ## properly set the GT and AD fields based on the population
    for (var in seq_along(colnames(vcf$DP))) {
        for (site in seq_along(rownames(vcf$DP))) {
            a <- pop[[var]]$cA[site]
            b <- pop[[var]]$cB[site]
            ## a and b are true or false
            if (a && b) {  # both true / 1
                vcf$GT[site, var, ] <- c(1, 1)
                vcf$AD[site, var, ] <- c(0, 2)
            } else if (!a && !b) {  # both false / 0
                vcf$GT[site, var, ] <- c(0, 0)
                vcf$AD[site, var, ] <- c(2, 0)
            } else {  # different from each other
                vcf$GT[site, var, ] <- c(0, 1)
                vcf$AD[site, var, ] <- c(1, 1)
            }
        }
    }

    vcf  # implicit return
}


SamplePopulationVCF <- function(sample.result, sites){
    if (! "sample" %in% class(sample.result)) {
        stop("sample.result must be of class sample")
    }
    require(abind)

    pop <- sample.result$population
    sample <- sample.result$sample
    vcf <- EmptyVCF(pop, sites)

    ## properly set the AD fields based on the population
    for (observation in sample) {
        var    <- observation[1]
        site   <- observation[2]
        chrom  <- observation[3]
        ## add 1 to allele to compensate for 1 based indexing so allele will
        ## either equal 1 if the site is FALSE or 2 if the site is TRUE
        if (chrom == 1) {
            allele <- as.numeric(pop[[var]]$cA[site]) + 1
        } else {
            allele <- as.numeric(pop[[var]]$cB[site]) + 1
        }
        vcf$AD[site, var, allele] <- 1 + vcf$AD[site, var, allele]
    }

    ## properly set the GT and DP fields based on the AD field
    vcf$DP <- apply(vcf$AD, 1:2, sum)  # sum both allele counts at each site

    gt.1 <- apply(vcf$AD, 1:2, function(depths) {
        if (depths[1] == 0 && depths[2] == 0) {
            NA_integer_
        } else if (depths[1] == 0) {
            1
        } else if (depths[2] == 0) {
            0
        } else {  # both non-zero
            0
        }
    })

    gt.2 <- apply(vcf$AD, 1:2, function(depths) {
        if (depths[1] == 0 && depths[2] == 0) {
            NA_integer_
        } else if (depths[1] == 0) {
            1
        } else if (depths[2] == 0) {
            0
        } else {  # both non-zero
            1
        }
    })

    vcf$GT <- abind(gt.1, gt.2, along=3)

    vcf  # implicit return
}

##' .. content for \description{} (no empty lines) ..
##'
##' It is assumed that the population parents are the first two entries in
##' population that the vcf is built from. The data will also only check for
##' correct imputation in progeny and ignore the parents.
##' @title 
##' @param correct.vcf 
##' @param test.vcf 
##' @param mask.list a list of 2 element vectors indicating the 
##' @return 
##' @author Jason Vander Woude
checkAccuracy <- function(correct.vcf, test.vcf) {
    if (! "vcf" %in% class(correct.vcf)) {
        stop("correct.vcf must be of class vcf")
    }
    if (! "vcf" %in% class(test.vcf)) {
        stop("test.vcf must be of class vcf")
    }

    n.variants <- ncol(correct.vcf$GT)

    gt.corr <- correct.vcf$GT[, 3:n.variants, ]
    gt.test <- test.vcf$GT[, 3:n.variants, ]

    ## bind together the genotype data structures to allow efficiently applying
    ## a function to them without having to index. Note that both gt.corr and
    ## gt.test are already 3D structure where the first dimension is the site
    ## index, the second is the variant index, and the third is the chromosome
    ## index.
    bound <- abind(gt.corr, gt.test, along=4)

    types <- c("ref.called.ref",
               "ref.called.het",
               "ref.called.alt",
               "het.called.ref",
               "het.called.het",
               "het.called.alt",
               "alt.called.ref",
               "alt.called.het",
               "alt.called.alt")
    enum <- 1:9  # note that m later assumes length 9
    names(enum) <- types

    ret <- list()

    ## apply over the first two dimensions to check each site, variant pair
    all.call.types <- apply(bound, 1:2, function(data) {
        ## sum the two genotype vector entries. If the sum is 0, then both
        ## entries are 0 and it is homozygous reference. If the sum is 1, then
        ## one of the entires is 0 and the other is 1, so it is heterozygous,
        ## and if the sum is 2, then both entries are 1 and it is homozygous
        ## alternate.
        corr <- sum(data[,1])
        test <- sum(data[,2])

        if (corr == 0) {
            if (test == 0) {
                enum["ref.called.ref"]
            } else if (test == 1) {
                enum["ref.called.het"]
            } else {  # test == 2
                enum["ref.called.alt"]
            }
        } else if (corr == 1) {
            if (test == 0) {
                enum["het.called.ref"]
            } else if (test == 1) {
                enum["het.called.het"]
            } else {  # test == 2
                enum["het.called.alt"]
            }
        } else {  # corr == 2
            if (test == 0) {
                enum["alt.called.ref"]
            } else if (test == 1) {
                enum["alt.called.het"]
            } else {  # test == 2
                enum["alt.called.alt"]
            }
        }
    })
    ret$all.call.types <- all.call.types

    ## summarize all calls into the nine types
    m <- sapply(enum, function(call.type) {
        sum(result == enum)
    })
    m <- matrix(m, nrow=3, ncol=3, byrow=TRUE)
    rownames(m) <- c("ref", "het", "alt")
    colnames(m) <- c("called.ref", "called.het", "called.alt")
    ## m :=         Call
    ##
    ##              r h a
    ##              e e l
    ##              f t t
    ##    A   ref [ *     ]
    ##    c   het [   *   ]
    ##    t   alt [     * ]
    ret$call.summary <- m

    #m.small <- matrix(nrow=2, ncol=2)

    ret  ## implicit return
}
