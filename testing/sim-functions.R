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
    gamete <- sample(1:2, size=1)  # uniform randomly choose 1 which gamete
    if (gamete==1) {
        base.gamete <- parent$cA
        recomb.gamete <-parent$cB
    } else {  # gamete==2
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
    p1$cA <- p1$cB <- rep(FALSE, times=n.sites)  # parent 1
    p2$cA <- p2$cB <- rep(TRUE, times=n.sites)  # parent 2
    class(p1) <- class(p2) <- c("genome", class(p1))

    ## not necessary because f1 will have either cA or cB be all TRUE and the
    ## other all FALSE, but this allows for introducing error later
    f1 <- cross(p1, p2, recomb.probs)

    f2s <- lapply(seq(n.mem), function(i){cross(f1, f1, recomb.probs)})
    fis <- f2s

    if (n.gen != 2) {
        for (i in 3:n.gen) {
            fis <- lapply(fis, function(member){cross(member, member, recomb.probs)})
        }
    }

    fis <- c(list(p1), list(p2), fis)  # prepend parents to list of progeny
    class(fis) <- c("population", class(fis))
    fis  # implicit return
}


## Choose the alleles that will be sampled
## This assumes that the parents are the first two entries in the population and
## will not be sampled
## TODO(Jason): allow for false reads
sample.population <- function(pop, rerr, coverage=1, type="uniform", parents=TRUE) {
    if (! "population" %in% class(pop)) {
        stop("pop must be of class population")
    }
    if (! type %in% c("uniform")) {
        stop("currently only type supported is uniform")
    }

    n.taxa <- length(pop)  # how many members
    n.loci <- length(pop[[1]]$cA)  # check length of one chrom of first variant
    ## 2 parents are not sampled and 2 chromosomes per member
    ## TODO(Jason): verify that this is how coverage is computed
    total <- (n.taxa - 2) * n.loci * 2
    n.coverage.loci <- total * coverage * 1/2
    taxa <- sample(3:n.taxa, size=n.coverage.loci, replace=TRUE)
    loci    <- sample(1:n.loci,    size=n.coverage.loci, replace=TRUE)
    chroms   <- sample(1:2,          size=n.coverage.loci, replace=TRUE)
    erroneous <- runif(n=n.coverage.loci) < rerr

    ## also sample both chromosomes of both parents at every site
    p.taxa <- rep(1:2, each=2*n.loci)  # 1,1,...,1,2,...,2,2
    p.chroms <- rep(rep(1:2, each=n.loci), times=2)
    p.loci <- rep(1:n.loci, times=2*2)  # 1,2,3,...,n,1,2,3,...
    p.erroneous <- rep(FALSE, 4*n.loci)

    taxa <- c(p.taxa, taxa)
    loci <- c(p.loci, loci)
    chroms <- c(p.chroms, chroms)
    erroneous <- c(p.erroneous, erroneous)

    res <- array(0, dim=c(n.taxa, n.loci, 2))

    for (read.index in seq_along(erroneous)) {
        which.taxa <- taxa[read.index]
        which.chrom <- chroms[read.index]
        which.loci <- loci[read.index]

        read <- pop[[which.taxa]][[which.chrom]][which.loci]

        if (erroneous[read.index])
            read <- !read

        which.allele <- ifelse(read, 1, 0)

        res[which.taxa, which.loci, which.allele+1] <-
            res[which.taxa, which.loci, which.allele+1] + 1
    }

    return(res)

    result <- list()
    class(result) <- c("sample", class(result))
    result$population <- pop
    result$sample <- list.of.cols(rbind(taxa, loci, chroms))

    result  # implicit return
}


population.to.vcf <- function(pop, n.gen, out.file, sites=NA) {
    pop <- pop.to.geno.matrix(pop)

    ## add generic row and column names if missing
    if (length(sites) == 1 && is.na(sites)) {
        sites <- 1:ncol(pop)
    } else if (length(sites) != ncol(pop))
        stop("sites has incompatible length")

    parent.names <- c("P_1", "P_2")
    progeny.names <- paste0("F", n.gen, "_", seq(nrow(pop)-2))
    names <- c(parent.names, progeny.names)

    lines <- c(
        "##fileformat=VCFv4.2",
        paste0(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
            paste0(names, collapse="\t"))
        )

    chrom <- "UN"
    ## make sure sci notation is not used in writing vcf
    pos <- sapply(sites, format, scientific=FALSE)
    id <- paste0("ID_", seq(ncol(pop)))
    ref <- "A"
    alt <- "G"
    qual <- "."
    filter <- "PASS"
    info <- "."
    format <- "GT:AD"

    gt <- apply(pop, 1:2, function(num) {
        if (num == 0)
            "0/0"
        else if (num == 1)
            "0/1"
        else if (num == 2)
            "1/1"
        else
            stop("invalid num; this should never happen")
    })

    data <- t(gt)

    data <- cbind(chrom, pos, id, ref, alt, qual, filter, info, format, data)
    data <- apply(data, 1, paste0, collapse="\t")

    data <- c(lines, data)

    con <- file(out.file, "w")
    writeLines(data, con)
    close(con)
}


sample.to.vcf <- function(sample, n.gen, out.file, sites=NA) {
    ## add generic row and column names if missing
    if (length(sites) == 1 && is.na(sites)) {
        sites <- 1:ncol(sample)
    } else if (length(sites) != ncol(sample))
        stop("sites has incompatible length")

    parent.names <- c("P_1", "P_2")
    progeny.names <- paste0("F", n.gen, "_", seq(nrow(sample)-2))
    names <- c(parent.names, progeny.names)

    ## if (is.null(colnames(sample))) {
    ##     if (length(sites) == 1 && is.na(sites)) {
    ##         colnames(sample) <- paste0("loc_", seq(ncol(sample)))
    ##     } else {
    ##         if (length(sites) == ncol(sample))
    ##             colnames(sample) <- paste0("loc_", sites)
    ##         else
    ##             stop("sites has incompatible length")
    ##     }
    ## }
    ## if (is.null(rownames(sample))) {
    ##     parent.names <- c("P_1", "P_2")
    ##     progeny.names <- paste0("F", n.gen, "_", seq(nrow(sample)-2))
    ##     rownames(sample) <- c(parent.names, progeny.names)
    ## }

    sample.to.ad <- function(sample) {
        ## collapse the 3rd dimension (which indicates number of reference reads
        ## and number of alternate reads) into text form for vcf
        apply(sample, 1:2, paste0, collapse=",")
    }

    sample.to.gt <- function(sample) {
        ## collapse the 3rd dimension (which indicates number of reference reads
        ## and number of alternate reads) into directly inferred genotype text
        ## form for vcf
        apply(sample, 1:2, read.to.gt)
    }

    read.to.gt <- function(read) {
        ## read is expected to be a numeric vector of length 2 where the first
        ## index is number of reference reads and the sencond is the number of
        ## alternate reads
        if (any(is.na(read)))
            stop("a read is na. this should never happen")
        if (all(read != 0))     # read both alleles
            "0/1"
        else if (read[1] != 0)  # read only reference
            "0/0"
        else if (read[2] != 0)  # read only alternate
            "1/1"
        else                    # no reads at locus
            "./."
    }

    lines <- c(
        "##fileformat=VCFv4.2",
        paste0(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
            paste0(names, collapse="\t"))
        )

    chrom <- "UN"
    ## make sure sci notation is not used in writing vcf
    pos <- sapply(sites, format, scientific=FALSE)
    id <- paste0("ID_", seq(ncol(sample)))
    ref <- "A"
    alt <- "G"
    qual <- "."
    filter <- "PASS"
    info <- "."
    format <- "GT:AD:DP"

    gt <- t(sample.to.gt(sample))
    ad <- t(sample.to.ad(sample))
    dp <- t(apply(sample, 1:2, sum))

    data <- matrix(paste0(gt, ":", ad, ":", dp), nrow=nrow(gt))

    data <- cbind(chrom, pos, id, ref, alt, qual, filter, info, format, data)
    data <- apply(data, 1, paste0, collapse="\t")

    data <- c(lines, data)

    con <- file(out.file, "w")
    writeLines(data, con)
    close(con)
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


## sample.population.vcf <- function(sample.result, sites){
##     if (! "sample" %in% class(sample.result)) {
##         stop("sample.result must be of class sample")
##     }
##     require(abind)

##     pop <- sample.result$population
##     sample <- sample.result$sample
##     vcf <- EmptyVCF(pop, sites)

##     ## properly set the AD fields based on the population
##     for (observation in sample) {
##         var    <- observation[1]
##         site   <- observation[2]
##         chrom  <- observation[3]
##         ## add 1 to allele to compensate for 1 based indexing so allele will
##         ## either equal 1 if the site is FALSE or 2 if the site is TRUE
##         if (chrom == 1) {
##             allele <- as.numeric(pop[[var]]$cA[site]) + 1
##         } else {
##             allele <- as.numeric(pop[[var]]$cB[site]) + 1
##         }
##         vcf$AD[site, var, allele] <- 1 + vcf$AD[site, var, allele]
##     }

##     ## properly set the GT and DP fields based on the AD field
##     vcf$DP <- apply(vcf$AD, 1:2, sum)  # sum both allele counts at each site

##     gt.1 <- apply(vcf$AD, 1:2, function(depths) {
##         if (depths[1] == 0 && depths[2] == 0) {
##             NA_integer_
##         } else if (depths[1] == 0) {
##             1
##         } else if (depths[2] == 0) {
##             0
##         } else {  # both non-zero
##             0
##         }
##     })

##     gt.2 <- apply(vcf$AD, 1:2, function(depths) {
##         if (depths[1] == 0 && depths[2] == 0) {
##             NA_integer_
##         } else if (depths[1] == 0) {
##             1
##         } else if (depths[2] == 0) {
##             0
##         } else {  # both non-zero
##             1
##         }
##     })

##     vcf$GT <- abind(gt.1, gt.2, along=3)

##     vcf  # implicit return
## }

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
checkAccuracy <- function(correct.vcf, call.vcf) {
    if (! "vcf" %in% class(correct.vcf)) {
        stop("correct.vcf must be of class vcf")
    }
    if (! "vcf" %in% class(call.vcf)) {
        stop("test.vcf must be of class vcf")
    }

    ## TODO(Jason): use factors here instead of strings to save space while
    ## maintaining clarity of code
    get.type <- function(gt) {
        ## sum the two genotype vector entries. If the sum is 0, then both
        ## entries are 0 and it is homozygous reference. If the sum is 1, then
        ## one of the entires is 0 and the other is 1, so it is heterozygous,
        ## and if the sum is 2, then both entries are 1 and it is homozygous
        ## alternate.
        n <- sum(gt)

        if (is.na(n)) {
            "unknown"
        } else if (n == 0) {
            "reference"
        } else if (n == 1) {
            "heterozygous"
        } else if (n == 2) {
            "alternate"
        } else {
            stop(paste("get.type cannot accept the value", n))
        }
    }

    generalize.type <- function(type) {
        ifelse(
            type %in% c("reference", "alternate"),
            "homozygous",
            as.character(type)
        )
    }

    get.quality <- function(corr, call) {
        f <- function(cor, cal) {
            if (cor == "unknown") {
                "unknown"
            } else if (cal == "unknown") {
                "skipped"
            } else if (as.character(cor) == as.character(cal)) {
                "correct"
            } else if (cor == "heterozygous") {
                "partial"
            } else {
                "wrong"
            }
        }
        mapply(f, corr, call)  # implicit return
    }

    variants        <- colnames(correct.vcf$GT)
    n.variants      <- length(variants)
    sites           <- rownames(correct.vcf$GT)
    n.sites         <- length(sites)

    progeny.indices <- 3:n.variants
    progeny         <- variants[progeny.indices]
    n.progeny       <- length(progeny)

    gt.corr         <- correct.vcf$GT[, progeny.indices, ]
    gt.call         <- call.vcf$GT[, progeny.indices, ]
    ad.corr         <- correct.vcf$AD[, progeny.indices, ]

    ## as.numeric and as.factor are used to convert the matrix returned by apply
    ## into a vector
    corr.types      <- as.factor(apply(gt.corr, 1:2, get.type))
    gen.corr.types  <- as.factor(generalize.type(corr.types))
    call.types      <- as.factor(apply(gt.call, 1:2, get.type))
    gen.call.types  <- as.factor(generalize.type(call.types))
    if (is.null(ad.corr)) {
        ref.depths <- alt.depths <- NA_integer_
    } else {
        ref.depths      <- as.numeric(apply(ad.corr, 1:2, function(ad) {ad[1]}))
        alt.depths      <- as.numeric(apply(ad.corr, 1:2, function(ad) {ad[2]}))
    }
    qualities       <- get.quality(corr.types, call.types)

    ## the use of as.numeric abover will turn the matrix into a vector in column
    ## major order, and because columns represent variants, every block of
    ## n.sites entries correspond to the same variant/progeny. This is why
    ## variants are repeated with "each" while sites are repeated with "times"
    variants       <- rep(progeny, each=n.sites)
    sites          <- rep(sites, times=n.progeny)

    ## implicitly return data.frame
    data.frame(
        variant       = variants,
        site          = sites,
        type          = corr.types,
        call.type     = call.types,
        gen.type      = gen.corr.types,
        gen.call.type = gen.call.types,
        ref.depth     = ref.depths,
        alt.depth     = alt.depths,
        depth         = ref.depths + alt.depths,
        quality       = qualities
    )
}


pop.to.geno.matrix <- function(pop) {
    reduce.member <- function(member) {
        member$cA + member$cB
    }
    do.call(rbind, lapply(pop, reduce.member))
}


simulate.ril.geno.matrix <- function(n.samples, genetic.dists, n.gen) {
    pop.to.geno.matrix(create.ril.pop(n.gen, n.samples, genetic.dists))
}


pheno.from.geno <- function(geno, betas, baseline, sd, n.sim, seed=NULL) {
    ## SNPs are columns in geno and plants are rows
    ## length(betas) should be length(SNPS)=ncol(geno) and are effects
    ## baseline

    if (length(baseline) != 1)
        stop("There should only be one baseline value")
    if (length(sd) != 1)
        stop("There should only be one standard deviation")
    if (length(betas) != ncol(geno))
        stop("length of betas and ncol of geno should be equal")

    pheno.means <- geno %*% (betas/2) + baseline

    observed.phenos <- function() {
        rnorm(n=length(pheno.means),
              mean=pheno.means,
              sd=sd)
    }

    ## allow setting random generator seed to get the same result every time for
    ## testing
    if (! is.null(seed))
        set.seed(seed)

    pheno.types <- replicate(n.sim, observed.phenos())

    colnames(pheno.types) <- paste0("sim_", seq(ncol(pheno.types)))
    rownames(pheno.types) <- rownames(geno)

    pheno.types
}


test.pheno.creation <- function(n.samples, genetic.dists, n.gen, betas, n.sims, baseline=0, sd=1) {
    pop <- simulate.ril.geno.matrix(n.samples, genetic.dists, n.gen)
    ## remove parents from population
    pop <- pop[3:nrow(pop), ]
    colnames(pop) <- paste0("loc_", seq(ncol(pop)))
    rownames(pop) <- paste0("F", n.gen, "_", seq(nrow(pop)))

    message("Here is the population:")
    print(pop)

    phenos <- pheno.from.geno(pop, betas, baseline, sd, n.sims)

    message("Here are the phenos:")
    print(phenos)
}


