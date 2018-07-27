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


## library for working with VCF files
require(vcfR, quietly=T)
require(abind, quietly=T)


LabyrinthImpute <- function (vcf, parents, generation, out.file,
                             use.fwd.bkwd=FALSE, use.viterbi=TRUE,
                             calc.posteriors=TRUE,
                             read.err=0.05, geno.err=0.05,
                             parallel=TRUE, cores=4) {
    ## begin timer
    total.timer <- new.timer()
    print.labyrinth.header()


    ## parameter verification
    if (use.fwd.bkwd && use.viterbi) {
        display(0, "Cannot use both fwd.bkwd and viterbi algorithms to impute\n")
        stop()
    } else if (!use.fwd.bkwd && !use.viterbi) {
        display(0, "Must select either fwd.bkwd or viterbi algorithm to impute\n")
        stop()
    }

    if (use.fwd.bkwd && !calc.posteriors) {
        display(0, "When using fwd.bkwd, posterior probabilities must be calculated")
        display(0, "These probabilities will be included in the output file\n")
    }

    ## vcf load code
    if (! inherits(vcf, "vcfR")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }


    ## vcf requirement verification code
    timer <- new.timer()
    display(0, "Verifying requirements of vcf")
    if (! all(parents %in% getSAMPLES(vcf)))
        stop("Some or all parents are not in the vcf")

    non.biallelic <- which(! is.biallelic(vcf))
    if (length(non.biallelic) != 0) {
        display(0, "The following sites are not biallelic:")
        chroms <- getCHROM(vcf)[non.biallelic]
        positions <- getPOS(vcf)[non.biallelic]
        for (i in seq_along(chroms)) {
            display(1, "\tCHR:", chroms[i], "\tPOS:", positions[i])
        }
    }
    non.hom.poly <- which( ! parent.hom.and.poly(vcf, parents))
    if (length(non.hom.poly) != 0) {
        display(0, "The parents are not homozygous within and ",
                "polymorphic between at the following sites:")
        chroms <- getCHROM(vcf)[non.hom.poly]
        positions <- getPOS(vcf)[non.hom.poly]
        ad <- getAD(vcf)[non.hom.poly, parents]
        for (i in seq_along(chroms)) {
            display(1, "\tCHR:", chroms[i],
                    "\tPOS:", positions[i],
                    paste0("\t", parents[1], ":"), ad[i, 1],
                    paste0("\t", parents[2], ":"), ad[i, 2])
        }
    }

    if (length(non.biallelic) != 0 ||
        length(non.hom.poly) != 0)
        stop("The vcf must be filtered before imputing; please run LabyrinthFilter first")
    display(1, "Completed in ", timer(), "\n")


    ## site.pair.transition.probs variable will be loaded
    timer <- new.timer()
    display(0, "Loading data specific to F", generation)
    source(paste0("./new-transition-probs/F", generation, ".R"))
    display(1, "Completed in ", timer(), "\n")


    ## parental imputation and recombination rate estimates
    timer <- new.timer()
    display(0, "Imputing ", parents[1], " and ", parents[2], " and determining recombination rates")
    parental.results <- determine.parents.and.recombs(vcf,
                                                      parents,
                                                      read.err,
                                                      generation,
                                                      parallel,
                                                      cores)
    display(1, "Completed in ", timer(), "\n")


    ## transition probability code
    timer <- new.timer()
    display(0, "Generating transition probabilities")
    transition.structures <- get.transition.structures(vcf,
                                                       generation,
                                                       recomb.probs,
                                                       parental.state,
                                                       parallel,
                                                       cores)
    display(1, "Completed in ", timer(), "\n")


    ## emission probability code
    timer <- new.timer()
    display(0, "Generating emission probabilities")
    emission.structures <- get.emission.structures(vcf,
                                                   parents,
                                                   read.err,
                                                   geno.err,
                                                   generation,
                                                   parallel,
                                                   cores)
    display(1, "Completed in ", timer(), "\n")


    ## imputation code
    timer <- new.timer()
    display(0, "Imputing missing sites")
    imp.res <- impute(vcf, parents, emission.structures, transition.structures,
                      parallel, cores, use.fwd.bkwd, calc.posteriors)
    display(1, "Completed in ", timer(), "\n")

    ## new vcf creation code
    timer <- new.timer()
    display(0, "Creating new vcf with imputed data")
    vcf <- update.vcf(vcf, imp.res)
    display(1, "Completed in ", timer(), "\n")

    write.vcf(vcf, paste0(out.file, ".vcf.gz"))
    display(0, "LaByRInth imputation completed in ", total.timer())

    invisible(vcf) ## implicit return
}


## remove all sites that are not homozygous within and polymorphic between for
## the parents and remove all sites that are not biallelic
LabyrinthFilter <- function(vcf, parents, out.file) {
    display(0, "Checking if vcf is an object or file")
    if (! inherits(vcf, "vcfR")) {
        display(0, "Loading vcf")
        vcf <- read.vcfR(vcf, verbose=F)
    }
    display(0, "Checking if parents are in the vcf")
    if (! all(parents %in% getSAMPLES(vcf)))
        stop("Some or all parents are not in the vcf")
    display(0, "Checking for sites that are not biallelic")
    non.biallelic <- ! is.biallelic(vcf)
    if (length(non.biallelic) != 0) {
        display(0, "The following sites are not biallelic and will be removed:")
        chroms <- getCHROM(vcf)[which(non.biallelic)]
        positions <- getPOS(vcf)[which(non.biallelic)]
        for (i in seq_along(chroms)) {
            display(0, "\tCHR:", chroms[i], "\tPOS:", positions[i])
        }
    }
    display(0, "Checking for sites where parents are not homozygous within and polymorphic between")
    non.hom.poly <- ! parent.hom.and.poly(vcf, parents)
    if (length(non.hom.poly) != 0) {
        display(0, "The parents are not homozygous within and polymorphic between at the following sites which will be removed:")
        chroms <- getCHROM(vcf)[which(non.hom.poly)]
        positions <- getPOS(vcf)[which(non.hom.poly)]
        ad <- getAD(vcf)[which(non.hom.poly), parents]
        for (i in seq_along(chroms)) {
            display(0, "\tCHR:", chroms[i],
                    "\tPOS:", positions[i],
                    paste0("\t", parents[1], ":"), ad[i, 1],
                    paste0("\t", parents[2], ":"), ad[i, 2])
        }
    }

    mask <- !non.biallelic & !non.hom.poly
    vcf@fix <- vcf@fix[mask, ]
    vcf@gt <- vcf@gt[mask, ]
    write.vcf(vcf, paste0(out.file, ".vcf.gz"), mask=TRUE)
    display(0, paste("\nFiltering is complete;", sum(!mask), "of", length(mask), "sites removed"))
    invisible(vcf)  # implicit return
}


## get a logical vector indicating at which sites the parents are homozygous
## within and polymorphic between.
parent.hom.and.poly <- function(vcf, parents) {
    ## TODO(Jason): ask Jesse what should be done about ones that fail. Maybe we
    ## can soften this in LaByRInth so that we use confidence that we have read
    ## an allele from parent 1 or parent 2

    ## count number of zeros in numeric vector
    n.zeros <- function(numeric)
        sum(numeric==0)

    ## return logical vector indicating positions of nonzero values
    nonzero <- function(numeric)
        numeric!=0

    if (length(parents) != 2) {
        stop("Length of parents must be 2")
    }
    par.mat <- getAD(vcf)[ , parents]
    apply(par.mat, 1, function(row) {
        (n.zeros(ad.to.num(row[1])) > 0          # at least one allele not read in P1
            && n.zeros(ad.to.num(row[2])) > 0    # at least one allele not read in P2
            && ! any(nonzero(ad.to.num(row[1])) & nonzero(ad.to.num(row[2]))) # no common reads
            && n.zeros(ad.to.num(row[1])) != 2    # TODO(Jason): remove condition and impute parents
            && n.zeros(ad.to.num(row[2])) != 2    # TODO(Jason): remove condition and impute parents
        )
    })
}


getSAMPLES <- function(vcf) {
    ## column names are the samples except the first column which is "FORMAT"
    colnames(vcf@gt)[-1]
}


## More convenient strsplit if length of vector is 1
str.split <- function(str, sep) {
    if (length(str) != 1) {
        warning("Only splitting the first element of the vector")
    }
    strsplit(str, sep)[[1]]
}


## return the total probabilities of each state at each site. Adapted from
## https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
## forward and backward probabilities are not technically correct because they
## are frequently being normalized. This is done to decrease the numerical
## instability that results when the probabilities become very low
fwd.bkwd <- function(emm, trans) {
    normalize <- function(mat, col) {
        s <- sum(mat[, col])
        if (s == 0){
            mat[, col] <- 1 / nrow(mat)
        } else {
            mat[, col] <- mat[, col] / s
        }
        mat
    }

    n.states <- nrow(emm)
    n.sites <- ncol(emm)

    f.probs <- b.probs <- matrix(data=0, nrow=n.states, ncol=n.sites)

    ## forward probabilities
    start.probs <- rep(1/n.states, n.states)
    f.probs[, 1] <- start.probs * emm[, 1]
    f.probs <- normalize(f.probs, 1)

    for (site in 2:n.sites) {
        t.index <- site - 1  # transition structure index
        prev.site <- site - 1
        for (to in 1:n.states) {
            f.probs[to, site] <-
                emm[to, site] * sum(trans[ , to, t.index] * f.probs[, prev.site])
         }
        f.probs <- normalize(f.probs, site)
    }

    ## backward probabilities
    end.probs <- rep(1, n.states)
    b.probs[, n.sites] <- end.probs
    b.probs <- normalize(b.probs, n.sites)

    for (site in (n.sites-1):1) {
        t.index <- site
        next.site <- site + 1
        for (from in 1:n.states) {
            b.probs[from, site] <-
                sum(trans[from, , t.index] * b.probs[, next.site] * emm[ , next.site])
        }
        b.probs <- normalize(b.probs, site)
    }

    res <- f.probs * b.probs
    for (site in 1:n.sites) {
        res <- normalize(res, site)
    }
    res
}


## If emm.log is TRUE, then the data passed in emm should already be
## log-scaled. Similarly with trans.log
viterbi <- function(emm, trans, emm.log=FALSE, trans.log=FALSE) {
    if (!emm.log)
        emm <- log(emm)
    if (!trans.log)
        trans <- log(trans)

    n.states <- nrow(emm)
    n.sites <- ncol(emm)

    path.tracker <- matrix(data=0, nrow=n.states, ncol=n.sites)

    start.probs <- log(rep(1/n.states, n.states))
    probs <- start.probs + emm[, 1]

    for (site in 2:n.sites) {
        t.index <- site - 1  # transition structure index
        prev.site <- site - 1
        new.probs <- probs
        for (to in 1:n.states) {
            x <- probs + trans[ , to, t.index]
            new.probs[to] <- max(x) + emm[to, site]
            path.tracker[to, site] <- which.max(x)
        }
        probs <- new.probs
    }

    ## reconstruct the path
    path <- rep(NA, n.sites)
    best.state <- which.max(probs)
    for (site in n.sites:1) {
        path[site] <- best.state
        best.state <- path.tracker[best.state, site]
    }

    list(path=path, prob=max(probs))
}


## viterbi <- function(emm, trans) {
##     n.states <- nrow(emm)
##     n.sites <- ncol(emm)

##     path.tracker <- matrix(data=0, nrow=n.states, ncol=n.sites)

##     start.probs <- log(rep(1/n.states, n.states))
##     probs <- start.probs + log(emm[, 1])

##     for (site in 2:n.sites) {
##         t.index <- site - 1  # transition structure index
##         prev.site <- site - 1
##         new.probs <- probs
##         for (to in 1:n.states) {
##             x <- probs + log(trans[ , to, t.index])
##             new.probs[to] <- max(x) + log(emm[to, site])
##             path.tracker[to, site] <- which.max(x)
##         }
##         probs <- new.probs
##     }

##     ## reconstruct the path
##     path <- rep(NA, n.sites)
##     best.state <- which.max(probs)
##     for (site in n.sites:1) {
##         path[site] <- best.state
##         best.state <- path.tracker[best.state, site]
##     }

##     list(path=path, prob=max(probs))
## }


## returns the emission probabilities for a given variant and chromosome. Note
## that the structure of the returned data is different than a traditional
## emission probability matrix. The first dimension of the data corresponds to
## the chromosome, and the second dimension corresponds to the sample. reads
## should be a vector of strings where each string is a numeric value followed
## by a comma and another numeric value. This is the format that the vcfR
## package stores the reads in. Computing all of the emission probabilities
## before they are needed may use more memory than if they were computed as
## needed; the reason for doing all of them before imputing is because in an
## earlier version of the code, the data was imputed under a variety of
## different models, and each model could utilize the same emission
## probabilities, so it was faster to compute all of the emission probabilities
## only once.
get.emission.structures <- function(vcf, parents, rerr, gerr, generation, parallel=F, cores=1) {

    ## determine at which sites parent 1 is reference
    str.ad <- getAD(vcf)[ , parents[1]]
    p1.is.ref <- sapply(str.ad, function(str) {
        ## split the string and check if reference read is nonzero
        ad.to.num(str)[1] != 0
    })

    perc.het <- 0.5^(generation - 1)
    perc.p1 <- perc.p2 <- (1 - perc.het) / 2

    states <- 1:4
    names(states) <- c("P1", "H1", "H2", "P2")

    ## probability of a read given a genotype/state
    prob <- function(state, read, p1.ref) {
        ad.num <- ad.to.num(read)
        n.ref <- ad.num[1]  # num reference reads
        n.alt <- ad.num[2]  # num alternate reads
        n <- n.ref + n.alt

        n1 <- ifelse(p1.ref, n.ref, n.alt)
        n2 <- ifelse(p1.ref, n.alt, n.ref)

        names(n1) <- NULL  # these pick up names from ifelse
        names(n2) <- NULL  # debugging is easier without them

	p1.read.emm <- choose(n, n1) * (1-rerr)^n1 * (rerr)^n2
	p2.read.emm <- choose(n, n2) * (1-rerr)^n2 * (rerr)^n1
	het.read.emm <- choose(n, n1) * (0.5)^n

	err.prob <- gerr * (perc.p1 * p1.read.emm + perc.p2 * p2.read.emm + perc.het * het.read.emm)

        if (state == states["P1"]) {
	    ret <- (1-gerr) * p1.read.emm + err.prob
            ## ret <- choose(n, n1) * (1-rerr)^n1 * (rerr)^n2
        } else if (state == states["P2"]) {
	    ret <- (1-gerr) * p2.read.emm + err.prob
            ## ret <- choose(n, n2) * (1-rerr)^n2 * (rerr)^n1
        } else if (state == states["H1"] || state == states["H2"]) {
	    ret <- (1-gerr) * het.read.emm + err.prob
            ## ret <- choose(n, n1) * (0.5)^n
        } else {
            stop("Invalid state; this should never happen")
        }
        ret  # implicit return
    }

    reads.emm.probs <- function(reads, p1.is.ref) {
        names(reads) <- NULL  # makes debugging easier
        ret.val <- do.call(rbind,
                lapply(states, function(state) {
                    ## mapply will repeat state as many times as necessary
                    mapply(FUN=prob, state, reads, p1.is.ref)
                }))

        ## normalize the emission probabilities so that the forward-backward
        ## algorithm is more stable

        ret.val <- ret.val / colSums(ret.val)[col(ret.val)]  # normalize each
                                        # column

        l <- length(states)

        ret.val[is.na(ret.val)] <- 1 / l  # default normalize na columns

        ret.val

    }

    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)
    samples <- getSAMPLES(vcf)
    ad <- getAD(vcf)

    listapply <- get.lapply(parallel, cores)


    ## Progress bar code
    ## -------------------------------------------------------------------------
    progress.env <- new.env()
    thefifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prog.env <- progress.env
    n.jobs <- length(u.chroms) * length(samples)
    ## -------------------------------------------------------------------------

    ## -------------------------------------------------------------------------
    writeBin(0, thefifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
    }  # else the forked process handles this
    ## -------------------------------------------------------------------------
    ## Progress bar code


    ret.val <- lapply(u.chroms, function(chrom) {
        ret.val.2  <- listapply(samples, function(sample) {
            ## emissions will be useless for parents, but that is fine
            ret.val.3 <- reads.emm.probs(ad[chroms==chrom,
                                            colnames(ad)==sample], p1.is.ref[chroms==chrom])

            ## Progress bar code
            ## -----------------------------------------------------------------
            writeBin(1/n.jobs, thefifo)  # update the progress bar info
            if (!parallel) {  # if running in serial mode
                prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
            }  # else the forked process handles this
            ## -----------------------------------------------------------------
            ## Progress bar code


            ret.val.3
        })
        names(ret.val.2) <- samples
        ret.val.2
    })


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code


    names(ret.val) <- u.chroms
    ret.val  # implicitly returned
}


get.transition.structures <- function(vcf, generation, recomb.probs,
                                      parental.state, parallel, cores) {

    states <- 1:4
    names(states) <- c("P1", "H1", "H2", "P2")

    listapply <- get.lapply(parallel, cores)

    trans.probs <- function(chrom) {
        trans.matrices <- lapply(recomb.probs[[chrom]], get.trans)

        result <- do.call(abind3, trans.matrices)
        result[result < 0] <- 0  # handle numerical errors
        result
    }

    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)


    ## Progress bar code
    ## -------------------------------------------------------------------------
    progress.env <- new.env()
    thefifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prog.env <- progress.env
    n.jobs <- length(u.chroms)
    ## -------------------------------------------------------------------------

    ## -------------------------------------------------------------------------
    writeBin(0, thefifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
    }  # else the forked process handles this
    ## -------------------------------------------------------------------------
    ## Progress bar code

    ret.val <- lapply(u.chroms, function(chrom) {

        ret.val.2 <- trans.probs(chrom)

        ## Progress bar code
        ## -----------------------------------------------------------------
        writeBin(1/n.jobs, thefifo)  # update the progress bar info
        if (!parallel) {  # if running in serial mode
            prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
        }  # else the forked process handles this
        ## -----------------------------------------------------------------
        ## Progress bar code

        ret.val.2

    })


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code


    names(ret.val) <- u.chroms
    ret.val  # implicitly returned
}


impute <- function(vcf, parents, emm.structures, trans.structures, parallel,
                   cores, use.fwd.bkwd, calc.posteriors) {

    listapply <- get.lapply(parallel, cores)

    p1.is.ref <- is.parent.ref(vcf, parents[1])

    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)
    samples <- getSAMPLES(vcf)

    impute.sample.chrom <- function(sample, chrom) {

        n.sites <- sum(chroms==chrom)  # boolean addition
        p1.ref <- p1.is.ref[chroms==chrom]

        res <- list()
        res$posteriors <- NULL

        if (use.fwd.bkwd || calc.posteriors) {
            if (sample == parents[1]) {
                posterior.mat <- rbind(rep(1, n.sites),
                                       0,
                                       0,
                                       0)
            } else if (sample == parents[2]) {
                posterior.mat <- rbind(0,
                                       0,
                                       0,
                                       rep(1, n.sites))
            } else {
                emm <- emm.structures[[chrom]][[sample]]
                trans <- trans.structures[[chrom]]

                posterior.mat <- fwd.bkwd(emm, trans)
            }

            ref <- ifelse(p1.ref, posterior.mat[1, ], posterior.mat[4, ])
            alt <- ifelse(p1.ref, posterior.mat[4, ], posterior.mat[1, ])
            het <- posterior.mat[2, ] + posterior.mat[3, ]

            posteriors <- rbind(ref, het, alt)

            phred.scaled <- apply(posteriors, 1:2, function(x) {
                ## 1-x is the probability the call is wrong
                ## use min to prevent infinity from getting throughk
                min(-10*log((1-x), base=10), 100)
            })

            res$posteriors <- apply(phred.scaled, 2, paste0, collapse=",")

            ## res$posteriors <- array(c(ref,   # homozygous reference allele prob
            ##                           het,   # heterozygous prob
            ##                           alt),  # homozygous alternate allele prob
            ##                         dim=c(n.sites, 1, 3))

            if (use.fwd.bkwd) {
                res$gt <- apply(posteriors, 2, function(states) {
                    c("0/0","0/1","1/1")[which.max(states)]
                })
            }
        }

        if (!use.fwd.bkwd) {
            if (sample == parents[1]) {
                best.path <- rep(1, n.sites)
            } else if (sample == parents[2]) {
                best.path <- rep(4, n.sites)
            } else {
                emm <- emm.structures[[chrom]][[sample]]
                trans <- trans.structures[[chrom]]

                path.and.prob <- viterbi(emm, trans)
                best.path <- path.and.prob$path
                ## In best.path,
                ##    1: parent 1
                ##    2: het type 1
                ##    3: het type 2
                ##    4: parent 2
            }

            res$gt <- sapply(seq_along(best.path), function(index) {
                state <- best.path[index]

                if (state == 1)
                    ifelse(p1.ref[index], "0/0", "1/1")
                else if (state == 2)
                    "0/1"
                else if (state == 3)
                    "1/0"
                else if (state == 4)
                    ifelse(p1.ref[index], "1/1", "0/0")
                else
                    stop("invalid state; this should never happen")
            })

        }

        res  # implicit return
    }


    get.posteriors <- function(list) {
        lapply(list, function(result) {
            result$posteriors
        })
    }


    get.listed.gt <- function(list) {
        lapply(list, function(result) {
            result$gt
        })
    }


    ## Progress bar code
    ## -------------------------------------------------------------------------
    progress.env <- new.env()
    thefifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prog.env <- progress.env
    n.jobs <- length(u.chroms) * length(samples)
    ## -------------------------------------------------------------------------

    ## -------------------------------------------------------------------------
    writeBin(0, thefifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
    }  # else the forked process handles this
    ## -------------------------------------------------------------------------
    ## Progress bar code


    ## TODO(Jason): run chroms w/in samples or even do conditional to run larger
    ## inside. This will have an advantage because for a given chromosome, the
    ## samples will all have a similar length and so the compute times for
    ## imputation should be fairly similar. When having the samples in the inner
    ## loop, the chromosomes for that sample have drastically different sizes
    ## and can have quite varying runtimes. Ideally,

    imputed.chroms <- lapply(u.chroms, function(chrom) {

        imputed.samples <- listapply(samples, function(sample) {

            ret.val <- impute.sample.chrom(sample, chrom)
            ## display(1, "Imputed chromosome ", chrom, " of sample ", sample)

            ## Progress bar code
            ## -----------------------------------------------------------------
            writeBin(1/n.jobs, thefifo)  # update the progress bar info
            if (!parallel) {  # if running in serial mode
                prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
            }  # else the forked process handles this
            ## -----------------------------------------------------------------
            ## Progress bar code


            ret.val  # implicit return
        })

        ret.val <- list()
        ret.val$posteriors <- do.call(cbind, get.posteriors(imputed.samples))
        ret.val$gt <- do.call(cbind, get.listed.gt(imputed.samples))

                                        #do.call(abind2, imputed.samples)  # implicit return
        ret.val  # implicit return

    })

    ret.val <- list()
    ret.val$posteriors <- do.call(rbind, get.posteriors(imputed.chroms))
    ret.val$gt <- do.call(rbind, get.listed.gt(imputed.chroms))

    ## combine columns of each imputed sample
    ## ret.val <- do.call(abind1, imputed.chroms)


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code


    ret.val  # implicit return
}


save.trans <- function(data, chrom, model) {
    dir <- "./.saved_trans"
    dir.create(dir, recursive=TRUE)
    saveRDS(data, paste0(dir, "/", paste0(model, collapse=""), ".rds"))
}


## get.trans <- function(data, chrom, model) {
##     dir <- "./.saved_trans"
##     readRDS(data, paste0(dir, "/", paste0(model, collapse=""), ".rds"))
## }


is.parent.ref <- function(vcf, parent) {
    str.ad <- getAD(vcf)[ , parent]
    sapply(str.ad, function(str) {
        ## split the string and check if reference read is nonzero
        ad.to.num(str)[1] != 0
    })
}


## modify the vcf in place to add imputed calls in gt
update.vcf <- function(vcf, impute.res) {
    ## phred.scaled.res <- apply(impute.res, 1:3, function(x) {
    ##     ## 1-x is the probability the call is wrong
    ##     ## use min to prevent infinity from getting throughk
    ##     min(-10*log((1-x), base=10), 100)
    ## })

    ## ## rounded.res <- apply(phred.scaled.res, 1:3, round)

    ## gp <- apply(phred.scaled.res, 1:2, paste0, collapse=",")

    ## gt <- apply(impute.res, 1:2, function(states) {
    ##     c("0/0","0/1","1/1")[which.max(states)]
    ## })

    gt <- impute.res$gt
    gp <- impute.res$posteriors
    ad <- getAD(vcf)

    n.rows <- nrow(ad)
    if (n.rows != nrow(gt))
        stop("conflict in number of rows; should never happen")

    if (is.null(gp)) {
        concat <- matrix(paste0(gt, ":", ad), nrow=n.rows)
        vcf@gt <- cbind("GT:AD", concat)
    } else {
        concat <- matrix(paste0(gt, ":", ad, ":", gp), nrow=n.rows)
        vcf@gt <- cbind("GT:AD:GP", concat)
    }

    colnames(vcf@gt) <- c("FORMAT", colnames(ad))
    rownames(vcf@gt) <- rownames(ad)

    vcf  # implicit return
}


getAD <- function(vcf) {
    ad <- extract.gt(vcf, "AD")
    ad[is.na(ad)] <- "0,0"  # replace NA entries
    ad
}


getGT <- function(vcf) {
    gt <- extract.gt(vcf, "GT")
    gt[is.na(gt)] <- "./."  # replace NA entries
    gt
}


getGP <- function(vcf) {
    extract.gt(vcf, "GP")
}


## convert string representation to numeric vector
ad.to.num <- function(str) {
    as.numeric(str.split(str, ","))

}


## convert string representation to numeric vector
gt.to.num <- function(str) {
    str <- gsub("\\|", "/", str)  # replace '|' with '/'
    suppressWarnings(as.numeric(str.split(str, "/")))

}


gp.to.num <- function(str) {
    ad.to.num(str)
}


all.bool.vec <- function(n) {
    if (n > 16)
        stop("all.bool.vec function not supported for n > 16")

    lapply(0:(2^n-1), function(config) {
        as.logical(intToBits(config)[1:n])
    })
}


## The abind{i} functions are useful for cleaning the code for do.call(abind,
## some.list) calls because otherwise the 'along' argument has to be appended to
## the 'some.list' argument
abind1 <- function(...) {
    abind(..., along=1)
}


abind2 <- function(...) {
    abind(..., along=2)
}


abind3 <- function(...) {
    abind(..., along=3)
}


## Progress monitor code from https://stackoverflow.com/questions/27726134/
## how-to-track-progress-in-mclapply-in-r-in-parallel-package
## TODO(Jason): don't use prefs$fifo, but instead try to use a fifo variable
## in the progress.env environment
ProgressMonitor <- function(env, parallel) {
    local({
        f <- fifo(tempfile(), open="w+b", blocking=T)
        if (!parallel) {  # don't fork if running serially
            return(f)
        }
        if (inherits(parallel:::mcfork(), "masterProcess")) {
            progress <- 0.0
            while (progress < 1 && !isIncomplete(f)) {
                progress <- PrintProgress(f, progress)
            }
            parallel:::mcexit()
        }
        f  # implicit return
    }, envir=env)
}


PrintProgress <- function(f, curr.prog) {
    msg <- readBin(f, "double")
    progress <- curr.prog + as.numeric(msg)
    cat(sprintf(paste0("    * ",  "Progress: %.2f%%\r"), progress * 100))
    progress  # implicit return
}


## TODO(Jason): Get this function to work; problem likely environments
InitiateProgress <- function(prefs, n.jobs) {
    ## progress.env <- new.env()
    ## prefs$fifo <- ProgressMonitor(progress.env)
    ## assign("progress", 0.0, envir=progress.env)
    ## prefs$prog.env <- progress.env
    ## prefs$n.jobs <- n.jobs
    ## assign("prefs$n.jobs", n.jobs, envir=.GlobalEnv)
}


MakeProgress <- function(fifo, n.jobs, prog.env, parallel) {
    writeBin(1/n.jobs, fifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(fifo, prog.env$progress)
    }  # else the forked process handles this
}


get.lapply <- function(parallel, cores=1) {
    serial.lapply <- function(..., mc.preschedule=F, mc.cores=cores) {
        lapply(...)
    }

    ## parallel.lapply <- function(..., mc.preschedule=F, mc.cores=cores) {
    ##     mclapply(..., mc.preschedule=mc.preschedule, mc.cores=mc.cores)
    ## }

    parallel.lapply <- function(...) {
        mclapply(..., mc.preschedule=F, mc.cores=cores)
    }

    ## implicit return
    if (parallel && cores != 1) {
        require(parallel, quietly=T)
        parallel.lapply
    } else {
        lapply
    }
}


new.timer <- function() {
    start <- Sys.time()
    function() {
        elapsed <- difftime(Sys.time(), start)
        runtime <- as.numeric(elapsed)
        units <- attr(elapsed, "units")
        paste(round(runtime, 2), units)
    }
}


display <- function(indent, ...) {
    if (! is.numeric(indent))
        stop("indent must be numeric")
    message(rep("   ", indent), " * ", ...)
}


print.labyrinth.header <- function() {

    ## the image looks funny because '\' in the displayed image must be '\\' in the code
    writeLines("")
    writeLines(" _____________________________________________________________________")
    writeLines("|          __          ____        ____  _____                        |")
    writeLines("|         / /         / __ \\      / __ \\/_  _/                        |")
    writeLines("|        / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __      |")
    writeLines("|       / /   / _  \\/ _  _/\\ \\/ / _  _/ / / /  | / /_  __/ /_/ /      |")
    writeLines("|      / /___/ /_/ / /_\\ \\  \\  / / \\ \\_/ /_/ /||/ / / / / __  /       |")
    writeLines("|     /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/        |")
    writeLines("|                                                                     |")
    writeLines("| LaByRInth: Low-coverage Biallelic R Imputation                      |")
    writeLines("| Copyright 2017 Jason Vander Woude and Nathan Ryder                  |")
    writeLines("| Licensed under the Apache License, Version 2.0                      |")
    writeLines("| Source code: github.com/Dordt-Statistics-Research/LaByRInth         |")
    writeLines("| Based on LB-Impute: github.com/dellaporta-laboratory/LB-Impute      |")
    writeLines("| Funding received from the National Science Foundation (IOS-1238187) |")
    writeLines("|_____________________________________________________________________|")
    writeLines("")

}


analyze <- function(orig, mask, imp) {

    load <- function(vcf) {
        if (! inherits(vcf, "vcfR"))
            read.vcfR(vcf, verbose=F)
        else
            vcf
    }

    check.correct <- function(orig.gt, imp.gt) {
        gsub("\\|", "/", orig.gt)  # replace '|' with '/'
        if (orig.gt == "1/0")
            orig.gt <- "0/1"
        orig.gt == imp.gt
    }

    check.masked <- function(orig.ad, mask.ad) {
        orig.ad != mask.ad
    }

    get.gt <- function(vcf)
        extract.gt(vcf, "GT")

    orig <- load(orig)
    mask <- load(mask)
    imp <- load(imp)

    depths <- sapply(getAD(orig), ad.to.num)

    post.probs <- getGP(imp)

    max.post.probs <- sapply(post.probs, function(gp) {
        max(gp.to.num(gp))
    })

    data <- data.frame(
        ref.depth  = depths[1, ],
        alt.depth  = depths[2, ],
        depth      = depths[1, ] + depths[2, ],
        correct    = mapply(FUN=check.correct, getGT(orig), getGT(imp)),
        masked     = mapply(FUN=check.masked, getAD(orig), getAD(mask)),
        max.prob   = max.post.probs
        ## post.probs = post.probs
    )

    data  # implicit return
}


quick.check <- function() {
    ## an <- analyze("./analysis/filtered_lakin_fuller.vcf", "./analysis/masked_files_LF_from_sandesh/masked_01_filtered_lakin_fuller.vcf", "./thresh-90.vcf.gz")
    ## an <- analyze("./analysis/first-25-orig.vcf.gz",
    ##               "./analysis/first-25-masked.vcf.gz", "./analysis/first-25-imp-f3.vcf.gz")
    ## interest <- an$correct[an$masked==T]; sum(interest) / length(interest)





    ## an <- analyze("./analysis/filtered_lakin_fuller.vcf", "./analysis/masked_files_LF_from_sandesh/masked_01_filtered_lakin_fuller.vcf", "~/Desktop/imputed/imputed-mask-01-f5-1e8.vcf.gz")

    for (qual in seq(from=0, to=100, by=5)) {
        interest <- an$correct[an$masked & an$max.prob >= qual]
        display(1, "thresh: ", qual, "\tn: ", length(interest), "\tquality: ", sum(interest) / length(interest))

    }
}


LabyrinthQC <- function(vcf) {
    ## vcf load code
    if (! inherits(vcf, "vcfR")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }

    gp <- getGP(vcf)
    gt <- getGT(vcf)
    ad <- getAD(vcf)

    highest.posterior <- apply(gp, 1:2, function(gp) {
        max(gp.to.num(gp))
    })

    ## y <- sapply(0:100, function(x) {sum(highest.posterior >= x) / length(highest.posterior)})
    ## plot((0:100)/100, y, type="l", main="Imputation Posterior Probabilities",
    ## xlab="Posterior Probabilities",
    ## ylab="Proportion of Sites with Higher Posterior",
    ## ylim=c(0,1), xlim=c(0,1))

    thresh.mask <- apply(highest.posterior, 1:2, `>=`, 90)

    gt <- ifelse(thresh.mask, gt, "./.")

    concat <- matrix(paste0(gt, ":", ad, ":", gp), nrow=nrow(ad))

    vcf@gt <- cbind("GT:AD:GP", concat)
    colnames(vcf@gt) <- c("FORMAT", colnames(ad))
    rownames(vcf@gt) <- rownames(ad)


    write.vcf(vcf, "thresh-90.vcf.gz")

    print("done")
}


get.ad.array <- function(vcf) {
    ads <- apply(getAD(vcf), 1:2, ad.to.num)
    abind(ads[1, , ], ads[2, , ], along=3)
}


determine.parents.and.recombs <- function(vcf, parents, rerr,
                                          generation, parallel=F, cores=1) {
    timer <- new.timer()
    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)
    all.ad <- get.ad.array(vcf)  # 3D array with 2 reads as the third dimension
    parent.indices <- which(colnames(all.ad) %in% parents)  # find parents
    ad <- all.ad[ , -parent.indices, ]  # remove parents

    rm(all.ad)
    listapply <- get.lapply(parallel, cores)


    hom.ref.read.emm <- (1-rerr)^all.ad[ , , 1] * (rerr)^all.ad[ , , 2]
    hom.alt.read.emm <- (1-rerr)^all.ad[ , , 2] * (rerr)^all.ad[ , , 1]
    het.read.emm     <-    (0.5)^all.ad[ , , 1] *  (0.5)^all.ad[ , , 2]

    rm(ad)

    rpgs <- abind(  # read probs given states
        hom.ref.read.emm,
        het.read.emm,
        het.read.emm,
        hom.alt.read.emm,

        along=3
    )

    rm(hom.ref.read.emm, hom.alt.read.emm, het.read.emm)

    parental.rpgs <- rpgs[ , parent.indices, ]
    rpgs <- rpgs[ , -parent.indices, ]

    get.trans.probs <- function(i) {
        j <- i+1

        ## snp.i will be constant along the third dimension and snp.j will be
        ## constant along the second dimension. By binding additional copies in
        ## this way we can replicate mathematical matrix multiplication with
        ## element-by-element multiplication allowing us to do all SNPs at once
        ## without needing an apply function which can be slower
        snp.i <- abind(rpgs[i, , ], rpgs[i, , ], rpgs[i, , ], rpgs[i, , ], along=3)
        snp.j <- abind(rpgs[j, , ], rpgs[j, , ], rpgs[j, , ], rpgs[j, , ], along=3)
        snp.j <- aperm(snp.j, c(1,3,2))

        rpgsp <- snp.i * snp.j  # read probs given state pair
        rm(snp.i, snp.j)

        get.objective.fun <- function(f) {
            function(r) {
                per.snp <- apply(rpgsp, 1, function(layer) {
                    sum(layer * f(r))
                })
                sum(log(per.snp))
            }
        }

        ## message("Starting double loop for ", i)
        log.liklihood.mat <- matrix(-Inf, nrow=16, ncol=16)
        recomb.val.mat <- matrix(0, nrow=16, ncol=16)
        for (x in 2:15) {      # 1 and 16 are not biallelic parental states and there
            for (y in 2:15) {  # is not optimal recombination probability
                site.pair.probs.fun <- site.pair.transition.probs[[x]][[y]]
                obj.fun <- get.objective.fun(site.pair.probs.fun)
                ## browser()

                init <- 0.1
                result <- optim(par=init,
                                obj.fun,
                                method="Brent",
                                lower=0,
                                upper=0.5,
                                control=list(ndeps=1e-2,  # step size
                                             fnscale=-1))
                recomb.val.mat[x,y] <- result$par
                log.liklihood.mat[x,y] <- result$value
                ## message(x, ",", y)
            }
        }

        list(logliklihoods = log.liklihood.mat,
             recombs = recomb.val.mat)
    }

    ## print("ready")
    ## browser()
    ## print("done")




    ## Progress bar code
    ## -------------------------------------------------------------------------
    progress.env <- new.env()
    thefifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prog.env <- progress.env
    n.jobs <- length(chroms) - length(u.chroms)
    ## -------------------------------------------------------------------------

    ## -------------------------------------------------------------------------
    writeBin(0, thefifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
    }  # else the forked process handles this
    ## -------------------------------------------------------------------------
    ## Progress bar code


    ret.val <- lapply(u.chroms, function(chrom) {
        indices <- which(chroms == chrom)  # which indices correspond with this chrom
        indices <- indices[-length(indices)]  # remove the last element

        ret.val.2  <- listapply(indices, function(index) {

            ret.val.3 <- get.trans.probs(index)

            ## Progress bar code
            ## -----------------------------------------------------------------
            writeBin(1/n.jobs, thefifo)  # update the progress bar info
            if (!parallel) {  # if running in serial mode
                prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
            }  # else the forked process handles this
            ## -----------------------------------------------------------------
            ## Progress bar code

            ret.val.3
        })

        ret.val.2

        ## do.call(abind3, ret.val.2)  # bind transmissions in 3rd dim and return
    })


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code


    names(ret.val) <- u.chroms
    transmissions <- ret.val



    ## ## Construct emission liklihood matrix. There are 16 possible parental states at each SNP
    ## p <- 0.99  # probability of a site in parents being heterozygous
    ## log.penalty <- log(c(p^2, p*(1-p), p*(1-p), p^2,
    ##                      p*(1-p), (1-p)^2, (1-p)^2, p*(1-p),
    ##                      p*(1-p), (1-p)^2, (1-p)^2, p*(1-p),
    ##                      p^2, p*(1-p), p*(1-p), p^2))

    emissions <- apply(parental.rpgs, 1, function(snp.depths) {
        result <- sapply(1:16, function(i) {
            i <- i-1            # deal with base 1 indexing
            p1 <- floor(i / 4)  # bit shift i right twice to get two most sig bits
            p2 <- i - p1*4      # two least significant bits
            log(snp.depths[1, (p1+1)]) +
                log(snp.depths[2, (p2+1)]) +
                log.penalty[i+1]
        })
        names(result) <- NULL
        result  # implicit return
        ## browser()
        ## print(round(result,2))
        ## print(order(result, decreasing=TRUE))
    })
    ## print(timer())
    ## browser()
    ## print(timer())

    ## TODO(Jason): verify that at least 2 sites are in the chromosome before runnin viterbi
    parental.models <- listapply(u.chroms, function(chrom) {
        viterbi(emissions[, chroms == chrom],
                do.call(abind3, lapply(transitions[[chrom]], function(elem) {elem$logliklihoods})),
                emm.log = TRUE,
                trans.log = TRUE)
    })

    recombs <- listapply(u.chroms, function(chrom) {
        parental.model <- parental.models[[chrom]]
        recomb.matrices <- lapply(transitions[[chrom]], function(elem) {elem$recombs})
        if (length(parental.model) != length(recomb.matrices) + 1)
            stop("This should never happen. Error in model lengths")

        sapply(seq_along(recomb.matrices), function(i) {
            recomb.matrices[parental.model[i], parental.model[i+1]]
        })
    })

    list(parental.models <- parental.models,
         recombs <- recombs)  # implicit return
}


log.lik.path <- function(states, em, tr) {
    if (length(states) != ncol(em))
        stop("length of states must match number of cols of em")
    if (length(states) != dim(tr)[3] + 1)
        stop("length of states must match number of cols of tr + 1")

    em.log.liks <- sum(sapply(seq_len(ncol(em)), function(i) {
        em[states[i], i]
    }))

    tr.log.liks <- sum(sapply(seq_len(dim(tr)[3]), function(i) {
        tr[states[i], states[i+1], i]
    }))

    log(1/16) + em.log.liks + tr.log.liks
    ## em[states[1],1] +
    ##     em[states[2],2] +
    ##     em[states[3],3] +
    ##     em[states[4],4] +
    ##     em[states[5],5] +
    ##     tr[states[1], states[2], 1] +
    ##     tr[states[2], states[3], 2] +
    ##     tr[states[3], states[4], 3] +
    ##     tr[states[4], states[5], 4]
}
