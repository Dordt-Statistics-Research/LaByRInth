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


LabyrinthImpute <- function (vcf, parents, generation, out.file, parental.imputation,
                             use.fwd.bkwd=FALSE, use.viterbi=TRUE,
                             calc.posteriors=TRUE, geno.err=0.05,
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

    if (cores == 1 && parallel) {
        display(0, "Cores cannot be 1 if parallel is true\n")
        stop()
    }

    if (cores > 1 && !parallel) {
        display(0, "Cores cannot be greater than 1 if parallel is false\n")
        stop()
    }

    if (cores < 1) {
        display(0, "Cores cannot be less than 1\n")
        stop()
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
    ## non.hom.poly <- which( ! parent.hom.and.poly(vcf, parents))
    ## if (length(non.hom.poly) != 0) {
    ##     display(0, "The parents are not homozygous within and ",
    ##             "polymorphic between at the following sites:")
    ##     chroms <- getCHROM(vcf)[non.hom.poly]
    ##     positions <- getPOS(vcf)[non.hom.poly]
    ##     ad <- getAD(vcf)[non.hom.poly, parents]
    ##     for (i in seq_along(chroms)) {
    ##         display(1, "\tCHR:", chroms[i],
    ##                 "\tPOS:", positions[i],
    ##                 paste0("\t", parents[1], ":"), ad[i, 1],
    ##                 paste0("\t", parents[2], ":"), ad[i, 2])
    ##     }
    ## }

    if (length(non.biallelic) != 0) {
        display(1, "The vcf must be filtered before imputing; please run LabyrinthFilter first\n")
        stop()
    }
    display(1, "Completed in ", timer(), "\n")




    ## Variable initialization and function loading
    timer <- new.timer()
    display(0, "Loading data specific to F", generation)
    ## site.pair.transition.probs variable will be loaded
    source(paste0("./new-transition-probs/F", generation, ".R"))
    display(1, "Completed in ", timer(), "\n")


    timer <- new.timer()
    display(0, "Storing markers, samples, and read depths", generation)
    snp.chroms <- getCHROM(vcf)
    u.chroms <- unique(snp.chroms)
    ad <- get.ad.array(vcf)
    sample.names <- getSAMPLES(vcf)
    marker.names <- getID(vcf)
    parental.results <- readRDS(parental.imputation)
    display(1, "Completed in ", timer(), "\n")


    timer <- new.timer()
    display(0, "Estimating erroneous read rate from heterozygosity")
    read.err <- estimate.read.err(ad, generation)
    display(1, "Error rate estimated to be ", round(read.err*100, 2), "%")
    display(1, "Completed in ", timer(), "\n")


    ## emission probability code
    timer <- new.timer()
    display(0, "Generating emission probabilities")
    emission.structure <- get.emission.structures(ad,
                                                  sample.names,
                                                  marker.names,
                                                  read.err)
    display(1, "Completed in ", timer(), "\n")




    ## ## parental imputation and recombination rate estimates
    ## timer <- new.timer()
    ## display(0, "Imputing parents and estimating recombination rates")
    ## ## parental.results <- determine.parents.and.recombs(emission.structure,
    ## ##                                                   parents,
    ## ##                                                   snp.chroms,
    ## ##                                                   marker.names,
    ## ##                                                   sample.names,
    ## ##                                                   parallel,
    ## ##                                                   cores)
    ## ## load("testing/temp_storage_6.RData")
    ## parental.results <- readRDS("parental-results-est-careful.rds")
    ## saveRDS(parental.results, "parental-results-est-careful.rds")
    ## display(1, "Completed in ", timer(), "\n")




    ## transition probability code
    timer <- new.timer()
    display(0, "Generating transition probabilities")
    transition.structures <- get.transition.structures(snp.chroms,
                                                       marker.names,
                                                       parental.results,
                                                       parallel,
                                                       cores)
    display(1, "Completed in ", timer(), "\n")




    ## imputation code
    timer <- new.timer()
    display(0, "Imputing missing sites")
    imp.res <- impute(parents,
                      emission.structure,
                      transition.structures,
                      parental.results,
                      sample.names,
                      snp.chroms,
                      parallel,
                      cores,
                      use.fwd.bkwd,
                      calc.posteriors)
    display(1, "Completed in ", timer(), "\n")

    ## new vcf creation code
    timer <- new.timer()
    display(0, "Creating new vcf with imputed data")
    vcf <- update.vcf(vcf, imp.res$gt, imp.res$posteriors)
    write.vcf(vcf, paste0(out.file, ".vcf.gz"))
    display(1, "Completed in ", timer(), "\n")


    display(0, "LaByRInth imputation completed in ", total.timer())
    invisible(vcf) ## implicit return
}


## remove all sites that are not homozygous within and polymorphic between for
## the parents and remove all sites that are not biallelic
LabyrinthFilter <- function(vcf, parents, out.file, hom.poly=FALSE) {
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
    ## display(0, "Checking for sites where parents are not homozygous within and polymorphic between")
    ## non.hom.poly <- ! parent.hom.and.poly(vcf, parents)
    ## if (length(non.hom.poly) != 0) {
    ##     display(0, "The parents are not homozygous within and polymorphic between at the following sites which will be removed:")
    ##     chroms <- getCHROM(vcf)[which(non.hom.poly)]
    ##     positions <- getPOS(vcf)[which(non.hom.poly)]
    ##     ad <- getAD(vcf)[which(non.hom.poly), parents]
    ##     for (i in seq_along(chroms)) {
    ##         display(0, "\tCHR:", chroms[i],
    ##                 "\tPOS:", positions[i],
    ##                 paste0("\t", parents[1], ":"), ad[i, 1],
    ##                 paste0("\t", parents[2], ":"), ad[i, 2])
    ##     }
    ## }

    mask <- !non.biallelic## & !non.hom.poly
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

    state.names <- rownames(emm)
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
    rownames(res) <- state.names
    res
}


## If emm.log is TRUE, then the data passed in emm should already be
## log-scaled. Similarly with trans.log
viterbi <- function(emm, trans, emm.log=FALSE, trans.log=FALSE) {
    if (!emm.log)
        emm <- log(emm)
    if (!trans.log)
        trans <- log(trans)

    state.names <- rownames(emm)
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
get.emission.structures <- function(ad, sample.names, marker.names, read.err, geno.err) {

    k <- choose(ad[ , , 1] + ad[ , , 2], ad[ , , 1])  # n choose k

    hom.ref.read.emm <- k * (1-read.err)^ad[ , , 1] * (read.err)^ad[ , , 2]
    hom.alt.read.emm <- k * (1-read.err)^ad[ , , 2] * (read.err)^ad[ , , 1]
    het.read.emm     <- k *        (0.5)^ad[ , , 1] *      (0.5)^ad[ , , 2]

    ## read probs given states
    ## first dimension is the markers/loci/SNPs
    ## second dimension is the members of the population
    ## third dimension is the state (hom ref, het type I, het type II, hom alt)
    rpgs <- abind(
        hom.ref.read.emm,
        het.read.emm,
        het.read.emm,
        hom.alt.read.emm,

        along=3
    )

    dimnames(rpgs) <- list(marker.names,
                           sample.names,
                           c("hom.ref", "het.I", "het.II", "hom.alt"))

    rpgs  # implicit return
}


get.transition.structures <- function(snp.chroms, marker.names,
                                      parental.results, parallel, cores) {
    u.chroms <- unique(snp.chroms)
    listapply <- get.lapply(parallel, cores)

    trans.probs <- function(chrom) {
        parental.model <- parental.results$parental.models[[chrom]]$path
        recombs <- parental.results$recombs[[chrom]]  # vector of functions

        ## site.pair.transition.probs was loaded when a file in
        ## new-transition-probs was sourced
        trans.matrices <- lapply(seq_along(recombs), function(i) {
            ## for marker i and i+1, get the transition matrix based on the
            ## parental states at both sites. This "matrix" is actually a
            ## function of the recombination probability which can be called
            ## with recombs[i] to get a 4x4 matrix
            mat <- site.pair.transition.probs[[parental.model[i]]][[parental.model[i+1]]](recombs[i])
            ## normalize so that every rowSum is 0 since these are probabilities
            mat / rep(rowSums(mat), 4)
        })

        result <- do.call(abind3, trans.matrices)
        result[result < 0] <- 0  # handle numerical errors

        marker.names <- marker.names[snp.chroms == chrom]
        trans.names <- paste0(marker.names[-length(marker.names)],
                              ".",
                              marker.names[-1])

        dimnames(result) <- list(c("hom.ref", "het.I", "het.II", "hom.alt"),
                                 c("hom.ref", "het.I", "het.II", "hom.alt"),
                                 trans.names)
        result
    }


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


## The following labeling schemes are used. Both the Viterbi algorithm and the
## forward backward algorithm will keep track of the following four state types
## when determining the best path
##
## hom.ref     het.I     het.II     hom.alt
## 0|0         0|1       1|0        1|1
##
## When the Viterbi algorithm is used, a single opimal sequence of states is
## determined which thus includes haplotypic information and can thus
## distinguish betwee het.I and het.II. Thus the calls are phased and the
## vertical bar is used in the VCF file as shown above. If the forward backward
## algorithm is used, then it makes more sense to combine the het.I and het.II
## types because the VCF format does not allow for specifying posterior
## probabilities for both het.I and het.II. Then the following scheme is used:
##
## hom.ref     het     hom.alt
## 0/0         0/1     1/1
impute <- function(parents, emm.structures, trans.structures, parental.results,
                   sample.names, snp.chroms, parallel, cores, use.fwd.bkwd,
                   calc.posteriors) {

    listapply <- get.lapply(parallel, cores)
    u.chroms <- unique(snp.chroms)

    ## Given two probabilites of odd recombinations (r1 between sites x and y;
    ## r2 between sites y and z) what is the probability of an odd number of
    ## recombinations between sites x and z
    combine.recombs <- function(r1, r2) {r1*(1-r2) + r2*(1-r1)}

    impute.sample.chrom <- function(sample, chrom) {

        n.sites <- sum(snp.chroms==chrom)  # boolean addition
        parent.paths <- parental.model.to.parents(
            parental.results$parental.models[[chrom]]$path
        )

        emm <- t(emm.structures[chrom==snp.chroms, sample, ])
        trans <- trans.structures[[chrom]]

        res <- list()
        res$posteriors <- NULL

        if (use.fwd.bkwd || calc.posteriors) {
            if (sample %in% parents) {
                posterior.mat <- rbind("hom.ref" = rep(1/3, n.sites),
                                       "het.I"   = rep(1/6, n.sites),
                                       "het.II"  = rep(1/6, n.sites),
                                       "hom.alt" = rep(1/3, n.sites))
            } else {
                posterior.mat <- fwd.bkwd(emm, trans)
            }

            posteriors <- rbind("hom.ref" = posterior.mat["hom.ref", ],
                                "het"     = posterior.mat["het.I", ] + posterior.mat["het.II", ],
                                "hom.alt" = posterior.mat["hom.alt", ])

            phred.scaled <- apply(posteriors, 1:2, function(x) {
                ## 1-x is the probability the call is wrong
                ## use min to prevent infinity from getting through
                min(-10*log((1-x), base=10), 100)
            })

            res$posteriors <- apply(phred.scaled, 2, paste0, collapse=",")

            if (use.fwd.bkwd) {
                if (sample == parents[1]) {
                    res$gt <- c("0|0", "0|1", "1|0", "1|1")[parent.paths[[1]]]
                } else if (sample == parents[2]) {
                    res$gt <- c("0|0", "0|1", "1|0", "1|1")[parent.paths[[2]]]
                } else {
                    res$gt <- apply(posteriors, 2, function(states) {
                        c("0/0","0/1","1/1")[which.max(states)]
                    })
                }
            }
        }

        if (!use.fwd.bkwd) {
            if (sample == parents[1]) {
                best.path <- parent.paths[[1]]
            } else if (sample == parents[2]) {
                best.path <- parent.paths[[2]]
            } else {
                relevant <- ! apply(emm, 2, function(state.probs) {
                    all(state.probs == state.probs[1])
                })

                which.relevant <- which(relevant)
                names(which.relevant) <- NULL

                if (length(which.relevant) < 2) {
                    best.path <- rep(5, length(relevant))
                } else {

                    path.and.prob <- viterbi(emm, trans)
                    best.path <- path.and.prob$path

                    inter.relevant.probs <- sapply.pairs(which.relevant, function(a, b) {exp(log.lik.path(best.path, emm, trans, a, b))})
                    names(inter.relevant.probs) <- NULL

                    THRESHOLD <- 1e-5
                    for (unlikely.index in which(inter.relevant.probs < THRESHOLD)) {
                        start <- which.relevant[unlikely.index] + 1
                        end <- which.relevant[unlikely.index + 1] - 1

                        best.path[start:end] <- 5
                    }

                    best.path
                }
            }

            res$gt <- c("0|0", "0|1", "1|0", "1|1", "./.")[best.path]
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
    n.jobs <- length(u.chroms) * length(sample.names)
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

        imputed.samples <- listapply(sample.names, function(sample) {

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
update.vcf <- function(vcf, gt, posteriors=NULL) {
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

    gp <- posteriors
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


ad.to.gt <- function(str.mat) {
    
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
    writeLines("|     ======  L O W - C O V E R A G E   B I A L L E L I C  ======     |")
    writeLines("|     =======   R - P A C K A G E   I M P U T A T I O N   =======     |")
    ## writeLines("|                                                                     |")
    ## writeLines("|         Copyright 2017 Jason Vander Woude and Nathan Ryder          |")
    ## writeLines("|           Licensed under the Apache License, Version 2.0            |")
    writeLines("|_____________________________________________________________________|")
    writeLines("")

    ## ## the image looks funny because '\' in the displayed image must be '\\' in the code
    ## writeLines("")
    ## writeLines(" _____________________________________________________________________")
    ## writeLines("|          __          ____        ____  _____                        |")
    ## writeLines("|         / /         / __ \\      / __ \\/_  _/                        |")
    ## writeLines("|        / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __      |")
    ## writeLines("|       / /   / _  \\/ _  _/\\ \\/ / _  _/ / / /  | / /_  __/ /_/ /      |")
    ## writeLines("|      / /___/ /_/ / /_\\ \\  \\  / / \\ \\_/ /_/ /||/ / / / / __  /       |")
    ## writeLines("|     /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/        |")
    ## writeLines("|                                                                     |")
    ## writeLines("| LaByRInth: Low-coverage Biallelic R Imputation                      |")
    ## writeLines("| Copyright 2017 Jason Vander Woude and Nathan Ryder                  |")
    ## writeLines("| Licensed under the Apache License, Version 2.0                      |")
    ## writeLines("| Source code: github.com/Dordt-Statistics-Research/LaByRInth         |")
    ## writeLines("| Based on LB-Impute: github.com/dellaporta-laboratory/LB-Impute      |")
    ## writeLines("| Funding received from the National Science Foundation (IOS-1238187) |")
    ## writeLines("|_____________________________________________________________________|")
    ## writeLines("")

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


determine.parents.and.recombs <- function(emm.structure, parents, snp.chroms,
                                          marker.names, sample.names, parallel,
                                          cores) {

    listapply <- get.lapply(parallel, cores)

    u.chroms <- unique(snp.chroms)
    which.parents <- sample.names %in% parents
    parental.rpgs <- emm.structure[ , which.parents, ]  # read probs given states
    rpgs <- emm.structure[ , !which.parents, ]


    ## Progress bar code
    ## -------------------------------------------------------------------------
    progress.env <- new.env()
    thefifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prog.env <- progress.env
    n.jobs <- (length(snp.chroms) - length(u.chroms))*14^2
    ## -------------------------------------------------------------------------

    ## -------------------------------------------------------------------------
    writeBin(0, thefifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
    }  # else the forked process handles this
    ## -------------------------------------------------------------------------
    ## Progress bar code



    transitions <- lapply(u.chroms, function(chrom) {
        indices <- which(snp.chroms == chrom)  # which indices correspond with this chrom
        indices <- indices[-length(indices)]  # remove the last element

        ret.val.2  <- listapply(indices, function(i) {

            ## ret.val.3 <- get.trans.probs(index)

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
            essential.layers <- rowSums(rpgs[i, , ]) != 4 & rowSums(rpgs[j, , ]) != 4

            dimnames(rpgsp) <- list(sample.names[ !which.parents ],
                                    c("hom.ref", "het.I", "het.II", "hom.alt"),
                                    c("hom.ref", "het.I", "het.II", "hom.alt"))


            rm(snp.i, snp.j)

            get.objective.fun <- function(f) {
                const.per.snp <- apply(rpgsp[!essential.layers,,], 1, function(layer) {
                    sum(layer * f(0))
                })
                k <- sum(log(const.per.snp))

                function(r) {
                    f.of.r <- f(r)
                    per.snp <- apply(rpgsp[essential.layers,,], 1, function(layer) {
                        sum(layer * f.of.r)
                    })
                    sum(log(per.snp)) + k  # implicit return from objective.fun
                }  # get.objective.fun returns this function
            }

            log.liklihood.mat <- matrix(-Inf, nrow=16, ncol=16)
            recomb.val.mat <- matrix(0, nrow=16, ncol=16)
            for (x in 2:15) {      # 1 and 16 are not biallelic parental states and there
                for (y in 2:15) {  # is not optimal recombination probability
                    site.pair.probs.fun <- site.pair.transition.probs[[x]][[y]]
                    obj.fun <- get.objective.fun(site.pair.probs.fun)

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

                    ## Progress bar code
                    ## -----------------------------------------------------------------
                    writeBin(1/n.jobs, thefifo)  # update the progress bar info
                    if (!parallel) {  # if running in serial mode
                        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
                    }  # else the forked process handles this
                    ## -----------------------------------------------------------------
                    ## Progress bar code
                }
            }

            mymat <- log(exp(log.liklihood.mat) / rep(rowSums(exp(log.liklihood.mat)), 16))
            mymat[is.nan(mymat)] <- -Inf
            log.liklihood.mat <- mymat

            list(logliklihoods = log.liklihood.mat,
                 recombs = recomb.val.mat)

            ## ret.val.3
        })

        ret.val.2

        ## do.call(abind3, ret.val.2)  # bind transmissions in 3rd dim and return
    })


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code

    names(transitions) <- u.chroms




    ## ## Construct emission liklihood matrix. There are 16 possible parental states at each SNP
    p <- 0.99  # probability of a site in parents being heterozygous
    log.penalty <- log(c(p^2, p*(1-p), p*(1-p), p^2,
                         p*(1-p), (1-p)^2, (1-p)^2, p*(1-p),
                         p*(1-p), (1-p)^2, (1-p)^2, p*(1-p),
                         p^2, p*(1-p), p*(1-p), p^2))

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

    ## TODO(Jason): verify that at least 2 sites are in the chromosome before
    ## runnin viterbi







    ## ###### REMOVE ######
    ## transitions <- readRDS("testing/transitions_1.rds")

    parental.models <- listapply(u.chroms, function(chrom) {
        viterbi(emissions[, snp.chroms == chrom],
                do.call(abind3, lapply(transitions[[chrom]], function(elem) {elem$logliklihoods})),
                emm.log = TRUE,
                trans.log = TRUE)
    })
    names(parental.models) <- u.chroms

    recombs <- listapply(u.chroms, function(chrom) {
        parental.model <- parental.models[[chrom]]$path
        recomb.matrices <- lapply(transitions[[chrom]], function(elem) {elem$recombs})
        if (length(parental.model) != length(recomb.matrices) + 1)
            stop("This should never happen. Error in model lengths")

        sapply(seq_along(recomb.matrices), function(i) {
            recomb.matrices[[i]][parental.model[i], parental.model[i+1]]
        })
    })
    names(recombs) <- u.chroms

    list(parental.models = parental.models,
         recombs = recombs)  # implicit return
}


## parental.model is a vector of values in range 1-16 inclusive. This indicates
## the state of the two parents for each site. The value is 4 times the value of
## the state of parent 1 plus the state of parent 2. This function returns a
## list where the first element is the vector of states of parent 1 and the
## second is the vector of states of parent 2
parental.model.to.parents <- function(parental.model) {
    parents <- sapply(parental.model, function(call) {
        call <- call-1            # deal with base 1 indexing
        p1 <- floor(call / 4)  # bit shift i right twice to get two most sig bits
        p2 <- call - p1*4      # two least significant bits
        c(p1+1, p2+1)
    })

    list(
        parents[1, ],
        parents[2, ]
    )
}


LabyrinthAnalyze <- function(orig, masked, imputed) {
    ## timer <- new.timer()
    ## display(0, "Reading original file")
    ## display(0, "Reading masked file")
    ## display(0, "Reading filtered file")
    gt.o <- getGT(orig)
    gt.m <- getGT(masked)
    gt.i <- getGT(imputed)

    ## vector.indices
    masked.sites <-
        (gt.o != "./." & gt.o != ".|.") &
        (gt.m == "./." | gt.m == ".|.")

    same <-
        ((gt.o == "0/0" | gt.o == "0|0") & (gt.i == "0/0" | gt.i == "0|0")) |
        ((gt.o == "0/1" | gt.o == "0|1" | gt.o == "1/0" | gt.o == "1|0") &
         (gt.i == "0/1" | gt.i == "0|1" | gt.i == "1/0" | gt.i == "1|0")) |
        ((gt.o == "1/1" | gt.o == "1|1") & (gt.i == "1/1" | gt.i == "1|1"))

    n.same <- sum(masked.sites & same)
    n.masked <- sum(masked.sites)
    accuracy <-  n.same / n.masked

    w <- which(masked.sites & !same, arr.ind = TRUE)
    ad <- getAD(orig)

    browser()
    display(0, "Accuracy: ", n.same, "/", n.masked, " = ", accuracy)
}


estimate.read.err <- function(ad, n.gen) {
    depths <- ad[ , , 1] + ad[ , , 2]
    relevant <- depths > 1  # depth of 0 or 1 gives no info on heterozygosity
    ref <- ad[ , , 1][relevant]
    alt <- ad[ , , 2][relevant]

    h <- 0.5^(n.gen - 1)  # expected proportion of sites that are heterozygous

    ## objective function of r, the probability of an erroneous read at the
    ## haplotype level
    obj.fun <- function(r) {
        sum(log((1-h)/2 * (r^ref * (1-r)^alt + r^alt * (1-r)^ref) + h*0.5^(ref + alt)))
            ## mapply(
            ##     function(n.ref, n.alt) {
            ##        (1-h)/2 * (r^n.ref * (1-r)^n.alt + r^n.alt * (1-r)^n.ref) +
            ##             h*0.5^(n.ref + n.alt)
            ##     }
            ## , refs, alts)
        ## ))
    }

    init <- 0.01
    result <- optim(par=init,
                    obj.fun,
                    method="Brent",
                    lower=0,
                    upper=0.5,
                    control=list(ndeps=1e-3,  # step size
                                 fnscale=-1))

    result$par
}


sapply.pairs <- function(vec, fun) {
    if (length(vec) < 2) {
        stop("first argument must have length >= 2")
    }
    v1 <- vec[-length(vec)]  # remove last entry
    v2 <- vec[-1]  # remove first entry

    mapply(fun, v1, v2)
}


log.lik.path <- function(states, em, tr, index1, index2) {
    if (length(states) != ncol(em))
        stop("length of states must match number of cols of em")
    if (length(states) != dim(tr)[3] + 1)
        stop("length of states must match number of cols of tr + 1")

    em <- log(em[ , index1:index2, drop=FALSE])
    tr <- log(tr[ , , index1:(index2 - 1), drop=FALSE])
    states <- states[index1:index2]

    em.log.liks <- sum(sapply(seq_len(ncol(em)), function(i) {
        em[states[i], i]
    }))

    tr.log.liks <- sum(sapply(seq_len(dim(tr)[3]), function(i) {
        tr[states[i], states[i+1], i]
    }))

    em.log.liks + tr.log.liks
}
