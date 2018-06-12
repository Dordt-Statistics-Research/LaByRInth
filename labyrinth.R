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
                             read.err=0.05, geno.err=0.05,
                             recomb.dist=1e6, parallel=TRUE,
                             cores=4, quiet=FALSE) {

    ## begin timer
    total.timer <- new.timer()
    print.labyrinth.header()


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


    ## transition probability code
    trans.file <- paste0("./transition-structures-f", generation, ".rds")
    ## if (file.exists(trans.file)) {
    ##     timer <- new.timer()
    ##     display(0, "Loading transition probabilities")
    ##     transition.structures <- readRDS(trans.file)
    ##     display(1, "Completed in ", timer(), "\n")
    ## } else {
        timer <- new.timer()
        display(0, "Generating transition probabilities")
        transition.structures <- get.transition.structures(vcf, generation, recomb.dist, parallel, cores)
        display(1, "Completed in ", timer(), "\n")
        transition.structures <<- transition.structures  # for debugging
    ##     saveRDS(transition.structures, trans.file)  # for debugging
    ## }

    ## emission probability code
    emm.file <- paste0("./emission-structures-f", generation, ".rds")
    ## if (file.exists(emm.file)) {
    ##     timer <- new.timer()
    ##     display(0, "Loading emission probabilities")
    ##     emission.structures <- readRDS(emm.file)
    ##     display(1, "Completed in ", timer(), "\n")
    ## } else {
        timer <- new.timer()
        display(0, "Generating emission probabilities")
        emission.structures <- get.emission.structures(vcf, parents, read.err, geno.err, generation, parallel, cores)
        display(1, "Completed in ", timer(), "\n")
        emission.structures <<- emission.structures  # for debugging
    ##     saveRDS(emission.structures, emm.file)  # for debugging
    ## }


    ## imputation code
    timer <- new.timer()
    display(0, "Imputing missing sites")
    imp.res <- impute(vcf, parents, emission.structures, transition.structures,
                      parallel, cores)
    display(1, "Completed in ", timer(), "\n")

    ## new vcf creation code
    timer <- new.timer()
    display(0, "Creating new vcf with imputed data")
    vcf <- update.vcf(vcf, imp.res)
    display(1, "Completed in ", timer(), "\n")

    write.vcf(vcf, out.file)
    display(0, "LaByRInth imputation completed in ", total.timer())

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
        stop("Lenght of parent must be 2")
    }
    par.mat <- getAD(vcf)[ , parents]
    apply(par.mat, 1, function(row) {
        (n.zeros(ad.to.num(row[1])) > 0
            && n.zeros(ad.to.num(row[2])) > 0
            && ! any(nonzero(ad.to.num(row[1]))
                     & nonzero(ad.to.num(row[2]))))
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

    ## forward probabilities
    f.probs <- matrix(data=0, nrow=n.states, ncol=n.sites)
    b.probs <- f.probs

    start.probs <- rep(1/n.states, n.states)
    f.probs[, 1] <- start.probs * emm[, 1]
    f.probs <- normalize(f.probs, 1)

    for (site in 2:n.sites) {
        t.index <- site - 1
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


## returns the emission probabilities for a given variant and chromosome. Note
## that the structure of the returned data is different than a traditional
## emission probability matrix. The first dimension of the data corresponds to
## the chromosome, and the second dimension corresponds to the sample. reads
## should be a vector of strings where each string is a numeric value followed
## by a comma and another numeric value. This is the format that the vcfR
## package stores the reads in.
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
        names(reads) <- NULL  # makes debuggin easier
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


get.transition.structures <- function(vcf, generation, recomb.dist, parallel, cores) {

    states <- 1:4
    names(states) <- c("P1", "H1", "H2", "P2")

    listapply <- get.lapply(parallel, cores)

    source(paste0("./new-transition-probs/F", generation, ".R"))
    all.model.trans <- do.call(paste0("get.all.model.trans.F", generation), list())

    phys.recomb.prob <- function(dist, recomb.dist) {
        ## this is twice the value that LB-Impute used
        (1 - exp(-1.0 * dist / recomb.dist))
    }


    trans.probs <- function(chrom, model) {
        positions <- getPOS(vcf)[getCHROM(vcf) == chrom]
        distances <- diff(positions)  # difference b/w successive positions

        phys.recomb.probs <- lapply(distances, phys.recomb.prob, recomb.dist)

        trans.matrices <- lapply(phys.recomb.probs, function(phys.r.prob) {
            ## each element of list all.model.trans is a function of
            ## phsyical recombination distance
            all.model.trans[[model]](phys.r.prob)
        })

        ## helper code for running do.call(abind, ...)
        abind.args <- c(
            list(along = 3),
            trans.matrices
        )

        do.call(abind, abind.args)
    }

    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)
    n.selfings <- generation - 1
    models <- 1:(4^n.selfings)


    ## Progress bar code
    ## -------------------------------------------------------------------------
    progress.env <- new.env()
    thefifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prog.env <- progress.env
    n.jobs <- length(u.chroms) * length(models)
    ## -------------------------------------------------------------------------

    ## -------------------------------------------------------------------------
    writeBin(0, thefifo)  # update the progress bar info
    if (!parallel) {  # if running in serial mode
        prog.env$progress <- PrintProgress(thefifo, prog.env$progress)
    }  # else the forked process handles this
    ## -------------------------------------------------------------------------
    ## Progress bar code

    ret.val <- lapply(u.chroms, function(chrom) {
        ret.val.2  <- listapply(models, function(model) {

            ret.val.3 <- trans.probs(chrom, model)


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
        names(ret.val.2) <- NULL  # just to be clear the models aren't named
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


impute <- function(vcf, parents, emm.structures, trans.structures, parallel, cores) {

    listapply <- get.lapply(parallel, cores)

    p1.is.ref <- is.parent.ref(vcf, parents[1])

    abind1 <- function(...) {
        abind(..., along=1)
    }

    abind2 <- function(...) {
        abind(..., along=2)
    }

    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)
    samples <- getSAMPLES(vcf)

    impute.sample.chrom <- function(sample, chrom) {

        n.sites <- sum(chroms==chrom)  # boolean addition
        p1.ref <- p1.is.ref[chroms==chrom]

        if (sample == parents[1]) {
            normalized <- rbind(rep(1, n.sites),
                                0,
                                0,
                                0)
        } else if (sample == parents[2]) {
            normalized <- rbind(0,
                                0,
                                0,
                                rep(1, n.sites))
        } else {
            emm <- emm.structures[[chrom]][[sample]]
            model.trans <- trans.structures[[chrom]]
            n.models <- length(model.trans)

            list.of.posterior.matrices <-
                lapply(model.trans, function(trans) {
                    fwd.bkwd(emm, trans)
                })

            summed.posteriors <- Reduce(`+`, list.of.posterior.matrices)
            normalized <- 1/n.models * summed.posteriors

        }

        ref <- ifelse(p1.ref, normalized[1, ], normalized[4, ])
        alt <- ifelse(p1.ref, normalized[4, ], normalized[1, ])
        het <- normalized[2, ] + normalized[3, ]

        array(c(ref,                     # homozygous reference allele prob
                alt,                     # homozygous alternate allele prob
                het),                    # heterozygous prob
              dim=c(n.sites, 1, 3))
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

        do.call(abind2, imputed.samples)  # implicit return

    })

    ## combine columns of each imputed sample
    ret.val <- do.call(abind1, imputed.chroms)


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


get.trans <- function(data, chrom, model) {
    dir <- "./.saved_trans"
    readRDS(data, paste0(dir, "/", paste0(model, collapse=""), ".rds"))
}


is.parent.ref <- function(vcf, parent) {
    str.ad <- getAD(vcf)[ , parent]
    sapply(str.ad, function(str) {
        ## split the string and check if reference read is nonzero
        ad.to.num(str)[1] != 0
    })
}


## modify the vcf in place to add imputed calls in gt
update.vcf <- function(vcf, impute.res) {
    rounded.res <- apply(impute.res, 1:3, round, digits=2)
    rounded.res <- apply(rounded.res, 1:3, `*`, 100)

    gp <- apply(rounded.res, 1:2, paste0, collapse=",")

    gt <- apply(impute.res, 1:2, function(states) {
        c("0/0","1/1","0/1")[which.max(states)]
    })

    ad <- getAD(vcf)

    n.rows <- nrow(ad)
    if (n.rows != nrow(gt))
        stop("conflict in number of rows; should never happen")

    concat <- matrix(paste0(gt, ":", ad, ":", gp), nrow=n.rows)

    vcf@gt <- cbind("GT:AD:GP", concat)
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
    gt[is.na(gt)] <- "0/0"  # replace NA entries
    gt
}


getGP <- function(vcf) {
    extract.gt(vcf, "GP")
}


## convert string representation to numeric vector
ad.to.num <- function(str) {
    as.numeric(str.split(str, ","))

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
    writeLines("|          __          ____        _____  _____                       |")
    writeLines("|         / /         / __ \\      / ___ \\/_  _/                       |")
    writeLines("|        / /   ____  / /_/ /_  __/ /__/ / / / __   __________  __     |")
    writeLines("|       / /   / _  \\/ _  _/\\ \\/ / _  __/ / / /  | / /_  __/ /_/ /     |")
    writeLines("|      / /___/ /_/ / /_\\ \\  \\  / / \\ \\ _/ /_/ /||/ / / / / __  /      |")
    writeLines("|     /_____/_/ /_/______/  /_/_/  /_/_____/_/ |__/ /_/ /_/ /_/       |")
    writeLines("|                                                                     |")
    writeLines("| LaByRInth: Low-coverage Biallelic R Imputation                      |")
    writeLines("| Copyright 2017 Jason Vander Woude and Nathan Ryder                  |")
    writeLines("| Licensed under the Apache License, Version 2.0                      |")
    writeLines("| Source code: github.com/Dordt-Statistics-Research/LaByRInth         |")
    writeLines("| Based on LB-Impute: github.com/dellaporta-laboratory/LB-Impute      |")
    writeLines("| Funding recieved from the National Science Foundation (IOS-1238187) |")
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

    browser()

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

    browser()

    thresh.mask <- apply(highest.posterior, 1:2, `>=`, 90)

    gt <- ifelse(thresh.mask, gt, "./.")

    concat <- matrix(paste0(gt, ":", ad, ":", gp), nrow=nrow(ad))

    vcf@gt <- cbind("GT:AD:GP", concat)
    colnames(vcf@gt) <- c("FORMAT", colnames(ad))
    rownames(vcf@gt) <- rownames(ad)


    write.vcf(vcf, "thresh-90.vcf.gz")

    print("done")
}

## LabyrinthImpute("./analysis/first-25-masked.vcf.gz", c("LAKIN", "FULLER"), 5, "./analysis/first-25-imp-f5.vcf.gz")

## INPROGRESS
## write code from emission to calls
## decide how to present the probabilities in the output vcf
## don't compute the emission probabilities for H2 because they are the same as H1
## cache the emission probabilities because they will often be the same
## estimate generation using equations for P(G|R) which are in the Libre Calc sheet
## fix the forward-backward algorithm from numerical instability
## check if extract.gt with na param as F is sufficient for both ad and gt fields

## check if cbind error is still occurring on the smaller vcf file
##   if so, normalize the emission probabilities

## look into whether it is worth changing the precomputed R files into C files

## do the incorrect calls typically happen at sites with low confidence in
## parental information?

## To tell Tintle
##  * scrapped almost all original code
##  * using vcf library
##  * made code more readible and more maintainable
##  * implemented desired improvement (multi-model)
##  * now using forward-backward instead of viterbi
##     * posterior probabilities
##  * logo!
##  * currently testing lakin-fuller imputed as f2
##  * currently debugging


## THIS SHOWS THAT THE POSTERIOR PROBABILITIES ARE INFLATED
## > interest <- an$correct[an$masked & an$max.prob >= 10 & an$max.prob <= 20]; sum(interest) / length(interest)
## [1] NaN
## > interest <- an$correct[an$masked & an$max.prob >= 20 & an$max.prob <= 30]; sum(interest) / length(interest)
## [1] NaN
## > interest <- an$correct[an$masked & an$max.prob >= 30 & an$max.prob <= 40]; sum(interest) / length(interest)
## [1] 0.5
## > length(interest)
## [1] 2
## > interest <- an$correct[an$masked & an$max.prob >= 40 & an$max.prob <= 50]; sum(interest) / length(interest)
## [1] 0.3142857
## > interest <- an$correct[an$masked & an$max.prob >= 50 & an$max.prob <= 60]; sum(interest) / length(interest)
## [1] 0.3191489
## > interest <- an$correct[an$masked & an$max.prob >= 60 & an$max.prob <= 70]; sum(interest) / length(interest)
## [1] 0.4369748
## > interest <- an$correct[an$masked & an$max.prob >= 70 & an$max.prob <= 80]; sum(interest) / length(interest)
## [1] 0.5182482
## > interest <- an$correct[an$masked & an$max.prob >= 80 & an$max.prob <= 90]; sum(interest) / length(interest)
## [1] 0.6733668
## > interest <- an$correct[an$masked & an$max.prob >= 90 & an$max.prob <= 100]; sum(interest) / length(interest)
## [1] 0.8978936


## > quick.check()
##     * thresh: 0	n: 6732	quality: 0.867944147355912
##     * thresh: 5	n: 6732	quality: 0.867944147355912
##     * thresh: 10	n: 6732	quality: 0.867944147355912
##     * thresh: 15	n: 6732	quality: 0.867944147355912
##     * thresh: 20	n: 6732	quality: 0.867944147355912
##     * thresh: 25	n: 6732	quality: 0.867944147355912
##     * thresh: 30	n: 6732	quality: 0.867944147355912
##     * thresh: 35	n: 6732	quality: 0.867944147355912
##     * thresh: 40	n: 6732	quality: 0.867944147355912
##     * thresh: 45	n: 6725	quality: 0.868550185873606
##     * thresh: 50	n: 6705	quality: 0.870096942580164
##     * thresh: 55	n: 6670	quality: 0.873013493253373
##     * thresh: 60	n: 6620	quality: 0.877190332326284
##     * thresh: 65	n: 6572	quality: 0.879945222154595
##     * thresh: 70	n: 6512	quality: 0.884674447174447
##     * thresh: 75	n: 6451	quality: 0.888389396992714
##     * thresh: 80	n: 6389	quality: 0.892158397245265
##     * thresh: 85	n: 6317	quality: 0.894253601393066
##     * thresh: 90	n: 6219	quality: 0.897893552018009
##     * thresh: 95	n: 6080	quality: 0.903125
##     * thresh: 100	n: 4174	quality: 0.914710110206037
