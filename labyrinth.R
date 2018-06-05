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
require(vcfR)
require(parallel)

LabyrinthImpute <- function (vcf, parents, out.file="", use.only.ad=TRUE,
                            leave.all.calls=TRUE, ref.alt.by.parent=TRUE,
                            recomb.double=TRUE, read.err=0.05,
                            genotype.err=0.05, recomb.dist=1e6,
                            write=TRUE, parallel=TRUE, cores=4,
                            quiet=FALSE) {

    ## If vcf is not already of class vcfR, assume it is a filename and try to
    ## read it. If not, read.vcfR will throw an error
    message("Checking if vcf is an object or file")
    if (! inherits(vcf, "vcfR")) {
        message("Loading vcf")
        vcf <- read.vcfR(vcf, verbose=F)
    }
    message("Checking if parents are in the vcf")
    if (! all(parents %in% getVARIANTS(vcf)))
        stop("Some or all parents are not in the vcf")
    message("Checking for sites that are not biallelic")
    non.biallelic <- which(! is.biallelic(vcf))
    if (length(non.biallelic) != 0) {
        message("The following sites are not biallelic:")
        chroms <- getCHROM(vcf)[non.biallelic]
        positions <- getPOS(vcf)[non.biallelic]
        for (i in seq_along(chroms)) {
            message("\tCHR:", chroms[i], "\tPOS:", positions[i])
        }
    }
    message("Checking for sites where parents are not homozygous within and polymorphic between")
    non.hom.poly <- which( ! parent.hom.and.poly(vcf, parents))
    if (length(non.hom.poly) != 0) {
        message("The parents are not homozygous within and polymorphic between at the following sites:")
        chroms <- getCHROM(vcf)[non.hom.poly]
        positions <- getPOS(vcf)[non.hom.poly]
        ad <- extract.gt(vcf, "AD")[non.hom.poly, parents]
        for (i in seq_along(chroms)) {
            message("\tCHR:", chroms[i],
                    "\tPOS:", positions[i],
                    paste0("\t", parents[1], ":"), ad[i, 1],
                    paste0("\t", parents[2], ":"), ad[i, 2])
        }
    }
    if (length(non.biallelic) != 0 ||
        length(non.hom.poly) != 0)
        stop("The vcf must be filtered before imputing; please run LabyrinthFilter first")

    message("Generating emission probabilities")
    ## INPROGRESS

    emission.structures <- get.emission.structures(vcf, parents, read.err)
    global.emission.structures <<- get.emission.structures(vcf, parents, read.err)
    browser()

    message("LaByRInth imputation complete")

}

## remove all sites that are not homozygous within and polymorphic between for
## the parents and remove all sites that are not biallelic
LabyrinthFilter <- function(vcf, parents, out.file) {
    message("Checking if vcf is an object or file")
    if (! inherits(vcf, "vcfR")) {
        message("Loading vcf")
        vcf <- read.vcfR(vcf, verbose=F)
    }
    message("Checking if parents are in the vcf")
    if (! all(parents %in% getVARIANTS(vcf)))
        stop("Some or all parents are not in the vcf")
    message("Checking for sites that are not biallelic")
    non.biallelic <- ! is.biallelic(vcf)
    if (length(non.biallelic) != 0) {
        message("The following sites are not biallelic and will be removed:")
        chroms <- getCHROM(vcf)[which(non.biallelic)]
        positions <- getPOS(vcf)[which(non.biallelic)]
        for (i in seq_along(chroms)) {
            message("\tCHR:", chroms[i], "\tPOS:", positions[i])
        }
    }
    message("Checking for sites where parents are not homozygous within and polymorphic between")
    non.hom.poly <- ! parent.hom.and.poly(vcf, parents)
    if (length(non.hom.poly) != 0) {
        message("The parents are not homozygous within and polymorphic between at the following sites which will be removed:")
        chroms <- getCHROM(vcf)[which(non.hom.poly)]
        positions <- getPOS(vcf)[which(non.hom.poly)]
        ad <- extract.gt(vcf, "AD")[which(non.hom.poly), parents]
        for (i in seq_along(chroms)) {
            message("\tCHR:", chroms[i],
                    "\tPOS:", positions[i],
                    paste0("\t", parents[1], ":"), ad[i, 1],
                    paste0("\t", parents[2], ":"), ad[i, 2])
        }
    }

    mask <- !non.biallelic & !non.hom.poly
    vcf@fix <- vcf@fix[mask, ]
    vcf@gt <- vcf@gt[mask, ]
    write.vcf(vcf, paste0(out.file, ".vcf.gz"), mask=TRUE)
    message(paste("\nFiltering is complete;", sum(!mask), "of", length(mask), "sites removed"))
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
    par.mat <- extract.gt(vcf, "AD")[ , parents]
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
fwd.bkwd <- function(emm, trans) {
    n.states <- nrow(emm)
    n.sites <- ncol(emm)

    ## forward probabilities
    f.probs <- matrix(data=0, nrow=n.states, ncol=n.sites)
    b.probs <- f.probs

    start.probs <- rep(1/n.states, n.states)
    f.probs[, 1] <- start.probs * emm[, 1]

    for (site in 2:n.sites) {
        t.index <- site - 1
        prev.site <- site - 1
        for (to in 1:n.states) {
            f.probs[to, site] <-
                emm[to, site] * sum(trans[ , to, t.index] * f.probs[, prev.site])
        }
    }
    fwd.prob <- sum(f.probs[, n.sites])

    ## backward probabilities
    end.probs <- rep(1, n.states)
    b.probs[, n.sites] <- end.probs

    for (site in (n.sites-1):1) {
        t.index <- site
        next.site <- site + 1
        for (from in 1:n.states) {
            b.probs[from, site] <-
                sum(trans[from, , t.index] * b.probs[, next.site] * emm[ , next.site])
        }
    }
    bkw.prob <- sum(start.probs * b.probs[, 1] * emm[ , 1])

    f.probs * b.probs / fwd.prob  # implicit return
}


## returns the emission probabilities for a given variant and chromosome. Note
## that the structure of the returned data is different than a traditional
## emission probability matrix. The first dimension of the data corresponds to
## the chromosome, and the second dimension corresponds to the sample. reads
## should be a vector of strings where each string is a numeric value followed
## by a comma and another numeric value. This is the format that the vcfR
## package stores the reads in.
get.emission.structures <- function(vcf, parents, rerr) {

    str.ad <- extract.gt(vcf, "AD")[ , parents[1]]
    p1.is.ref <- sapply(str.ad, function(str) {
        ## split the string and check if reference read is nonzero
        ad.to.num(str)[1] != 0
    })

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

        if (state == states["P1"])
            ret <- choose(n, n1) * (1-rerr)^n1 * (rerr)^n2
        else if (state == states["P2"])
            ret <- choose(n, n2) * (1-rerr)^n2 * (rerr)^n1
        else if (state == states["H1"] || state == states["H2"])
            ret <- choose(n, n1) * (0.5)^n
        else
            stop("Invalid state; this should never happen")


        ## ret <- switch(state,
        ##        P1 = choose(n, n1) * (1-rerr)^n1 * (rerr)^n2,
        ##        P2 = choose(n, n2) * (1-rerr)^n2 * (rerr)^n1,
        ##        H1 = choose(n, n1) * (0.5)^n
        ##        H2 = choose(n, n1) * (0.5)^n
        ##        )  # implicit return
        ret
    }

    reads.emm.probs <- function(reads, p1.is.ref) {
        names(reads) <- NULL  # makes debuggin easier
        ## browser()
        do.call(rbind,
                lapply(states, function(state) {
                    ## mapply will repeat state as many times as necessary
                    mapply(FUN=prob, state, reads, p1.is.ref)
                }))
    }

    chroms <- getCHROM(vcf)
    u.chroms <- unique(chroms)
    samples <- getSAMPLES(vcf)
    ad <- extract.gt(vcf, "AD")

    listapply <- function(..., mc.preschedule=F, mc.cores=4) {
        mclapply(..., mc.preschedule=mc.preschedule, mc.cores=mc.cores)
    }

    ret.val <- listapply(u.chroms, function(chrom) {
        ret.val.2  <- listapply(samples, function(sample) {
            ## emissions will be useless for parents, but that is fine
            ## if (sample == parents[1])
            ##     return rep(states["P1"],
            reads.emm.probs(ad[chroms==chrom, colnames(ad)==sample], p1.is.ref[chroms==chrom])
        })
        names(ret.val.2) <- samples
        ret.val.2
    })

}


## convert string representation to numeric vector
ad.to.num <- function(str)
    as.numeric(str.split(str, ","))



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
    cat(sprintf(paste0(" *  ",  "Progress: %.2f%%\r"), progress * 100))
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
