## TODO(Jason): Allow user to spefify threshold of probability to call a
## genotype without haplotype information. E.g. if there is not good

## TODO(Jason): Try writing an empty file immediately in case they give a bad
## destination and provide nice directory DNE message or file exists message
LabyrinthImpute <- function(vcf, parents, out.file, generation, use.only.ad=TRUE,
                            leave.all.calls=TRUE, ref.alt.by.parent=TRUE,
                            recomb.double=TRUE, read.err=0.05,
                            genotype.err=0.05, recomb.dist=1e6,
                            write=TRUE, parallel=TRUE, cores=4,
                            quiet=FALSE) {

    start.time <- Sys.time()
    pseudo.start.time <- start.time

    ## Create a preferences objects containing all preferences
    prefs <- list()
    class(prefs)            <- "prefs"

    ## Algorithm parameters
    prefs$recomb.double     <- recomb.double
    prefs$read.err          <- read.err
    prefs$genotype.err      <- genotype.err
    prefs$recomb.dist       <- recomb.dist
    prefs$generation        <- generation
    ## should the GT info be inferred from the AD info
    prefs$use.only.ad       <- use.only.ad
    ## Should non-imputed sites be in the output VCF file
    prefs$leave.all.calls   <- leave.all.calls
    prefs$parents           <- parents

    prefs$states            <- 3     # currently only support for 2 parents
    ## TODO(Jason): implement this feature
    prefs$ref.alt.by.parent <- FALSE # Should the reference and alternate be
                                     # switched in the output so that parent 1
                                     # is always reference and parent 2 is
                                     # always alternate

    ## Logistic parameters
    prefs$quiet             <- quiet
    prefs$cores             <- cores
    prefs$parallel          <- parallel
    prefs$write             <- write
    prefs$out.file          <- out.file

    ValidatePreferences(prefs)

    prefs$lapply <- get.lapply(prefs)

    print.labyrinth.header()

    writeLines(paste0(" *  Running imputation in ", ifelse(prefs$parallel,
               paste0("parallel (", prefs$cores, " cores)\n"), "serial\n")))


    ## if the vcf is actually a filename, load the vcf object
    if (! inherits(vcf, "vcf")) {
        ## if the vcf is not even a string, stop
        if (! inherits(vcf, "character")) {
            stop("vcf must either be a filename or an object of class vcf")
        }
        writeLines(paste0(" *  Loading VCF file ", file))

        vcf <- VCF(file, prefs)

        end.time <- Sys.time()
        time <- difftime(end.time, pseudo.start.time)
        pseudo.start.time <- end.time
        runtime <- as.numeric(time)
        units <- attr(time, "units")
        writeLines(paste0(" *  VCF loaded in ", round(runtime, 2), " ", units, "\n"))
    }


    ## vcf$gen.dists <- get.gen.dists(vcf)


    chroms <- vcf$chrom.names
    variants <- vcf$variant.names
    parent.geno <- ResolveHomozygotes(vcf, prefs$parents)
    n.chrom <- length(chroms)
    n.variants <- length(variants)
    n.sites <- nrow(parent.geno)
    prefs$n.jobs <- n.variants * n.chrom

    writeLines(paste0(" *  Imputing ",
                      (n.variants - 2), " variants at ",  # parents aren't imputed
                      n.chrom, " chromosomes (",
                      n.sites, " sites)"))
    writeLines(paste0(" *  ", prefs$n.jobs,
                      " imputations will run (",
                      (n.variants - 2), " x ", n.chrom, ")"))  # no parents

    ## Progress monitor code from https://stackoverflow.com/questions/27726134/
    ## how-to-track-progress-in-mclapply-in-r-in-parallel-package
    ## TODO(Jason): don't use prefs$fifo, but instead try to use a fifo variable
    ## in the progress.env environment
    progress.env <- new.env()
    prefs$fifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prefs$prog.env <- progress.env

    ## Actually run the imputation
    result <- do.call(cbind, prefs$lapply(variants, function(variant) {
        LabyrinthImputeSample(vcf, variant, parent.geno, prefs)
    }, mc.preschedule=FALSE, mc.cores=prefs$cores))

    close(prefs$fifo)

    ## Assign the GT of the vcf object according to the results
    for (site in 1:n.sites) {
        for (var in 1:n.variants) {
            call <- result[site, var]
                if (is.na(call)) {
                    vcf$GT[site, var, ] <- NA_integer_
                } else if (call %in% 1:2) {
                    ## set GT as homozygous same allele as corresponding parent
                    vcf$GT[site, var, ] <- parent.geno[site, call]
                } else if (call == 4 || call == 8 || call == 4+8) {
                    vcf$GT[site, var, ] <- 0:1  # heterozygous
                } else {
                    ## partial imputation will be ignored
                    vcf$GT[site, var, ] <- NA_integer_
                }
        }
    }

    vcf$header.lines <- get.labyrinth.vcf.header()

    end.time <- Sys.time()
    time <- difftime(end.time, pseudo.start.time)
    pseudo.start.time <- end.time
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste0(" *  Imputation completed in ", round(runtime, 2), " ", units, "\n"))

    if (prefs$write) {
        ## Replace spaces and colons in the date with dashes
        if (prefs$out.file == "") {
            ## Write to the same directory as before, but prepend "LaByRInth" to the name
            ## TODO(Jason): problem: this won't work because file dne
            qualified.name <- str.split(file, "/")
            name <- qualified.name[length(qualified.name)]
            if (grepl("/", file)) {
                ## replace everything after and including the last '/' with '/'
                ## then concatonate that directory with "LaByRInth_" and name
                prefs$out.file <- paste0(gsub("[/][^/]*$", "/", file), "LaByRInth_", name)
            } else {
                prefs$out.file <- paste0("LaByRInth_", name)
            }
        }

        ## TODO(Jason): prefs$out.file is not updated if it is ""
        writeLines(paste0(" *  Writing results to ", prefs$out.file))

        WriteVCF(vcf, prefs$out.file)

        end.time <- Sys.time()
        time <- difftime(end.time, pseudo.start.time)
        pseudo.start.time <- end.time
        runtime <- as.numeric(time)
        units <- attr(time, "units")
        writeLines(paste0(" *  File write completed in ", round(runtime, 2), " ", units, "\n"))
    }
    end.time <- Sys.time()
    time <- difftime(end.time, start.time)
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste0(" *  LaByRInth completed in ", round(runtime, 2), " ", units, "\n\n"))
    invisible(vcf)  # implicit invisible return
}


LabyrinthImputeSample <- function(vcf, sample, parent.geno, prefs) {

    if (sample %in% prefs$parents) {
        ## Find which parent number it is and repeat that number as many times
        ## as necessary. i.e. parent 1 will retrun 11111...11111
        return(rep(match(sample, prefs$parents), nrow(parent.geno)))
    }

    chroms <- vcf$chrom.names

    do.call(c, prefs$lapply(chroms, function(chrom) {  # c is the concatenate R function
        result <- LabyrinthImputeChrom(vcf, sample, chrom, parent.geno, prefs)
        writeBin(1/prefs$n.jobs, prefs$fifo)  # update the progress bar info
        if (!prefs$parallel) {  # if running in serial mode
            prefs$prog.env$progress <- PrintProgress(prefs$fifo, prefs$prog.env$progress)
        }  # else the forked process handles this
        result  # implicit return
    }, mc.preschedule=FALSE, mc.cores=prefs$cores))
}


LabyrinthImputeChrom <- function(vcf, sample, chrom, parent.geno, prefs) {

    if (length(sample) != 1) {
        stop("Length of sample must be 1")
    }
    if (length(chrom) != 1) {
        stop("Length of chrom must be 1")
    }

    print("imputing chrom")
    browser()
    emission.probs <- GetProbabilities(vcf, sample, chrom, parent.geno, prefs)

    site.pos <- sapply(rownames(emission.probs), function(name) {
        as.numeric(str.split(name, ":")[2])
    })

    ## This is the sites where both parents were called (not NA) and where they
    ## are different from each other. It is a boolean vector indicating whether
    ## the site is relevant, thus the length is the same as the length of the
    ## final imputation for this sample and chromosome

    relevant.sites <- GetRelevantProbabiltiesIndex(vcf, chrom, parent.geno, prefs)
    #writeLines(paste0(chrom, " ", sample, ": ", paste0(relevant.sites, collapse="")))

    informative.sites <- apply(emission.probs, 1, function(row) {
        !all(row == row[1])
    })
    ## relevant.sites <- informative.sites

    if (length(relevant.sites) != length(informative.sites)) {
        stop("Site index arrays differ")
    }
    relevant.sites <- relevant.sites & informative.sites

    ## If there are not enough markers according to user preference (or if there
    ## are 0), then do not do the imputation and return a path of NA's of the
    ## correct length
    n.relevant.sites <- sum(relevant.sites)  # boolean addition
    if (n.relevant.sites < 1) {
        full.path <- rep(NA_integer_, length(relevant.sites))
    } else {
        names(relevant.sites) <- NULL  # Makes debugging easier

        ## distances between relevant sites
        dists <- diff(site.pos[relevant.sites])

        relevant.probs <- emission.probs[relevant.sites, , drop=F]
        class(relevant.probs) <- "probs"

        path <- multiviterbi(relevant.probs, dists, prefs)

        full.path <- relevant.sites

        path.index <- 1
        filler <- NA_integer_
        ## The missing calls that were not relevant will be filled back in
        ## to create the full path from the relevant part of the path
        for (i in seq_along(relevant.sites)) {
            if (relevant.sites[i]) {  # if the site was relevant
                full.path[i] <- path[path.index]  # set to next call
                path.index <- path.index + 1  # increment call index
            } else {
                ## If we can safely decrement the index and elements are same
                if (path.index > 1 &&
                    path.index <= length(path) &&
                    !is.na(path[path.index - 1]) &&
                    !is.na(path[path.index]) &&
                    path[path.index - 1] == path[path.index]) {
                    filler <- path[path.index]
                } else {
                    filler <- NA_integer_
                }
                full.path[i] <- filler
            }
        }
    }

    ## At this stage full.path has entries of 1 through 7, or NA which indicates
    ## the call at that site according to the following table
    ## 1: Homozygous and the allele is the same as parent 1
    ## 2: Homozygous and the allele is the same as parent 2
    ## 4: Heterozygous
    ## -------------------------------------------------------------------------
    ## In the same way that binary counting works, we can use these 3 'basis'
    ## values to explain what other numbers indicate. This is shown below
    ## -------------------------------------------------------------------------
    ## 3 (1+2): Homozygous, but the actual allele is unknown
    ## 5 (1+4): One of the alleles matches parent 1, but the other is unknown
    ## 6 (2+4): One of the alleles matches parent 2, but the other is unknown
    ## 7 (1+2+4): Nothing is known about the alleles
    ## NA: The site was not imputed

    full.path  # implicit return
}
