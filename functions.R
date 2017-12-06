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


library(abind)
library(parallel)

ensure_writability <- function(filename) {
  # The '/' is literal.
  # The '[^/]*' will match any character other than '/' any number of time.
  # The '$' only matches the end of the line
  # Basically, this strips any trailing characters and the last '/' to get
  # the directory of the file
  dir <- gsub("/[^/]*$", "", filename)
  if (!dir.exists(dir))
    dir.create(dir, recursive=T)
}


# Function that will call the standard writeLines
# function after first ensuring that the directory
# being saved to exists. If the directory does not
# exist it is created.
force_write_lines <- function(data, filename) {
  ensure_writability(filename)
  writeLines(data, filename)
}


translate <- function(sample) {
    ## make the following replacements
    ## 1 -> 0
    ## 2 -> 1
    ## 4 -> 2
    ## anything else -> 5
    ifelse(is.na(sample), 5, ifelse(sample==1, 0, ifelse(sample==2, 1, ifelse(sample==4, 2, 5))))
}


##' Produce an image showing the difference between two others
##'
##' Input consists of two pgm type images representing the genotype calls of the
##'      vcf. The output is a red and green image where a red pixel indicates
##'      that the two input ima
##' @title
##' @param image1 path to first pgm image
##' @param image2 path to second pgm image
##' @param out.image path to ppm image representing the difference
##' @return
##' @author Jason Vander Woude
make.diff.ppm <- function(image1, image2, out.image) {j
    im1 <- data.table:::fread(image1)
    im2 <- data.table:::fread(image2)

    if (nrow(im1) != nrow(im2)) {
        stop("heights of images do not match")
    }
    if (ncol(im1) != ncol(im2)) {
        stop("widths of images do not match")
    }

    diff.image <- im1==im2
    sink(out.image)
    writeLines("P3")
    writeLines(paste(ncol(im1), nrow(im1), "1"))  # the last one indicates that RGB values will be 0 or 1
    invisible(apply(ifelse(diff.image, "0 1 0", "1 0 0"), 1, function(row) {writeLines(paste0(row, collapse=" "))}))
    sink()

    perc.diff <- 100 * sum(!diff.image) / (ncol(diff.image) * nrow(diff.image))
    writeLines(sprintf("There is a %.2f%% difference in the files", perc.diff))
}


strip.vcf <- function(file) {
    lines <- readLines(file)

    name <- str.split(file, "/")
    name <- name[length(name)]
    out.file <- paste0(gsub("[/][^/]*$", "/", file), "stripped_", name)

    sink(out.file)
    sapply(lines, function(line) {
        components <- str.split(line, "\t")
        writeLines(paste0(sapply(components, function(component){
            str.split(component, ":")[1]
        }), collapse="\t"))
    })
    sink()
}


temp <- function() {
    chroms <-  c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "UN")
    for (chrom in chroms) {
        make.diff(paste0("LaByRInth_", chrom, ".pgm"),
                  paste0("LB-Impute_", chrom, ".pgm"),
                  paste0("Diff_",      chrom, ".ppm"))
    }
}

vcf.to.pgm <- function(vcf.file, parents, out.file) {
    data <- data.table:::fread(vcf.file, blank.lines.skip=TRUE)
    parent.data <- data[, parents, with=FALSE]
    format.col <- which(colnames(data)=="FORMAT")
    data <- data[, (format.col+1):ncol(data)]
    data <- data[, !colnames(data) %in% parents, with=FALSE]
    ncol <- nrow(data)
    nrow <- ncol(data)
    data <- apply(data, 1:2, function(x){str.split(as.character(x), ":")[1]})
    parent.data <- apply(parent.data, 1:2, function(x){str.split(as.character(x), ":")[1]})

    none <- 5

    data2 <- sapply(1:nrow(data), function(row) {
        zero <- which(as.character(parent.data[row, ])=="0/0")
        one <- which(as.character(parent.data[row, ])=="1/1")
        if (length(zero)==0 || length(one)==0) {
            res <- rep(none, length(data[row, ]))
        } else {
            res <- mclapply(data[row, ], function(call) {
                if (call=="0/0") {
                    val <- zero-1
                } else if (call=="1/1") {
                    val <- one-1
                } else if (call=="0/1") {
                    val <- 2
                } else if (call=="1/0") {
                    val <- 2
                } else {
                    val <- none
                }
                val  # implicit return
            }, mc.cores=4)
        }
        res  # implicit return
    })
    sink(out.file)
    writeLines("P2")
    writeLines(paste(ncol, nrow, none))
    tmp <- apply(data2, 1, function(row){writeLines(paste0(row, collapse=" "))})
    sink()
}


##' Create an image representing an imputation result
##'
##' Generate an image where a black pixel represents a call for parent 1, a dark
##'      grey pixel represents a call for parent 2, a light grey pixel
##'      represents a call for heterozygous, and a white pixel represents a call
##'      that was not made.
##' @title
##' @param result result of LabyrinthImpute
##' @param image.dir location to save the images
##' @return
##' @author Jason Vander Woude
make.image.dir <- function(result, image.dir, parents) {
    sites <- rownames(result)
    chroms <- unique(sapply(sites, function(site) {str.split(site, ":")[1]}))
    variants <- colnames(result)[! colnames(result) %in% parents]
    nrow <- length(variants)

    for (chrom in chroms) {
        chrom.result <- result[grepl(chrom, sites), ]
        ncol <- nrow(chrom.result)
        sink(paste0(image.dir, "/LaByRInth_", chrom, ".pgm"))
        writeLines("P2")
        writeLines(paste(ncol, nrow, "5"))
        for (variant in variants) {
            calls <- chrom.result[, variant]
            writeLines(paste0(translate(calls), collapse=" "))
        }
        sink()
    }

    sink(paste0(image.dir, "/LaByRInth_all.pgm"))
    ncol <- nrow(result)
    writeLines("P2")
    writeLines(paste(ncol, nrow, "5"))
    for (variant in variants) {
        calls <- result[, variant]
        writeLines(paste0(translate(calls), collapse=" "))
    }
    sink()
}


result.to.pgm <- function(labyrinth.result, parents, out.file) {

    sites <- rownames(labyrinth.result)
    chroms <- unique(sapply(sites, function(site) {str.split(site, ":")[1]}))
    variants <- colnames(labyrinth.result)[! colnames(labyrinth.result) %in% parents]
    nrow <- length(variants)

    sink(out.file)
    ncol <- nrow(labyrinth.result)
    writeLines("P2")
    writeLines(paste(ncol, nrow, "5"))
    for (variant in variants) {
        calls <- labyrinth.result[, variant]
        writeLines(paste0(translate(calls), collapse=" "))
    }
    sink()
}


str.to.num <- function(str, sep) {
    as.numeric(str.split(str, sep))
}


## More convenient strsplit if length of vector is 1
str.split <- function(str, sep) {
    if (length(str) != 1) {
        warning("Only splitting the first element of the vector")
    }
    strsplit(str, sep)[[1]]
}


nsites.probs <- function(probs) {
    nrow(probs)
}


nstates.probs <- function(probs) {
    ncol(probs)
}


## Used for creating a 3-D structure
reorder <- function(x, n) {
    as.vector(sapply(1:n, function(i) {
        v <- rep(FALSE, n)
        v[i] <- TRUE
        x[v]
    }))
}


vec.or <- function(vec) {
    Reduce(`|`, vec)
}


vcf.to.numeric <- function(str.vec) {
    sapply(str.vec, function(str) {
        if (check.int(str)) {
            as.numeric(str)
        } else if (str == ".") {
            NA_integer_
        } else {
            stop("Attempted to convert invalid VCF string")
        }
    })
}


print.labyrinth.header <- function() {
    writeLines("\n")
    writeLines("+---------------------------------------------------------------------+")
    writeLines("| LaByRInth: Low-coverage Biallelic R Imputation                      |")
    writeLines("| Copyright 2017 Jason Vander Woude and Nathan Ryder                  |")
    writeLines("| Licensed under the Apache License, Version 2.0                      |")
    writeLines("| Source code: github.com/Dordt-Statistics-Research/LaByRInth         |")
    writeLines("| Based on LB-Impute: github.com/dellaporta-laboratory/LB-Impute      |")
    writeLines("| Funding recieved from the National Science Foundation (IOS-1238187) |")
    writeLines("+---------------------------------------------------------------------+")
    writeLines("")
}


print.labyrinth <- function(...) {
    writeLines(paste0(" *  ", ...))
}


## A way to check if a string contains only numeric characters. Code comes from
## stackoverflow.com/questions/13301437/how-to-check-if-the-value-is-numeric
check.int <- function(vec){
    sapply(vec, function(N) {
        !length(grep("[^[:digit:]]", as.character(N)))})
}


##' Generate a path from a path tracker
##'
##' Given a path tracker (which is a specific type of array produced by viterbi
##'     algorithm) and indices representing the final hidden states of the path,
##'     compute the paths that ended in those states. The returned path will be
##'     a vector of integers where integers which are powers of 2 represent the
##'     the state with the 0-based index which is the log_2 of the number.
##'     (e.g. 1,2,4,8,16,... represent the states associated with the 1st, 2nd,
##'     3rd, 4th, 5th, etc row)
##' @title
##' @param path.tracker the boolean 3D array in which the path info is embedded
##' @param boolean.indices which rows do optimal paths end at
##' @return the path
##' @author Jason Vander Woude
generatePath <- function(path.tracker, boolean.indices) {
    if (length(boolean.indices) != dim(path.tracker)[1]) {
        stop("boolean.indices has wrong length")
    }
    path <- rep(NA, dim(path.tracker)[2])
    base <- 0:(dim(path.tracker)[1] - 1)
    powers <- 2**base
    for (i in length(path):1) {
        path[i] <- sum(powers[boolean.indices])
        slice <- path.tracker[, i, ]
        slice.optimal <- slice[boolean.indices, , drop=F]
        boolean.indices <- as.logical(apply(slice.optimal, 2, vec.or))
    }
    path  # implicitly return path
}


##' Find the most probable paths using the viterbi algorithm
##'
##' See http://homepages.ulb.ac.be/~dgonze/TEACHING/viterbi.pdf for details
##'     about the viterbi algorithm and understanding this implementation. This
##'     imlementation is designed to handle multiple paths with identical
##'     probabilities.
##' @title
##' @param probs emission probabilities
##' @param dists vector of distances between sites
##' @param prefs imputation preferences object
##' @return the most probable sequence of states
##' @author Jason Vander Woude
viterbi <- function(probs, dists, prefs) {
    nstates <- nstates.probs(probs)
    path.size <- nsites.probs(probs)

    ## 3D array
    paths.tracker <- array(NA, dim=c(nstates, path.size, nstates))

    ## This will keep track of the overall probabilities of each of the {nstates}
    ## final paths, so that the probability of the final paths do not have to be
    ## computed again. The probabilities are initialized to the emission
    ## probabilities of the first site of the probs matrix which is the first
    ## row
    if (path.size < 1) {
        stop("viterbi requires that probs have > 0 rows")
    }
    probs.tracker <- log(probs[1, ])

    ## Hard code the first column to the vector 1,2,...,nstates as an index
    ## This is what the generatePath function will need
    paths.tracker[, 1, ] <- diag(TRUE, nstates)

    if (path.size != 1) {  # if the path size is 1 just use the emission probs
        for (site in 2:path.size) {  # for each site in the path

            dist <- dists[site - 1]

            ## log of the probability for each possible hidden state at this site
            probs.tracker <- sapply(1:nstates, function(state) {

                extension.probs <- sapply(1:nstates, function(i) {
                    ## log of the probability of being at state i before and
                    ## transitioning to the 'state' state
                    probs.tracker[i] + log(transProb(i, state, dist, prefs))
                })

                optimal.indices <- extension.probs==max(extension.probs)
                ## use <<- to assign to a variable outside the current scope
                paths.tracker[state, site, ] <<- optimal.indices
                max(extension.probs) + log(probs[site, state])  # return prob
            })
        }
    }

    ## The code above has already computed the optimal path, but that
    ## information is encoded within the paths.tracker matrix and needs to be
    ## extracted. That is what generatePath will do when passed the path.tracker
    ## matrix and the indices of the optimal path.
    generatePath(paths.tracker, probs.tracker==max(probs.tracker))  # return best path
}


##' Find the transission probability between hidden states
##'
##' Use the equations specified in the LB-Impute paper to compute the
##'     transmission probability between two sites.
##' @title
##' @param a first site
##' @param b second site
##' @param dist the distance between sites a and b
##' @param prefs a preferences object
##' @return the transmission probability between these hidden states
##' @author Jason Vander Woude
transProb <- function(a, b, dist, prefs) {
    if (a == b) {
        ## If there is no recombination
        0.5 * (1 + exp(-dist/prefs$recomb.dist))
    } else if (a %in% 1:2 && b %in% 1:2 && prefs$recomb.double) {
        ## Double recomb occurred and has square of probability of single recomb
        ## This only works for biparental
        (0.5 * (1 - exp(-dist/prefs$recomb.dist)))**2
    } else {
        ## Some type of recombination occurred that has regular probability
        0.5 * (1 - exp(-dist/prefs$recomb.dist))
    }
}


##' Extract the information from a vcf file and save it as a vcf object
##'
##' The returned vcf object will have the following: variants, header.lines,
##'     variant.names, chrom.names, GT, AD.
##' @title
##' @param file the path to the vcf file
##' @return the vcf object created from the file
##' @author Jason Vander Woude
VCF <- function(file, prefs, required=c("AD","GT")) {
    ## TODO(Jason): Save rds version of vcd to impute again with different prefs

    ## TODO(Jason): add filtering step to remove non-biallelic calls

    vcf <- list()
    class(vcf) <- "vcf"

    vcf$variants <- readLines(file)
    isComment <- sapply(vcf$variants, function(line){substr(line,1,1) == "#"})

    ## Remove the header, but save so it can be restored later if desired. In
    ## the current implementation, the output vcf file does not contain the
    ## original header except the first line which should specify the version
    vcf$header.lines <- vcf$variants[isComment]
    vcf$variant.lines <- vcf$variants[!isComment]
    vcf$variants <- vcf$variant.lines

    ## Make table
    vcf$variants <- do.call(rbind,
                            lapply(vcf$variants,
                                   function(line){str.split(line, "\t")}))

    header <- vcf$header.lines[length(vcf$header.lines)]  #get column heading
    header <- substr(header, 2, nchar(header)) #remove leading '#'
    colnames(vcf$variants) <- str.split(header, "\t")

    formatExample <- vcf$variants[1, "FORMAT"]
    field.names <- str.split(formatExample, ":")

    ## ## Verify that required fields are in the VCF
    ## if (prefs$use.only.ad) {
    ##     required.fields <- "AD"  # could also use GQ
    ## } else {
    ##     required.fields <- c("GT", "AD")  # could also use GQ
    ## }
    ## if (!all(required.fields %in% field.names)) {
    ##     stop(paste("VCF file does not contain all required fields.",
    ##                "Required fields are", toString(required.fields)))
    ## }
    ## field.indices <- match(required.fields, field.names)
    ## names(field.indices) <- required.fields

    ## The "FORMAT" column is the last one before the variants start
    format.col <- match("FORMAT", colnames(vcf$variants))
    n.variants <- ncol(vcf$variants) - format.col
    n.sites <- nrow(vcf$variants)

    samples <- vcf$variants[ , (format.col + 1):ncol(vcf$variants)]
    rownames(samples) <- paste0(vcf$variants[, "CHROM"], ":", vcf$variants[, "POS"])

    vcf$variant.names <- colnames(samples)
    vcf$chrom.names <- unique(vcf$variants[, "CHROM"])

    ## AD and GT section
    ## Third dim = 2 because all sites must be biallelic
    available.fields <- c()
    if ("AD" %in% field.names) {
        vcf$AD <- array(NA_integer_, dim=c(n.sites, n.variants, 2))
        vcf$DP <- matrix(NA_integer_, nrow=n.sites, ncol=n.variants)
        colnames(vcf$AD) <- colnames(samples)
        rownames(vcf$AD) <- rownames(samples)
        colnames(vcf$DP) <- colnames(samples)
        rownames(vcf$DP) <- rownames(samples)
        available.fields <- c(available.fields, "AD")
    }
    if ("GT" %in% field.names) {
        vcf$GT <- array(NA_integer_, dim=c(n.sites, n.variants, 2))
        ## vcf$GT.coded <- matrix(NA_integer_, nrow=n.sites, ncol=n.variants)
        colnames(vcf$GT) <- colnames(samples)
        rownames(vcf$GT) <- rownames(samples)
        ## colnames(vcf$GT.coded) <- colnames(samples)
        ## rownames(vcf$GT.coded) <- rownames(samples)
        available.fields <- c(available.fields, "GT")
    } else {
        stop("VCF format requires that a GT (genotype) field be included")
    }
    ## Set up index array to get field position by name
    field.indices <- match(available.fields, field.names)
    names(field.indices) <- available.fields

    progress.env <- new.env()
    prefs$fifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prefs$prog.env <- progress.env
    prefs$n.jobs <- n.sites
    ## TODO(Jason): clean up the progress code
    ##InitiateProgress(prefs, n.jobs=n.sites*n.variants)

    for (r in 1:n.sites) {
        ## make sure site is biallelic
        if (length(str.split(vcf$variants[r, "ALT"], ",")) != 1 ) {
            vcf$AD[r, , ] <- NA_integer_
            vcf$GT[r, , ] <- NA_integer_
            writeLines(paste(" *  Removing", rownames(samples)[r],
        "from sites as it is not biallelic"))
            next()
        }
        for (c in 1:n.variants) {
            if ("AD" %in% field.names) {
                ## Separate out AD field and split it into two read depths
                ret.val.ad <- vcf.to.numeric(str.split(
                    str.split(samples[r, c], ":")[field.indices["AD"]], ","))
                if (length(ret.val.ad) > 2) {
                    stop(paste0("There are more than two observed alleles at ",
                                rownames(samples)[r], " for variant ", colnames(samples)[c]))
                }
                if (length(ret.val.ad) < 2 && !is.na(ret.val.ad)) {
                    stop(paste0("There are less than two observed alleles at ",
                                rownames(samples)[r], " for variant ", colnames(samples)[c]))
                }
                vcf$AD[r, c, ] <- ret.val.ad
                vcf$DG <- sum(ret.val.ad)
            }

            if ("GT" %in% field.names) {
                if (prefs$use.only.ad && "AD" %in% field.names) {
                    ret.val.gt <- GenotypeFromDepth(ret.val.ad)
                    vcf$GT[r, c, ] <- ret.val.gt
                    vcf$GT.coded <- CodifyGT(ret.val.gt)
                } else { ## TODO(Jason): warn of else condition outside loop
                    ## Separate out GT field and split it into two calls
                    ret.val.gt <- str.split(samples[r,c], ":")[field.indices["GT"]]
                    ret.val.gt <- gsub("\\|", "/", ret.val.gt)  # replace '|' with '/'
                    ret.val.gt <- vcf.to.numeric(str.split(ret.val.gt, "/"))
                    if (length(ret.val.gt) > 2) {
                        stop(paste0("There are more than two genotype calls at ",
                                    rownames(samples)[r], " for variant ", colnames(samples)[c]))
                    }
                    if (length(ret.val.gt) < 2) {
                        stop(paste0("There are less than two genotype calls at ",
                                    rownames(samples)[r], " for variant ", colnames(samples)[c]))
                    }
                    vcf$GT[r, c, ] <- ret.val.gt
                }
            }
        }

        MakeProgress(prefs)

    }

    vcf  # implicit return
}


##' Get a subset of the data from the vcf object
##'
##' Get data pertaining to the specified field and subset it by samples and
##'     chromosomes.
##' @title
##' @param vcf an object of class vcf
##' @param field the name of a single field ("GT" or "AD")
##' @param samples a vector of indices or names of variants
##' @param chromosomes a vector of chromosome names to subset by
##' @return a 3-dimensional array representing the subset of the data
##' @author Jason Vander Woude
Get <- function(vcf, field, samples, chromosomes=NULL) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (length(field) != 1) {
        stop("Length of field must be 1")
    }
    if (is.null(chromosomes)) {
        chromosomes <- vcf$chrom.names
    }
    rows <- vcf$variants[, "CHROM"] %in% chromosomes
    vcf[[field]][rows, samples, , drop=F]
}


##' Resolve heterozygous sites in the samples
##'
##' Returns a matrix of genotype calls for the samples such that the entry is 0
##'     if all calls were for the reference allele; 1 if all calls were for the
##'     alternate allele; and NA if there were calls for both the reference and
##'     alternate allele.
##' @title
##' @param vcf an object of class vcf
##' @param samples a vector of indices or names of variants
##' @param prefs a preferences object
##' @return a matrix of genotypes (0, 1, NA)
##' @author Jason Vander Woude
ResolveHomozygotes <- function(vcf, samples) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    ## TODO(Jason): We can set the parents to NA now, but they should be set
    ## back to the original call before turning the children imputation states
    ## back into genotype calls (I think). That is we don't want to call all
    ## children as NA just because one of the parents was unknown. This may not
    ## be an issue since we are filtering the parents, but if we upgrade the
    ## code to handle non-HWPB (homozygous within and polymorphic between) this
    ## may become important

    ## TODO(Jason): This is the step where we are making an absolute claim on
    ## the true genotype of the parents which means we should be accounting for
    ## a possible read error here. If the allelic depths are 1 and 100 for
    ## example, we are currently claiming that the site is heterozygous when in
    ## reality it is probably homozygous with a read error. In addition, if we
    ## are more sure about the parents at a site, we should be more heavily
    ## weighting the emission probabilities at that site in the viterbi. Right
    ## now all sites have equal weight, but it is possible that there are times
    ## we get the parents wrong. Again this is probably not a huge concern with
    ## the code as is, but especially if we begin utilizing non-HWPB data, this
    ## is a feature that should be incorporated.

    ## TODO(Jason): Again, if utilizing non-HWPB data, firs remove false
    ## homozygosity because otherwise variants will only be caled homozygous of
    ## the variant seen or heterozygous even if they were actually called
    ## homozygous of the other.
    genotype <- Get(vcf, "GT", samples)
    allele.counts <- Get(vcf, "AD", samples, vcf$chrom.names)
    ret.val <- genotype[, , 1]  # initilize to first slice in 3rd dim
    for (r in 1:nrow(genotype)) {  # by chromosome site
        for (c in 1:ncol(genotype)) {  # by variant/sample
            ## Observed alleles from a specific site of a variant chromosome
            alleles <- genotype[r, c, ]
            if (any(is.na(alleles))) {  # One of the alleles is NA
                ret.val[r, c] <- NA
            } else if (all(alleles == alleles[1])) {  # Check if all are same
                ret.val[r, c] <- alleles[1]
            } else if (sum(allele.counts[r, c, ] != 0) == 1) {  # Only 1 counted
                ret.val[r, c] <- alleles[as.logical(allele.counts[r, c, ])]
            } else {  # Contradictory calls
                ret.val[r, c] <- NA
            }
        }
    }
    ret.val  # implicit return
}


##' Get matrix of emission probabilities
##'
##' Computes the emission probabilities for each state based on allelic depth of
##'     coverage (using the binomial assumption). Sample must be of length 1 and
##'     the entries of the parent.geno matrix must be either 0, 1, or NA.
##' @title
##' @param vcf an object of class vcf
##' @param sample an index or name of a variant
##' @param parent.geno a matrix of parental genotypes
##' @param prefs a preferences object
##' @return a matrix of posterior probabilities
##' @author Jason Vander Woude
GetProbabilities <- function(vcf, sample, chrom, parent.geno, prefs) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (length(sample) != 1) {
        stop("Length of sample must be 1")
    }
    if (length(chrom) != 1) {
        stop("Length of sample must be 1")
    }

    chrom.indices <- grepl(paste0(chrom, ":"), rownames(parent.geno))
    parent.geno <- parent.geno[chrom.indices, ]

    gt <- Get(vcf, "GT", sample, chrom)
    ad <- Get(vcf, "AD", sample, chrom)

    ret.val <- matrix(NA_integer_, nrow = nrow(gt), ncol = prefs$states)
    rownames(ret.val) <- rownames(gt)
    colnames(ret.val) <- c(colnames(parent.geno), "HETEROZYGOUS")
    class(ret.val) <- "prob"

    ## probability constraints from original LB-Impute code
    max.allowed <- 1 - (2 * prefs$genotype.err)
    min.allowed <- prefs$genotype.err
    rerr <- prefs$read.err

    for (row in 1:nrow(ret.val)) {
        ## gt and ad are both 3D structures where the first index is the site,
        ## the second index is the sample/variety, and the third index allows
        ## getting the various genotypes (e.g. a 0 in both indices of the last
        ## dimension represents 0/0). The '1' in [row, 1, ] is because there is
        ## guaranteed to be only one sample
        geno.calls <- gt[row, 1, ]
        allele.counts <- ad[row, 1, ]

        if (all(is.na(geno.calls))) {
            ret.val[row, ] <- rep(1, prefs$states)
        } else if (any(is.na(geno.calls))) {
            stop("Some but not all genotype calls are NA")  # TODO(Jason): remove
        } else {
            ref.calls <- allele.counts[1]
            alt.calls <- allele.counts[2]

            ## I think this is outdated code, but I'll leave it for awhile just
            ## in case
            ## if (is.na(ref.calls)) {
            ##     ref.calls <- 0
            ## }
            ## if (is.na(alt.calls)) {
            ##     alt.calls <- 0
            ## }

            ## Calculate the emission probabilities for this site
            ref.prob <- (1 - rerr)**ref.calls * (rerr)**alt.calls
            alt.prob <- (1 - rerr)**alt.calls * (rerr)**ref.calls
            hom.prob <- (0.5)**(ref.calls + alt.calls)  # homozygous
            max.prob <- max(ref.prob, alt.prob, hom.prob)

            normalize <- function(x) {
                ## TODO(Jason): Correction: max.allowed should be swapped with
                ## (max.allowed - min.allowed)
                x / max.prob * max.allowed + min.allowed
            }

            ## TODO(Jason): This seems incorrect. If the parent is unknown (NA)
            ## then why should the probability of being that parent be the max
            ## of alt.prob and ref.prob?
            for (state in 1:(prefs$states - 1)) {
                if (is.na(parent.geno[row, state])) {
                    ret.val[row, state] <- normalize(max(alt.prob, ref.prob))
                } else if (parent.geno[row, state] == 0) {
                    ret.val[row, state] <- normalize(ref.prob)
                } else if (parent.geno[row, state] == 1) {
                    ret.val[row, state] <- normalize(alt.prob)
                } else {
                    stop("Parental genotype was not NA, 0, or 1")
                }
            }
            ret.val[row, prefs$states] <- normalize(hom.prob)
        }
    }
    ret.val  # implicit return
}


## Determine which rows/sites have parents that are HWPB (homozygous within and
## polymorphic between)
GetRelevantProbabiltiesIndex <- function(vcf, chromosomes, parent.geno, prefs) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    rows <- vcf$variants[, "CHROM"] %in% chromosomes
    apply(parent.geno[rows, , drop=F], 1, function(calls) {
        ## parents should be distinct and non-NA
        !anyNA(calls) && length(calls) == length(unique(calls))
    })
}


## TODO(Jason): Allow user to spefify threshold of probability to call a
## genotype without haplotype information. E.g. if there is not good

## TODO(Jason): Try writing an empty file immediately in case they give a bad
## destination and provide nice directory DNE message or file exists message
LabyrinthImpute <- function(file, parents, out.file="", use.only.ad=TRUE,
                            leave.all.calls=TRUE, ref.alt.by.parent=TRUE,
                            recomb.double=TRUE, read.err=0.05,
                            genotype.err=0.05, recomb.dist=1e6,
                            write=TRUE, parallel=TRUE, cores=4,
                            quiet=FALSE) {

    ## Create a preferences objects containing all preferences
    prefs <- list()
    class(prefs)            <- "prefs"

    ## Algorithm parameters
    prefs$recomb.double     <- recomb.double
    prefs$read.err          <- read.err
    prefs$genotype.err      <- genotype.err
    prefs$recomb.dist       <- recomb.dist
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

    prefs  # implicit return

    start.time <- Sys.time()
    pseudo.start.time <- start.time

    ## if (is.null(prefs)) {
    ##     prefs <- GetDefaultPreferences()
    ## }

    ValidatePreferences(prefs)

    ## Determine whether to run in parallel and how many cores to use
    if (prefs$parallel) {
        require(parallel)
        prefs$lapply <- function(...,
                                 mc.preschedule=F,
                                 mc.cores=prefs$cores) {
            mclapply(..., mc.preschedule=mc.preschedule, mc.cores=mc.cores)
        }
    } else {
        prefs$lapply <- function(...,
                                 mc.preschedule=F,
                                 mc.cores=prefs$cores) {
            lapply(...)
        }
    }

    print.labyrinth.header()

    writeLines(paste0(" *  Running imputation in ", ifelse(prefs$parallel,
               paste0("parallel (", prefs$cores, " cores)\n"), "serial\n")))
    writeLines(paste0(" *  Loading VCF file ", file))

    vcf <- VCF(file, prefs)

    end.time <- Sys.time()
    time <- difftime(end.time, pseudo.start.time)
    pseudo.start.time <- end.time
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste0(" *  VCF loaded in ", round(runtime, 2), " ", units, "\n"))

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

    colnames(result) <- variants
    rownames(result) <- rownames(parent.geno)

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

        writeLines(paste0(" *  Writing results to ", prefs$out.file))

        ## TODO(Jason): don't use sink()
        sink(prefs$out.file)
        writeLines(vcf$header[1])  # Add header
        writeLines("##LaByRInth=<ID=Imputation,Version=1.0,Description=\"Code can be found at github.com/Dordt-Statistics-Research/LaByRInth\">")
        writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
        names <- vcf$header[length(vcf$header)]
        writeLines(names)  # Add header
        names <- str.split(names, "\t")

        ## The text on all the lines before the actual genotype call data
        prefix <- cbind(vcf$variants[, 1:which(names=="INFO")], "GT")
        ## TODO(Jason): this probably isn't necessary if I use the vcf$samples variable
        ref.geno <- Get(vcf, "GT", prefs$parents)
        ref.geno <- apply(ref.geno, 1:2, function(genotypes) {
            ## Concatenate the genotypes with a slash between
            text <- paste0(genotypes, collapse="/")
            text <- gsub("NA", ".", text)
        })

        ## TODO(Jason): this won't work unless I find a way around sink()
        ## progress.env <- new.env()
        ## prefs$fifo <- ProgressMonitor(progress.env)
        ## assign("progress", 0.0, envir=progress.env)
        ## prefs$prog.env <- progress.env
        ## prefs$n.jobs <- n.sites

        ## TODO(Jason): this should be cleaned up and there should be an option
        ## to use or not use partion imputation
        for (i in 1:n.sites) {
            ## Write prefix columns with tab seperation
            cat(paste0(prefix[i, ], collapse="\t"))
            cat("\t")
            for (j in 1:n.variants) {
                call <- result[i, j]
                if (is.na(call)) {
                    text <- "./."
                } else if (call %in% 1:2) {
                    text <- ref.geno[i, call]
                } else if (call == 4) {
                    text <- "0/1"
                } else if (call %in% c(3,5,6,7)) {
                    if (call == 3) {
                        text <- "./." # used to be X/X
                    } else if (call == 5) {
                        parent.call <- parent.geno[i, 1]  # 6 indicates parent 1
                        if (is.na(parent.call)) {
                            ## TODO(Jason): I'm not sure what to do here
                            text <- "./."
                        } else if (parent.call == 0) {
                            text <- "0/."
                        } else if (parent.call == 1) {
                            text <- "./1"
                        } else {
                            stop("Unexpected imputation result for value 5")
                        }
                    } else if (call == 6) {
                        parent.call <- parent.geno[i, 2]  # 6 indicates parent 2
                        if (is.na(parent.call)) {
                            ## TODO(Jason): I'm not sure what to do here
                            text <- "./."
                        } else if (parent.call == 0) {
                            text <- "0/."
                        } else if (parent.call == 1) {
                            text <- "./1"
                        } else {
                            stop("Unexpected imputation result for value 6")
                        }
                    } else if (call == 7) {
                        text <- "?/?"
                    }
                } else {
                    stop(paste("Invalid genotype call:", call))
                }
                cat(text)
                if (j != n.variants) {  # add tab if not last column
                    cat("\t")
                }
            }
            if (i != n.sites) {  # add newline if not last row
                cat("\n")
            }

            ## MakeProgress(prefs)
        }
        sink()  # turn off sink
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
    invisible(result)  # implicit return
}


## TODO(Jason): add option to switch ref and alt when filtering
LabyrinthFilter <- function(file, parents, out.file="", use.only.ad=TRUE,
                            leave.all.calls=TRUE, ref.alt.by.parent=TRUE,
                            recomb.double=TRUE, read.err=0.05,
                            genotype.err=0.05, recomb.dist=1e6,
                            write=TRUE, parallel=TRUE, cores=4,
                            quiet=FALSE) {

    ## Create a preferences objects containing all preferences
    prefs <- list()
    class(prefs)            <- "prefs"

    ## Algorithm parameters
    prefs$recomb.double     <- recomb.double
    prefs$read.err          <- read.err
    prefs$genotype.err      <- genotype.err
    prefs$recomb.dist       <- recomb.dist
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

    start.time <- Sys.time()
    pseudo.start.time <- start.time

    ## if (is.null(prefs)) {
    ##     prefs <- GetDefaultPreferences()
    ## }

    ## Determine whether to run in parallel and how many cores to use
    if (prefs$parallel) {
        require(parallel)
        prefs$lapply <- function(...,
                                 mc.preschedule=F,
                                 mc.cores=prefs$cores) {
            mclapply(..., mc.preschedule=mc.preschedule, mc.cores=mc.cores)
        }
    } else {
        prefs$lapply <- function(...,
                                 mc.preschedule=F,
                                 mc.cores=prefs$cores) {
            lapply(...)
        }
    }

    print.labyrinth.header()

    writeLines(paste0(" *  Running filter in ", ifelse(prefs$parallel,
               paste0("parallel (", prefs$cores, " cores)\n"), "serial\n")))
    writeLines(paste0(" *  Loading VCF file ", file))

    vcf <- VCF(file, prefs)
    end.time <- Sys.time()
    time <- difftime(end.time, pseudo.start.time)
    pseudo.start.time <- end.time
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste0(" *  VCF loaded in ", round(runtime, 2), " ", units, "\n"))

    parent.geno <- ResolveHomozygotes(vcf, prefs$parents)
    chroms <- vcf$chrom.names
    relevant.sites <- GetRelevantProbabiltiesIndex(vcf, chroms, parent.geno, prefs)
    n.sites <- length(relevant.sites)
    prefs$n.jobs <- n.sites

    ## Progress monitor code from https://stackoverflow.com/questions/27726134/
    ## how-to-track-progress-in-mclapply-in-r-in-parallel-package
    ## TODO(Jason): don't use prefs$fifo, but instead try to use a fifo variable
    ## in the progress.env environment
    progress.env <- new.env()
    prefs$fifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prefs$prog.env <- progress.env



    if (prefs$out.file == "") {
        ## Write to the same directory as before, but prepend "LaByRInth" to the name
        qualified.name <- str.split(file, "/")
        name <- qualified.name[length(qualified.name)]
        if (grepl("/", file)) {
            ## replace everything after and including the last '/' with '/'
            ## then concatonate that directory with "LaByRInth_" and name
            prefs$out.file <- paste0(gsub("[/][^/]*$", "/", file), "filtered_", name)
        } else {
            prefs$out.file <- paste0("filtered_", name)
        }
    }
    writeLines(paste0(" *  Writing results to ", prefs$out.file))

    sink(prefs$out.file)
    writeLines(vcf$header.lines)
    sapply(1:n.sites, function(site.num) {
        if (relevant.sites[site.num]) {
            writeLines(vcf$variant.lines[site.num])
        }
        writeBin(1/prefs$n.jobs, prefs$fifo)  # update the progress bar info
    })
    close(prefs$fifo)
    sink()  # turn off sink

    end.time <- Sys.time()
    time <- difftime(end.time, pseudo.start.time)
    pseudo.start.time <- end.time
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste0(" *  Filter and file write completed in ", round(runtime, 2), " ", units, "\n"))

    end.time <- Sys.time()
    time <- difftime(end.time, start.time)
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste0(" *  LaByRInth completed in ", round(runtime, 2), " ", units, "\n\n"))
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

    ## informative.sites <- apply(emission.probs, 1, function(row) {
    ##     !all(row == row[1])
    ## })
    ## relevant.sites <- informative.sites

    ## if (length(relevant.sites) != length(informative.sites)) {
    ##     stop("Site index arrays differ")
    ## }
    ## relevant.sites <- relevant.sites & informative.sites

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

        path <- viterbi(relevant.probs, dists, prefs)

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


GenotypeFromDepth <-  function(allelic.depths) {
    ad <- allelic.depths

    if (length(ad) != 2 && !all(is.na(ad))) {
        stop("GenotypeFromDepth does not support non-biallelic reads")
    }

    if (all(is.na(ad))) {
        genotype <- NA  # Why am I not using c(NA, NA) instead of NA
    } else if (all(ad == 0)) {
        genotype <- NA  # Why am I not using c(NA, NA) instead of NA
    } else if (all(ad != 0)) {
        genotype <- c(0, 1)
    } else if (ad[1] == 0) {
        genotype <- c(1, 1)
    } else {  # ad[2] == 0
        genotype <- c(0, 0)
    }
    genotype  # implicit return
}


## GetDefaultPreferences <- function() {
##     prefs <- list()
##     class(prefs)            <- "prefs"

##     ## Algorithm parameters
##     prefs$recomb.double     <- TRUE
##     prefs$read.err          <- 0.05
##     prefs$genotype.err      <- 0.05
##     prefs$recomb.dist       <- 1000000
##     prefs$states            <- 3     # currently only support for 2 parents
##     prefs$use.only.ad       <- TRUE  # should the GT info be inferred from the
##                                      # AD info
##     prefs$leave.all.calls   <- TRUE  # Should non-imputed sites be in the output
##                                      # VCF file
##     prefs$ref.alt.by.parent <- FALSE # Should the reference and alternate be
##                                      # switched in the output so that parent 1
##                                      # is always reference and parent 2 is
##                                      # always alternate

##     ## Logistic parameters
##     prefs$quiet             <- FALSE
##     prefs$cores             <- 4
##     prefs$parallel          <- TRUE
##     prefs$write             <- TRUE
##     prefs$out.file          <- ""

##     prefs  # implicit return
## }


##' Checking preferences for the correct variable type
##'
##'
##' @title
##' @param prefs
##' @return
##' @author Jason Vander Woude and Nathan Ryder
ValidatePreferences <- function(prefs) {
    ## TODO(Jason): check for NA logicals
    if (!inherits(prefs, "prefs")) {
        stop("prefs must be of class 'prefs'")
    }
    if (length(prefs$parents) != 2) {
        stop("exactly 2 parents must be specified")
    }
    if (!is.logical(prefs$recomb.double)) {
        stop("'recomb.double' must be of type logical")
    }
    if (!(0 <= prefs$read.err && prefs$read.err < 1) ||
        !(0 <= prefs$genotype.err && prefs$genotype.err < 1)) {
        stop("'error' values must be between 0 and 1")
    }
    if (!is.numeric(prefs$recomb.dist) || !(prefs$recomb.dist > 0)) {
        stop("recombination distance ('recomb.dist') must be a number greater than 0")
    }
    if (prefs$states != length(prefs$parents) + 1) {
        stop("illegal number of states")
    }
    if (!is.logical(prefs$quiet)) {
        stop("'quiet' must be of type logical")
    }
    if (!is.logical(prefs$parallel)) {
        stop("'parallel' must be of type logical")
    }
    if (!is.numeric(prefs$cores) ||
        !(prefs$cores >= 1) ||
        ceiling(prefs$cores) != (prefs$cores)) {
        stop("'cores' should be an integer greater than or equal to 1")
    }
    if (!is.logical(prefs$use.only.ad) || is.na(prefs$use.only.ad)) {
        stop("'use.only.ad' must be a non-NA logical")
    }
    if (!is.logical(prefs$write) || is.na(prefs$write)) {
        stop("'write' must be a non-NA logical")
    }
}


## Progress monitor code from https://stackoverflow.com/questions/27726134/
## how-to-track-progress-in-mclapply-in-r-in-parallel-package
## TODO(Jason): don't use prefs$fifo, but instead try to use a fifo variable
## in the progress.env environment
ProgressMonitor <- function(env) {
    local({
        f <- fifo(tempfile(), open="w+b", blocking=T)
        if (!prefs$parallel) {  # don't fork if running serially
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


MakeProgress <- function(prefs) {
    writeBin(1/prefs$n.jobs, prefs$fifo)  # update the progress bar info
    if (!prefs$parallel) {  # if running in serial mode
        prefs$prog.env$progress <- PrintProgress(prefs$fifo, prefs$prog.env$progress)
    }  # else the forked process handles this
}


WriteVCF <- function(vcf, f) {
    con <- file(f)
    writeLines(make.vcf.lines(vcf), con)
    close(con)
}


ad.to.str <- function(ad.layer) {
    text <- paste0(ad.layer, collapse=",")
    #text <- gsub("NA", ".", text)
}

gt.to.str <- function(gt.layer) {
    text <- paste0(gt.layer, collapse="/")
    text <- gsub("NA", ".", text)
}

ad.to.dp.str <- function(ad.layer) {
    text <- as.character(sum(ad.layer))
    #text <- gsub("NA", ".", text)
}

make.vcf.lines <- function(vcf) {
    ad.str <- apply(vcf$AD, 1:2, ad.to.str)
    gt.str <- apply(vcf$GT, 1:2, gt.to.str)
    dp.str <- apply(vcf$AD, 1:2, ad.to.dp.str)

    combined <- abind(gt.str, ad.str, dp.str, along=3)
    data.strings <- apply(combined, 1:2, paste0, collapse=":")

    names <- str.split(vcf$header[length(vcf$header)], "\t")
    prefix.strings <- vcf$variants[, 1:which(names=="INFO")]

    content.strings <- cbind(prefix.strings, "GT:AD:DP", data.strings)
    all.strings <- rbind(names, content.strings)

    lines <- apply(all.strings, 1, paste0, collapse="\t")

    c(vcf$header[1],
               "##LaByRInth=<ID=Imputation,\"Version=1.0,Description=\"Code can be found at github.com/Dordt-Statistics-Research/LaByRInth\">",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               lines) # implicit return
}


CodifyGT <- function(gt) {
    ## gt is a vector of genotypes
    switch(type,
           mean = 1,
           median = 2,
           trimmed = 3)
}


AnalyzeImputationsRDS <- function(imp, orig, mask, N=30000) {
    four.d.data <- abind(imp$GT, orig$GT, mask$GT, along=4)

    ## df <- data.frame(num=rep(NA_integer_, N), txt=rep("", N),  # as many cols as you need
    ##              stringsAsFactors=FALSE)          # you don't know levels yet

    find.masked.sites <- function(slice) {
        correct <- 1
        partial <- 2
        skipped <- 3
        wrong   <- 4
        ## Assume orignal VCF file had no partial calls (i.e. "./1" or "0/.")
        i <- slice[,1]  # imp
        o <- slice[,2]  # orig
        m <- slice[,3]  # mask
        if (all(is.na(m)) && !any(is.na(o))) {
            ## i <- sort(slice$imp[is.na(slice$imp)]
            if (all(is.na(i))) {
                skipped
            } else if (any(is.na(i))) {
                non.na.i <- i[!is.na(i)]  # will be length 1 since i is length 2
                if (non.na.i %in% o) {
                    partial
                } else {
                    wrong
                }
            } else if (all(sort(i) == sort(o))) {
                correct
            } else {
                wrong
            }
        } else {
            0
        }
    }

    find.depth <- function(strip) {
        if (strip[1]) {
            sum(strip[2:length(strip)])
        } else {
            0
        }
    }

    masked.sites <- apply(four.d.data, 1:2, find.masked.sites)
    print("found masked sites")

    mask.w.ad <- abind(masked.sites, orig$AD, along=3)

    depths <- as.numeric(apply(mask.w.ad, 1:2, find.depth))  # matrix to vector
    qualities <- as.numeric(masked.sites)  # matrix to vector

    relevant <- which(qualities!=0)

    depths <- depths[relevant]
    qualities <- qualities[relevant]
    browser()
    df <- data.frame(depth=depths,
                     quality=factor(qualities,
                                    levels=1:4,
                                    labels = c("correct",
                                               "partial",
                                               "skipped",
                                               "wrong")))

    df
}
