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



## This file contains a number of functions that are used throughout the rest of
## the code


translate <- function(sample) {
    ## make the following replacements
    ## 1 -> 0
    ## 2 -> 1
    ## 4 -> 2
    ## anything else -> 5
    ifelse(is.na(sample), 5, ifelse(sample==1, 0, ifelse(sample==2, 1, ifelse(sample==4, 2, 5))))
}


## TODO(Jason): don't use sink
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


nsites.probs <- function(probs) {
    nrow(probs)
}


nstates.probs <- function(probs) {
    ncol(probs)
}


vcf.to.numeric <- function(str.vec) {
    sapply(str.vec, function(str) {
        if (check.int(str)) {
            as.numeric(str)
        } else if (str == ".") {
            NA_integer_
        } else {
            stop(paste("Attempted to convert invalid VCF string:", str))
        }
    })
}


get.lapply <- function(prefs) {
    serial.lapply <- function(..., mc.preschedule=F, mc.cores=prefs$cores) {
        lapply(...)
    }

    parallel.lapply <- function(..., mc.preschedule=F, mc.cores=prefs$cores) {
        mclapply(..., mc.preschedule=mc.preschedule, mc.cores=mc.cores)
    }

    ## implicit return
    if (prefs$parallel) {
        require(parallel)
        parallel.lapply
    } else {
        serial.lapply
    }
}


get.labyrinth.vcf.header <- function() {
    "##fileformat=VCFv4.0"
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
##    prefix.strings <- vcf$field.data[, 1:which(names=="INFO")]

    content.strings <- cbind(prefix.strings, "GT:AD:DP", data.strings)
    all.strings <- rbind(names, content.strings)

    lines <- apply(all.strings, 1, paste0, collapse="\t")

    c(vcf$header[1],
               "##LaByRInth=<ID=Imputation,\"Version=1.0,Description=\"Code can be found at github.com/Dordt-Statistics-Research/LaByRInth\">",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               lines) # implicit return
}


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


AnalyzeImputationsRDS <- function(imp, orig, mask, N=30000) {
    ## TODO(Jason): without an original mask matrix it is impossible to tell
    ## what percent of masked sites were imputed correctly. But it is possible
    ## to tell the percent on non-NA masked sites which is really what we care about
    four.d.data <- abind(imp$GT, orig$GT, mask$GT, along=4)

    ## df <- data.frame(num=rep(NA_integer_, N), txt=rep("", N),  # as many cols as you need
    ##              stringsAsFactors=FALSE)          # you don't know levels yet

    find.masked.sites <- function(slice) {
        correct  <- 1
        partial  <- 2
        skipped  <- 3
        wrong    <- 4
        unmasked <- 5
        ## Assume orignal VCF file had no partial calls (i.e. "./1" or "0/.")
        i <- slice[,1]  # imp
        o <- slice[,2]  # orig
        m <- slice[,3]  # mask
        if (all(is.na(m)) && !all(is.na(o))) { # all masked are na, no original are na
            ##browser()
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
            ##print(slice)
            unmasked
        }
    }

    find.depth <- function(strip) {
        ##browser()
        if (strip[1]) {
            sum(strip[2:length(strip)])
        } else {
            0
        }
    }

#    browser()
    masked.sites <- apply(four.d.data, 1:2, find.masked.sites)
    print("found masked sites")

    mask.w.ad <- abind(masked.sites, orig$AD, along=3)

    depths <- as.numeric(apply(mask.w.ad, 1:2, find.depth))  # matrix to vector
    qualities <- as.numeric(masked.sites)  # matrix to vector

    relevant <- which(qualities!=5)
#    relevant <- which(qualities!=-1)  # all sites

    depths <- depths[relevant]
    qualities <- qualities[relevant]
    df <- data.frame(depth=depths,
#                     quality=qualities)
                     quality=factor(qualities,
                                    levels=1:4,
                                    labels = c("correct",
                                               "partial",
                                               "skipped",
                                               "wrong")))

    df
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
