## Copyright 2018 Jason Vander Woude
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


##               __          ____        ____  _____
##              / /         / __ \      / __ \/_  _/
##             / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __
##            / /   / _  \/ _  _/\ \/ / _  _/ / / /  | / /_  __/ /_/ /
##           / /___/ /_/ / /_\ \  \  / / \ \_/ /_/ /||/ / / / / __  /
##          /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/
##
##                  L O W - C O V E R A G E   B I A L L E L I C
##                    R - P A C K A G E   I M P U T A T I O N
##

version <- function() {"3.2.14"}

################################################################################
###################### PRIMARY TOP-LEVEL FUNCTIONS FOR USER ####################
################################################################################


##' Remove bad data prior to imputing
##'
##' Remove all sites that are not biallelic and possibly sites which are not
##' homozygous within and polymorphic between.
##'
##' @param vcf File path or vcfR object to impute.
##' @param out.file File path of output file.
##' @param parents Character vector with names of the two parents.
##' @param require.hom.poly Logical indicating if the sites that are kept are
##'        required to by homozygous within and polymorphic between (i.e. should
##'        the parents be homozygous for different alleles at every kept site
##'        according to the genotype (gt) field in the vcf). The implementation
##'        of LB-Impute that LaByRInth is based on required this, but this
##'        version of LaByRInth does not. This functionality is mostly included
##'        for legacy purposes and should generally be set to false.
##' @return A vcfR object with all sites removed that don't meet the specified
##'         criteria.
##' @examples
##' input <- system.file(
##'     "extdata",
##'     "vcf-files",
##'     "original-lakin-fuller-sample.vcf",
##'     package = "labyrinth",
##'     mustWork = TRUE)
##' output <- tempfile(fileext="-filtered-result.vcf.gz")
##' print(output)
##' result <- LabyrinthFilter(
##'     vcf = input,
##'     out.file = output,
##'     parents = c("LAKIN", "FULLER"),
##'     require.hom.poly = TRUE)
##' @author Jason Vander Woude
##' @export
LabyrinthFilter <- function(vcf, out.file, parents, require.hom.poly=FALSE) {

    ## begin timer
    total.timer <- new.timer()


    ## check classes of arguments
    if (class(vcf) != "character" && class(vcf) != "vcfR")
        stop("vcf must be of class character or vcfR\n")
    if (class(parents) != "character")
        stop("parents must be of class character\n")
    if (class(out.file) != "character")
        stop("out.file must be of class character\n")
    if (class(require.hom.poly) != "logical")
        stop("require.hom.poly must be of class logical\n")


    ## file verification
    if (file.exists(out.file))
        stop("Output file already exists; please choose another name\n")
    if (! verify.file.extension(out.file, ".vcf.gz"))
        stop("Output file name must end with '.vcf.gz'\n")
    if (! verify.file.dir.exists(out.file))
        stop("Directory of the output file does not exist; please create it\n")



    ## load vcf if needed
    if (inherits(vcf, "character")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- vcfR::read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }


    display(0, "Checking if parents are in the vcf\n")
    if (! all(parents %in% getSAMPLES(vcf))) {
        text <- "Some or all parents are not in the vcf"
        if (all(toupper(parents) %in% toupper(getSAMPLES(vcf))))
            text <- paste0(text, "; try checking the capitalization of the names")
        stop(text)
    }
    display(0, "Checking for sites that are not biallelic")
    biallelic <- vcfR::is.biallelic(vcf)
    if (any(!biallelic)) {
        display(1, "The following sites are not biallelic and will be removed:")
        chroms <- vcfR::getCHROM(vcf)[which(! biallelic)]
        positions <- vcfR::getPOS(vcf)[which(! biallelic)]
        refs <- vcf@fix[! biallelic, "REF"]
        alts <- vcf@fix[! biallelic, "ALT"]
        for (i in seq_along(chroms)) {
            display(2, "CHR:", chroms[i],
                    "\tPOS:", positions[i],
                    "\tREF:", refs[i],
                    "\tALT:", alts[i])
        }
    }
    message("")
    mask <- biallelic

    if (require.hom.poly) {
        display(0, "Checking for sites where parents are not ",
                   "homozygous within and polymorphic between")
        hom.poly <- parent.hom.and.poly(vcf, parents)
        if (any(! hom.poly)) {
            display(1, "The parents are not homozygous within and polymorphic ",
                       "between at the following sites which will be removed:")
            chroms <- vcfR::getCHROM(vcf)[! hom.poly]
            positions <- vcfR::getPOS(vcf)[! hom.poly]
            gt <- getGT(vcf)[! hom.poly, parents, drop=FALSE]
            for (i in seq_along(chroms)) {
                display(2, "CHR:", chroms[i],
                        "\tPOS:", positions[i],
                        paste0("\t", parents[1], ":"), gt[i, 1],
                        paste0("\t", parents[2], ":"), gt[i, 2])
            }
        }
        message("")
        mask <- mask & hom.poly
    }

    any.removed <- any(mask==FALSE)

    if (any.removed) {
        display(0, "Removing sites and saving vcf\n")
        vcf@fix <- vcf@fix[mask, ]
        vcf@gt <- vcf@gt[mask, ]
        vcfR::write.vcf(vcf, out.file)
    } else {
        display(0, "All sites okay; copying original vcf file\n")
        vcfR::write.vcf(vcf, out.file)
    }


    display(0, paste("Summary:", sum(!mask),
                     "of", length(mask), "sites removed\n"))


    display(0, "LaByRInth filtering completed in ", total.timer(), "\n")
    invisible(vcf)
}


##' Parental imputation
##'
##' Imputes both of the parents by maximizing likelihood of data across
##' all progeny at each site.
##'
##' @param vcf File path or vcfR object to impute.
##' @param out.file File path of output file.
##' @param parents Character vector with names of the two parents.
##' @param breed.scheme Population type (e.g. F2, F5, F1BC1).
##' @param geno.err Estimate of proportion of calls genotyped incorrectly.
##' @param parent.het Estimate of proportion of sites in each parent that are
##'        truly heterozygous.
##' @param parallel Logical indicating if imputation should be run in parallel
##'        or serial.
##' @param cores Numeric indicating how many sub-processes should be spawned if
##'        running in parallel.
##' @return A vcfR object with parental data replaced by imputation results.
##' @examples
##' input <- system.file(
##'     "extdata",
##'     "vcf-files",
##'     "filtered-lakin-fuller-sample.vcf",
##'     package = "labyrinth",
##'     mustWork = TRUE)
##' output <- tempfile(fileext="-parental-result.rds")
##' print(output)
##' result <- LabyrinthImputeParents(
##'     vcf = input,
##'     out.file = output,
##'     parents = c("LAKIN", "FULLER"),
##'     breed.scheme = "F5",
##'     progeny.het = 0.5^4, # expected heterozygosity of F5 population
##'     geno.err = 0.015,
##'     parent.het = 0.005,
##'     parallel = FALSE,
##'     cores = 1)
##' @author Jason Vander Woude
##' @export
LabyrinthImputeParents <- function (vcf, out.file, parents, breed.scheme,
                                    progeny.het=NULL, read.err=NULL,
                                    geno.err=0.01, parent.het=0.01,
                                    parallel=TRUE, cores=4) {
    ## begin timer
    total.timer <- new.timer()


    ## check classes of arguments
    if (class(vcf) != "character" && class(vcf) != "vcfR")
        stop("vcf must be of class character or vcfR\n")
    if (class(parents) != "character")
        stop("parents must be of class character\n")
    if (class(breed.scheme) != "character")
        stop("breed.scheme must be of class character\n")
    if (class(progeny.het) != "numeric" && !is.null(progeny.het))
        stop("progeny.het must be of class numeric (or NULL)\n")
    if (class(read.err) != "numeric" && !is.null(read.err))
        stop("read.err must be of class numeric (or NULL)\n")
    if (class(out.file) != "character")
        stop("out.file must be of class character\n")
    if (class(geno.err) != "numeric")
        stop("geno.err must be of class numeric\n")
    if (class(parent.het) != "numeric")
        stop("parent.het must be of class numeric\n")
    if (class(parallel) != "logical")
        stop("parallel must be of class logical\n")
    if (class(cores) != "numeric")
        stop("cores must be of class numeric\n")


    ## file verification
    if (file.exists(out.file))
        stop("Output file already exists; please choose another name\n")
    if (! verify.file.extension(out.file, ".rds"))
        stop("Output file name must end with '.rds'\n")
    if (! verify.file.dir.exists(out.file))
        stop("Directory of the output file does not exist; please create it\n")


    ## parameter verification
    if (cores == 1 && parallel)
        stop("Cores cannot be 1 if parallel is true\n")
    if (cores > 1 && !parallel)
        stop("Cores cannot be greater than 1 if parallel is false\n")
    if (cores < 1)
        stop("Cores cannot be less than 1\n")
    if (length(parents) != 2)
        stop("parents must have length 2\n")
    if ((is.null(progeny.het) && is.null(read.err)) ||
        (!is.null(progeny.het) && !is.null(read.err)))
        stop("Exactly one of progeny.het and read.err must be specified. ",
             "If read.err is specified, that will be used as the probability ",
             "of a sequencing error. If progeny.het is specified, the ",
             "probability of a sequencing error will be estimated from the ",
             "population characteristics.")
    if (!is.null(read.err) && read.err > 0.5)
        stop("read.err must be between 0 and 0.5 (inclusive)")
    if (!is.null(read.err) && read.err > 0.1)
        warning("read.err set to ", round(read.err*100, 2), "%. ",
                "This is unusually high. Proceeding anyway.")


    ## Variable initialization and function loading
    timer <- new.timer()
    display(0, "Loading data specific to ", breed.scheme, " populations")
    ## site.pair.transition.probs variable will be loaded
    tryCatch({
        trans.file <- system.file("extdata",
                                  "transition-probs",
                                  paste0(breed.scheme, ".R"),
                                  package = "labyrinth",
                                  mustWork = TRUE)
        source(trans.file)
    }, error = function(e) {
        stop(paste("Cannot find file ", trans.file, "; the specified ",
                   "breeding scheme may not be supported by default. ",
                   "Additional breeding schemes can be added; for instructions, ",
                   "see the inst/extdata/multi-model-symbolics.sage file in the ",
                   "labyrinth package directory."))
    })
    display(1, "Completed in ", timer(), "\n")



    ## load vcf if needed
    if (inherits(vcf, "character")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- vcfR::read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }



    ## vcf requirement verification code
    timer <- new.timer()
    display(0, "Verifying requirements of vcf")
    if (! all(parents %in% getSAMPLES(vcf))) {
        text <- "Some or all parents are not in the vcf"
        if (all(toupper(parents) %in% toupper(getSAMPLES(vcf))))
            text <- paste0(text, "; try checking the capitalization of the names")
        stop(text)
    }

    non.biallelic <- which(! vcfR::is.biallelic(vcf))
    if (length(non.biallelic) != 0) {
        display(0, "The following sites are not biallelic:")
        chroms <- vcfR::getCHROM(vcf)[non.biallelic]
        positions <- vcfR::getPOS(vcf)[non.biallelic]
        for (i in seq_along(chroms)) {
            display(1, "\tCHR:", chroms[i], "\tPOS:", positions[i])
        }
    }

    if (length(non.biallelic) != 0) {
        stop("The vcf must be filtered before imputing; please run LabyrinthFilter first\n")
    }
    display(1, "Completed in ", timer(), "\n")


    timer <- new.timer()
    display(0, "Generating and caching variables")
    result <- list()

    ## save parameters in result
    result$vcf <- vcf
    result$parents <- parents
    result$breed.scheme <- breed.scheme
    result$geno.err <- geno.err
    result$parent.het <- parent.het
    result$site.pair.transition.probs <- site.pair.transition.probs

    ## save additional information so that it doesn't have to be computed again
    result$snp.chroms   <- snp.chroms   <- vcfR::getCHROM(vcf)
    result$u.chroms     <- u.chroms     <- unique(snp.chroms)
    result$ad           <- ad           <- get.ad.array(vcf)
    result$sample.names <- sample.names <- getSAMPLES(vcf)
    result$marker.names <- marker.names <- vcfR::getID(vcf)

    display(1, "Completed in ", timer(), "\n")



    if (is.null(read.err)) {
        timer <- new.timer()
        display(0, "Estimating erroneous read rate from heterozygosity")
        read.err <- estimate.read.err(ad, progeny.het)
        display(1, "Error rate estimated to be ", round(read.err*100, 2), "%")
        display(1, "Completed in ", timer(), "\n")
    }
    result$read.err <- read.err



    ## emission probability code
    timer <- new.timer()
    display(0, "Generating emission probabilities")
    emission.structure <- get.emission.structures(ad,
                                                  sample.names,
                                                  marker.names,
                                                  progeny.het,
                                                  read.err,
                                                  geno.err)
    result$emission.structure <- emission.structure
    display(1, "Completed in ", timer(), "\n")



    ## parental imputation and recombination rate estimates
    timer <- new.timer()
    display(0, "Imputing parents and estimating recombination rates")
    parental.results <- determine.parents.and.recombs(emission.structure,
                                                      parents,
                                                      snp.chroms,
                                                      marker.names,
                                                      sample.names,
                                                      parent.het,
                                                      parallel,
                                                      cores)
    result$parent.models <- parental.results
    class(result) <- c("parental.imputation", class(result))
    saveRDS(result, out.file)
    display(1, "Completed in ", timer(), "\n")


    display(0, "LaByRInth parental imputation completed in ", total.timer(), "\n")
    invisible(result)
}


##' Offspring imputation
##'
##' Imputes all of the progeny (offspring).
##'
##' @param parental An object of class parental.imputation (a result of the
##'        LabyrinthParentalImputation function) which will contain information
##'        about the population including the parental imputation information.
##' @param out.file File path of output file.
##' @param parallel Logical indicating if imputation should be run in parallel
##'        or serial.
##' @param cores Numeric indicating how many sub-processes should be spawned if
##'        running in parallel.
##' @return A vcfR object with both parents and progeny imputed.
##' @examples
##' input <- system.file(
##'     "extdata",
##'     "vcf-files",
##'     "parental-lakin-fuller-sample.rds",
##'     package = "labyrinth",
##'     mustWork = TRUE)
##' output <- tempfile(fileext="-progeny-result.vcf.gz")
##' print(output)
##' result <- LabyrinthImputeProgeny(
##'     parental = input,
##'     out.file = output,
##'     parallel = FALSE,
##'     cores = 1)
##' @author Jason Vander Woude
##' @export
LabyrinthImputeProgeny <- function (parental, out.file,
                                    parallel=TRUE, cores=4) {
    ## begin timer
    total.timer <- new.timer()


    ## check classes of arguments
    if (class(parental) != "character" && class(parental) != "parental.imputation")
        stop("vcf must be of class character or parental.imputation\n")
    if (class(out.file) != "character")
        stop("out.file must be of class character\n")
    if (class(parallel) != "logical")
        stop("parallel must be of class logical\n")
    if (class(cores) != "numeric")
        stop("cores must be of class numeric\n")


    ## load parental file if needed
    if (inherits(parental, "character")) {
        if (! file.exists(parental)) {
            stop(paste("file", parental, "does not exist"))
        }

        timer <- new.timer()
        display(0, "Loading parental imputation file")

        tryCatch({
            parental <- readRDS(parental)
        }, error = function(e) {
            stop(paste("If specifying a parental file, the file must be",
                       "a file generated by the function",
                       "LabyrinthImputeProgeny"))
        })
        display(1, "Completed in ", timer(), "\n")
    }


    ## file verification
    if (file.exists(out.file))
        stop("Output file already exists; please choose another name\n")
    if (! verify.file.extension(out.file, ".vcf.gz"))
        stop("Output file name must end with '.vcf.gz'\n")
    if (! verify.file.dir.exists(out.file))
        stop("Directory of the output file does not exist; please create it\n")
    if (cores == 1 && parallel)
        stop("Cores cannot be 1 if parallel is true\n")
    if (cores > 1 && !parallel)
        stop("Cores cannot be greater than 1 if parallel is false\n")
    if (cores < 1)
        stop("Cores cannot be less than 1\n")


    timer <- new.timer()
    display(0, "Restoring variables from parental imputation result")
    fwd.bkwd.threshold <- 0  # use LabyrinthQualityControl for other thresholds
    vcf <- parental$vcf
    parents <- parental$parents
    parent.models <- parental$parent.models
    breed.scheme <- parental$breed.scheme
    read.err <- parental$read.err
    geno.err <- parental$geno.err
    site.pair.transition.probs <- parental$site.pair.transition.probs
    snp.chroms <- parental$snp.chroms
    u.chroms <- parental$u.chroms
    ad <- parental$ad
    sample.names <- parental$sample.names
    marker.names <- parental$marker.names
    emission.structure <- parental$emission.structure
    display(1, "Completed in ", timer(), "\n")



    ## transition probability code
    timer <- new.timer()
    display(0, "Computing progeny transition probabilities")
    transition.structures <- get.transition.structures(snp.chroms,
                                                       marker.names,
                                                       parent.models,
                                                       site.pair.transition.probs,
                                                       parallel,
                                                       cores)
    gc()
    display(1, "Completed in ", timer(), "\n")



    ## imputation code
    timer <- new.timer()
    display(0, "Imputing missing sites")
    imp.res <- impute(parents,
                      emission.structure,
                      transition.structures,
                      parent.models,
                      sample.names,
                      snp.chroms,
                      parallel,
                      cores,
                      fwd.bkwd.threshold)
    gc()
    display(1, "Completed in ", timer(), "\n")

    ## new vcf creation code
    timer <- new.timer()
    display(0, "Creating new vcf with imputed data")
    vcf <- update.vcf(vcf, imp.res$gt, imp.res$posteriors)
    vcfR::write.vcf(vcf, out.file)
    display(1, "Completed in ", timer(), "\n")


    display(0, "LaByRInth progeny imputation completed in ", total.timer(), "\n")
    invisible(vcf) ## implicit return
}


##' Ensure quality of imputed calls
##'
##' LaByRInth will impute every site for every member of the population and
##' include a marginal posterior, so any calls with too low of a posterior
##' probability will be removed (i.e. a site for a population member that was
##' called with low probability of correctness will be "uncalled").
##'
##' @param vcf File path or vcfR object to impute.
##' @param out.file File path of output file.
##' @param min.posterior Numeric specifying the minimum probability for any call
##'        that is kept.
##' @param parallel Logical indicating if imputation should be run in
##'        parallel or serial.
##' @param cores Numeric indicating how many sub-processes should be
##'        spawned if running in parallel.
##' @return A vcfR object with all low probability sites removed.
##' @examples
##' input <- system.file(
##'     "extdata",
##'     "vcf-files",
##'     "progeny-lakin-fuller-sample.vcf",
##'     package = "labyrinth",
##'     mustWork = TRUE)
##' output <- tempfile(fileext="-quality-result.vcf.gz")
##' print(output)
##' result <- LabyrinthQualityControl(
##'     vcf = input,
##'     out.file = output,
##'     min.posterior = 0.8,
##'     parallel = FALSE,
##'     cores = 1)
##' @author Jason Vander Woude
##' @export
LabyrinthQualityControl <- function(vcf, out.file, min.posterior,
                                      parallel=TRUE, cores=4) {
    ## begin timer
    total.timer <- new.timer()


    ## check classes of arguments
    if (class(vcf) != "character" && class(vcf) != "vcfR")
        stop("vcf must be of class character or vcfR\n")
    if (class(out.file) != "character")
        stop("out.file must be of class character\n")
    if (class(parallel) != "logical")
        stop("parallel must be of class logical\n")
    if (class(cores) != "numeric")
        stop("cores must be of class numeric\n")
    if (class(min.posterior) != "numeric")
        stop("min.posterior must be of class numeric\n")
    if (class(vcf) != "character" && class(vcf) != "vcfR")
        stop("vcf must be of class character or vcfR\n")


    ## file verification
    if (file.exists(out.file))
        stop("Output file already exists; please choose another name\n")
    if (! verify.file.extension(out.file, ".vcf.gz"))
        stop("Output file name must end with '.vcf.gz'\n")
    if (! verify.file.dir.exists(out.file))
        stop("Directory of the output file does not exist; please create it\n")


    ## load vcf if needed
    if (inherits(vcf, "character")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- vcfR::read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }


    timer <- new.timer()
    display(0, "Masking sites with posterior probability below ", min.posterior)

    listapply <- get.lapply(parallel, cores)
    min.phred <- prob.to.phred(min.posterior)

    formats <- vcf@gt[ , "FORMAT"]

    new.data <- t(sapply(seq_along(formats), function(i) {
        format.components <- str.split(formats[i], ":")
        gp.index <- format.components == "GP"
        gt.index <- format.components == "GT"

        ## The '-1' is to ignore the 'FORMAT' column of the vcfR@gt
        ## matrix. vcf.entry will be something such as "0/0:1,0:8.6,0.6,0.0"
        ## which is GT:AD:GP.
        sapply(vcf@gt[i, -1], function(vcf.entry) {

            components <- str.split(vcf.entry, ":")
            posteriors <- as.numeric(str.split(components[gp.index], ","))

            if (max(posteriors) < min.phred)
                components[gt.index] <- "./."

            paste0(components, collapse=":")
        })
    }))

    rownames <- rownames(vcf@gt)
    colnames <- colnames(vcf@gt)

    vcf@gt <- cbind(formats, new.data)
    rownames(vcf@gt) <- rownames
    colnames(vcf@gt) <- colnames

    vcfR::write.vcf(vcf, out.file)
    display(1, "Completed in ", timer(), "\n")


    display(0, "LaByRInth quality control completed in ", total.timer(), "\n")
    invisible(vcf) ## implicit return
}


##' Run entire imputation including filtering and quality control
##'
##' Run all four other functions (LabyrinthFilter, LabyrinthImputeParents,
##' LabyrinthImputeProgeny, and LabyrinthQualityControl) in sequence and save
##' the intermediate files to a temporary directory. Running this function will
##' produce equivalent results as running all four functions separately; the
##' only difference will be where the files are saved. This is provided simply
##' for user convenience.
##'
##' @param vcf File path or vcfR object to impute.
##' @param out.file File path of output file.
##' @param parents Character vector with names of the two parents.
##' @param breed.scheme Population type (e.g. F2, F5, F1BC1).
##' @param min.posterior Numeric specifying the minimum probability for any call
##'        that is kept.
##' @param geno.err Estimate of proportion of calls genotyped incorrectly.
##' @param parent.het Estimate of proportion of sites in each parent that are
##'        truly heterozygous.
##' @param require.hom.poly Logical indicating if the sites that are kept are
##'        required to by homozygous within and polymorphic between (i.e. should
##'        the parents be homozygous for different alleles at every kept site
##'        according to the genotype (gt) field in the vcf). The implementation
##'        of LB-Impute that LaByRInth is based on required this, but this
##'        version of LaByRInth does not. This functionality is mostly included
##'        for legacy purposes and should generally be set to false.
##' @param parallel Logical indicating if imputation should be run in parallel
##'        or serial.
##' @param cores Numeric indicating how many sub-processes should be spawned if
##'        running in parallel.
##' @return A vcfR object with all sites removed that don't meet the filtering
##'         criteria, both parents and progeny imputed and low probability sites
##'         removed.
##' @examples
##' input <- system.file(
##'     "extdata",
##'     "vcf-files",
##'     "original-lakin-fuller-sample.vcf",
##'     package = "labyrinth",
##'     mustWork = TRUE)
##' output <- tempfile(fileext="-result.vcf.gz")
##' result <- LabyrinthImpute(
##'     vcf = input,
##'     out.file = output,
##'     parents = c("LAKIN", "FULLER"),
##'     breed.scheme = "F5",
##'     progeny.het = 0.5^4, # expected heterozygosity of F5 population
##'     min.posterior = 0.8,
##'     geno.err = 0.015,
##'     parent.het = 0.005,
##'     require.hom.poly = TRUE,
##'     parallel = FALSE,
##'     cores = 1)
##' @author Jason Vander Woude
##' @export
LabyrinthImpute <- function(vcf, out.file, parents, breed.scheme, progeny.het=NULL,
                            read.err=NULL, min.posterior, geno.err=0.01, parent.het=0.01,
                            require.hom.poly=FALSE, parallel=TRUE, cores=4) {

    total.timer <- new.timer()

    dir <- tempdir()  # create temporary directory
    filtered.file    <- paste0(dir, "/filtered.vcf.gz")
    parental.file    <- paste0(dir, "/parental.rds")
    all.imputed.file <- paste0(dir, "/all-imputed.vcf.gz")

    filtered.result <- LabyrinthFilter(
        vcf = vcf,
        out.file = filtered.file,
        parents = parents,
        require.hom.poly = require.hom.poly
    )

    parental.result <- LabyrinthImputeParents(
        vcf = filtered.result,
        out.file = parental.file,
        parents = parents,
        breed.scheme = breed.scheme,
        geno.err = geno.err,
        parent.het = parent.het,
        parallel = parallel,
        cores = cores
    )

    all.imputed.result <- LabyrinthImputeProgeny(
        parental = parental.result,
        out.file = all.imputed.file,
        parallel = parallel,
        cores = cores
    )

    quality.result <- LabyrinthQualityControl(
        vcf = all.imputed.result,
        out.file = out.file,
        min.posterior = min.posterior,
        parallel = parallel,
        cores = cores
    )

    display(0, "The intermediate files generated were saved as follows:")
    display(1, "Filtered population: ", filtered.file)
    display(1, "Parental imputation: ", parental.file)
    display(1, "Full imputation: ", all.imputed.file, "\n")

    display(0, "LaByRInth full process completed in ", total.timer(), "\n")
    invisible(quality.result) ## implicit return
}


##' Calculate expected heterozygosity in a population
##'
##' Calculate expected heterozygosity in a population. Supported breeding
##' schemes are "F1BC1" and "F<x>" where <x> is replaced by any (possibly
##' multi-digit) number.
##'
##' @param breed.scheme Supported breeding scheme
##' @return Expected heterozygousity proportion
##' @author Jason Vander Woude
##' @export
LabyrinthCalcProgenyHet <- function(breed.scheme) {
    RIL.pattern <- "^F[0-9]*$" # ^ and $ match beginning and end
                               # F is literal
                               # [0-9]* matches any number of numeric characters
    F1BC1.str <- "F1BC1"

    if (class(breed.scheme) != "character") {
        stop("breed.scheme must be of class character\n")
    }
    if (length(breed.scheme) != 1) {
        stop("LabyrinthCalcProgenyHet is not vectorized. breed.scheme must have length 1\n")
    }

    if (breed.scheme == F1BC1.str) {
        0.5 # F1BC1 populations should have 50% heterozygosity
    } else if (length(grep(RIL.pattern, breed.scheme))!=0) {
        gen <- as.integer(substr(breed.scheme, 2, nchar(breed.scheme)))
        # F1 has 100% heterozygosity, and the heterozygosity should be cut in
        # half each generation after
        0.5^(gen - 1)
    } else {
        stop("Unsupported breed.scheme")
    }
}





################################################################################
######################## SECONDARY LABYRINTH FUNCTIONS #########################
################################################################################


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
get.emission.structures <- function(ad, sample.names, marker.names,
                                    progeny.het, read.err, geno.err) {

    perc.het <- progeny.het
    perc.hom.ref <- perc.hom.alt <- (1 - perc.het) / 2

    k <- choose(ad[ , , 1] + ad[ , , 2], ad[ , , 1])  # n choose k

    hom.ref.read.emm <- k * (1-read.err)^ad[ , , 1] * (read.err)^ad[ , , 2]
    hom.alt.read.emm <- k * (1-read.err)^ad[ , , 2] * (read.err)^ad[ , , 1]
    het.read.emm     <- k *        (0.5)^ad[ , , 1] *      (0.5)^ad[ , , 2]


    ## by accounting for the probability that a given site was erroneously
    ## assembled we are in essence more heavily weighting the heterozygous
    ## emission probabilities
    err.emm <- perc.hom.ref * hom.ref.read.emm +
               perc.hom.alt * hom.alt.read.emm +
               perc.het * het.read.emm

    err.hom.ref.read.emm <- (1 - geno.err) * hom.ref.read.emm + (geno.err) * err.emm
    err.hom.alt.read.emm <- (1 - geno.err) * hom.alt.read.emm + (geno.err) * err.emm
    err.het.read.emm     <- (1 - geno.err) * het.read.emm     + (geno.err) * err.emm

    ## read probs given states
    ## first dimension is the markers/loci/SNPs
    ## second dimension is the members of the population
    ## third dimension is the state (hom ref, het type I, het type II, hom alt)
    rpgs <- abind::abind(
        err.hom.ref.read.emm,
        err.het.read.emm,
        err.het.read.emm,
        err.hom.alt.read.emm,

        along=3
    )

    dimnames(rpgs) <- list(marker.names,
                           sample.names,
                           c("hom.ref", "het.I", "het.II", "hom.alt"))

    rpgs  # implicit return
}


get.transition.structures <- function(snp.chroms, marker.names,
                                      parent.models,
                                      site.pair.transition.probs,
                                      parallel, cores) {

    u.chroms <- unique(snp.chroms)
    listapply <- get.lapply(parallel, cores)

    trans.probs <- function(chrom) {

        ## parental.model <- parental.results$parental.models[[chrom]]$path
        ## recombs <- parental.results$recombs[[chrom]]  # vector of functions
        model <- parent.models[[chrom]]$model
        recombs <- parent.models[[chrom]]$recombs

        ## site.pair.transition.probs was loaded when a file in
        ## new-transition-probs was sourced
        trans.matrices <- lapply(seq_along(recombs), function(i) {
            ## for marker i and i+1, get the transition matrix based on the
            ## parental states at both sites. This "matrix" is actually a
            ## function of the recombination probability which can be called
            ## with recombs[i] to get a 4x4 matrix
            mat <- site.pair.transition.probs[[model[i]]][[model[i+1]]](recombs[i])
            ## normalize so that every rowSum is 0 since these are probabilities
            mat <- mat / rep(rowSums(mat), 4)
            mat[is.nan(mat)] <- 0.25
            mat
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


## return the total probabilities of each state at each site. Adapted from
## https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
## forward and backward probabilities are not technically correct because they
## are frequently being normalized. This is done to decrease the numerical
## instability that results when the probabilities become very low
fwd.bkwd <- function(emm, trans, temp) {
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


## TODO(Jason): the initial probabilities are not correct here or in the
## fwd.bkwd

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


## The following labeling schemes are used. Both the Viterbi algorithm and the
## forward backward algorithm will keep track of the following four state types
## when determining the best path
##
## hom.ref     het.I     het.II     hom.alt
## 0|0         0|1       1|0        1|1
##
## When the Viterbi algorithm is used, a single optimal sequence of states is
## determined which thus includes haplotypic information and can thus
## distinguish between het.I and het.II. Thus the calls are phased and the
## vertical bar is used in the VCF file as shown above. If the forward backward
## algorithm is used, then it makes more sense to combine the het.I and het.II
## types because the VCF format does not allow for specifying posterior
## probabilities for both het.I and het.II. Then the following scheme is used:
##
## hom.ref     het     hom.alt
## 0/0         0/1     1/1
impute <- function(parents, emm.structures, trans.structures, parent.models,
                   sample.names, snp.chroms, parallel, cores,
                   fwd.bkwd.threshold) {

    listapply <- get.lapply(parallel, cores)
    u.chroms <- unique(snp.chroms)

    ## Given two probabilities of odd recombinations (r1 between sites x and y;
    ## r2 between sites y and z) what is the probability of an odd number of
    ## recombinations between sites x and z
    combine.recombs <- function(r1, r2) {r1*(1-r2) + r2*(1-r1)}

    impute.sample.chrom <- function(sample, chrom) {

        n.sites <- sum(snp.chroms==chrom)  # Boolean addition
        parent.paths <- extract.each.parent(
            parent.models[[chrom]]$model
        )

        emm <- t(emm.structures[chrom==snp.chroms, sample, ])
        trans <- trans.structures[[chrom]]

        res <- list()
        res$posteriors <- NULL

        if (sample %in% parents) {
            posterior.mat <- rbind("hom.ref" = rep(1/3, n.sites),
                                   "het.I"   = rep(1/6, n.sites),
                                   "het.II"  = rep(1/6, n.sites),
                                   "hom.alt" = rep(1/3, n.sites))
        } else {
            posterior.mat <- fwd.bkwd(emm, trans, sample=="HincII_F2-01" && as.character(chrom)=="2")
        }

        posteriors <- rbind("hom.ref" = posterior.mat["hom.ref", ],
                            "het"     = posterior.mat["het.I", ] + posterior.mat["het.II", ],
                            "hom.alt" = posterior.mat["hom.alt", ])

        phred.scaled <- apply(posteriors, 1:2,  prob.to.phred)

        res$posteriors <- apply(phred.scaled, 2, paste0, collapse=",")

        if (sample == parents[1]) {
            res$gt <- c("0|0", "0|1", "1|0", "1|1")[parent.paths[[1]]]
        } else if (sample == parents[2]) {
            res$gt <- c("0|0", "0|1", "1|0", "1|1")[parent.paths[[2]]]
        } else {
            res$gt <- apply(posteriors, 2, function(states) {
                c("0/0","0/1","1/1")[which.max(states)]
            })
        }

        which.remove <- apply(posteriors, 2, function(states) {
            max(states) < fwd.bkwd.threshold
        })

        res$gt[which.remove] <- "./."

        gc()
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

        ret.val  # implicit return

    })

    ret.val <- list()
    ret.val$posteriors <- do.call(rbind, get.posteriors(imputed.chroms))
    ret.val$gt <- do.call(rbind, get.listed.gt(imputed.chroms))


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code


    ret.val  # implicit return
}


determine.parents.and.recombs <- function(emm.structure, parents, snp.chroms,
                                          marker.names, sample.names,
                                          parent.het, parallel, cores) {

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

    ## transitions will be a list with an entry for every chromosome. Each
    ## chromosome entry will be a list with an entry for every adjacent pair of
    ## markers in the chromosome (one fewer entries than the number of markers
    ## in the chromosome). Each chromosome-markerpair entry is a list
    ## containing two elements: each a 16x16 matrix indicating pairs of parental
    ## states, the first containing the maximized log-liklihoods, and the second
    ## containing the recombination values at which they were maximized.
    transitions <- lapply(u.chroms, function(chrom) {
        indices <- which(snp.chroms == chrom)  # which indices correspond with this chrom
        indices <- indices[-length(indices)]  # remove the last element

        ret.val.2  <- listapply(indices, function(i) {

            j <- i+1

            ## snp.i will be constant along the third dimension and snp.j will be
            ## constant along the second dimension. By binding additional copies in
            ## this way we can replicate mathematical matrix multiplication with
            ## element-by-element multiplication allowing us to do all SNPs at once
            ## without needing an apply function which can be slower
            snp.i <- abind::abind(rpgs[i, , ], rpgs[i, , ], rpgs[i, , ], rpgs[i, , ], along=3)
            snp.j <- abind::abind(rpgs[j, , ], rpgs[j, , ], rpgs[j, , ], rpgs[j, , ], along=3)
            snp.j <- aperm(snp.j, c(1,3,2))

            rpgsp <- snp.i * snp.j  # read probs given state pair
            essential.layers <- rowSums(rpgs[i, , ]) != 4 & rowSums(rpgs[j, , ]) != 4

            dimnames(rpgsp) <- list(sample.names[ !which.parents ],
                                    c("hom.ref", "het.I", "het.II", "hom.alt"),
                                    c("hom.ref", "het.I", "het.II", "hom.alt"))


            rm(snp.i, snp.j)
            gc()

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
                                    upper=1,
                                    control=list(ndeps=1e-2,  # step size
                                                 fnscale=-1))
                    recomb.val.mat[x,y] <- result$par
                    log.liklihood.mat[x,y] <- result$value

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

            if (all(log.liklihood.mat == -Inf)) {
                ## TODO(Jason): This is a temporary fix due to numerical
                ## instabilities. If the probabilities get too small, taking the
                ## exponent above returns 0 and then taking the log again
                ## returns negative infinity, and if all the probabilities are
                ## negative infinity, then all paths will have -Inf log
                ## probability and the wrong one will likely be returned. There
                ## should really be a check added so that if the probability of
                ## the path is -Inf, it is computed again another way
                log.liklihood.mat[2:15, 2:15] <- 0
            }

            list(logliklihoods = log.liklihood.mat,
                 recombs = recomb.val.mat)

        })

        ret.val.2

    })


    ## Progress bar code
    ## -------------------------------------------------------------------------
    close(thefifo)
    ## -------------------------------------------------------------------------
    ## Progress bar code

    names(transitions) <- u.chroms




    ## Construct emission likelihood matrix. There are 16 possible parental states at each SNP
    p <- 1 - parent.het  # probability of a site in parents being homozygous
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
                log(snp.depths[2, (p2+1)])# +
               log.penalty[i+1]
        })
        names(result) <- NULL
        result  # implicit return
    })

    ## TODO(Jason): verify that at least 2 sites are in the chromosome before
    ## running viterbi


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

    result <- lapply(u.chroms, function(chrom) {
        list(model   = parental.models[[chrom]]$path,
             recombs = recombs[[chrom]])
    })
    names(result) <- u.chroms
    result

}





################################################################################
######################### TERTIARY LABYRINTH FUNCTIONS #########################
################################################################################


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
    par.mat <- getGT(vcf)[ , parents]
    apply(par.mat, 1, function(row) {
        row <- gsub("\\|", "/", row)  # replace '|' with '/'
        (all(row==c("0/0", "1/1")) ||
         all(row==c("1/1", "0/0")) ||
         all(row==c("./.", "0/0")) ||
         all(row==c("./.", "1/1")) ||
         all(row==c("0/0", "./.")) ||
         all(row==c("1/1", "./.")))
    })
}


## modify the vcf in place to add imputed calls in gt
update.vcf <- function(vcf, gt, posteriors=NULL) {
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


get.lapply <- function(parallel, cores=1) {
    if (parallel && cores != 1) {
        if (require(parallel, quietly=T)) {
            function(...) {
                parallel::mclapply(..., mc.preschedule=F, mc.cores=cores)
            }  # implicit return
        } else {
            stop(paste0("In order to utilize parallel computation, the ",
                        "package 'parallel' must be installed"))
        }
    } else {
        lapply  # implicit return
    }
}


## parental.model is a vector of values in range 1-16 inclusive. This indicates
## the state of the two parents for each site. The value is 4 times the value of
## the state of parent 1 plus the state of parent 2. This function returns a
## list where the first element is the vector of states of parent 1 and the
## second is the vector of states of parent 2
extract.each.parent <- function(parental.model) {
    parents <- sapply(parental.model, function(call) {
        call <- call-1         # deal with base 1 indexing
        p1 <- floor(call / 4)  # bit shift i right twice to get two most sig bits
        p2 <- call - p1*4      # two least significant bits
        c(p1+1, p2+1)
    })

    list(
        parents[1, ],
        parents[2, ]
    )
}


estimate.read.err <- function(ad, expected.het) {
    depths <- ad[ , , 1] + ad[ , , 2]
    relevant <- depths > 1  # depth of 0 or 1 gives no info on heterozygosity
    ref <- ad[ , , 1][relevant]
    alt <- ad[ , , 2][relevant]

    h <- expected.het

    ## objective function of r, the probability of an erroneous read at the
    ## haplotype level
    obj.fun <- function(r) {
        sum(log((1-h)/2 * (r^ref * (1-r)^alt + r^alt * (1-r)^ref) + h*0.5^(ref + alt)))
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





################################################################################
############################## VCF DATA ACCESS CODE ############################
################################################################################


getSAMPLES <- function(vcf) {
    ## column names are the samples except the first column which is "FORMAT"
    colnames(vcf@gt)[-1]
}


getAD <- function(vcf) {
    ad <- vcfR::extract.gt(vcf, "AD")
    ad[is.na(ad)] <- "0,0"  # replace NA entries
    ad
}


getGT <- function(vcf) {
    gt <- vcfR::extract.gt(vcf, "GT")
    gt[is.na(gt)] <- "./."  # replace NA entries
    gt
}


getGP <- function(vcf) {
    vcfR::extract.gt(vcf, "GP")
}





################################################################################
############################# PHRED CONVERSION CODE ############################
################################################################################


prob.to.phred <- function(x) {
    ## 1-x is the probability the call is wrong
    ## use min to prevent infinity from getting through
    min(-10*log((1-x), base=10), 100)
}


phred.to.prob <- function(y) {
    -(10^(y / -10) - 1)
}





################################################################################
############################# STRING TO NUMERIC CODE ###########################
################################################################################


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


get.ad.array <- function(vcf) {
    ads <- apply(getAD(vcf), 1:2, ad.to.num)
    abind::abind(ads[1, , ], ads[2, , ], along=3)
}





################################################################################
########################## ABIND HELPER FUNCTIONS CODE #########################
################################################################################


## The abind{i} functions are useful for cleaning the code for do.call(abind,
## some.list) calls because otherwise the 'along' argument has to be appended to
## the 'some.list' argument
abind1 <- function(...) {
    abind::abind(..., along=1)
}


abind2 <- function(...) {
    abind::abind(..., along=2)
}


abind3 <- function(...) {
    abind::abind(..., along=3)
}





################################################################################
########################### GENERAL UTILITY FUNCTIONS ##########################
################################################################################


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


verify.file.extension <- function(file, extension) {
    ## split the extension by periods and remove the leading empty character
    extension.parts <- str.split(extension, "\\.")[-1]
    n.parts <- length(extension.parts)
    file.parts <- str.split(file, "\\.")

    check <- rev(extension.parts) == rev(file.parts)[1:n.parts]
    !any(is.na(check)) && all(check)
}


verify.file.dir.exists <- function(file) {
    dir.exists(dirname(file))
}


## More convenient strsplit if length of vector is 1
str.split <- function(str, sep) {
    if (length(str) != 1) {
        warning("Only splitting the first element of the vector")
    }
    strsplit(str, sep)[[1]]
}





################################################################################
############################# PROGRESS MONITOR CODE ############################
################################################################################


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
