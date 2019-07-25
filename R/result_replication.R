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


get.datasets <- function() {c("Lakin-Fuller", "HincII", "RsaI", "IBM-RIL")}

get.all.parents <- function() {
    list("Lakin-Fuller"  = c("LAKIN",      "FULLER"),
         "HincII"        = c("HincII_B73", "HincII_CG"),
         "RsaI"          = c("RsaI_B73",   "RsaI_CG"),
         "IBM-RIL"       = c("B73",        "M17"))
}

## I don't actually know what generation IBM-RIL is, so I'm just saying 7
get.breed.schemes <- function() {
    schemes <- c("F5", "F2", "F2", "F7")
    names(schemes) <- get.datasets()
    schemes
}

## minimum call depths found with trial and error to get ~1% called sites masked
get.min.mask.depths <- function() {
    min.mask.depths <- c(7, 90, 12, 11)
    names(min.mask.depths) <- get.datasets()
    min.mask.depths
}

## directory associated with extra data files for package
extdata.dir <- function() {system.file("extdata", package = "labyrinth", mustWork = TRUE)}


original.file <- function(dataset.name) {
    file.name <- c("Lakin-Fuller" = "LakinFuller_GBSv2_20170509.vcf.gz",
                   "HincII"       = "HincII_ordered.vcf.gz",
                   "RsaI"         = "RsaI_ordered.vcf.gz",
                   "IBM-RIL"      = "IBM-RIL_ordered.vcf.gz"
      )[dataset.name]
    paste0(extdata.dir(), "/datasets/original_files/", file.name)
}

## for the masked and filtered datasets
dataset.dir <- function(dataset.name) {
    paste0(extdata.dir(), "/datasets/", dataset.name, "/")
}

prep.dir <- function(dataset.name, base.dir) {
    paste0(base.dir, "/", dataset.name, "/")
    ## base.dir
}

imputation.dir <- function(dataset.name, base.dir, algorithm) {
    if (algorithm %in% c("LaByRInth", "LB-Impute")) {
        paste0(prep.dir(dataset.name, base.dir), "/", algorithm, "/")
        ## prep.dir(dataset.name, base.dir)
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}

filtered.file <- function(dataset.name, base.dir) {
    paste0(prep.dir(dataset.name, base.dir), "/", dataset.name, "_filtered.vcf.gz")
}

masked.file <- function(dataset.name, base.dir, mask.ID = 1) {
    paste0(prep.dir(dataset.name, base.dir), "/", dataset.name, "_masked_", mask.ID, ".vcf.gz")
}

uncompressed.masked.file <- function(dataset.name, base.dir, mask.ID = 1) {
    paste0(prep.dir(dataset.name, base.dir), "/", dataset.name, "_masked_", mask.ID, ".vcf")
}

imputed.parents.file <- function(dataset.name, base.dir, algorithm, config = 1) {
    path <- paste0(imputation.dir(dataset.name, base.dir, algorithm), "/",
                   dataset.name, "_", algorithm, "_")
    if (algorithm=="LaByRInth") {
        paste0(path, "parental-imputation-result_", config, ".rds")
    } else if (algorithm=="LB-Impute") {
        paste0(path, "imputed-parents_", config, ".vcf")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}

imputed.progeny.file <- function(dataset.name, base.dir, algorithm,  parental.config = 1, progeny.config = 1) {
    if (algorithm=="LaByRInth") {
        paste0(imputation.dir(dataset.name, base.dir, algorithm), "/imputed-progeny_", parental.config, "_", progeny.config, ".vcf.gz")
    } else if (algorithm=="LB-Impute") {
        paste0(imputation.dir(dataset.name, base.dir, algorithm), "/imputed-progeny_", parental.config, "_", progeny.config, ".vcf")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}


##' Reproduce the prepared data of the LaByRInth publication
##'
##' This function will filter and mask one of the four data files from the
##' LaByRInth publication so that the results can be verified. The filtered and
##' masked result files will be placed in the specified output directory.
##'
##' @param dataset one of the following strings: "Lakin-Fuller", "HincII",
##'        "RsaI", or "IBM-RIL"
##' @param output.dir a directory where all subdirectories and filtered and
##'        masked files will be placed.
##' @return NULL
##' @author Jason Vander Woude
##' @export
LabyrinthPreparePublicationData <- function(dataset, output.dir) {

    ## Ensure dataset is a valid option
    if (! dataset %in% get.datasets()) {
        stop("dataset must be one of these strings: '", paste0(get.datasets(),
        collapse="', '"), "'")
    }

    ## Ensure the original dataset file stil exists in the package
    if (! file.exists(original.file(dataset))) {
        stop("The ", dataset,
             " dataset is missing from the labyrinth package. Try re-installing the package.")
    }

    ## Create the outpur directory if needed
    if (! dir.exists(output.dir)) {
        dir.create(prep.dir(dataset, output.dir), recursive=TRUE)
    }

    ## Filter the datasets
    f.file <- filtered.file(dataset, output.dir)
    if (! dir.exists(dirname(f.file))) {
        dir.create(dirname(f.file), recursive=TRUE)
    }
    display(0, "Filtering the ", dataset, " dataset\n")
    if (file.exists(f.file)) {
        display(1, f.file, " already exists, so ", dataset,
                  " will not be filtered again.\n")
    } else {
        LabyrinthFilter(vcf      = original.file(dataset),
                        out.file = f.file,
                        parents  = get.all.parents()[[dataset]],
                        require.hom.poly = FALSE)
    }

    ## Mask the datasets
    m.file <- masked.file(dataset, output.dir)
    if (! dir.exists(dirname(m.file))) {
        dir.create(dirname(m.file), recursive=TRUE)
    }
    display(0, "Masking the ", dataset, " dataset\n")
    if (file.exists(m.file)) {
        display(1, m.file, " already exists, so ", dataset,
                " will not be masked again.\n")
    } else {
        set.seed(0)  # for repeatability
        LabyrinthMask(vcf       = filtered.file(dataset, output.dir),
                      parents   = get.all.parents()[[dataset]],
                      out.file  = m.file,
                      depth     = get.min.mask.depths()[dataset],
                      lik.ratio = 100,
                      rerr      = 0.05)
    }

    NULL
}


##' Impute the masked publication data
##'
##' First run LabyrinthPreparePublicationData,
##'
##' @param dataset one of the following strings: "Lakin-Fuller", "HincII",
##'        "RsaI", or "IBM-RIL"
##' @param output.dir a directory where all subdirectories and filtered and
##'        masked files will be placed. MUST BE THE SAME AS THE output.dir
##'        PARAMETER USED WHEN RUNNING LabyrinthPreparePublicationData.
##' @param parallel Logical indicating if imputation should be run in parallel
##'        or serial.
##' @param cores Numeric indicating how many sub-processes should be spawned if
##'        running in parallel.
##' @return NULL
##' @author Jason Vander Woude
##' @export
LabyrinthImputePublicationData <- function(dataset, output.dir, parallel=FALSE, cores=1) {

    geno.errs       <- c(0.001, 0.01, 0.1)
    parent.het      <- 0.01
    parents         <- get.all.parents()[[dataset]]
    breed.scheme    <- get.breed.schemes()[dataset]

    display(0, "The ", dataset, " dataset will be imputed with ", length(geno.errs), " different parameter configurations. This could take a few hours. If you are not using Windows, you can set the parallel argument to TRUE and and the cores argument to the number CPUs on your machine to run this in parallel.")

    ## Ensure dataset is a valid option
    if (! dataset %in% get.datasets()) {
        stop("dataset must be one of these strings: '", paste0(get.datasets(),
        collapse="', '"), "'")
    }

    ## Ensure the original dataset file stil exists in the package
    if (! file.exists(original.file(dataset))) {
        stop("The ", dataset,
             " dataset is missing from the labyrinth package. Try re-installing the package.")
    }

    ## Create the outpur directory if needed
    if (! dir.exists(output.dir)) {
        dir.create(imputed.dir(dataset), recursive=TRUE)
    }

    ## Check if masked datasets exist
    m.file <- masked.file(dataset, output.dir)
    if (! file.exists(m.file)) {
        stop("Masked file is missing. Run LabyrinthPreparePublicationData first.")
    }

    ## Ensure imputed datasets can be saved, creating the directory if needed
    example.parent.file <- imputed.parents.file(dataset, output.dir, "LaByRInth")
    example.progeny.file <- imputed.progeny.file(dataset, output.dir, "LaByRInth")
    if (! dir.exists(dirname(example.parent.file))) {
        dir.create(dirname(example.parent.file), recursive=TRUE)
    }
    if (! dir.exists(dirname(example.progeny.file))) {
        dir.create(dirname(example.progeny.file), recursive=TRUE)
    }

    ## Impute the datasets
    for (i in seq_along(geno.errs)) {
        display(0, "Beginning full imputation of the dataset ", dataset,
                " with genotype error ", geno.errs[i])

        par.file        <- imputed.parents.file(dataset,
                                                output.dir,
                                                "LaByRInth",
                                                config=i)
        out.file        <- imputed.progeny.file(dataset,
                                                output.dir,
                                                "LaByRInth",
                                                parental.config=i,
                                                progeny.config=1)

        geno.err        <- geno.errs[i]

        if (! file.exists(par.file)) {
            LabyrinthImputeParents(vcf               = m.file,
                                   parents           = parents,
                                   breed.scheme      = breed.scheme,
                                   out.file          = par.file,
                                   geno.err          = geno.err,
                                   parent.het        = parent.het,
                                   parallel          = parallel,
                                   cores             = cores)
        } else {
            display(1, "Parental file ", par.file, " already exists and will be used")
        }

        if (! file.exists(out.file)) {
            LabyrinthImputeProgeny(parental          = readRDS(par.file),
                                   out.file          = out.file,
                                   use.fwd.bkwd      = use.fwd.bkwd,
                                   calc.posteriors   = use.fwd.bkwd,
                                   viterbi.threshold = NA,  # irrelevant,
                                   parallel          = parallel,
                                   cores             = cores)
        } else {
            display(1, "Imputed file ", out.file, " already exists and will be used")
        }
    }

    NULL
}


##' Mask the most confident calls in a VCF file for validation
##'
##' Compute the liklihood of the the reads for each of the three genotypes
##' (homozygous reference, heterozygous, and homozygous alternate) and return
##' the ratio of the the greatest liklihood over the sum of the other
##' two. This gives an indication of likly it is that the inferred call is
##' the true genotype.

##' @param vcf File path or vcfR object to impute.
##' @param parents Character vector with names of the two parents.
##' @param out.file File path of output file.
##' @param depth Minimum read depth of masked sites.
##' @param lik.ratio Minimum likelihood ratio of masked sits.
##' @param rerr Estimated read error rate of the sequencing process (used in
##'        determining the likelihood of the reads).
##' @return A vcfR object with all sites removed that meet the masking criteria.
##' @author Jason Vander Woude
LabyrinthMask <- function(vcf, parents, out.file, depth=0, lik.ratio=100, rerr=0.05) {

    ## file verification
    if (file.exists(out.file)) {
        stop("Output file already exists; please choose another name\n")
    }

    if (! verify.file.extension(out.file, ".vcf.gz")) {
        stop("Output file name must end with '.vcf.gz'\n")
    }

    if (! verify.file.dir.exists(out.file)) {
        stop("Directory of the outuput file does not exist; please create it\n")
    }

    if (rerr > 0.5) {
        stop("The value of rerr should be in the range [0, 0.5]")
    }

    ## Compute the liklihood of the the reads for each of the three genotypes
    ## (homozygous reference, heterozygous, and homozygous alternate) and return
    ## the ratio of the the greatest liklihood over the sum of the other
    ## two. This gives an indication of likly it is that the inferred call is
    ## the true genotype.
    get.liklihood.ratio <- function(depth.a, depth.b, rerr) {
        if (length(depth.a) > 1)
            stop("length of depth.a must be 1")
        if (length(depth.b) > 1)
            stop("length of depth.b must be 1")

        pa <- (1-rerr)^depth.a * (rerr)^depth.b
        pb <- (rerr)^depth.a * (1-rerr)^depth.b
        ph <- (0.5)^depth.a * (0.5)^depth.b

        if (depth.a == 0) {
            if (depth.b == 0)
                0.5
            else
                pb / (pa + ph)
        } else if (depth.b == 0) {
            pa / (pb + ph)
        } else {
            ph / (pa + pb)
        }
    }

    ## vcf load code
    if (! inherits(vcf, "vcfR")) {
        vcf <- vcfR::read.vcfR(vcf, verbose=F)
    }

    ad.arr <- get.ad.array(vcf)

    ## the third dimension of ad.arr has length two. The first value is the
    ## number of reference reads, and the second is the number of alternate
    ## reads. These are summed to get a total number of reads.
    depths <- apply(ad.arr, 1:2, sum)
    liklihood.ratios <- apply(ad.arr, 1:2, function(reads) {
        get.liklihood.ratio(reads[1], reads[2], rerr)
    })

    parental         <- matrix(rep(getSAMPLES(vcf) %in% parents, nrow(ad.arr)),
                               nrow=nrow(ad.arr),
                               byrow=TRUE)
    sufficient.depth <- depths >= depth
    sufficient.ratio <- liklihood.ratios >= lik.ratio
    mask <- sufficient.depth & sufficient.ratio & !parental

    n.parent.reads <- sum(depths[parental])
    n.reads <- sum(depths) - n.parent.reads
    n.total <- nrow(depths) * (ncol(depths) - 2) # total number of non-parent sites
    n.called <- sum(depths != 0 & !parental)
    n.masked <- sum(mask)

    message(" * Sites are considered progeny/marker pairs")

    message(" * ", n.called, " sites in the progeny are called of ", n.total, " total sites (",
            round(n.called / n.total, 3)*100, "%)")

    message(" * There were ", n.reads,
            " total reads in the progeny across all sites: an average depth of ",
            round(n.reads / n.total, 3), " reads per progeny site")

    message(" * ", n.masked, " sites in the progeny will be masked of ", n.total, " total sites (",
            round(n.masked / n.total, 3)*100, "%)")

    message(" * ", n.masked, " sites in the progeny will be masked of ", n.called, " called sites (",
            round(n.masked / n.called, 3)*100, "%)")
    message("")


    new.ad <- getAD(vcf)
    new.ad[mask] <- "0,0"

    new.gt <- getGT(vcf)
    new.gt[mask] <- "./."

    new.dp <- depths
    new.dp[mask] <- "0"


    concat <- matrix(paste0(new.gt, ":", new.ad, ":", new.dp), nrow=nrow(new.ad))
    vcf@gt <- cbind("GT:AD:DP", concat)

    colnames(vcf@gt) <- c("FORMAT", colnames(new.ad))
    rownames(vcf@gt) <- rownames(new.ad)

    write.vcf(vcf, out.file)
    invisible(vcf)  # implicit return
}


## original <- read.vcfR(filtered.file(dataset))
## masked <- read.vcfR(in.file)

## analysis.df <- do.call(rbind, lapply(seq_along(geno.errs), function(i) {
##     parental.result <- readRDS(parental.file(dataset, config=i))
##     imputed <- read.vcfR(imputed.file(dataset, parental.config=i, progeny.config=1))
##     good.imputed <- LabyrinthUncall(vcf = imputed, min.posterior = 0.9)
##     really.good.imputed <- LabyrinthUncall(vcf = imputed, min.posterior = 0.99)

##     rbind(
##         LBImputeAnalyze(original,
##                         masked,
##                         imputed,
##                         parental.result$read.err,
##                         parental.result$geno.err,
##                         NA,
##                         "everything"),

##         LBImputeAnalyze(original,
##                         masked,
##                         good.imputed,
##                         parental.result$read.err,
##                         parental.result$geno.err,
##                         NA,
##                         "> 90% posterior"),

##         LBImputeAnalyze(original,
##                         masked,
##                         really.good.imputed,
##                         parental.result$read.err,
##                         parental.result$geno.err,
##                         NA,
##                         "> 99% posterior")
##     )
## }))

## write.csv(analysis.df, analysis.file(dataset), row.names = FALSE)
