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



get.datasets <- function() {
    c("Lakin-Fuller",
      "HincII",
      "RsaI",
      "IBM-RIL",
      "Sim-A-F1BC1")
}

## get.sim.datasets <- function() {
##     c("Lakin-Fuller",
##       "HincII",
##       "RsaI",
##       "IBM-RIL")
## }

## get.real.datasets <- function() {
##     c("Sim-A-F1BC1")
## }

get.masks <- function(dataset) {
    if (dataset %in% c("Lakin-Fuller",
                       "HincII",
                       "RsaI",
                       "IBM-RIL"))
    {
        list(list(ID="all_of_top_1_percent",   mask.prop=0.01,  top.prop=0.01),
             list(ID="all_of_top_3_percent",   mask.prop=0.03,  top.prop=0.03),
             list(ID="all_of_top_5_percent",   mask.prop=0.05,  top.prop=0.05),
             list(ID="all_of_top_10_percent",  mask.prop=0.10,  top.prop=0.10),
             list(ID="all_of_top_50_percent",  mask.prop=0.50,  top.prop=0.50),
             list(ID="all_of_top_75_percent",  mask.prop=0.75,  top.prop=0.75))

    } else if (dataset == "Sim-A-F1BC1") {
        list(list(ID="sampled"))
    } else {
        stop("did not find masks for dataset ", dataset)
    }
}

get.configs <- function() {
    list(list(ID="geno_0.05", geno.err=0.05))
}

get.parents <- function(dataset) {
    list("Lakin-Fuller" = c("LAKIN",      "FULLER"),
         "HincII"       = c("HincII_B73", "HincII_CG"),
         "RsaI"         = c("RsaI_B73",   "RsaI_CG"),
         "IBM-RIL"      = c("B73",        "M17"),
         "Sim-A-F1BC1"  = c("par1",       "par2")
         )[[dataset]]
}

get.fsfhap.contributions <- function(dataset) {
    list("Lakin-Fuller" = c(0.5, 0.5),
         "HincII"       = c(0.5, 0.5),
         "RsaI"         = c(0.5, 0.5),
         "IBM-RIL"      = c(0.5, 0.5),
         "Sim-A-F1BC1"  = c(0.75, 0.25)
         )[[dataset]]
}

get.fsfhap.inbreeding.coef <- function(dataset) {
    1 - get.progeny.het(dataset)
}

## I don't actually know what generation IBM-RIL is, so I'm just saying 7
get.breed.scheme <- function(dataset) {
    c("Lakin-Fuller" = "F5",
      "HincII"       = "F2",
      "RsaI"         = "F2",
      "IBM-RIL"      = "F7",
      "Sim-A-F1BC1"  = "F1BC1"
      )[dataset]
}

get.read.err.rate <- function(dataset) {
    c("Lakin-Fuller" = 0.0138, # 1.38%
      "HincII"       = 0.0033, # 0.33%
      "RsaI"         = 0.0009, # 0.09%
      "IBM-RIL"      = 0.0008, # 0.08%
      "Sim-A-F1BC1"  = 0       # 0%
      )[dataset]
}

get.progeny.het <- function(dataset) {
    LabyrinthCalcProgenyHet(get.breed.scheme(dataset))
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
    paste0(prep.dir(dataset.name, base.dir), "/",
           dataset.name, "_filtered.vcf.gz")
}

masked.file <- function(dataset.name, base.dir, mask.ID) {
    paste0(prep.dir(dataset.name, base.dir), "/",
           dataset.name, "_masked_", mask.ID, ".vcf.gz")
}

fsfhap.dir <- function(dataset.name, base.dir) {
    paste0(prep.dir(dataset.name, base.dir), "/FSF-Hap/")
}

fsfhap.filtered.vcf.file <- function(dataset.name, base.dir) {
    paste0(fsfhap.dir(dataset.name, base.dir), "/",
           dataset.name, "_filtered_no_parents.vcf.gz")
}

fsfhap.masked.vcf.file <- function(dataset.name, base.dir, mask.ID) {
    paste0(fsfhap.dir(dataset.name, base.dir), "/",
           dataset.name, "_masked_", mask.ID, "_no_parents.vcf.gz")
}

fsfhap.pedigree.file <- function(dataset.name, base.dir, mask.ID) {
    paste0(fsfhap.dir(dataset.name, base.dir), "/",
           dataset.name, "_filtered_pedigree.txt")
}

uncompressed.masked.file <- function(dataset.name, base.dir, mask.ID) {
    paste0(prep.dir(dataset.name, base.dir), "/", dataset.name, "_masked_", mask.ID, ".vcf")
}

imputed.parents.file <- function(dataset.name, base.dir, algorithm, mask.ID, config.ID) {
    path <- paste0(imputation.dir(dataset.name, base.dir, algorithm), "/",
                   dataset.name, "_", algorithm, "_")
    if (algorithm=="LaByRInth") {
        paste0(path, "parental-imputation-result_masked_", mask.ID, "_config_", config.ID, ".rds")
    } else if (algorithm=="LB-Impute") {
        paste0(path, "imputed-parents_masked_", mask.ID, "_config_", config.ID, ".vcf")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}

imputed.progeny.file <- function(dataset.name, base.dir, algorithm,
                                 mask.ID, progeny.config.ID) {
        path <- paste0(imputation.dir(dataset.name, base.dir, algorithm), "/",
                       dataset.name, "_", algorithm, "_imputed-progeny_masked_",
                       mask.ID, "_config_", progeny.config.ID)
    if (algorithm=="LaByRInth") {
        paste0(path, ".vcf.gz")
    } else if (algorithm=="LB-Impute") {
        paste0(path, ".vcf")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}


##' Reproduce the prepared data of the LaByRInth publication
##'
##' This function will filter and mask one of the five data files from the
##' LaByRInth publication so that the results can be verified. The filtered and
##' masked result files will be placed in the specified output directory.
##'
##' @param dataset one of the following strings: "Lakin-Fuller", "HincII",
##'        "RsaI", "IBM-RIL", "Sim-A-F1BC1
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
                        parents  = get.parents(dataset),
                        require.hom.poly = FALSE)
    }

    for (index in seq_along(get.masks(dataset))) {
        mask <- get.masks(dataset)[[index]]

        m.file <- masked.file(dataset, output.dir, mask$ID)
        if (! dir.exists(dirname(m.file))) {
            dir.create(dirname(m.file), recursive=TRUE)
        }
        display(0, "Masking the ", dataset, " dataset with ", mask$ID, "\n")
        if (file.exists(m.file)) {
            display(1, m.file, " already exists, so ", dataset,
                    " will not be masked in this way again.\n")
        } else {
            set.seed(index)  # for repeatability
            LabyrinthMask(vcf       = filtered.file(dataset, output.dir),
                          parents   = get.parents(dataset),
                          out.file  = m.file,
                          mask.prop = mask$mask.prop,
                          top.prop = mask$top.prop,
                          rerr      = get.read.err.rate(dataset))
        }
    }

    invisible(NULL)
}


##' Impute the masked publication data
##'
##' First run LabyrinthPreparePublicationData,
##'
##' @param dataset one of the following strings: "Lakin-Fuller", "HincII",
##'        "RsaI", "IBM-RIL", "Sim-A-F1BC1
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
LabyrinthImputePublicationData <- function(dataset, output.dir, parallel=FALSE, cores=1, mask.IDs=NULL) {

    parent.het      <- 0.01
    parents         <- get.parents(dataset)
    breed.scheme    <- get.breed.scheme(dataset)
    progeny.het     <- get.progeny.het(dataset)

    display(0, "The ", dataset, " masked datasets will be imputed with various parameter configurations. This could take a few hours. If you are not using Windows, you can set the parallel argument to TRUE and and the cores argument to the number CPUs on your machine to run this in parallel.")

    ## Ensure dataset is a valid option
    if (! dataset %in% get.datasets()) {
        stop("dataset must be one of these strings: '", paste0(get.datasets(),
        collapse="', '"), "'")
    }

    ## Create the outpur directory if needed
    if (! dir.exists(output.dir)) {
        dir.create(imputed.dir(dataset), recursive=TRUE)
    }

    for (mask in get.masks(dataset)) {

        if (! is.null(mask.IDs)) {
            if (! mask$ID %in% mask.IDs) {
                display(1, "skipping mask ", mask$ID)
                next
            }
        }

        ## Check if masked datasets exist
        m.file <- masked.file(dataset, output.dir, mask$ID)
        if (! file.exists(m.file)) {
            stop("Masked file is missing. Run LabyrinthPreparePublicationData first.")
        }

        ## Ensure imputed datasets can be saved, creating the directory if needed
        example.parent.file <- imputed.parents.file(
            dataset, output.dir, "LaByRInth", mask$ID, "TEMP-CONFIG")
        example.progeny.file <- imputed.progeny.file(
            dataset, output.dir, "LaByRInth", mask$ID, "TEMP-CONFIG")
        if (! dir.exists(dirname(example.parent.file))) {
            dir.create(dirname(example.parent.file), recursive=TRUE)
        }
        if (! dir.exists(dirname(example.progeny.file))) {
            dir.create(dirname(example.progeny.file), recursive=TRUE)
        }

        ## Impute the datasets
        configs <- get.configs()
        for (i in seq_along(configs)) {
            display(0, "Beginning full imputation of the dataset ", dataset,
                    " with genotype error ", configs[[i]]$geno.err, "\n")

            par.file        <- imputed.parents.file(dataset,
                                                    output.dir,
                                                    "LaByRInth",
                                                    mask = mask$ID,
                                                    config = configs[[i]]$ID)
            out.file        <- imputed.progeny.file(dataset,
                                                    output.dir,
                                                    "LaByRInth",
                                                    mask = mask$ID,
                                                    progeny.config = configs[[i]]$ID)

            if (! file.exists(par.file)) {
                LabyrinthImputeParents(vcf               = m.file,
                                       out.file          = par.file,
                                       parents           = parents,
                                       breed.scheme      = breed.scheme,
                                       progeny.het       = progeny.het,
                                       geno.err          = configs[[i]]$geno.err,
                                       parent.het        = parent.het,
                                       parallel          = parallel,
                                       cores             = cores)
            } else {
                display(1, "Parental file ", par.file, " already exists and will be used")
            }

            if (! file.exists(out.file)) {
                LabyrinthImputeProgeny(parental          = readRDS(par.file),
                                       out.file          = out.file,
                                       parallel          = parallel,
                                       cores             = cores)
            } else {
                display(1, "Imputed file ", out.file, " already exists and will be used")
            }
        }
    }

    invisible(NULL)
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
LabyrinthMask <- function(vcf, parents, out.file, mask.prop, top.prop, rerr=0.05) {

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
    if (mask.prop > top.prop) {
        stop("mask.prop must be less than or equal to top.prop")
    }
    if (mask.prop > 1 || top.prop > 1) {
        stop("mask.prop and top.prop must be in the range [0,1]")
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

    ## get array of indices indicating a decreasing ordering where progeny are
    ## indexed first (! parental) and then liklihood ratios are used.
    call.confidence.ordering <- order(!parental,
                                      liklihood.ratios,
                                      decreasing=TRUE)
    ## how many sites must be masked
    progeny.depths <- depths[!parental]
    n.called <- sum(progeny.depths != 0)
    ## how many called sites constitute the top 'top.prop'
    n.top <- ceiling(n.called * top.prop)
    n.masked <- ceiling(n.called * mask.prop)
    top.indices <- call.confidence.ordering[1:n.top]
    ## randomly select mask indices
    if (n.top == n.masked) {
        mask.indices <- top.indices
    } else {
        mask.indices <- sample(top.indices, size=n.masked)
    }

    ## the mask to use
    mask <- sapply(liklihood.ratios, function(x) {FALSE})
    mask[mask.indices] <- TRUE

    if (any(mask & parental)) {
        stop("Error in code")
    }

    ## other metrics
    n.total <- length(liklihood.ratios[!parental]) # total progeny sites in vcf
    masked.depths <- depths[mask]
    min.depth <- min(masked.depths)
    n.progeny.reads <- sum(progeny.depths)
    n.progeny.reads.masked <- sum(masked.depths)
    min.liklihood.ratio <- min(liklihood.ratios[mask])


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

    vcfR::write.vcf(vcf, out.file)

    message(" * Sites are considered taxa/marker pairs")
    message(" * All metrics are with respect to the progeny (parents are ignored)")

    message(" * ", n.called, " sites in the progeny are called of ", n.total, " total sites (",
            round(n.called / n.total, 3)*100, "%)")

    message(" * There were ", n.progeny.reads,
            " total reads in the progeny across all sites: an average depth of ",
            round(n.progeny.reads / n.total, 3), " reads per site")
    message("")
    message(" * ", n.masked, " sites in the progeny will be masked of ", n.total, " total sites (",
            round(n.masked / n.total, 3)*100, "%)")

    message(" * ", n.masked, " sites in the progeny will be masked of ", n.called, " called sites (",
            round(n.masked / n.called, 3)*100, "%)")
    message(" * ", n.progeny.reads.masked, " reads in the progeny will be masked of ",
            n.progeny.reads, " reads in the progeny (",
            round(n.progeny.reads.masked / n.progeny.reads, 3)*100, "%)")
    message(" * There masked file will have ", n.progeny.reads - n.progeny.reads.masked,
            " remaining reads in the progeny across all sites: an average depth of ",
            round((n.progeny.reads - n.progeny.reads.masked) / n.total, 3), " reads per site")
    message(" * All masked sites in the progeny have a call depth of at least ", min.depth)
    message(" * All masked sites in the progeny have a liklihood ratio of at least ", min.liklihood.ratio)
    message("\n")

    invisible(vcf)  # implicit return
}


LabyrinthPrepareFSFHap <- function(dataset, output.dir) {
    ## Ensure dataset is a valid option
    if (! dataset %in% get.datasets()) {
        stop("dataset must be one of these strings: '", paste0(get.datasets(),
        collapse="', '"), "'")
    }

    for (mask in get.masks(dataset)) {

        ## Check if masked datasets exist
        m.file <- masked.file(dataset, output.dir, mask$ID)
        if (! file.exists(m.file)) {
            stop("Masked file is missing. Run LabyrinthPreparePublicationData first.")
        }

        ## Create the filtered output directory if needed
        filtered.vcf.file <- fsfhap.filtered.vcf.file(dataset, output.dir)
        if (! dir.exists(dirname(filtered.vcf.file))) {
            dir.create(dirname(filtered.vcf.file), recursive=TRUE)
        }

        ## Create the masked vcf output directory if needed
        masked.vcf.file <- fsfhap.masked.vcf.file(dataset, output.dir, mask.ID)
        if (! dir.exists(dirname(masked.vcf.file))) {
            dir.create(dirname(masked.vcf.file), recursive=TRUE)
        }

        ## Create the pedigree output directory if needed
        ped.file <- fsfhap.pedigree.file(dataset, output.dir, mask.ID)
        if (! dir.exists(dirname(ped.file))) {
            dir.create(dirname(ped.file), recursive=TRUE)
        }

        ## Check if filtered datasets exist
        f.file <- filtered.file(dataset, output.dir)
        if (! file.exists(m.file)) {
            stop("Filtered file is missing. Run LabyrinthPreparePublicationData first.")
        }

        ## Check if masked datasets exist
        m.file <- masked.file(dataset, output.dir, mask.ID)
        if (! file.exists(m.file)) {
            stop("Masked file is missing. Run LabyrinthPreparePublicationData first.")
        }

        ## save new vcf files with parents removed
        f.vcf <- vcfR::read.vcfR(f.file)
        m.vcf <- vcfR::read.vcfR(m.file)
        parents <- get.parents(dataset)
        taxa <- getSAMPLES(f.vcf)
        is.progeny <- !(taxa %in% parents)
        cols.to.keep <- c(TRUE, is.progeny) # first column is 'FORMAT'
        out.f.vcf <- f.vcf                  # copy vcf
        out.m.vcf <- m.vcf                  # copy vcf

        out.f.vcf@gt <- f.vcf@gt[ , cols.to.keep] # remove parents
        out.m.vcf@gt <- m.vcf@gt[ , cols.to.keep] # remove parents
        vcfR::write.vcf(out.f.vcf, filtered.vcf.file)
        vcfR::write.vcf(out.m.vcf, masked.vcf.file)

        ## save pedigree file
        family <- dataset                   # name of the family/dataset
        progeny <- taxa[is.progeny]         # all progeny, no parents
        contributions <- get.fsfhap.contributions(dataset)
        coef.of.inbreeding <- get.fsfhap.inbreeding.coef(dataset)
        ## all areguments except `progeny` have length 1, so are repeated for every
        ## row of this matrix
        data <- cbind("Family" = family,
                      "Taxa" = progeny,
                      "Parent1" = parents[1],
                      "Parent2" = parents[2],
                      "Contribution1" = contributions[1],
                      "Contribution2" = contributions[2],
                      "F" = coef.of.inbreeding)
        write.table(data, ped.file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
}


##' Analyze all masks and imputation configurations of a dataset
##'
##' Analyze all masks and imputation configurations of a dataset
##'
##' @param dataset one of the following strings: "Lakin-Fuller", "HincII",
##'        "RsaI", "IBM-RIL", "Sim-A-F1BC1"
##' @param output.dir a directory where all subdirectories and filtered and
##'        masked files will be placed. MUST BE THE SAME AS THE output.dir
##'        PARAMETER USED WHEN RUNNING LabyrinthPreparePublicationData.
##' @param parallel Logical indicating if imputation should be run in parallel
##'        or serial.
##' @param cores Numeric indicating how many sub-processes should be spawned if
##'        running in parallel.
##' @return list containing analysis results
##' @author Jason Vander Woude
##' @export
LabyrinthAnalyzePublicationData <- function(dataset, output.dir, algorithm) {
    masks <- get.masks()
    names(masks) <- sapply(masks, `[`, "ID") # apply the index function to get mask IDs
    lapply(masks, function(mask) {

        configs <- get.configs()
        names(configs) <- sapply(configs, `[`, "ID")
        lapply(configs, function(config) {

            tryCatch({
                files <- list(orig    = filtered.file(dataset, output.dir),
                              masked  = masked.file(dataset, output.dir, mask$ID)
                              imputed = imputed.progeny.file(dataset, output.dir,
                                                             algorithm, mask$ID,
                                                             config$ID)
                              )

                vcfs <- lapply(files, vcfR::read.vcfR)

                do.call(LabyrinthAnalyze, vcfs)
            }, warning = function(warning_condition) {
                NULL
            }, error = function(error_condition) {
                NULL
            })
        })
    })
}
