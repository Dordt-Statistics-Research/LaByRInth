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

#    vcf$field.data <- vcf$variants[, 1:which(colnames(vcf$variants)=="FORMAT")]
#    colnames(vcf$field.data) <- str.split(header, "\t")


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
                    #vcf$GT.coded <- CodifyGT(ret.val.gt)
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
