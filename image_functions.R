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



## This file contains functions related to generating images of the vcf files
## and imputation results



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


vcf.to.pgm <- function(vcf.file, parents, out.file) {
    require(parallel)
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


temp <- function() {
    chroms <-  c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "UN")
    for (chrom in chroms) {
        make.diff(paste0("LaByRInth_", chrom, ".pgm"),
                  paste0("LB-Impute_", chrom, ".pgm"),
                  paste0("Diff_",      chrom, ".ppm"))
    }
}
