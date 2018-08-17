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


################################################################################
############################# LIBRARY REQUIREMENTS #############################
################################################################################


require(vcfR, quietly=T)
require(abind, quietly=T)
require(digest, quietly=T)





LabyrinthAnalyze <- function(orig, masked, imputed) {

    gt.o <- getGT(orig)
    gt.m <- getGT(masked)
    gt.i <- getGT(imputed)

    ## vector.indices
    masked.sites <-
        (gt.o != "./." & gt.o != ".|.") &
        (gt.m == "./." | gt.m == ".|.")

    same <-
        ((gt.o == "0/0" | gt.o == "0|0") & (gt.i == "0/0" | gt.i == "0|0")) |
        ((gt.o == "0/1" | gt.o == "0|1" | gt.o == "1/0" | gt.o == "1|0") &
         (gt.i == "0/1" | gt.i == "0|1" | gt.i == "1/0" | gt.i == "1|0")) |
        ((gt.o == "1/1" | gt.o == "1|1") & (gt.i == "1/1" | gt.i == "1|1"))

    skipped <-
        (gt.i == "./.")

    data.frame(
        n.sites   = prod(dim(gt.o)),
        n.masked  = n.masked <- sum(masked.sites),
        n.same    = n.same <- sum(masked.sites & same),
        n.skipped = n.skipped <- sum(masked.sites & skipped),
        n.wrong   = n.masked - n.same - n.skipped,
        accuracy  = n.same / (n.masked - n.skipped),
        quality   = n.same / n.masked
    )
}


LabyrinthAnalyzePerSNP <- function(orig, masked, imputed) {

    gt.o <- getGT(orig)
    gt.m <- getGT(masked)
    gt.i <- getGT(imputed)

    ## vector.indices
    masked.sites <-
        (gt.o != "./." & gt.o != ".|.") &
        (gt.m == "./." | gt.m == ".|.")

    same <-
        ((gt.o == "0/0" | gt.o == "0|0") & (gt.i == "0/0" | gt.i == "0|0")) |
        ((gt.o == "0/1" | gt.o == "0|1" | gt.o == "1/0" | gt.o == "1|0") &
         (gt.i == "0/1" | gt.i == "0|1" | gt.i == "1/0" | gt.i == "1|0")) |
        ((gt.o == "1/1" | gt.o == "1|1") & (gt.i == "1/1" | gt.i == "1|1"))

    skipped <-
        (gt.i == "./.")

    snp.ids <- getID(orig)
    marker.chroms <- getCHROM(orig)
    per.snp.df <- do.call(rbind, lapply(seq_along(snp.ids), function(i) {
        id        <- snp.ids[i]
        chrom     <- marker.chroms[i]

        n.masked  <- sum(masked.sites[id, ])
        n.same    <- sum(masked.sites[id, ] & same[id, ])
        n.skipped <- sum(masked.sites[id, ] & skipped[id, ])

        data.frame(
            snp       = id,
            chrom     = chrom,
            n.masked  = n.masked,
            n.same    = n.same,
            n.skipped = n.skipped,
            n.wrong   = n.masked - n.same - n.skipped,
            accuracy  = n.same / (n.masked - n.skipped),
            quality   = n.same / n.masked
        )
    }))

    class(per.snp.df) <- c("LaByRInthAnalysis", class(per.snp.df))
    per.snp.df
}


per.chrom.df <- function(df) {
    if (! inherits(df, "LaByRInthAnalysis"))
        stop("df must be of class LaByRInthAnalysis")

    do.call(rbind, lapply(unique(df$chrom), function(chrom) {

        rows <- df$chrom == chrom

        data.frame(
            chrom     = chrom,
            n.snps    = sum(rows),
            n.masked  = n.masked <- sum(df$n.masked[rows]),
            n.same    = n.same <- sum(df$n.same[rows]),
            n.skipped = n.skipped <- sum(df$n.skipped[rows]),
            n.wrong   = sum(df$n.wrong[rows]),
            accuracy  = n.same / (n.masked - n.skipped),
            quality   = n.same / n.masked
        )
    }))
}


LabyrinthAnalyzePosteriors <- function(orig, masked, imputed, resolution=10) {

    gt.o <- getGT(orig)
    gt.m <- getGT(masked)
    gt.i <- getGT(imputed)
    gp   <- getGP(imputed)

    ## vector.indices
    masked.sites <-
        (gt.o != "./." & gt.o != ".|.") &
        (gt.m == "./." | gt.m == ".|.")

    same <-
        ((gt.o == "0/0" | gt.o == "0|0") & (gt.i == "0/0" | gt.i == "0|0")) |
        ((gt.o == "0/1" | gt.o == "0|1" | gt.o == "1/0" | gt.o == "1|0") &
         (gt.i == "0/1" | gt.i == "0|1" | gt.i == "1/0" | gt.i == "1|0")) |
        ((gt.o == "1/1" | gt.o == "1|1") & (gt.i == "1/1" | gt.i == "1|1"))

    skipped <-
        (gt.i == "./.")

    wrong <- !same & !skipped

    ## this is not efficient as most of the posteriors don't need to
    ## be computed
    posteriors <- apply(gp, 1:2, function(gp.string) {
        phred.to.prob(max(as.numeric(str.split(gp.string, ","))))
    })

    regions <- seq(0, 1, length.out = resolution+1)
    res <- sapply.pairs(regions, function(low, high) {
        which.sites <- masked.sites & posteriors > low & posteriors <= high & !skipped
        c(n.sites   = n.sites   <- sum(which.sites),
          ## n.correct = n.correct <- sum(which.sites & same),
          ## n.wrong   = sum(which.sites & wrong),
          prop.correct = sum(which.sites & same) / n.sites)
    })
    colnames(res) <- sapply.pairs(regions, function(low, high) {paste0(low, "-", high)})

    df <- data.frame(n.sites       = as.integer(res[1, ]),
                     min.posterior = regions[-length(regions)],
                     prop.correct  = res[2, ],
                     max.posterior = regions[-1])

    rownames(df) <- NULL

    df  # implicit return
}


LabyrinthPlotPosteriors <- function(orig, masked, imputed, resolution=10) {

    gt.o <- getGT(orig)
    gt.m <- getGT(masked)
    gt.i <- getGT(imputed)
    gp   <- getGP(imputed)

    ## vector.indices
    masked.sites <-
        (gt.o != "./." & gt.o != ".|.") &
        (gt.m == "./." | gt.m == ".|.")

    same <-
        ((gt.o == "0/0" | gt.o == "0|0") & (gt.i == "0/0" | gt.i == "0|0")) |
        ((gt.o == "0/1" | gt.o == "0|1" | gt.o == "1/0" | gt.o == "1|0") &
         (gt.i == "0/1" | gt.i == "0|1" | gt.i == "1/0" | gt.i == "1|0")) |
        ((gt.o == "1/1" | gt.o == "1|1") & (gt.i == "1/1" | gt.i == "1|1"))

    skipped <-
        (gt.i == "./.")

    wrong <- !same & !skipped


    temp <- gp[masked.sites]
    temp <- sapply(temp, function(gp.string) {
        phred.to.prob(max(as.numeric(str.split(gp.string, ","))))
    })
    posteriors <- rep(0, length(masked.sites))
    posteriors[masked.sites] <- temp


    regions <- seq(0, 1, length.out = resolution+1)
    res <- sapply.pairs(regions, function(low, high) {
        which.sites <- masked.sites & posteriors > low & posteriors <= high & !skipped
        c(n.sites   = n.sites   <- sum(which.sites),
          ## n.correct = n.correct <- sum(which.sites & same),
          ## n.wrong   = sum(which.sites & wrong),
          prop.correct = sum(which.sites & same) / n.sites)
    })
    colnames(res) <- sapply.pairs(regions, function(low, high) {paste0(low, "-", high)})

    df <- data.frame(n.sites       = as.integer(res[1, ]),
                     min.posterior = regions[-length(regions)],
                     prop.correct  = res[2, ],
                     max.posterior = regions[-1],
                     call.coverage = 1 - cumsum(as.integer(res[1, ])) / sum(masked.sites))

    rownames(df) <- NULL

    df  # implicit return
}


## LBImputeAnalyze <- function(orig, masked, imputed, out.file,
##                             read.err, geno.err, recomb.dist) {

##     gt.o <- getGT(orig)
##     gt.m <- getGT(masked)
##     gt.i <- getGT(imputed)

##     ## vector.indices
##     masked.sites <-
##         (gt.o != "./." & gt.o != ".|.") &
##         (gt.m == "./." | gt.m == ".|.")

##     same <-
##         ((gt.o == "0/0" | gt.o == "0|0") & (gt.i == "0/0" | gt.i == "0|0")) |
##         ((gt.o == "0/1" | gt.o == "0|1" | gt.o == "1/0" | gt.o == "1|0") &
##          (gt.i == "0/1" | gt.i == "0|1" | gt.i == "1/0" | gt.i == "1|0")) |
##         ((gt.o == "1/1" | gt.o == "1|1") & (gt.i == "1/1" | gt.i == "1|1"))

##     skipped <-
##         (gt.i == "./.")


##     n.masked <- sum(masked.sites)
##     n.same <- sum(masked.sites & same)
##     n.skipped <- sum(masked.sites & skipped)

##     result <- data.frame(
##         read.err = read.err,
##         geno.err = geno.err,
##         recomb.dist = recomb.dist,
##         n.masked = n.masked,
##         n.correct = n.same,
##         n.skipped = n.skipped,
##         n.wrong = n.masked - n.same - n.skipped
##     )

##     saveRDS(result, out.file)
##     result  # implicit return
## }
