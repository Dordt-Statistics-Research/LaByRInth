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





mimic <- function(vcf, full.out, sample.out, parents, n.gen, read.err,
                  parental.results=NULL, use.50=FALSE, parallel=TRUE, cores=4) {

    total.timer <- new.timer()
    print.labyrinth.sim.header()

    ## vcf load code
    if (! inherits(vcf, "vcfR")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }

    display(5, "BIALLELIC")

    timer <- new.timer()
    display(0, "Generating necessary variables")

    ad <- get.ad.array(vcf)
    depths <- ad[ , , 1] + ad[ , , 2]  # total read depth at each marker/member pair

    snp.chroms <- getCHROM(vcf)
    u.chroms <- unique(snp.chroms)
    ## ad <- get.ad.array(vcf)
    sample.names <- getSAMPLES(vcf)
    parent.indices <- which(sample.names %in% parents)
    marker.names <- getID(vcf)
    n.mem <- length(sample.names)
    display(1, "Completed in ", timer(), "\n")



    if (is.null(parental.results)) {
        ## emission probability code
        timer <- new.timer()
        display(0, "Generating emission probabilities")
        emission.structure <- get.emission.structures(ad,
                                                      sample.names,
                                                      marker.names,
                                                      read.err)
        display(1, "Completed in ", timer(), "\n")




        ## parental imputation and recombination rate estimates
        timer <- new.timer()
        display(0, "Imputing parents and estimating recombination rates")
        parental.results <- determine.parents.and.recombs(emission.structure,
                                                          parents,
                                                          snp.chroms,
                                                          marker.names,
                                                          sample.names,
                                                          parallel,
                                                          cores)
        display(1, "Completed in ", timer(), "\n")
    }



    parental.models <- lapply(parental.results$parental.models,
                              function(chrom.model){
                                  parental.model.to.parents(chrom.model$path)
                              })
    recomb.probs <- parental.results$recombs
    if (! use.50) {
        recomb.probs <- lapply(recomb.probs, function(chrom.recombs) {
            result <- chrom.recombs
            result[result > 0.49] <- 0.05
            result
        })
    }

    gt <- do.call(rbind, lapply(seq_along(parental.models), function(i) {
        p1 <- sequence.to.genome(parental.models[[i]][[1]])
        p2 <- sequence.to.genome(parental.models[[i]][[2]])
        f1 <- cross(p1, p2, recomb.probs[[i]])

        mimic.chrom.common.f1(p1,
                              p2,
                              f1,
                              recomb.probs[[i]],
                              parent.indices,
                              n.mem,
                              n.gen)
    }))

    cols <- colnames(vcf@gt)
    rows <- rownames(vcf@gt)
    vcf@gt <- cbind("GT", gt)

    colnames(vcf@gt) <- cols
    rownames(vcf@gt) <- rows

    write.vcf(vcf, full.out)

    if (! all(dim(depths) == dim(gt)))
        stop("This should never happen. Error in depth and gt dimensions")

    numeric.gt <- matrix(c("0/0"=0, "0/1"=1, "1/1"=2)[gt], nrow=nrow(gt))
    data.str <- abind(depths, numeric.gt, along=3)
    probs.mat <- matrix(c(1-read.err,    read.err,
                                 0.5,    0.5,
                            read.err,    1-read.err), nrow=3, byrow=T)

    sample.ad <- apply(data.str, 1:2, function(depth.and.gt) {
        depth <- depth.and.gt[1]
        gt <- depth.and.gt[2]

        if (depth == 0)
            return(c(0, 0))

        reads <- sample(0:1, size=depth, replace=TRUE, prob=probs.mat[gt+1, ])

        n.alt <- sum(reads)
        n.ref <- depth - n.alt

        c(n.ref, n.alt)
    })
    sample.ad <- aperm(sample.ad, c(2,3,1))

    sample.gt <- apply(sample.ad, 1:2, function(depths) {
        if (depths[1] == 0) {
            if (depths[2] == 0) {
                "./."
            } else {
                "1/1"
            }
        } else {
            if (depths[2] == 0) {
                "0/0"
            } else {
                "0/1"
            }
        }
    })

    sample.ad <- apply(sample.ad, 1:2, paste0, collapse=",")

    ## ## TODO(Jason): This code is ad-hoc and could be cleaned
    ## sample.ad <- matrix(mapply(function(depth, gt) {
    ##     if (depth == 0)
    ##         return("0,0")

    ##     if (gt == "0/0")
    ##         probs <- c(1-read.err, read.err)
    ##     else if (gt == "0/1")
    ##         probs <- c(0.5, 0.5)
    ##     else if (gt == "1/1")
    ##         probs <- c(read.err, 1-read.err)
    ##     else
    ##         stop("This should never happen. Error in gt type")

    ##     reads <- sample(0:1, size=depth, replace=TRUE, prob=probs)
    ##     n.alt <- sum(reads)
    ##     paste0((depth - n.alt), ",", n.alt)
    ## }, depths, gt), nrow=nrow(depths), ncol=ncol(depths))

    ## sample.gt <- apply(sample.ad, 1:2, function(ad) {
    ##     depths <- str.split(ad)
    ##     if (all(depths=="0"))
    ##         gt <- "./."
    ##     else if (depths[1] == "0")
    ##         gt <- "1/1"
    ##     else if (depths[2] == "0")
    ##         gt <- "0/0"
    ##     else
    ##         gt <- "0/1"

    ##     gt  # implicit return
    ## })

    concat <- matrix(paste0(sample.gt, ":", sample.ad, ":", depths), nrow=nrow(sample.gt))
    cols <- colnames(vcf@gt)
    rows <- rownames(vcf@gt)
    vcf@gt <- cbind("GT:AD:DP", concat)

    colnames(vcf@gt) <- cols
    rownames(vcf@gt) <- rows

    write.vcf(vcf, sample.out)
}


mimic.chrom.common.f1 <- function(p1, p2, f1, recomb.probs, parent.indices, n.mem, n.gen) {
    gts <- lapply(seq_len(n.mem),
                  function(i) {
                      if (i == parent.indices[1]) {
                          member.to.gt.str(p1)
                      } else if (i == parent.indices[2]) {
                          member.to.gt.str(p2)
                      } else {
                          member.to.gt.str(
                              generate.selfed.member.common.f1(f1, n.gen, recomb.probs)
                          )
                      }
                  })

    do.call(cbind, gts)
}


mimic.chrom <- function(p1, p2, recomb.probs, parent.indices, n.mem, n.gen) {
    gts <- lapply(seq_len(n.mem),
                  function(i) {
                      if (i == parent.indices[1]) {
                          member.to.gt.str(p1)
                      } else if (i == parent.indices[2]) {
                          member.to.gt.str(p2)
                      } else {
                          member.to.gt.str(
                              generate.selfed.member(p1, p2, n.gen, recomb.probs)
                          )
                      }
                  })

    do.call(cbind, gts)
}


sequence.to.genome <- function(seq) {
    individual <- list()
    class(individual) <- c("genome", class(individual))

    individual$cA <- as.logical(sapply(seq, function(state) {floor((state - 1) / 2)}))
    individual$cB <- as.logical(seq - 1 - individual$cA*2)

    individual
}


cross <- function(p1, p2, recomb.probs) {
    if (! "genome" %in% class(p1)) {
        stop("p1 must be of type genome")
    }
    if (! "genome" %in% class(p2)) {
        stop("p2 must be of type genome")
    }
    if (length(p1) != length(p2)) {
        stop("p1 and p2 have different lengths")
    }
    if (FALSE) {
        ## TODO(Jason): check that cA and cB for p1 and p2 are all same length
        ## and length of recomb.probs is one less
    }
    if (FALSE) {
        ## TODO(Jason): warn if probability of recombination is greater than 50% across
    }

    offspring <- list()
    class(offspring) <- c("genome", class(offspring))

    offspring$cA <- gamete(p1, recomb.probs)
    offspring$cB <- gamete(p2, recomb.probs)

    offspring  # implicit return
}


gamete <- function(parent, recomb.probs) {
    if (! "genome" %in% class(parent)) {
        stop("parent must be of type genome")
    }
    gamete <- sample(1:2, size=1)  # uniform randomly choose 1 which gamete
    if (gamete==1) {
        base.gamete <- parent$cA
        recomb.gamete <-parent$cB
    } else {
        base.gamete <- parent$cB
        recomb.gamete <-parent$cA
    }

    ## vector indicating if a physical recombination will occur between each
    ## pair of adjacent sites
    recomb.bools <- runif(length(recomb.probs)) < recomb.probs

    gamete.result <- base.gamete  # initialization
    odd.recombs <- FALSE  # indicate even or odd num of recombs up to point

    for (i in seq_along(recomb.bools)) {
        if (recomb.bools[i]) {
                                        # a recombination has occurred between index i and i+1 of
                                        # gamete.result
            odd.recombs <- !odd.recombs
        }
        if (odd.recombs) {
            ## the +1 is necessary because recomb.bools is one element
            ## smaller than base.gamete, recomb.gamete, and gamete.result
            ## and because index 1 of gamete.result is implied to be the
            ## same as index 1 of base.gamete by design
            gamete.result[i+1] <- recomb.gamete[i+1]
        }
    }
    gamete.result  # implicit return
}


generate.selfed.member.common.f1 <- function(f1, n.gen, recomb.probs) {
    member <- f1
    for (i in seq_len(n.gen - 1)) {
        member <- cross(member, member, recomb.probs)
    }
    member
}



generate.selfed.member <- function(p1, p2, n.gen, recomb.probs) {
    member <- cross(p1, p2, recomb.probs)
    for (i in seq_len(n.gen - 1)) {
        member <- cross(member, member, recomb.probs)
    }
    member
}


member.to.gt.str <- function(member) {
    gt <- mapply(sum, member$cA, member$cB)  # boolean addition will be 0, 1, or 2
    c("0/0", "0/1", "1/1")[gt + 1]
}


print.labyrinth.sim.header <- function() {

    ## the image looks funny because '\' in the displayed image must be '\\' in the code
    writeLines("")
    writeLines(" _____________________________________________________________________")
    writeLines("|          __          ____        ____  _____                        |")
    writeLines("|         / /         / __ \\      / __ \\/_  _/                        |")
    writeLines("|        / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __      |")
    writeLines("|       / /   / _  \\/ _  _/\\ \\/ / _  _/ / / /  | / /_  __/ /_/ /      |")
    writeLines("|      / /___/ /_/ / /_\\ \\  \\  / / \\ \\_/ /_/ /||/ / / / / __  /       |")
    writeLines("|     /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/        |")
    writeLines("|                                                                     |")
    writeLines("|     ========  P O P U L A T I O N   S I M U L A T O R  ========     |")
    writeLines("|_____________________________________________________________________|")
    writeLines("")
}






