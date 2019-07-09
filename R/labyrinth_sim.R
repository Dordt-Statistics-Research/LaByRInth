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


################################################################################
############################# LIBRARY REQUIREMENTS #############################
################################################################################


require(vcfR, quietly=T)
require(abind, quietly=T)
require(digest, quietly=T)





LabyrinthMask <- function(vcf, parents, out.file, depth=0, lik.ratio=100, rerr=0.01) {

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
        vcf <- read.vcfR(vcf, verbose=F)
    }

    ad.arr <- get.ad.array(vcf)
    depths <- apply(ad.arr, 1:2, sum)  # sum individual allelic reads
    liklihood.ratios <- apply(ad.arr, 1:2, function(reads) {
        get.liklihood.ratio(reads[1], reads[2], rerr)
    })

    parental         <- matrix(rep(getSAMPLES(vcf) %in% parents, nrow(ad.arr)),
                               nrow=nrow(ad.arr),
                               byrow=TRUE)
    sufficient.depth <- depths >= depth
    sufficient.ratio <- liklihood.ratios >= lik.ratio
    mask <- sufficient.depth & sufficient.ratio & !parental

    n.total <- prod(dim(depths))
    n.called <- sum(depths != 0)
    n.masked <- sum(mask)

    message(" * Sites are considered progeny/marker pairs")
    message(" * ", sum(mask), " sites will be masked of ", n.total, " total sites (",
            round(n.masked / n.total, 3)*100, "%)")
    message(" * ", sum(mask), " sites will be masked of ", n.called, " called sites (",
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


LabyrinthMimic <- function(parental, full.out, sample.out, parallel=TRUE, cores=4) {

    total.timer <- new.timer()
    print.labyrinth.header()


    if (! verify.file.extension(full.out, ".vcf.gz")) {
        stop("Output file (full.out) name must end with '.vcf.gz'\n")
    }

    if (! verify.file.extension(sample.out, ".vcf.gz")) {
        stop("Output file (sample.out) name must end with '.vcf.gz'\n")
    }


    ## Load variables from parental
    timer <- new.timer()
    display(0, "Restoring variables from parental imputation result")

    vcf <- parental$vcf
    parents <- parental$parents
    parent.models <- parental$parent.models
    generation <- parental$generation
    read.err <- parental$read.err
    geno.err <- parental$geno.err
    site.pair.transition.probs <- parental$site.pair.transition.probs
    snp.chroms <- parental$snp.chroms
    u.chroms <- parental$u.chroms
    ad <- parental$ad
    sample.names <- parental$sample.names
    marker.names <- parental$marker.names
    emission.structure <- parental$emission.structure

    depths <- ad[ , , 1] + ad[ , , 2]  # total read depth at each marker/member pair
    parent.indices <- which(sample.names %in% parents)
    n.mem <- length(sample.names)

    display(1, "Completed in ", timer(), "\n")


    ## if (! use.50) {
    ##     recomb.probs <- lapply(recomb.probs, function(chrom.recombs) {
    ##         result <- chrom.recombs
    ##         result[result > 0.49] <- 0.05
    ##         result
    ##     })
    ## }

    gt <- do.call(rbind, lapply(seq_along(parent.models), function(i) {
        indiv.parent.models <- extract.each.parent(parent.models$model)
        p1 <- sequence.to.genome(indiv.parent.models[[1]])
        p2 <- sequence.to.genome(indiv.parenta.models[[2]])
        f1 <- cross(p1, p2, parent.models$recombs[[i]])

        mimic.chrom.common.f1(p1,
                              p2,
                              f1,
                              recomb.probs[[i]],
                              parent.indices,
                              n.mem,
                              generation)
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

    concat <- matrix(paste0(sample.gt, ":", sample.ad, ":", depths), nrow=nrow(sample.gt))
    cols <- colnames(vcf@gt)
    rows <- rownames(vcf@gt)
    vcf@gt <- cbind("GT:AD:DP", concat)

    colnames(vcf@gt) <- cols
    rownames(vcf@gt) <- rows

    write.vcf(vcf, sample.out)
}


## When the parents are not totally homozygous for different alleles, then it
## matters whether a single F1 plant is selfed to create the population of F2s
## or if the parents are crossed multiple times to create a generation of
## different F1 plants
mimic.chrom.common.f1 <- function(p1, p2, f1, recomb.probs, parent.indices,
                                  n.mem, n.gen) {

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
