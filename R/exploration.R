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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param model a string such as "F2", or "F7"
##' @return 
##' @author Jason Vander Woude
hom.poly.vcf.to.common.ref.alt.parents <- function(vcf, parents) {

    ## load vcf if needed
    if (inherits(vcf, "character")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- vcfR::read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }

    if (!all(parent.hom.and.poly(vcf, parents))) {
        stop("Parents are not homozygous within and polymorphic between")
    }

    p1.gt.str <- getGT(vcf)[, parents[1]]
    p2.gt.str <- getGT(vcf)[, parents[2]]

    is.ref <- function(str) {str == "0/0"}
    is.alt <- function(str) {str == "1/1"}

    p1.is.ref <- is.ref(p1.gt.str) | is.alt(p2.gt.str)

    fix <- vcf@fix
    fix[!p1.is.ref, "REF"] <- vcf@fix[!p1.is.ref, "ALT"] # swap alt and ref
    fix[!p1.is.ref, "ALT"] <- vcf@fix[!p1.is.ref, "REF"] # swap alt and ref
    fix[ , "INFO"] <- "none" # remove info which may have dealt with ref/alt

    swap.geno <- function(strs) {
        ## replace 0's with _'s then
        ## replace 1's with 0's then
        ## replace _'s with 1's
        gsub("_","1",gsub("1","0",gsub("0", "_", strs)))
    }

    swap.reads <- function(strs) {
        sapply(strs, function(str) {
            reads <- ad.to.num(str)
            paste0(reads[2], ",", reads[1])
        })
    }

    gt <- getGT(vcf)
    ad <- getAD(vcf)
    for (row in which(!p1.is.ref)) {
        gt[row, ] <- swap.geno(gt[row, ])
        ad[row, ] <- swap.reads(ad[row, ])
    }

    concat <- matrix(paste0(gt, ":", ad), nrow=nrow(gt))
    vcf@gt <- cbind("GT:AD", concat)
    colnames(vcf@gt) <- c("FORMAT", colnames(gt))
    rownames(vcf@gt) <- rownames(gt)
    vcf@fix <- fix
    vcf@meta <- c("##fileformat=VCFv4.2",
                  paste0("##filedate=",format(Sys.Date(),"%Y%m%d")),
                  paste0("##source=LaByRInth_version_", version()))

    vcf
}




calculate.transitional.characteristics <- function(vcf, parents, model) {
    if (model %in% c("F2","F3","F4","F5","F6","F7")) {
        tryCatch({
            trans.file <- system.file("extdata",
                                      "transition-probs",
                                      paste0(model, ".R"),
                                      package = "LaByRInth",
                                      mustWork = TRUE)
            source(trans.file)
        }, error = function(e) {
            stop(paste("Invalid model"))
        })

        ## load vcf if needed
        if (inherits(vcf, "character")) {
            timer <- new.timer()
            display(0, "Loading vcf")
            vcf <- vcfR::read.vcfR(vcf, verbose=F)
            display(1, "Completed in ", timer(), "\n")
        }

        ## convert vcf and change the reference allele such that parent 1 is
        ## always homozygous reference and parent 2 is always homozygous
        ## alternate
        vcf <- hom.poly.vcf.to.common.ref.alt.parents(vcf, parents)

        ## site.pair.transition.probs is loaded when trans.file is sourced
        ## above. It is a list containing 16 entries. Each entry is itself a
        ## list containing 16 functions. Each of these 16*16=256 functions takes
        ## a single parameter r. The previous line instantiates a symbolic
        ## variable which applies symmetrically to both parents because this is
        ## a parameter of the species. An organism has 2 homologous chromosomes
        ## and we examine two markers on this homologous pair.  Each of those
        ## two homologous chromosomes has two sister chromatids. When a cell
        ## undergoes meiosis, four gametes will be produced, and for a fixed
        ## gamete, the allele at the first marker came from one of the four
        ## sister chromatids. Call this fixed gamete G and the corresponding
        ## sister chromatid S. The paremeter r represents the probability that
        ## there is some kind of crossover event between the two markers such
        ## that allele at the second marker of gamete G belongs to a chromatid
        ## which is not S or the sister of S. The parameter is specified this
        ## way because if there is a crossover event such that the allele of G
        ## at the second marker belongs to S or the sister of S, the type of
        ## allele in G is the same as if no crossover event had occurred.
        ##
        ## This 2D list can be used to determine what the transition
        ## probabilites are between two markers if the value r is known and the
        ## allele type (reference or alternate) of both homologs of both parents
        ## are known at both markers. The first index of the 2D list corresponds
        ## to marker 1 and the second index corresponds to marker 2. Indexing
        ## should be viewed in terms of binary representations. Specifically,
        ## the binary number obtained where the four bits are as follows
        ##
        ##     Most significant/left bit   (bit 3): allele of parent 1 homolog 1
        ##                                 (bit 2): allele of parent 1 homolog 2
        ##                                 (bit 1): allele of parent 2 homolog 1
        ##     least significant/right bit (bit 0): allele of parent 2 homolog 2
        ##
        ## would correspond to the the first or second index (in a 0-based
        ## indexing system) of the 2D list when the alleles are for the first or
        ## second marker respectively. The allele should be 0 if it is the
        ## reference allele and should be 1 if it is the alternate
        ## allele. Because R uses a 1-based indexing system, it is necessary to
        ## add 1 to the above number to obtain the correct index.
        ## The file inst/extdata/multi-model-symbolics.sage should be viewed for
        ## more details.
        ##
        ## Because the vcf file has been modified such that parent 1 is always
        ## homozygous reference and parent 2 is always homozygous alternate,
        ## then the appropriate binary value is 0011 for both markers. 0011 is
        ## binary for 3, so we use index 4 in both dimensions of this matrix
        ## because of the addition of 1 to compensate for the indexing
        ## type. Thus, transition.prob.matrix.generator will be a function that
        ## takes a single parameter r and returns a transition matrix. The
        ## returned matrix (call it T) is also indexed by binary
        ## representations.
        ##
        ## To use the transition matrix T, for a taxa/progeny/member of the
        ## population, encode a possible genetic state of the first marker as a
        ## binary representation of a number, i, as follows:
        ##
        ##     Most significant/left bit   (bit 1): allele of progeny homolog 1
        ##     Least significant/right bit (bit 0): allele of progeny homolog 2
        ##
        ## Encode the genetic state of a possible second marker as the number j
        ## in the same way. Then in a 0-based indexing system, row i and column
        ## j of T is the probability that the chosen progeny would be produced
        ## by the speficied parents (those used to index into the 2D
        ## list). Because R uses 1-based indexing the probability is actually
        ## obtained using row i+1 and column j+1. This probability will be
        ## written as T[i+1,j+1]. To obtain a transition probability to use in a
        ## hidden Markov Model (HMM), this probability should be conditioned on
        ## the first marker of the progeny having a genetic state as defined by
        ## the number i. That is, the probability T[i+1,j+1] must be divided by
        ## the sum of row i+1 of T (this sum is T[i+1,1] + T[i+1,2] + T[i+1,3] +
        ## T[i+1,4]).
        transition.prob.matrix.generator <- site.pair.transition.probs[[4]][[4]]


        ## A useful feature of these matrices is that while conditioning on a
        ## row gives transition probabilities, if no conditioning is used
        ## (i.e. the matrix is used as is) then the matrix encodes the
        ## probabilities of each possible marker pair genetic makeup. Thus in
        ## any population, and for any i,j if a pair of markers has an
        ## associated transition parameter of r, then the value T(r)[i,j] is the
        ## expected proportion of the population in which the first marker in the
        ## pair has configuration i and the second marker has configuration
        ## j. In a large population, the expected value can be approximated by
        ## the proportion of the population having such configurations. This
        ## idea works in reverse as well meaning that by first checking for each
        ## configuration pair <i,j> the proportion of the population matching that
        ## configuration, the value of r can be empirically estimated by
        ## choosing the value r such that the matrix T(r) most closely matches
        ## what is observed in the population. The matrix indices can be
        ## summarized as follows where ref is reference, alt is alternate, hom
        ## is homozygous, het t.1 is heterozygous type 1 (homolog 1 is reference
        ## and homolog 2 is alternate), and het t.2 is heterozygous type 2
        ## (homolog 1 is alternate and homolog 2 is reference).
        ##
        ##            | marker1 -> marker2
        ##    --------+--------------------
        ##     T[1,1] | hom ref -> hom ref
        ##     T[4,4] | hom alt -> hom alt
        ##    --------+--------------------
        ##     T[1,4] | hom ref -> hom alt
        ##     T[4,1] | hom alt -> hom ref
        ##    --------+--------------------
        ##     T[2,1] | het t.1 -> hom ref
        ##     T[3,1] | het t.2 -> hom ref
        ##     T[1,2] | hom ref -> het t.1
        ##     T[1,3] | hom ref -> het t.2
        ##     T[4,2] | hom alt -> het t.1
        ##     T[4,3] | hom alt -> het t.2
        ##     T[2,3] | het t.1 -> hom alt
        ##     T[2,4] | het t.2 -> hom alt
        ##    --------+--------------------
        ##     T[2,2] | het t.1 -> het t.1
        ##     T[3,3] | het t.2 -> het t.2
        ##    --------+--------------------
        ##     T[2,3] | het t.1 -> het t.2
        ##     T[3,2] | het t.2 -> het t.1
        ##
        ## Further, it turns out that when both parents are homozygous for
        ## opposite alleles at adjacent markers, and the progeny come from an
        ## F{n} population where n is a positive integer, then the following
        ## hold for the matrix T returned by transition.prob.matrix.generator
        ## regardless of the parameter r. The following equations are in 1-based
        ## indexing.
        ##
        ## T[1,1] = T[4,4]
        ## T[1,4] = T[4,1]
        ## T[2,1] = T[3,1] = T[1,2] = T[1,3] = T[4,2] = T[4,3] = T[2,4] = T[3,4]
        ## T[2,2] = T[3,3]
        ## T[3,2] = T[2,3]
        ##
        ## Thus, rather than categorizing each progeny into one of 16 categories
        ## they can be placed in one of 5 categories which should make the
        ## proportions more accurately reflect the expected values because the
        ## sample size will be larger. Further, because configurations <2,2> and
        ## <2,3> and <3,2> and <3,3> all represent progeny that are heterozygous
        ## at both markers, it is not possible to distinguish between them when
        ## looking at reads (because it is unknown which homolog a read came
        ## from), so these states will be combined as well. Thus each progeny in
        ## the population (that has reads) will be counted toward exactly one of the
        ## following 4 categories:
        ##
        ##         Category     |         Expected Proportion
        ##    ------------------|------------------------------------
        ##     hom X -> hom X   |T[1,1] + T[4,4]
        ##    ------------------|------------------------------------
        ##     hom X -> hom !X  |T[1,4] + T[4,1]
        ##    ------------------|------------------------------------
        ##     het -> hom OR    |T[2,1] + T[3,1] + T[1,2] + T[1,3]
        ##     hom -> het       |+ T[4,2] + T[4,3] + T[2,4] + T[3,4]
        ##    ------------------|------------------------------------
        ##     het -> het       |T[2,2] + T[3,3] + T[3,2] + T[2,3]
        ##
        types <- c("X -> X", "X -> !X", "het -> hom", "het -> het")
        condensed.trans.characteristics.generator <- function(r) {
            T <- transition.prob.matrix.generator(r)
            return.val <- c(T[1,1] + T[4,4],
                            T[1,4] + T[4,1],
                            T[2,1] + T[3,1] + T[1,2] + T[1,3]
                            + T[4,2] + T[4,3] + T[2,4] + T[3,4],
                            T[2,2] + T[3,3] + T[3,2] + T[2,3])

            names(return.val) <- c("X -> X","X -> !X","het -> hom","het -> het")
            return.val
        }


        ## Convert a genotype string to a value 0,1,2 indicating the number of
        ## alternate alleles. May result in NA's
        numeric.gt <- function(gt.str) {sum(gt.to.num(gt.str))}

        estimate.r <- function(marker1, marker2) {
            numeric.gt.1 <- sapply(gt[marker1, ], numeric.gt)
            numeric.gt.2 <- sapply(gt[marker2, ], numeric.gt)

            ## determine which progeny have a called genotype at both markers
            ## and remove those that don't
            valid <- !is.na(numeric.gt.1) & !is.na(numeric.gt.2)
            num.valid <- sum(valid)
            numeric.gt.1 <- numeric.gt.1[valid]
            numeric.gt.2 <- numeric.gt.2[valid]

            ## Rows correspond to marker1 being 0,1,2 respectively and columns
            ## correspond to marker2 being 0,1,2 respectively
            type.layout <- matrix(
                c("X -> X"    ,    "het -> hom",    "X -> !X"   ,
                  "het -> hom",    "het -> het",    "het -> hom",
                  "X -> !X"   ,    "het -> hom",    "X -> X"    ),
                byrow = TRUE, nrow=3)

            check.type <- function(gt.1, gt.2) {
                ## the `+1` is required to change from 0,1,2 into R's 1-based
                ## indexing so the indices will be from 1,2,3
                type.layout[gt.1 + 1, gt.2 + 1] # return value
            }

            observ.recomb.types <- mapply(check.type, numeric.gt.1, numeric.gt.2)

            ## compute proportion of TRUE in a vector
            prop <- function(vec) {
                sum(vec) / length(vec)
            }

            ## This is a vector of length 4 indicating the proportion of each
            ## transition type actually occurring in the population
            props <- sapply(types, function(t) {prop(observ.recomb.types == t)})

            objective.function <- function(r) {
                expected.props <- condensed.trans.characteristics.generator(r)
                ## Because this is a proof of concept method and not what is
                ## actually used, we will want to minimize the visual distance
                ## between the data points and theory lines, so we minimize the
                ## sum of absolute values of differences
                sum((expected.props - props)^2)
            }

            ## DON'T USE THIS NUMERICAL OPTIMIZATION BECAUSE THE FUNCTIONS ARE
            ## NOT CONVEX
            ##
            ## result <- optim(par=1/2,    # initial value
            ##                 fn=objective.function,
            ##                 method="Brent",
            ##                 lower=0,
            ##                 upper=1)
            ## result <- result$par

            possible.r <- seq(from=0, to=1, by=1/1000)
            objective.vals <- sapply(possible.r, objective.function)

            ## return estimated r value and proportions and total number
            c(possible.r[which.min(objective.vals)], props, num.valid)
        }


        ## Next, the markers will be broken down by chromosome because no
        ## markers from separate chromosomes should be compared. Then, in each
        ## chromosome the appove categorization will take place.
        chroms <- unique(vcf@fix[ , "CHROM"])
        all.markers <- vcf@fix[ , "ID"]
        # get all member names by removing "FORMAT" column which is first
        population.names <- colnames(vcf@gt)[2:ncol(vcf@gt)]
        gt <- getGT(vcf)[ , !(population.names %in% parents)] # non-parental genotypes

        r.estimate.data.frames <- lapply(chroms, function(chrom) {
            print(chrom)
            correct.chrom <- vcf@fix[ , "CHROM"] == chrom # find rows of same chrom
            markers <- all.markers[correct.chrom] # names of markers in chrom
            positions <- as.numeric(vcf@fix[all.markers %in% markers, "POS"])
            names(positions) <- markers
            if (length(markers) > 2) {
                ## estimate r for each adjacent pair of markers
                r.ests.and.props.and.n.valid <-
                    matrix(sapply.pairs(markers, estimate.r),
                           byrow=TRUE, ncol=6)

                ## skip r estimate and number of valid
                all.props <- r.ests.and.props.and.n.valid[ , 2:5]
                colnames(all.props) <- types

                phys.dists <-
                    sapply.pairs(markers,
                                 function(m1, m2) {positions[m2]-positions[m1]})

                df1 <- data.frame(marker1=markers[1:(length(markers)-1)],
                                  marker2=markers[2:length(markers)],
                                  phys.dist=phys.dists,
                                  r=r.ests.and.props.and.n.valid[ , 1],
                                  n=r.ests.and.props.and.n.valid[ , 6])
                df2 <- as.data.frame(all.props)

                ## return a data frame for the chromosome
                cbind(df1, df2)

            } else {
                data.frame()  # return empty data frame
            }
        })

        ## bind all results together
        result <- do.call(rbind, r.estimate.data.frames)

    } else {
        stop("Invalid model")
    }
}


## rows <- clean.vcf@fix[ , "CHROM"] %in% c("1A","1B","1C","1D"); small.vcf@fix <- clean.vcf@fix[rows, ]; small.vcf@gt <- clean.vcf@gt[rows, ]


characteristics.to.ggplot.format <- function(plot.data) {
    trans.types <- c("X -> X", "X -> !X", "het -> hom", "het -> het")
    data.frames <- lapply(trans.types, function(type) {
        inner.df <- plot.data
        inner.df[ , "proportion"] <- inner.df[ , type]
        inner.df[ , "type"] <- type
        for (t in trans.types) {
            inner.df[ , t] <- NULL
        }
        inner.df
    })

    ggplot.data <- do.call(rbind, data.frames)
    rownames(ggplot.data) <- NULL       # remove row names
    ggplot.data                         # return value
}



## plot.data.to.ggplot.format <- function(plot.data) {
##     trans.types <- c("X -> X", "X -> !X", "het -> hom", "het -> het")
##     data.frames <- apply(plot.data, 1, function(row.of.plot.data) {
##         ## inner.df <- data.frame()
##         browser()
##         inner.df <- data.frame(marker1=I(row.of.plot.data["marker1"]),
##                                marker2=I(row.of.plot.data["marker2"]),
##                                r=I(row.of.plot.data["r"]),
##                                type=I(trans.types),
##                                proportion=I(row.of.plot.data[trans.types]),

##                                row.names=NULL)

##         ## inner.df$marker1 <- as.character(inner.df$marker1) # R thinks it is a factor
##         ## inner.df$marker2 <- as.character(inner.df$marker2) # R thinks it is a
##         ##                                 # factor
##         ## browser()
##         ## inner.df$r <- as.numeric(inner.df$r) # R thinks it is a factor
##         ## inner.df$type <- as.character(inner.df$type) # R thinks it is a factor
##         ## inner.df$proportion <- as.numeric(inner.df$proportion) # R thinks it is a factor

##         ## df.columns.1 <- rbind(row.of.plot.data[c("marker1", "marker2", "r", "X -> X")],
##         ##                       row.of.plot.data[c("marker1", "marker2", "r", "X -> !X")],
##         ##                       row.of.plot.data[c("marker1", "marker2", "r", "het -> hom")],
##         ##                       row.of.plot.data[c("marker1", "marker2", "r", "het -> het")])
##         ## df.columns.2 <- data.frame(type=c("X -> X", "X -> !X", "het -> hom",
##         ##                                   "het -> het"))
##         ## cbind(df.columns.1, df.columns.2)

##         inner.df
##     })

##     ggplot.data <- do.call(rbind, data.frames)
##     rownames(ggplot.data) <- NULL       # remove row names
##     ggplot.data                         # return value
## }



generate.plot <- function(ggplot.data, transition.proportions.generator) {
    require(ggplot2)
    get.theory.function <- function(type) {
        function(r) {
            sapply(r, function(each) {
                transition.proportions.generator(each)[type]
            })
        }
    }
    x.max <- 1
    color1 <- "#37474F"
    color2 <- '#CACACA'
    xlims <- c(0, x.max)
    ylims <- c(0, 1)
    ## gentext <- ifelse(gen==2, "Generation", "Generations")
    thickness <- 0.5
    mod.data <- ggplot.data[ggplot.data$r <= x.max, ]
    ggplot(mod.data, aes(x=r, y=proportion, color=type)) + geom_point(size=1.75) +
        stat_function(aes(color="X -> X"), fun=get.theory.function("X -> X"),
                      xlim=xlims, size=thickness) +
        stat_function(aes(color="X -> !X"),
                      fun=get.theory.function("X -> !X"), xlim=xlims, size=thickness) +
        stat_function(aes(color="het -> hom"),
                      fun=get.theory.function("het -> hom"), xlim=xlims, size=thickness) +
        stat_function(aes(color="het -> het"),
                      fun=get.theory.function("het -> het"), xlim=xlims, size=thickness) +
        ## scale_color_manual(values = c("X -> X" = 1,
        ##                               "X -> !X" = 2,
        ##                               "het -> hom" = 3,
        ##                               "het -> het" = 4)) +
        coord_cartesian(xlim = xlims) +
        coord_cartesian(ylim = ylims) +
        scale_y_continuous(breaks=seq(0, 1, 0.2)) +
        scale_x_continuous(breaks=seq(0, x.max, 0.2)) +
        ggtitle(paste0("Transitional Characteristics")) +
        xlab("Recombination Parameter (r)") +
        ylab("Proportion of Population") +
        scale_color_discrete(name="Type") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(fill = '#37474F', colour = '#37474F'),
              panel.background = element_rect(fill = color2, color = color2),
              text = element_text(size=20),
              legend.key = element_rect(fill = color1, color=color1),
              legend.background = element_rect(fill = color1, color=color1),
              axis.text.x=element_text(color = color2),
              axis.text.y=element_text(color = color2),
              axis.title.x=element_text(color = color2, vjust = -3),
              axis.title.y=element_text(color = color2, vjust = 3),
              legend.title=element_text(color = color2),
              legend.text=element_text(color = color2),
              title=element_text(color = color2),
              plot.margin = unit(c(1,1,1,1), "cm"),
              legend.key.size = unit(3, 'lines'))
}



get.condensed.characeristic.generator <- function(vcf, parents, model) {
    if (model %in% c("F2","F3","F4","F5","F6","F7")) {
        tryCatch({
            trans.file <- system.file("extdata",
                                      "transition-probs",
                                      paste0(model, ".R"),
                                      package = "LaByRInth",
                                      mustWork = TRUE)
            source(trans.file)
        }, error = function(e) {
            stop(paste("Invalid model"))
        })

        ## load vcf if needed
        if (inherits(vcf, "character")) {
            timer <- new.timer()
            display(0, "Loading vcf")
            vcf <- vcfR::read.vcfR(vcf, verbose=F)
            display(1, "Completed in ", timer(), "\n")
        }

        vcf <- hom.poly.vcf.to.common.ref.alt.parents(vcf, parents)

        transition.prob.matrix.generator <- site.pair.transition.probs[[4]][[4]]

        types <- c("X -> X", "X -> !X", "het -> hom", "het -> het")
        condensed.trans.characteristics.generator <- function(r) {
            T <- transition.prob.matrix.generator(r)
            return.val <- c(T[1,1] + T[4,4],
                            T[1,4] + T[4,1],
                            T[2,1] + T[3,1] + T[1,2] + T[1,3]
                            + T[4,2] + T[4,3] + T[2,4] + T[3,4],
                            T[2,2] + T[3,3] + T[3,2] + T[2,3])

            names(return.val) <- types
            return.val
        }
        condensed.trans.characteristics.generator # return value
    } else {
        stop("Invalid model")
    }
}


## return vcf object with data only from specified chromosomes
filter.vcf.by.chrom <- function(vcf, chroms) {

    ## load vcf if needed
    if (inherits(vcf, "character")) {
        timer <- new.timer()
        display(0, "Loading vcf")
        vcf <- vcfR::read.vcfR(vcf, verbose=F)
        display(1, "Completed in ", timer(), "\n")
    }

    keep <- vcf@fix[ , "CHROM"] %in% chroms

    fix <- vcf@fix[keep, ]
    fix[ , "INFO"] <- "none" # remove info which may have dealt with ref/alt
    vcf@fix <- fix
    vcf@gt <- vcf@gt[keep, ]
    vcf@meta <- c("##fileformat=VCFv4.2",
                  paste0("##filedate=",format(Sys.Date(),"%Y%m%d")),
                  paste0("##source=LaByRInth_version_", version()))

    vcf
}



do.it <- function() {
    parents <- c("LAKIN", "FULLER")
    hom.poly.vcf <- vcfR::read.vcfR("~/Desktop/LF-big-hom-poly.vcf.gz")
    nice.vcf <- hom.poly.vcf.to.common.ref.alt.parents(hom.poly.vcf, parents)
    characteristics.sq <-
        calculate.transitional.characteristics(nice.vcf, parents, "F5")
    saveRDS(characteristics.sq, "~/Desktop/LF-characteristics-sq.rds")
    gg.data <- characteristics.to.ggplot.format(characteristics.sq)
    generate.plot(gg.data, get.condensed.characeristic.generator(nice.vcf, parents, "F5"))


    parents <- c("RsaI_B73", "RsaI_CG")
    vcf <- vcfR::read.vcfR("~/Desktop/RsaI-filtered.vcf.gz")
    hom.poly.vcf <- LabyrinthFilter(vcf, "~/Desktop/RsaI-hom-poly.vcf.gz",
                                    parents, require.hom.poly=TRUE)
    nice.vcf <- hom.poly.vcf.to.common.ref.alt.parents(hom.poly.vcf, parents)
    RsaI.characteristics <-
        calculate.transitional.characteristics(nice.vcf, parents, "F2")
    saveRDS(RsaI.characteristics, "~/Desktop/RsaI-characteristics.rds")
    gg.data <- characteristics.to.ggplot.format(RsaI.characteristics)
    generate.plot(gg.data, get.condensed.characeristic.generator(nice.vcf, parents, "F2"))


    HincII.parents <- c("HincII_B73", "HincII_CG")
    HincII.vcf <- vcfR::read.vcfR("~/Desktop/HincII-masked_1.vcf.gz")
    HincII.hom.poly.vcf <- LabyrinthFilter(HincII.vcf, "~/Desktop/HincII-hom-poly.vcf.gz",
                                           HincII.parents, require.hom.poly=TRUE)
    HincII.hom.poly.vcf <- vcfR::read.vcfR("~/Desktop/HincII-hom-poly.vcf.gz")
    HincII.nice.vcf <- hom.poly.vcf.to.common.ref.alt.parents(HincII.hom.poly.vcf, HincII.parents)
    HincII.characteristics.sq <-
        calculate.transitional.characteristics(HincII.nice.vcf, HincII.parents, "F2")
    saveRDS(HincII.characteristics.sq, "~/Desktop/HincII-characteristics-sq.rds")
    HincII.characteristics.sq <- readRDS("~/Desktop/HincII-characteristics-sq.rds")
    HincII.gg.data <- characteristics.to.ggplot.format(HincII.characteristics.sq)
    generate.plot(HincII.gg.data, get.condensed.characeristic.generator(HincII.nice.vcf, HincII.parents, "F2"))



    


    LF.parents <- c("LAKIN", "FULLER")
    LF.vcf <- vcfR::read.vcfR("~/Desktop/LF-masked_1-chr_1A.vcf.gz")
    LF.hom.poly.vcf <- LabyrinthFilter(LF.vcf, "~/Desktop/LF-masked_1-chr_1A-hom-poly.vcf.gz",
                                       LF.parents, require.hom.poly=TRUE)
    LF.nice.vcf <- hom.poly.vcf.to.common.ref.alt.parents(LF.hom.poly.vcf, parents)
    LF.characteristics <-
        calculate.transitional.characteristics(LF.nice.vcf, LF.parents, "F5")

    LabyrinthImputeParents(vcf="~/Desktop/LF-masked_1-chr_1A.vcf.gz", out.file="~/Desktop/LF-parents-chr_1A.rds", parents=c("LAKIN", "FULLER"), generation=5)
    parental <- readRDS("~/Desktop/LF-parents-chr_1A.rds")
    markers <- parental$marker.names
    hom.poly.markers <- parental$parent.models[["1A"]]$model %in% c(4,13)
    ## ensure both surrounding markers were hom poly
    valid.transition.regions <- sapply.pairs(hom.poly.markers, function(m1, m2) {m1 && m2})
    LF.sub.chraractersitics <- LF.characteristics[valid.transition.regions, ]
    LF.sub.chraractersitics$r <- parental$parent.models[["1A"]]$recombs[valid.transition.regions, ]



    saveRDS(characteristics.sq, "~/Desktop/LF-characteristics-sq.rds")
    gg.data <- characteristics.to.ggplot.format(characteristics.sq)
    generate.plot(gg.data, get.condensed.characeristic.generator(nice.vcf, parents, "F5"))
    non.hom.poly.sites.in.parents <- !(parental$parent.models[["1A"]]$model %in% c(4,13))
    bad.marker.indices <- which(non.hom.poly.sites.in.parents)
    bad.pair.indices <- c(bad.marker.indices, bad.marker.indices-1)
    bad.marker.indices <- bad.marker.indices[bad.marker.indices >= 1 & bad.marker.indices <= length(markers) - 1]
    

}
