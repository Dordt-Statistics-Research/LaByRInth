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


##' Get matrix of emission probabilities
##'
##' Computes the emission probabilities for each state based on allelic depth of
##'     coverage (using the binomial assumption). Sample must be of length 1 and
##'     the entries of the parent.geno matrix must be either 0, 1, or NA.
##' @title
##' @param vcf an object of class vcf
##' @param sample an index or name of a variant
##' @param parent.geno a matrix of parental genotypes
##' @param prefs a preferences object
##' @return a matrix of posterior probabilities
##' @author Jason Vander Woude
GetProbabilities <- function(vcf, sample, chrom, parent.geno, prefs, pseudo.site.pos) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (length(sample) != 1) {
        stop("Length of sample must be 1")
    }
    if (length(chrom) != 1) {
        stop("Length of sample must be 1")
    }

    chrom.indices <- grepl(paste0(chrom, ":"), rownames(parent.geno))
    parent.geno <- parent.geno[chrom.indices, ]

    gt <- Get(vcf, "GT", sample, chrom)
    ad <- Get(vcf, "AD", sample, chrom)

    n.superstates <- 3 * (2* prefs$generation + 1)
    n.boundaries <- 2 * prefs$generation  # number total gamete crossovers
    parents <- colnames(parent.geno)
    states <- c(parents, "HET")
    superstates <- c(sapply(0:n.boundaries, function(boundary) {
        paste0(states, "-", boundary)
    }))
    ret.val <- matrix(NA_integer_, nrow = nrow(gt), ncol = n.superstates)
    rownames(ret.val) <- rownames(gt)
    colnames(ret.val) <- superstates
    class(ret.val) <- "prob"

    rerr <- prefs$read.err

    for (row in 1:nrow(ret.val)) {
        ## gt and ad are both 3D structures where the first index is the site,
        ## the second index is the sample/variety, and the third index allows
        ## getting the various genotypes (e.g. a 0 in both indices of the last
        ## dimension represents 0/0). The '1' in [row, 1, ] is because there is
        ## guaranteed to be only one sample
        geno.calls <- gt[row, 1, ]
        allele.counts <- ad[row, 1, ]

        if (any(is.na(geno.calls))) {
            ## This shouldn't happen, so it is a sanity check
            stop("Some but not all genotype calls are NA")
        }

        if (all(is.na(geno.calls))) {
            ## TODO(Jason): I don't think this is necessary
            ref.calls <- 0  # number of reference reads
            alt.calls <- 0  # number of alternate reads
        } else {
            ref.calls <- allele.counts[1]  # number of reference reads
            alt.calls <- allele.counts[2]  # number of alternate reads
        }

        ## theoretical fraction of sites that should be heterozygous
        perc.het <- 0.5^(prefs$generation - 1)
        ## theoretical fraction of sites that should be homozygous
        perc.hom <- 1 - perc.het
        ## total probability of the sequence of reads occurring
        p <- perc.hom/2 * (1-rerr)^ref.calls * (rerr)^alt.calls +
             perc.hom/2 * (1-rerr)^alt.calls * (rerr)^ref.calls +
             perc.het   * (1/2)^(ref.calls + alt.calls)

        ## Calculate the emission probabilities for this site
        ref.prob <- perc.hom/2 * (1-rerr)**ref.calls * (rerr)**alt.calls / p
        alt.prob <- perc.hom/2 * (1-rerr)**alt.calls * (rerr)**ref.calls / p
        hom.prob <- perc.het   * (1/2)**(ref.calls + alt.calls) / p


        r <- pseudo.site.pos[row]  # pseudo position of site in range [0, 0.5]
        n <- 2 * prefs$generation  # number of boundaries
        for (block in 0:n) {
            k <- block
            for (state in 1:2) {
                superstate <- 3*block + state - 1
                p <- choose(n, k) * r^k * (1-r)^(n-k)

                ## TODO(Jason): This seems incorrect. If the parent is unknown (NA)
                ## then why should the probability of being that parent be the max
                ## of alt.prob and ref.prob?
                if (is.na(parent.geno[row, state])) {
                    print("ID:12345")
                    browser()
                    ret.val[row, superstate] <- max(alt.prob, ref.prob)
                } else if (parent.geno[row, state] == 0) {
                    ret.val[row, superstate] <- ref.prob * p
                } else if (parent.geno[row, state] == 1) {
                    ret.val[row, superstate] <- alt.prob * p
                } else {
                    stop("Parental genotype was not NA, 0, or 1")
                }

            }
            p <- choose(n, k) * r^k * (1-r)^(n-k)
            superstate <- 3*block + 2
            ret.val[row, superstate] <- hom.prob * p

        }
    }
    ret.val  # implicit return
}
