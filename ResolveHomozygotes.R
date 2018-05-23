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


##' Resolve heterozygous sites in the samples
##'
##' Returns a matrix of genotype calls for the samples such that the entry is 0
##'     if all calls were for the reference allele; 1 if all calls were for the
##'     alternate allele; and NA if there were calls for both the reference and
##'     alternate allele.
##' @title
##' @param vcf an object of class vcf
##' @param samples a vector of indices or names of variants
##' @param prefs a preferences object
##' @return a matrix of genotypes (0, 1, NA)
##' @author Jason Vander Woude
ResolveHomozygotes <- function(vcf, samples) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    ## TODO(Jason): We can set the parents to NA now, but they should be set
    ## back to the original call before turning the children imputation states
    ## back into genotype calls (I think). That is we don't want to call all
    ## children as NA just because one of the parents was unknown. This may not
    ## be an issue since we are filtering the parents, but if we upgrade the
    ## code to handle non-HWPB (homozygous within and polymorphic between) this
    ## may become important

    ## TODO(Jason): This is the step where we are making an absolute claim on
    ## the true genotype of the parents which means we should be accounting for
    ## a possible read error here. If the allelic depths are 1 and 100 for
    ## example, we are currently claiming that the site is heterozygous when in
    ## reality it is probably homozygous with a read error. In addition, if we
    ## are more sure about the parents at a site, we should be more heavily
    ## weighting the emission probabilities at that site in the viterbi. Right
    ## now all sites have equal weight, but it is possible that there are times
    ## we get the parents wrong. Again this is probably not a huge concern with
    ## the code as is, but especially if we begin utilizing non-HWPB data, this
    ## is a feature that should be incorporated.

    ## TODO(Jason): Again, if utilizing non-HWPB data, first remove false
    ## homozygosity because otherwise variants will only be caled homozygous of
    ## the variant seen or heterozygous even if they were actually called
    ## homozygous of the other.
    genotype <- Get(vcf, "GT", samples)
    allele.counts <- Get(vcf, "AD", samples, vcf$chrom.names)
    ret.val <- genotype[, , 1]  # initilize to first slice in 3rd dim
    for (r in 1:nrow(genotype)) {  # by chromosome site
        for (c in 1:ncol(genotype)) {  # by variant/sample
            ## Observed alleles from a specific site of a variant chromosome
            alleles <- genotype[r, c, ]
            if (any(is.na(alleles))) {  # One of the alleles is NA
                ret.val[r, c] <- NA
            } else if (all(alleles == alleles[1])) {  # Check if all are same
                ret.val[r, c] <- alleles[1]
            } else if (sum(allele.counts[r, c, ] != 0) == 1) {  # Only 1 counted
                ret.val[r, c] <- alleles[as.logical(allele.counts[r, c, ])]
            } else {  # Contradictory calls
                ret.val[r, c] <- NA
            }
        }
    }
    ret.val  # implicit return
}
