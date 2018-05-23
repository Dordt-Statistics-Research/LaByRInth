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


## Determine which rows/sites have parents that are HWPB (homozygous within and
## polymorphic between)
GetRelevantProbabiltiesIndex <- function(vcf, chromosomes, parent.geno, prefs) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    rows <- vcf$variants[, "CHROM"] %in% chromosomes
    apply(parent.geno[rows, , drop=F], 1, function(calls) {
        ## parents should be distinct and non-NA
        !anyNA(calls) && length(calls) == length(unique(calls))
    })
}
