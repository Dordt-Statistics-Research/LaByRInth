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


## TODO(Jason): only require parallel on functions that use it (use lapply if not available)
require(abind)
source("image_functions.R")
source("utility_functions.R")
source("labyrinth_utility_functions.R")
source("viterbi.R")
source("generatePath.R")
source("VCF.R")
source("Get.R")
source("ResolveHomozygotes.R")
source("GetProbabilities.R")
source("LabyrinthImpute.R")
source("GenotypeFromDepth.R")



adjusted.ad <- function(){}


## gt should be length 2 array from gt field
prob.reads.given.geno <- function(n.ref, n.alt, gt, rerr) {
    if (length(gt) != 2)
        stop("length of gt must be 2")
    if (! all( gt %in% 0:1))
        stop("entries of gt must be 0 or 1")

    ## sum of genotypes, so 0 if hom ref, 1 if het, 2 if hom alt
    gt <- sum(gt)
    n <- n.ref + n.alt

    if (gt==0) {  # homozygous reference
        choose(n, n.ref) * (1-rerr)^n.ref * rerr^n.alt
    } else if (gt==1) {  # heterozygous
        choose(n, n.alt) * (1-rerr)^n.alt * rerr^n.ref
    } else {  # if (gt==2); homozygous alternate
        ## note that choose(n, n.ref) == choose(n, n.alt) by definition
        choose(n, n.ref) * (0.5)^(n.alt + n.ref)
    }
}


corrected.parents <- function(vcf, prefs) {
    
}





## outline
 take vcf and create adjustedAD
