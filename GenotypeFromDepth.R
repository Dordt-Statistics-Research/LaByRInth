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


GenotypeFromDepth <-  function(allelic.depths) {
    ad <- allelic.depths

    if (length(ad) != 2 && !all(is.na(ad))) {
        stop("GenotypeFromDepth does not support non-biallelic reads")
    }

    if (all(is.na(ad))) {
        genotype <- NA  # Why am I not using c(NA, NA) instead of NA
    } else if (all(ad == 0)) {
        genotype <- NA  # Why am I not using c(NA, NA) instead of NA
    } else if (all(ad != 0)) {
        genotype <- c(0, 1)
    } else if (ad[1] == 0) {
        genotype <- c(1, 1)
    } else {  # ad[2] == 0
        genotype <- c(0, 0)
    }
    genotype  # implicit return
}
