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


# There is a seperate file for each number of sites because the equations are
# quite large and make editing the files difficult as is, so further editing
# burden was diminished by seperating the equations into seperate files. Thus
# the code here is necessary to call the correct function, though in general
# this would be considered bad design.

source("precomputed/precomputed_2_site.R")
source("precomputed/precomputed_3_site.R")
source("precomputed/precomputed_4_site.R")
source("precomputed/precomputed_5_site.R")
source("precomputed/precomputed_6_site.R")

f1.gamete.probs <- function(recomb.probs) {
    n.sites = length(recomb.probs) + 1

    if (n.sites == 2) {
        initial.2.site(recomb.probs)
    } else if (n.sites == 3) {
        initial.3.site(recomb.probs)
    } else if (n.sites == 4) {
        initial.4.site(recomb.probs)
    } else if (n.sites == 5) {
        initial.5.site(recomb.probs)
    } else if (n.sites == 6) {
        initial.6.site(recomb.probs)
    } else {
        stop("Unsupported number of recomb.probs")
    }
}

next.gamete.probs <- function(prev.gamete.probs, recomb.probs) {
    n.sites = length(recomb.probs) + 1

    if (n.sites == 2) {
        recursive.2.site(prev.gamete.probs, recomb.probs)
    } else if (n.sites == 3) {
        recursive.3.site(prev.gamete.probs, recomb.probs)
    } else if (n.sites == 4) {
        recursive.4.site(prev.gamete.probs, recomb.probs)
    } else if (n.sites == 5) {
        recursive.5.site(prev.gamete.probs, recomb.probs)
    } else if (n.sites == 6) {
        recursive.6.site(prev.gamete.probs, recomb.probs)
    } else {
        stop("Unsupported number of recomb.probs")
    }
}
