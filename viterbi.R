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


viterbi <- function(emission.probs, qs, prefs) {
    nstates <- nstates.probs(emission.probs)
    path.size <- nsites.probs(emission.probs)

    ## 3D array
    paths.tracker <- array(NA, dim=c(nstates, path.size, nstates))

    ## This will keep track of the overall probabilities of each of the {nstates}
    ## final paths, so that the probability of the final paths do not have to be
    ## computed again. The probabilities are initialized to the emission
    ## probabilities of the first site of the probs matrix which is the first
    ## row
    if (path.size < 1) {
        stop("viterbi requires that probs have > 0 rows")
    }
    probs.tracker <- log(emission.probs[1, ])

    ## Hard code the first column to the vector 1,2,...,nstates as an index
    ## This is what the generatePath function will need
    paths.tracker[, 1, ] <- diag(TRUE, nstates)

    if (path.size != 1) {  # if the path size is 1 just use the emission probs
        for (site in 2:path.size) {  # for each site in the path

            q <- qs[site - 1]

            ## log of the probability for each possible hidden state at this site
            probs.tracker <- sapply(1:nstates, function(state) {

                extension.probs <- sapply(1:nstates, function(i) {
                    ## log of the probability of being at state i before and
                    ## transitioning to the 'state' state
                    probs.tracker[i] + log(TransProb(i, state, q))
                })

                optimal.indices <- extension.probs==max(extension.probs)
                ## use <<- to assign to a variable outside the current scope
                paths.tracker[state, site, ] <<- optimal.indices
                max(extension.probs) + log(emission.probs[site, state])  # return prob
            })
        }
    }

    ## The code above has already computed the optimal path, but that
    ## information is encoded within the paths.tracker matrix and needs to be
    ## extracted. That is what generatePath will do when passed the path.tracker
    ## matrix and the indices of the optimal path.
    generatePath(paths.tracker, probs.tracker==max(probs.tracker))  # return best path
}
