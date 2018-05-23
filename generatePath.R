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


##' Generate a path from a path tracker
##'
##' Given a path tracker (which is a specific type of array produced by viterbi
##'     algorithm) and indices representing the final hidden states of the path,
##'     compute the paths that ended in those states. The returned path will be
##'     a vector of integers where integers which are powers of 2 represent the
##'     the state with the 0-based index which is the log_2 of the number.
##'     (e.g. 1,2,4,8,16,... represent the states associated with the 1st, 2nd,
##'     3rd, 4th, 5th, etc row)
##' @title
##' @param path.tracker the boolean 3D array in which the path info is embedded
##' @param boolean.indices which rows do optimal paths end at
##' @return the path
##' @author Jason Vander Woude
generatePath <- function(path.tracker, boolean.indices) {
    if (length(boolean.indices) != dim(path.tracker)[1]) {
        stop("boolean.indices has wrong length")
    }
    path <- rep(NA, dim(path.tracker)[2])
    base <- 0:(dim(path.tracker)[1] - 1)
    powers <- 2**base
    for (i in length(path):1) {
        path[i] <- sum(powers[boolean.indices])
        slice <- path.tracker[, i, ]
        slice.optimal <- slice[boolean.indices, , drop=F]
        boolean.indices <- as.logical(apply(slice.optimal, 2, vec.or))
    }
    path  # implicitly return path
}
