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


initial.2.site <- function(recomb.probs) {
    if (length(recomb.probs) != 1) {
        stop("recomb.probs is incorrect length")
    }

    p_0 <- recomb.probs[1]

    ## implicit return
    c(-1/4*p_0 + 1/2,
      1/4*p_0,
      1/4*p_0,
      -1/4*p_0 + 1/2)
}


recursive.2.site <- function(prev.gamete.probs, recomb.probs) {
    if (length(prev.gamete.probs) != 2**2) {
        stop("prev.gamete.probs is incorrect length")
    }
    if (length(recomb.probs) != 1) {
        stop("recomb.probs is incorrect length")
    }

    p_0 <- recomb.probs[1]

    q_0 <- prev.gamete.probs[1]
    q_1 <- prev.gamete.probs[2]
    q_2 <- prev.gamete.probs[3]
    q_3 <- prev.gamete.probs[4]

    ## implicit return
    c(1/2*p_0*q_1*q_2 - 1/2*p_0*q_0*q_3 + q_0^2 + q_0*q_1 + q_0*q_2 + q_0*q_3,
      -1/2*p_0*q_1*q_2 + 1/2*p_0*q_0*q_3 + q_0*q_1 + q_1^2 + q_1*q_2 + q_1*q_3,
      -1/2*p_0*q_1*q_2 + 1/2*p_0*q_0*q_3 + q_0*q_2 + q_1*q_2 + q_2^2 + q_2*q_3,
      1/2*p_0*q_1*q_2 - 1/2*p_0*q_0*q_3 + q_0*q_3 + q_1*q_3 + q_2*q_3 + q_3^2)
}
