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


initial.3.site <- function(recomb.probs) {
    if (length(recomb.probs) != 2) {
        stop("recomb.probs is incorrect length")
    }

    p_0 <- recomb.probs[1]
    p_1 <- recomb.probs[2]

    ## implicit return
    c(1/4*(p_0 - 1)*(p_1 - 1) + 1/4,
      -1/4*p_0*(p_1 - 1),
      1/4*p_0*p_1,
      -1/4*(p_0 - 1)*p_1,
      -1/4*(p_0 - 1)*p_1,
      1/4*p_0*p_1,
      -1/4*p_0*(p_1 - 1),
      1/4*(p_0 - 1)*(p_1 - 1) + 1/4)
}


recursive.3.site <- function(prev.gamete.probs, recomb.probs) {
    if (length(prev.gamete.probs) != 2**3) {
        stop("prev.gamete.probs is incorrect length")
    }
    if (length(recomb.probs) != 2) {
        stop("recomb.probs is incorrect length")
    }

    p_0 <- recomb.probs[1]
    p_1 <- recomb.probs[2]

    q_0 <- prev.gamete.probs[1]
    q_1 <- prev.gamete.probs[2]
    q_2 <- prev.gamete.probs[3]
    q_3 <- prev.gamete.probs[4]
    q_4 <- prev.gamete.probs[5]
    q_5 <- prev.gamete.probs[6]
    q_6 <- prev.gamete.probs[7]
    q_7 <- prev.gamete.probs[8]

    ## implicit return
    c(-p_0*p_1*q_1*q_4 - 1/2*p_0*p_1*q_3*q_4 + p_0*p_1*q_0*q_5 + 1/2*p_0*p_1*q_2*q_5 - 1/2*p_0*p_1*q_1*q_6 + 1/2*p_0*p_1*q_0*q_7 + 1/2*p_0*q_1*q_2 - 1/2*p_0*q_0*q_3 + 1/2*p_0*q_1*q_4 + 1/2*p_1*q_1*q_4 + 1/2*p_1*q_2*q_4 + 1/2*p_1*q_3*q_4 - 1/2*p_0*q_0*q_5 - 1/2*p_1*q_0*q_5 - 1/2*p_1*q_0*q_6 + 1/2*p_0*q_1*q_6 - 1/2*p_0*q_0*q_7 - 1/2*p_1*q_0*q_7 + q_0^2 + q_0*q_1 + q_0*q_2 + q_0*q_3 + q_0*q_4 + q_0*q_5 + q_0*q_6 + q_0*q_7,
      p_0*p_1*q_1*q_4 + 1/2*p_0*p_1*q_3*q_4 - p_0*p_1*q_0*q_5 - 1/2*p_0*p_1*q_2*q_5 + 1/2*p_0*p_1*q_1*q_6 - 1/2*p_0*p_1*q_0*q_7 - 1/2*p_0*q_1*q_2 + 1/2*p_0*q_0*q_3 - 1/2*p_0*q_1*q_4 - 1/2*p_1*q_1*q_4 + 1/2*p_0*q_0*q_5 + 1/2*p_1*q_0*q_5 + 1/2*p_1*q_2*q_5 + 1/2*p_1*q_3*q_5 - 1/2*p_0*q_1*q_6 - 1/2*p_1*q_1*q_6 + 1/2*p_0*q_0*q_7 - 1/2*p_1*q_1*q_7 + q_0*q_1 + q_1^2 + q_1*q_2 + q_1*q_3 + q_1*q_4 + q_1*q_5 + q_1*q_6 + q_1*q_7,
      -1/2*p_0*p_1*q_3*q_4 + 1/2*p_0*p_1*q_2*q_5 - 1/2*p_0*p_1*q_1*q_6 - p_0*p_1*q_3*q_6 + 1/2*p_0*p_1*q_0*q_7 + p_0*p_1*q_2*q_7 - 1/2*p_0*q_1*q_2 + 1/2*p_0*q_0*q_3 - 1/2*p_1*q_2*q_4 + 1/2*p_0*q_3*q_4 - 1/2*p_0*q_2*q_5 - 1/2*p_1*q_2*q_5 + 1/2*p_1*q_0*q_6 + 1/2*p_1*q_1*q_6 + 1/2*p_0*q_3*q_6 + 1/2*p_1*q_3*q_6 - 1/2*p_0*q_2*q_7 - 1/2*p_1*q_2*q_7 + q_0*q_2 + q_1*q_2 + q_2^2 + q_2*q_3 + q_2*q_4 + q_2*q_5 + q_2*q_6 + q_2*q_7,
      1/2*p_0*p_1*q_3*q_4 - 1/2*p_0*p_1*q_2*q_5 + 1/2*p_0*p_1*q_1*q_6 + p_0*p_1*q_3*q_6 - 1/2*p_0*p_1*q_0*q_7 - p_0*p_1*q_2*q_7 + 1/2*p_0*q_1*q_2 - 1/2*p_0*q_0*q_3 - 1/2*p_0*q_3*q_4 - 1/2*p_1*q_3*q_4 + 1/2*p_0*q_2*q_5 - 1/2*p_1*q_3*q_5 - 1/2*p_0*q_3*q_6 - 1/2*p_1*q_3*q_6 + 1/2*p_1*q_0*q_7 + 1/2*p_1*q_1*q_7 + 1/2*p_0*q_2*q_7 + 1/2*p_1*q_2*q_7 + q_0*q_3 + q_1*q_3 + q_2*q_3 + q_3^2 + q_3*q_4 + q_3*q_5 + q_3*q_6 + q_3*q_7,
      p_0*p_1*q_1*q_4 + 1/2*p_0*p_1*q_3*q_4 - p_0*p_1*q_0*q_5 - 1/2*p_0*p_1*q_2*q_5 + 1/2*p_0*p_1*q_1*q_6 - 1/2*p_0*p_1*q_0*q_7 - 1/2*p_0*q_1*q_4 - 1/2*p_1*q_1*q_4 - 1/2*p_1*q_2*q_4 - 1/2*p_0*q_3*q_4 - 1/2*p_1*q_3*q_4 + 1/2*p_0*q_0*q_5 + 1/2*p_1*q_0*q_5 + 1/2*p_0*q_2*q_5 + 1/2*p_1*q_0*q_6 + 1/2*p_0*q_5*q_6 + 1/2*p_1*q_0*q_7 - 1/2*p_0*q_4*q_7 + q_0*q_4 + q_1*q_4 + q_2*q_4 + q_3*q_4 + q_4^2 + q_4*q_5 + q_4*q_6 + q_4*q_7,
      -p_0*p_1*q_1*q_4 - 1/2*p_0*p_1*q_3*q_4 + p_0*p_1*q_0*q_5 + 1/2*p_0*p_1*q_2*q_5 - 1/2*p_0*p_1*q_1*q_6 + 1/2*p_0*p_1*q_0*q_7 + 1/2*p_0*q_1*q_4 + 1/2*p_1*q_1*q_4 + 1/2*p_0*q_3*q_4 - 1/2*p_0*q_0*q_5 - 1/2*p_1*q_0*q_5 - 1/2*p_0*q_2*q_5 - 1/2*p_1*q_2*q_5 - 1/2*p_1*q_3*q_5 + 1/2*p_1*q_1*q_6 - 1/2*p_0*q_5*q_6 + 1/2*p_1*q_1*q_7 + 1/2*p_0*q_4*q_7 + q_0*q_5 + q_1*q_5 + q_2*q_5 + q_3*q_5 + q_4*q_5 + q_5^2 + q_5*q_6 + q_5*q_7,
      1/2*p_0*p_1*q_3*q_4 - 1/2*p_0*p_1*q_2*q_5 + 1/2*p_0*p_1*q_1*q_6 + p_0*p_1*q_3*q_6 - 1/2*p_0*p_1*q_0*q_7 - p_0*p_1*q_2*q_7 + 1/2*p_1*q_2*q_4 + 1/2*p_1*q_2*q_5 - 1/2*p_1*q_0*q_6 - 1/2*p_0*q_1*q_6 - 1/2*p_1*q_1*q_6 - 1/2*p_0*q_3*q_6 - 1/2*p_1*q_3*q_6 - 1/2*p_0*q_5*q_6 + 1/2*p_0*q_0*q_7 + 1/2*p_0*q_2*q_7 + 1/2*p_1*q_2*q_7 + 1/2*p_0*q_4*q_7 + q_0*q_6 + q_1*q_6 + q_2*q_6 + q_3*q_6 + q_4*q_6 + q_5*q_6 + q_6^2 + q_6*q_7,
      -1/2*p_0*p_1*q_3*q_4 + 1/2*p_0*p_1*q_2*q_5 - 1/2*p_0*p_1*q_1*q_6 - p_0*p_1*q_3*q_6 + 1/2*p_0*p_1*q_0*q_7 + p_0*p_1*q_2*q_7 + 1/2*p_1*q_3*q_4 + 1/2*p_1*q_3*q_5 + 1/2*p_0*q_1*q_6 + 1/2*p_0*q_3*q_6 + 1/2*p_1*q_3*q_6 + 1/2*p_0*q_5*q_6 - 1/2*p_0*q_0*q_7 - 1/2*p_1*q_0*q_7 - 1/2*p_1*q_1*q_7 - 1/2*p_0*q_2*q_7 - 1/2*p_1*q_2*q_7 - 1/2*p_0*q_4*q_7 + q_0*q_7 + q_1*q_7 + q_2*q_7 + q_3*q_7 + q_4*q_7 + q_5*q_7 + q_6*q_7 + q_7^2)
}
