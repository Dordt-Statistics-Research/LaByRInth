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


transProb <- function(a, b, q) {
    if (a == 1) {
        if (b == 1) {
            q^2
        } else if (b == 2) {
            (1/2 - q)^2
        } else if (b == 3) {
            q*(1/2 - q)
        } else if (b == 4) {
            q*(1/2 - q)
        } else {
            stop(paste("TransProb failed. a:", a, "b:", b))
        }
    } else if (a == 2) {
        if (b == 1) {
            (1/2 - q)^2
        } else if (b == 2) {
            q^2
        } else if (b == 3) {
            q*(1/2 - q)
        } else if (b == 4) {
            q*(1/2 - q)
        } else {
            stop(paste("TransProb failed. a:", a, "b:", b))
        }
    } else if (a == 3) {
        if (b == 1) {
            q*(1/2 - q)
        } else if (b == 2) {
            q*(1/2 - q)
        } else if (b == 3) {
            q^2
        } else if (b == 4) {
            (1/2 - q)^2
        } else {
            stop(paste("TransProb failed. a:", a, "b:", b))
        }
    } else if (a == 4) {
        if (b == 1) {
            q*(1/2 - q)
        } else if (b == 2) {
            q*(1/2 - q)
        } else if (b == 3) {
            (1/2 - q)^2
        } else if (b == 4) {
            q^2
        } else {
            stop(paste("TransProb failed. a:", a, "b:", b))
        }
    } else {
        stop(paste("TransProb failed. a:", a, "b:", b))
    }
}
