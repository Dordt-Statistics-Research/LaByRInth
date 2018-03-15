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


sites <- c(1:16)
n.sites <- length(sites)

p1 <- p2 <- list()  # empty object initialization
class(p1) <- c("genome", class(p1))
class(p2) <- c("genome", class(p2))
p1$cA <- p2$cA <- rep(1, times=n.sites)
p1$cB <- p2$cB <- rep(6, times=n.sites)
p1$sites <- p2$sites <- sites
fun <- function(site1, site2) {0.25}  # 50% chance of crossover between
                                      # every pair of sites
print(p1)
print(p2)
print(cross(p1, p2, fun))

## DONE: looks good up to here


## get sites
### sites <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
### ril.pop <- create.ril.pop(5, sites)
