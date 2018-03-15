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

## TEST 1:

sites <- c(1:16)
n.sites <- length(sites)

p1 <- p2 <- list()  # empty object initialization
class(p1) <- c("genome", class(p1))
class(p2) <- c("genome", class(p2))

p1$cA <- p2$cA <- cA <- rep(0, times=n.sites)
p1$cB <- p2$cB <- cB <- rep(1, times=n.sites)
p1$sites <- p2$sites <- sites
fun <- function(site1, site2) {0.25}  # 50% chance of crossover between
                                      # every pair of sites
print(p1)
print(p2)
print(cross(p1, p2, fun))

## DONE: looks good up to here


## TEST 2:

## Try to average 1 physical recombination
## In actuallity, the probability of a crossover happening between the inner
## chromatids is 1 - (1 - (1/(n.sites-1)))^(n.sites-1) = 0.6447356 so the
## probability of not having a crossover between these two is 0.3552644. Given
## that there is a 50% chance of selecting one of the unmodified parent
## chromosomes, there is in total a 0.5 + 0.5*0.3552644 = 0.6776322 chance that
## a given homologous chromosome in the offspring will be identical to one of
## the parental homologous chromosomes.
fun <- function(site1, site2) {1 / (n.sites - 1) / 2}

## Check this assumption over 10000 trials
## Each of the homologs of each cross should be identical to one of the parents
## about 67% of the time
n <- 10000
res <- lapply(1:n, function(i) {
    cross(p1, p2, fun)
})

test.cAs <- lapply(res, function(r) {r$cA})
test.cBs <- lapply(res, function(r) {r$cB})

perc <- function(vec) {sum(vec) / length(vec)}
## Both of these values should be about 67%
print(perc(sapply(test.cAs, function(vec) {identical(vec, cA) || identical(vec, cB)})))
print(perc(sapply(test.cBs, function(vec) {identical(vec, cA) || identical(vec, cB)})))

## DONE: percentages come out as expected










## get sites
### sites <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
### ril.pop <- create.ril.pop(5, sites)
