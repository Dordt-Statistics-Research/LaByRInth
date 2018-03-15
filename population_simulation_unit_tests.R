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

## TEST 1: 'cross' and 'gamete' functions

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

## DONE 1: looks good up to here


## TEST 2: 'cross' and 'gamete' functions

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

## DONE 2: percentages come out as expected


## TEST 3: 'get.recomb.profile.fun'

## Check the appearance of the curve to verify that the functions are specified
## as intended
sites <- seq(from=1000, to=2000, by=5)
profile.fun <- get.recomb.profile.fun(peak=5, trough=1, d=300, sites=sites)
plot(sites, profile.fun(sites), type="l")

## DONE 3: plot shows the expected curve


## TEST 4: 'get.recomb.profile.fun'

## Check the appearance of the curve on an example of real site values
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
profile.fun <- get.recomb.profile.fun(peak=5, trough=0.1, d=5e7, sites=sites)
plot(sites, profile.fun(sites), type="p")

## DONE 4: plot shows the expected curve


## TEST 5: 'get.recomb.prob.fun'

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
peak <- 5
trough <- 0.0001
d <- 5e5

prob.fun <- get.recomb.prob.fun(peak=peak, trough=trough, d=d, sites=sites)
profile.fun <- get.recomb.profile.fun(peak=peak, trough=trough, d=d, sites=sites)

## find midpoint between every pair of sites
mid.sites <- sapply.pairs(sites, function(s1, s2) {mean(c(s1, s2))})

## find probability of observed recombination between these pairs
mid.site.probs <- sapply.pairs(sites, function(s1, s2) {prob.fun(s1, s2)})

## plot these probabilities; note that some of them could be above 1 if sites
## are close together and and in peak regions of the profile
plot(mid.sites, mid.site.probs, type="p")

## plot the profile over top; where sites are farther apart, the probability
## should tend to be higher; where the profile of the site is larger, the
## probability should be larger
points(sites, profile.fun(sites), type="p", col="green")

## plot histogram and density of the probabilities
hist(mid.site.probs)
plot(density(mid.site.probs))

## DONE 5: everything looks as it should; quite dependent on peak, trough, and d



## TEST 6: 

### ril.pop <- create.ril.pop(5, sites)
