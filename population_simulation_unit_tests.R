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

recomb.probs <- rep(0.25, n.sites-1)

print(p1)
print(p2)
print(cross(p1, p2, recomb.probs))

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
recomb.probs <- rep(1 / (n.sites - 1) / 2, n.sites - 1)

## Check this assumption over 10000 trials
## Each of the homologs of each cross should be identical to one of the parents
## about 67% of the time
n <- 10000
res <- lapply(1:n, function(i) {
    cross(p1, p2, recomb.probs)
})

test.cAs <- lapply(res, function(r) {r$cA})
test.cBs <- lapply(res, function(r) {r$cB})

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

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
profile.fun <- get.recomb.profile.fun(peak=peak, trough=trough, d=d, sites=sites)

## find midpoint between every pair of sites
mid.sites <- sapply.pairs(sites, function(s1, s2) {mean(c(s1, s2))})

## plot these probabilities; note that some of them could be above 1 if sites
## are close together and and in peak regions of the profile
plot(mid.sites, recomb.probs, type="p")

## plot the profile over top; where sites are farther apart, the probability
## should tend to be higher; where the profile of the site is larger, the
## probability should be larger
points(sites, profile.fun(sites), type="p", col="green")

## plot histogram and density of the probabilities
hist(recomb.probs)
plot(density(recomb.probs))

## DONE 5: everything looks as it should; quite dependent on peak, trough, and d



## TEST 6: 'create.ril.pop'

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5

recomb.probs <- rep(0.5, n.sites-1)  # pick a number
n.tests <- 7
n.gen.to.test <- 5
n.mem <- 1000

set.seed(1)
## all of the following values printed should be very close to 0
for (gen in 2:n.gen.to.test) {
    homozygosities <- sapply(1:n.tests, function(i) {
        ril.pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, sites, recomb.probs)
        get.pop.homozygosity(ril.pop)  # implicit return
    })
    print(paste("Generation", gen))
    print(homozygosities - (1 - 0.5^(gen-1)))
}

## DONE 6: populations seem to be generating correctly


## TEST 7: simulate realistic population and save as vcf file

## note: peak=5, trough=0.1, d=5e5 looks similar to current lakin fuller imputation
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.1
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
n.gen <- 5
n.mem <- 100
cov <- 0.25

set.seed(1)
ril.pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(ril.pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(ril.pop, coverage=cov), sites)

source("functions.R")  # need vcf functions

WriteVCF(full.vcf, "./analysis/simulated_population/full_f5_test_1.vcf")
WriteVCF(sample.vcf, "./analysis/simulated_population/sample_f5_test_1.vcf")

vcf.to.pgm("./analysis/simulated_population/full_f5_test_1.vcf", c("P1", "P2"),
           "./analysis/simulated_population/full_f5_test_1.pgm")
vcf.to.pgm("./analysis/simulated_population/sample_f5_test_1.vcf", c("P1", "P2"),
           "./analysis/simulated_population/sample_f5_test_1.pgm")

## DONE 7: trying various parameters not shown here, everything seems to
## consistently look as expected


## TEST 8: try imputing sample data

source("functions.R")  # need vcf functions

## note: peak=5, trough=0.1, d=5e5 looks similar to current lakin fuller imputation
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.1
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
n.gen <- 5
n.mem <- 100
cov <- 0.25

set.seed(1)
ril.pop <- create.ril.pop(n.gen=n.gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(ril.pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(ril.pop, coverage=cov), sites)

WriteVCF(sample.vcf, "./analysis/simulated_population/sample_f5_test_1.vcf")

## Create a preferences objects containing all preferences
prefs <- list()
class(prefs)            <- "prefs"

## Algorithm parameters
prefs$recomb.double     <- TRUE
prefs$read.err          <- 0.05
prefs$genotype.err      <- 0.05
prefs$recomb.dist       <- 1e6
## should the GT info be inferred from the AD info
prefs$use.only.ad       <- FALSE # THIS DIFFERS FROM STANDARD PREFS
## Should non-imputed sites be in the output VCF file
prefs$leave.all.calls   <- TRUE
prefs$parents           <- c("P1", "P2")

prefs$states            <- 3     # currently only support for 2 parents
## TODO(Jason): implement this feature
prefs$ref.alt.by.parent <- FALSE # Should the reference and alternate be
                                        # switched in the output so that parent 1
                                        # is always reference and parent 2 is
                                        # always alternate

## Logistic parameters
prefs$quiet             <- FALSE
prefs$cores             <- 4
prefs$parallel          <- TRUE
prefs$write             <- TRUE
prefs$out.file          <- "./analysis/simulated_population/sample_f5_test_1_imputed.vcf"

LabyrinthImpute("./analysis/simulated_population/sample_f5_test_1.vcf",
    parents = c("P1", "P2"),
    out.file = "./analysis/simulated_population/sample_f5_test_1_imputed.vcf")
imputed.vcf <- VCF("./analysis/simulated_population/sample_f5_test_1_imputed.vcf", prefs)

LabyrinthImpute("./analysis/simulated_population/sample_f5_test_1.vcf",
    recomb.dist = 500e6,
    parents = c("P1", "P2"),
    out.file = "./analysis/simulated_population/sample_f5_test_1_imputed_2.vcf")
prefs$out.file          <- "./analysis/simulated_population/sample_f5_test_1_imputed.vcf"
imputed.2.vcf <- VCF(prefs$out.file, prefs)

LabyrinthImpute("./analysis/simulated_population/sample_f5_test_1.vcf",
    parents = c("P1", "P2"),
    out.file = "./analysis/simulated_population/sample_f5_test_1_imputed_3.vcf")
prefs$out.file <- "./analysis/simulated_population/black_hole"
imputed.3.vcf <- VCF("./analysis/simulated_population/sample_f5_test_1_imputed_3.vcf", prefs)


res <- checkAccuracy(full.vcf, imputed.vcf)
res.2 <- checkAccuracy(full.vcf, imputed.2.vcf)
res.3 <- checkAccuracy(full.vcf, imputed.3.vcf)

with(res, sum(quality=="correct") / length(quality))
with(res, sum(quality=="correct") / (length(quality) - sum(quality=="skipped")))

with(res.2, sum(quality=="correct") / length(quality))
with(res.2, sum(quality=="correct") / (length(quality) - sum(quality=="skipped")))

vcf.to.pgm("./analysis/simulated_population/sample_f5_test_1_imputed.vcf", c("P1", "P2"),
           "./analysis/simulated_population/sample_f5_test_1_imputed.pgm")


vcf.to.pgm("./analysis/simulated_population/sample_f5_test_1_imputed_3.vcf", c("P1", "P2"),
           "./analysis/simulated_population/sample_f5_test_1_imputed_3.pgm")

## using 500 times the recombination distance resulted in exactly the same imputation?
## only 7 heterozygous calls were made, so it is probably being penalized too
## much

## what is the true type when we get it wrong
summary(with(res, gen.type[quality=="wrong"]))
summary(with(res, gen.type[quality=="skipped"]))  ## 10% of skips are heterozygous

## recursive imputation??? need depth data though

## interesting that when imputing without a priori multiplying the heterozygous
## probabilities by 0.05, there are more het calls which makes sense, but there
## are fewer correct calls and few partial calls (which also makes sense) and
## fewer correct calls and twice as many wrong calls

## NOT DONE 8: more to do here


## TEST 9: sequential generations

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 2
n.mem <- 100

set.seed(1)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)

dir <- "analysis/simulated_population/breeding"
basename <- "test_00004"

full.vcf <- FullPopulationVCF(pop, sites)
filename <- paste0(basename, "_f", gen, ".ppm")
ppm.from.vcf(full.vcf, paste0(dir, filename))
for (gen in 3:6) {
    pop <- self.population(pop, recomb.probs)
    full.vcf <- FullPopulationVCF(pop, sites)
    filename <- paste0(basename, "_f", gen, ".ppm")
    ppm.from.vcf(full.vcf, paste0(dir, filename))
}

## NOT DONE 9


## TEST 10: imputing an F2

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 2  # F2
n.mem <- 100
cov = 1/16
qs <- compute.qs(recomb.probs, gen)

set.seed(1)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00005"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), file, qs=qs, parallel=TRUE)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)

## TEST 11

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 5  # F2
n.mem <- 100
cov = 0.25
qs <- compute.qs(recomb.probs, gen)

set.seed(1)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00006"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), file, qs=qs, parallel=TRUE)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)
with(acc, summary(call.type[quality=="wrong" | quality=="partial"]))

## TEST 12

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.1
d <- 5e6

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 2  # F2
n.mem <- 100
cov = 0.125
qs <- compute.qs(recomb.probs, gen)

set.seed(1)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00007"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), file, qs=qs, parallel=TRUE)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)
with(acc, summary(call.type[quality=="wrong" | quality=="partial"]))

## TEST 13

sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.01
d <- 5e6

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 2  # F2
n.mem <- 100
cov = 0.125
qs <- compute.qs(recomb.probs, gen)

set.seed(1)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00008"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), file, qs=qs, parallel=TRUE)
file <- paste0(dir, "/", basename, "_imputed_mod_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)
with(acc, summary(call.type[quality=="wrong" | quality=="partial"]))

## TEST 14:
p <- 0.1
q <- (2-p)/4
for (i in 1:5) {
    print(q)
    q <- -p*q/2 + p/8 + q
}
print(c(q^2, 2*q*(1/2-q), (1/2-q)^2, 2*q^2+2*(1/2-q)^2))



## TEST 15: F5 population
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 5  # F5
n.mem <- 100
cov = 1/16

set.seed(1)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00009"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), gen, file, parallel=TRUE)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)


## TEST 15: F5 population
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 5  # F5
n.mem <- 100
cov = 1/4

set.seed(2)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00010"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), gen, file, parallel=TRUE)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)
## These results don't look good


## TEST 15: F5 population
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5
recomb.dist <- 50 / trough * 1e6

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 2  # F3
n.mem <- 100
cov = 1/4

set.seed(2)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00011"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), gen, file,
                               parallel=FALSE, use.uncalled.sites=TRUE,
                               write=FALSE, recomb.dist=recomb.dist)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)
## 99.68% accuracy


## TEST 16: F3 population
sites  <- get.sites("analysis/simulated_population/lakin_fuller_sites/1A.csv")
n.sites <- length(sites)
peak <- 5
trough <- 0.0001
d <- 5e5
recomb.dist <- 50 / trough * 1e6

recomb.probs <- get.recomb.probs(peak=peak, trough=trough, d=d, sites=sites)
gen <- 3  # F3
n.mem <- 100
cov = 1/8

set.seed(2)
pop <- create.ril.pop(n.gen=gen, n.mem=n.mem, recomb.probs)
full.vcf <- FullPopulationVCF(pop, sites)
sample.vcf <- SamplePopulationVCF(SamplePopulation(pop, coverage=cov), sites)

dir <- "analysis/simulated_population/unit_tests"
basename <- "test_00012"
ensure_writability(paste0(dir, "/", basename))

file <- paste0(dir, "/", basename, "_full_f", gen, ".ppm")
ppm.from.vcf(full.vcf, file)
file <- paste0(dir, "/", basename, "_sample_f", gen, ".ppm")
ppm.from.vcf(sample.vcf, file)
file <- paste0(dir, "/", basename, "_imputed_f", gen, ".vcf")
imputed.vcf <- LabyrinthImpute(sample.vcf, c("P1", "P2"), gen, file,
                               parallel=FALSE, write=FALSE, recomb.dist=recomb.dist,
                               use.uncalled.sites=TRUE, weight.uncalled=TRUE, weight.called=TRUE, seperate.het=TRUE)

file <- paste0(dir, "/", basename, "_imputed_f", gen, ".ppm")
ppm.from.vcf(imputed.vcf, file)

acc <- CheckAccuracy(full.vcf, imputed.vcf)
Summarize(acc)
## use.uncalled.sites          weight    _uncalled_      _called_        seperate.het
##               TRUE            TRUE                                            TRUE       78% (no almost no het calls)
##              FALSE            TRUE                                            TRUE       78% (no almost no het calls)
##              FALSE           FALSE                                            TRUE       78% (no almost no het calls)
##               TRUE           FALSE                                            TRUE       78% (no almost no het calls)
##               TRUE           FALSE                                           FALSE       22% (almost all het)
##               TRUE            TRUE                                           FALSE       59% (63% het calls)
##              FALSE            TRUE                                           FALSE       ERROR
##                                                                                          99% (splotchy though)
##               TRUE                           TRUE        FALSE               FALSE       22% (almost all het)
##               TRUE                          FALSE         TRUE               FALSE       22% (almost all het)
##               TRUE                           TRUE         TRUE               FALSE       22% (almost all het)
##               TRUE                          FALSE        FALSE               FALSE       22% (almost all het)

##               TRUE                           TRUE        FALSE                TRUE       77.71565% (almost no het)
##               TRUE                          FALSE         TRUE                TRUE       77.65267% (almost no het)
##               TRUE                           TRUE         TRUE                TRUE       77.65267% (almost no het)
##               TRUE                          FALSE        FALSE                TRUE       


