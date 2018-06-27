
 ## set.seed(0)

## recomb.probs <- runif(9, 0, 0.1)
## n.gen <- 2

## genos <- simulate.ril.geno.matrix(300, recomb.probs, n.gen)


## pop <- create.ril.pop(3, 5, recomb.probs)

## samples <- lapply(1:10, function(i) {SamplePopulation(pop, 0, 4)})


## pop.to.geno.matrix(pop)
## for (sample in samples) {
##     browser()
##     print(apply(sample, 1:2, paste0, collapse="/"))
##     writeLines(sample.to.vcf(sample, n.gen))
## }

## set.seed(0)
## base.name <- "10000skip512"
## ## sites <- get.sites("./analysis/1A.csv")  # sites from Lakin-Fuller chrom 1A
## sites <- seq(from=0, by=1e4, length.out=512)
## n.gen <- 2                               # F2
## recomb.dist <- 1e6                       # default recomb dist in the algorithm
## n.mem <- 312                             # 300 progeny
## read.err <- 0                            # prob of erroneous read
## geno.err <- 0                            # prob of assembly error
## coverage <- 1                            # coverage level on per-locus basis



## dists <-  diff(sites)
## recomb.probs <- sapply(dists, function(dist) {(1 - exp(-1.0 * dist / recomb.dist)) / 2})

## pop <- create.ril.pop(n.gen, n.mem, recomb.probs)
## population.to.vcf(pop, n.gen, pop.file, sites)
## sample <- sample.population(pop, read.err, coverage)
## sample.to.vcf(sample, n.gen, sample.file, sites)

source("population_simulation.R")

generate.pop.and.sample.vcfs <- function(base.name, sites, n.gen, recomb.dist,
                                         n.mem, read.err, geno.err, coverage) {
    dir <- paste0("datasets/sim-",
                  base.name, "-",
                  "F", n.gen, "-",
                  n.mem, "-",
                  length(sites), "-",
                  format(recomb.dist, scientific=FALSE), "-",
                  coverage, "-",
                  read.err, "-",
                  geno.err)

    if (dir.exists(dir))
        warning(paste("Directory", dir, "already exists."))

    dir.create(dir, recursive=TRUE)

    pop.file <- paste0(dir, "/sim-pop.vcf")
    sample.file <- paste0(dir, "/sim-sample.vcf")

    dists <-  diff(sites)
    recomb.probs <- sapply(dists,
                           function(dist) {
                               (1 - exp(-1.0 * dist / recomb.dist)) / 2
                           })

    pop <- create.ril.pop(n.gen, n.mem, recomb.probs)     # create population
    population.to.vcf(pop, n.gen, pop.file, sites)        # save population
    sample <- sample.population(pop, read.err, coverage)  # sample population
    sample.to.vcf(sample, n.gen, sample.file, sites)      # save sampled population

    c("pop" = pop.file, "sample" = sample.file)
}


prep.for.analysis <- function(orig, sampled, imputed, perfect) {
    load.vcf <- function(vcf, name) {
        if (! inherits(vcf, "vcfR")) {
            timer <- new.timer()
            display(0, "Loading ", name, " vcf")
            vcf <- read.vcfR(vcf, verbose=F)
            display(1, "Completed in ", timer(), "\n")
        }
        vcf
    }

    obj <- list()

    obj$orig <- load(orig, "original")
    obj$sampled <- load(sampled, "sampled")
    obj$imputed <- load(imputed, "imputed")

}


sites.from.dists <- function(dists) {
    sites <- rep(0, length(dists))

    for (i in seq_along(dists)) {
        sites[i+1] <- sites[i] + dists[i]
    }

    sites
}


run.1 <- function() {
    ## An F2
    set.seed(0)

    base.name <- "1-normal-dist"

    n.sites <- 500
    n.gen <- 2                               # F2
    recomb.dist <- 1e6                       # default recomb dist in the algorithm
    n.mem <- 312                             # 300 progeny
    read.err <- 0                            # prob of erroneous read
    geno.err <- 0                            # prob of assembly error
    coverage <- 1                            # coverage level on per-locus basis

    mean.dist <- 10000
    sd.dist <- 3000
    dists <- round(rnorm(n.sites-1, mean.dist, sd.dist))
    sites <- sites.from.dists(dists)

    generate.pop.and.sample.vcfs(base.name, sites, n.gen, recomb.dist,
                                 n.mem, read.err, geno.err, coverage)
}


run.2 <- function() {
    set.seed(0)

    base.name <- "2-normal-dist"

    n.sites <- 500
    n.gen <- 2                               # F2
    recomb.dist <- 1e6                       # default recomb dist in the algorithm
    n.mem <- 312                             # 300 progeny
    read.err <- 0                            # prob of erroneous read
    geno.err <- 0                            # prob of assembly error
    coverage <- 1                            # coverage level on per-locus basis

    mean.dist <- 100000
    sd.dist <- 30000
    dists <- round(rnorm(n.sites-1, mean.dist, sd.dist))
    sites <- sites.from.dists(dists)

    generate.pop.and.sample.vcfs(base.name, sites, n.gen, recomb.dist,
                                 n.mem, read.err, geno.err, coverage)
}


run.3 <- function() {
    ## An F5 population
    set.seed(0)

    base.name <- "3-normal-dist"

    n.sites <- 500
    n.gen <- 5                               # F2
    recomb.dist <- 1e6                       # default recomb dist in the algorithm
    n.mem <- 300                             # 300 progeny
    read.err <- 0                            # prob of erroneous read
    geno.err <- 0                            # prob of assembly error
    coverage <- 1                            # coverage level on per-locus basis

    mean.dist <- 10000
    sd.dist <- 3000
    dists <- round(rnorm(n.sites-1, mean.dist, sd.dist))
    sites <- sites.from.dists(dists)

    generate.pop.and.sample.vcfs(base.name, sites, n.gen, recomb.dist,
                                 n.mem, read.err, geno.err, coverage)
}



run.4 <- function() {
    ## An F5 population
    set.seed(0)

    base.name <- "4-normal-dist"

    n.sites <- 500
    n.gen <- 5                               # F2
    recomb.dist <- 1e6                       # default recomb dist in the algorithm
    n.mem <- 300                             # 300 progeny
    read.err <- 0                            # prob of erroneous read
    geno.err <- 0                            # prob of assembly error
    coverage <- 0.125                        # coverage level on per-locus basis

    mean.dist <- 10000
    sd.dist <- 3000
    dists <- round(rnorm(n.sites-1, mean.dist, sd.dist))
    sites <- sites.from.dists(dists)

    generate.pop.and.sample.vcfs(base.name, sites, n.gen, recomb.dist,
                                 n.mem, read.err, geno.err, coverage)
}



run.bulk.1 <- function() {
    ## An F5 population
    set.seed(0)

    base.name <- "bulk-run-1"

    params <- list(generations = c(2, 5),
                   distances = c(1e3, 1e5),
                   coverages = c(1, 4),
                   read.errs = c(0, 0.05),
                   geno.errs = c(0, 0.05))

    configs <- do.call(expand.grid, params)

    ## Actually run the things
    sapply(seq(nrow(configs)), function(config_index) {

        print(paste0(config_index))

        base.name <- paste0(base.name, "-", configs$distances[config_index])

        n.gen <- configs$generations[config_index]
        mean.dist <- configs$distances[config_index]
        coverage <- configs$coverages[config_index]
        read.err <- configs$read.errs[config_index]
        geno.err <- configs$geno.errs[config_index]

        n.sites <- 500                           # number of SNP markers
        recomb.dist <- 1e6                       # default recomb dist in the algorithm
        n.mem <- 100                             # number of progeny

        sd.dist <- mean.dist / 3
        dists <- round(rnorm(n.sites-1, mean.dist, sd.dist))
        dists[dists < 0] <- mean.dist  # prevent negative distances

        sites <- sites.from.dists(dists)

        files <- generate.pop.and.sample.vcfs(base.name, sites, n.gen, recomb.dist,
                                              n.mem, read.err, geno.err, coverage)

        out.file <- sub(".vcf", "-imputed.vcf.gz", files["sample"])
        if (file.exists(out.file)) {
            print(paste("File", out.file, "already exists"))
            return(NULL)
        }

        LabyrinthImpute(vcf = files["sample"],
                        parents = c("P_1", "P_2"),
                        generation = n.gen,
                        out.file = out.file,
                        use.fwd.bkwd = FALSE,
                        use.viterbi = TRUE,
                        calc.posteriors = TRUE,
                        read.err = read.err,
                        geno.err = geno.err,
                        recomb.dist = recomb.dist,
                        parallel = TRUE,
                        cores = 20)
    })
}



run.bulk.2 <- function() {
    ## An F5 population
    set.seed(0)

    base.name <- "bulk-run-2"

    params <- list(generations = c(2, 5),
                   distances = c(1e3, 1e5),
                   coverages = c(1),
                   read.errs = c(0),
                   geno.errs = c(0))

    configs <- do.call(expand.grid, params)

    ## Actually run the things
    sapply(seq(nrow(configs)), function(config_index) {

        print(paste0(config_index))

        base.name <- paste0(base.name, "-", configs$distances[config_index])

        n.gen <- configs$generations[config_index]
        mean.dist <- configs$distances[config_index]
        coverage <- configs$coverages[config_index]
        read.err <- configs$read.errs[config_index]
        geno.err <- configs$geno.errs[config_index]

        n.sites <- 500                           # number of SNP markers
        recomb.dist <- 1e6                       # default recomb dist in the algorithm
        n.mem <- 100                             # number of progeny

        sd.dist <- mean.dist / 3
        dists <- round(rnorm(n.sites-1, mean.dist, sd.dist))
        dists[dists < 0] <- mean.dist  # prevent negative distances

        sites <- sites.from.dists(dists)

        files <- generate.pop.and.sample.vcfs(base.name, sites, n.gen, recomb.dist,
                                              n.mem, read.err, geno.err, coverage)

        out.file <- sub(".vcf", "-imputed.vcf.gz", files["sample"])
        if (file.exists(out.file)) {
            print(paste("File", out.file, "already exists"))
            return(NULL)
        }

        LabyrinthImpute(vcf = files["sample"],
                        parents = c("P_1", "P_2"),
                        generation = n.gen,
                        out.file = out.file,
                        use.fwd.bkwd = FALSE,
                        use.viterbi = TRUE,
                        calc.posteriors = TRUE,
                        read.err = read.err,
                        geno.err = geno.err,
                        recomb.dist = recomb.dist,
                        parallel = TRUE,
                        cores = 20)
    })
}
