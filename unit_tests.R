## prefs
    prefs <- list()
    class(prefs)            <- "prefs"

    ## Algorithm parameters
    prefs$parents           <- c("LAKIN", "FULLER")
    prefs$resolve.conflicts <- FALSE
    prefs$recomb.double     <- FALSE
    prefs$read.err          <- 0.05
    prefs$genotype.err      <- 0.05
    prefs$recomb.err        <- 0.05
    prefs$recomb.dist       <- 100000
    prefs$min.markers       <- 1
    prefs$states            <- length(prefs$parents) + 1

    ## Logistics
    prefs$quiet             <- FALSE
    prefs$cores             <- 4
    prefs$parallel          <- TRUE
    prefs$write             <- TRUE

## end prefs



## Check the generatePath algorithm
path.tracker.a <- matrix(c(3,1,1,1,1,1,3,2,2,2,3,3),byrow=T,nrow=3)
path.tracker.a
generatePath(path.tracker.a, 1)
generatePath(path.tracker.a, 2)
generatePath(path.tracker.a, 3)

probs <- matrix(c(0.5,0.4,0.1,0.3,0.3,0.4,0.7,0.1,0.2,0.5,0.4,0.1),nrow=4,byrow=T)
probs
transProb <- function(a, b, dists, prefs) {}


## Check
entry <- "0/0:1,0:1:66:0,3,36"




## Check vcf file creation
alphabet <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")
reorder(alphabet, 3)
array(reorder(alphabet, 3), dim=c(2,2,3))
mat <- matrix(c("a|b|c", "d/e/f", "g/h|i", "j|k/l"), nrow=2)
mat
submat <- gsub("\\|", "/", mat)
submat
alphabet <- unlist(strsplit(submat, "/"))
alphabet
array(reorder(alphabet, 3), dim=c(2,2,3))




## vcf
#vcfobj <- VCF(file="../LakinFuller_GBSv2_20170509.vcf")
vcfobj <- VCF(file="../vcf_files/vcf25x2.vcf")
resolved.parents <- ResolveHomozygotes(vcfobj, test.prefs$parents)
parent.6b.depths <- Get(vcfobj, "AD", test.prefs$parents, "6B")
lakin.6b.depths <- Get(vcfobj, "AD", "LAKIN", "6B")
parent.map <- GetProbabilities(vcfobj, "LAKIN", resolved.parents, test.prefs)



## imputation logistics
vcf.file <- "../LakinFuller_004.vcf"
vcf.file.big <- "../LakinFuller_001.vcf"
vcf.obj <- VCF(vcf.file)
vcf.obj.big <- VCF(vcf.file.big)

## Full imputation
prefs$parallel <- FALSE
prefs$cores <- 40
prefs$quiet <- TRUE
prefs$write <- FALSE
prefs$viterbi.testing <- TRUE
## TODO(Jason): sometimes the path is all 1's. Why should this be the case if it is never 2's or 3's for instance in "U6202-140"
test.impute <- LabyrinthImputeHelper(vcf.obj, prefs)


## Chromosome imputation to test viterbi correctness
parent.geno <- ResolveHomozygotes(vcf.obj, prefs$parents)
sample <- "U6202-080"
chrom <- "1B"

data <- Get(vcf.obj, "GT", sample, chrom)


test.chrom.impute <- LabyrinthImputeChrom(vcf.obj, sample, chrom, parent.geno, prefs)
## test.sample <- LabyrinthImputeSample(vcf.obj, sample, parent.geno, prefs)


Display <- function(vecs) {
    invisible(sapply(vecs, function(vec) {
        vec <- sapply(vec, function(elem) {
            if (is.na(elem)) {
                "-"
            } else if (is.logical(elem)) {
                if (elem) {
                    "#"
                } else {
                    "."
                }
            } else {
                elem
            }
        })
        writeLines(paste0(vec, collapse=""))
    }))
}

lakin <- Get(vcf.obj, "GT", "LAKIN", chrom)
fuller <- Get(vcf.obj, "GT", "FULLER", chrom)
Display(list(lakin,
             fuller,
             data[, , 1],
             test.chrom.impute))

StatVerify <- function(vcf.obj, parent.geno, prefs) {
    variants <- vcf.obj$variant.names[!(vcf.obj$variant.names %in% prefs$parents)]
    print("A")
    lapply(vcf.obj$chrom.names, function(ch) {
        print("B")
        unlist(mclapply(variants, function(va) {
            res <- LabyrinthImputeChrom(vcf.obj, va, ch, parent.geno, prefs)
            print(res)
            res
        }, mc.cores=40))
    })
    ## keep track of the percentage of the time that the forward and the reverse
    ## viterbi algorithm agree
}

prefs$parallel <- TRUE
prefs$viterbi.testing <- TRUE
library(parallel)
parent.geno <- ResolveHomozygotes(vcf.obj.big, prefs$parents)
stats.verify2 <- StatVerify(vcf.obj.big, parent.geno, prefs)


## Updated generatePath
slices <- list()
slices[[1]] <- matrix(c(T,F,F,
                        F,T,F,
                        F,F,T), 3, byrow=T)
slices[[2]] <- matrix(c(F,T,T,
                        F,F,T,
                        F,F,T), 3, byrow=T)
slices[[3]] <- matrix(c(T,F,F,
                        T,F,F,
                        F,T,F), 3, byrow=T)
slices[[4]] <- matrix(c(T,F,F,
                        F,T,T,
                        F,T,F), 3, byrow=T)
slices[[5]] <- matrix(c(T,F,F,
                        T,T,F,
                        F,F,T), 3, byrow=T)
slices[[6]] <- matrix(c(T,F,F,
                        F,T,F,
                        F,F,T), 3, byrow=T)
slices[[7]] <- matrix(c(F,T,F,
                        F,F,T,
                        F,F,T), 3, byrow=T)
slices[[8]] <- matrix(c(F,T,F,
                        F,F,T,
                        T,F,F), 3, byrow=T)
path.tracker <- array(NA, dim=c(3,8,3))
for (i in 1:8) {
    path.tracker[, i, ] <- slices[[i]]
}

generatePath(path.tracker, c(T,T,F))
generatePath(path.tracker, c(F,T,T))
generatePath(path.tracker, c(F,F,T))


## Checking imputation results for soft imputation (imputation that uses
## 1: homozygous parent 1
## 2: homozygous parent 2
## 3: homozygous unknown parent
## 4: heterozygous
## 5: one of the alleles is from parent 1 and the other is unknown
## 6: one of the alleles is form parent 2 and the other is unknown
## 7: unknown (./.)
Display(list( test.impute[, 1]))
writeLines("\n")
Display(list(test.impute[, 2]))
writeLines("\n")
Display(list(test.impute[, 3]))




## Small Test Imputation
options(width = 160)

#### Prefs initialization
prefs <- list()
class(prefs)            <- "prefs"

#### Algorithm parameters
prefs$parents           <- c("LAKIN", "FULLER")
prefs$resolve.conflicts <- FALSE
prefs$recomb.double     <- FALSE
prefs$read.err          <- 0.05
prefs$genotype.err      <- 0.05
prefs$recomb.err        <- 0.05
prefs$recomb.dist       <- 100000
prefs$min.markers       <- 1
prefs$states            <- length(prefs$parents) + 1

#### Logistics
prefs$quiet             <- FALSE
prefs$cores             <- 4
prefs$parallel          <- TRUE
prefs$write             <- TRUE

#### Run the imputation
#vcf.obj <- VCF(file="../vcf_files/vcf_5x1.vcf")
#vcf.obj <- VCF(file="../vcf_files/vcf_25x1_simple_sequential.vcf")
#vcf.obj <- VCF(file="../vcf_files/vcf_25x1_simple_magnitude_sequential.vcf")
#vcf.obj <- VCF(file="../vcf_files/vcf_25x1_simple_magnitude_sequential_lakin_misread.vcf")
#vcf.obj <- VCF(file="../vcf_files/sandesh_lakin_fuller_clean.vcf")
#vcf.obj <- VCF(file="../vcf_files/vcf_25x2.vcf")
vcf.obj <- VCF(file="../vcf_files/sandesh_lakin_fuller_clean_5_variant.vcf")
test.impute <- LabyrinthImputeHelper(vcf.obj, prefs)

parent.geno <- ResolveHomozygotes(vcf.obj, prefs$parents)
sample <- "U6202-080"
chrom <- "1A"
data <- Get(vcf.obj, "GT", sample, chrom)

lakin <- Get(vcf.obj, "GT", "LAKIN", chrom)
fuller <- Get(vcf.obj, "GT", "FULLER", chrom)

test.chrom.impute <- LabyrinthImputeChrom(vcf.obj, sample, chrom, parent.geno, prefs)

Display(list(parent.geno[,"LAKIN"]))
Display(list(parent.geno[,"FULLER"]))
Display(list(test.chrom.impute))
