# What is LaByRInth?

LaByRInth stands for <b>L</b>ow-coverage <b>B</b>iallelic <b>R</b> <b>I</b>mputation. It is a genetic imputation algorithm based heavily on the <a href="https://github.com/dellaporta-laboratory/LB-Impute">LB-Impute algorithm</a>). It is designed to impute the missing sites in a low-coverage biallelic population where the parental genotypes are known and primarily homozygous.

The original purpose of the project was to port LB-Impute from Java to R to make it more practical to use in genetics research. Making the algorithm available in R will hopefully increase the ease of its use in larger projects. We also wanted to implement some of the improvements that the authors of LB-Impute mentioned in their <a href="https://www.ncbi.nlm.nih.gov/pubmed/26715670">paper</a>,so we have designed LaByRInth with these in mind. One of these was computational performance, so we have designed LaByRInth to be able to run in parallel on multiple cores to increase the speed of an imputation. A second improvement that they offered was to use the entire genome of a sample when finding the optimal imputation instead of a sliding window across the genome which we have done by using an implementation of the forward-backward algorithm.

LaByRInth version 3 includes several significant methodological improvements as well. It is tailored to inbred populations from two common ancestral, homozygous parents. While LB-Impute used a general model of the expected genetic structure of inbred lines, LaByRInth is able to impute different generations of inbred lines differently to improve accuracy (i.e. an F2 is imputed differently than an F5). Using a standard biological model, genetic recombination probabilities have been pre-computed symbolically and are included with LaByRInth; the probabilities appropriate to the specified generation of inbreeding are instantiated and used at run-time. This results in more accurate imputation which is still relatively fast.



# Installation

LaByRInth can be installed directly as an R package on your system by using the `devtools` package to install it from here as shown below:
```r
## Install devtools from CRAN if it is not installed already
if (!require("devtools")) install.packages("devtools")
## Use devtools to install LaByRInth
devtools::install_github('Dordt-Statistics-Research/LaByRInth')
```
It can then be used by running the following command as usual:
```r
## Attach LaByRInth
library(LaByRInth)
```



# Usage

LaByRInth provides the following four functions:
1. `LabyrinthFilter`
2. `LabyrinthImputeParents`
3. `LabyrinthImputeProgeny`
4. `LabyrinthQualityAssurance`
For any imputation, these four functions should all be used (and in the order listed). Each function requires an input file and will generate an output file. Each output file should be used as the input file for the next function on the above list. The reason that these four functions were not combined into a single function is that having intermediate result files may be helpful if any part of the imputation is to be rerun on different parameters. LabyrinthFilter actually may not be necessary for every dataset (although running it won't hurt). The function LabyrinthImputeParents will check that the VCF file has only biallelic sites, and if any non-biallelic sites are detected, it will exit with a message stating that LabyrinthFilter must be run. LabyrinthImputeParents will impute the parents of the dataset so that the imputation of the progeny (offspring) will be more accurate. This is typically the most time consuming step of the process. The result file will contain all information about the population that is necessary for LabyrinthImputeProgeny. LabyrinthImputeProgeny will use this file and generate a VCF file with all members of the population (parents and progeny) imputed at every site. However, because every site is imputed, there may be some that are of low quality, so LabyrinthQualityAssurance should be used to remove any calls with low expected quality (the expected quality information is obtained from the forward-backward algorithm).



# Help

The LaByRInth package includes more detailed help information on each of these four functions. In an R session, try any of the following commands:
```r
help(LabyrinthFilter)
help(LabyrinthImputeParents)
help(LabyrinthImputeProgeny)
help(LabyrinthQualityAssurance)
```



# Examples

The LaByRInth package also includes a small sample dataset which is part of the F<sub>5</sub> Lakin-Fuller population obtained from the Department of Agronomy at Kansas State University (this dataset is available in full [here](./paper_analysis/data/original_files/LakinFuller_GBSv2_20170509.vcf.gz)). In addition, there are examples included of how to run each function. You can view these examples in the help files as shown above, or run them directly by trying any of the following commands:
```r
example(LabyrinthFilter)
example(LabyrinthImputeParents)
example(LabyrinthImputeProgeny)
example(LabyrinthQualityAssurance)
```
Below is an example of how a dataset should be processed from start to finish using LaByRInth.

## Example 1
```r
library(LaByRInth)

dir.create("LaByRInth_example")
original.file <- system.file(
                     "extdata",
                     "vcf-files",
                     "original-lakin-fuller-sample.vcf",
                     package = "LaByRInth",
                     mustWork = TRUE
		 )
filtered.file    <- "./LaByRInth_example/filtered.vcf.gz"
parental.file    <- "./LaByRInth_example/parental.rds"
all.imputed.file <- "./LaByRInth_example/all-imputed.vcf.gz"
quality.file     <- "./LaByRInth_example/quality-imputed.vcf.gz"

parents <- c("LAKIN", "FULLER")

filtered.result <- LabyrinthFilter(
                       vcf = original.file,
                       out.file = filtered.file,
                       parents = parents,
                       require.hom.poly = FALSE  # should generally be false
		   )

parental.result <- LabyrinthImputeParents(
                       vcf = filtered.file,      # or vcf = filtered.result
                       out.file = parental.file,
                       parents = parents,
                       generation = 5,           # Lakin-Fuller is F5
                       geno.err = 0.015,         # estimate
                       parent.het = 0.005,       # estimate
                       parallel = FALSE,         # if able, set TRUE
                       cores = 1                 # more cores is faster
		   )

all.imputed.result <- LabyrinthImputeProgeny(
                          parental = parental.file,    # or parental = parental.result
                          out.file = all.imputed.file,
                          parallel = FALSE,            # if able, set TRUE
                          cores = 1                    # more cores is faster
                      )

quality.result <- LabyrinthQualityAssurance(
                      vcf = all.imputed.file,  # or vcf = all.imputed.result
                      out.file = quality.file,
                      min.posterior = 0.8,     # require at least 80% probabilty of correctness
                      parallel = FALSE,        # if able, set TRUE
                      cores = 1                # more cores is faster
                  )
```
