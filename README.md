# What is LaByRInth?

LaByRInth stands for <b>L</b>ow-coverage <b>B</b>iallelic <b>R</b> <b>I</b>mputation. It is a genetic imputation algorithm based heavily on the <a href="https://github.com/dellaporta-laboratory/LB-Impute">LB-Impute algorithm</a>). It is designed to impute the missing sites in a low-coverage biallelic population where the parental genotypes are known.

The original purpose of the project was to port LB-Impute from Java to R to make it more practical to use in genetics research. Making the algorithm available in R will hopefully increase the ease of its use in larger projects. We also wanted to implement some of the improvements that the authors of LB-Impute mentioned in their <a href="https://www.ncbi.nlm.nih.gov/pubmed/26715670">paper</a>,so we have designed LaByRInth with these in mind. One of these was computational performance, so we have designed LaByRInth to be able to run in parallel on multiple cores to increase the speed of an imputation. A second improvement that they offered was to use the entire genome of a sample when finding the optimal imputation instead of a sliding window across the genome which we have done by using the standard dynamic programming implementation of the Viterbi algorithm.

LaByRInth version 2 includes a significant methodological change as well. It is tailored to inbred populations from two common ancestral, homozygous parents. While LB-Impute used a general model of the expected genetic structure of inbred lines, LaByRInth is able to impute different generations of inbred lines differently to improve accuracy (i.e. an F2 is imputed differently than an F5). Using a standard biological model, genetic recombination probabilities were pre-computed symbolically and are included with LaByRInth; the probabilities appropriate to the specified generation of inbreeding are instantiated at run-time.



# Install

The code for the latest release can be downloaded <a href="https://github.com/Dordt-Statistics-Research/LaByRInth/releases">here</a>. Unzip the files and add them into your project directory and load all of the LaByRInth functions into an R script by including this line in your R script:
```r
source("labyrinth.R")
```



# Usage

The following is the list of LaByRInth functions that are designed for you to use. Any of the other functions in `labyrinth.R` are designed for use within the file and should not need to be used.



## LabyrinthFilter

This function will filter a VCF file. That is, it will remove all of the sites from the file in which the parents are not homozygous within and polymorphic between. This can be used before running LabyrinthImpute; if you try to run LabyrinthImpute on a file that does not meet these conditions, the function will exit and request that you run LabyrinthFilter before proceeding. Its function definition is as follows:

```r
LabyrinthFilter <- function(vcf, parents, out.file)
```

### vcf
This is the input file to be filtered. It must be a either a VCF format file that can be read by the vcfR library or a vcfR object.

### parents
This is a vector containing strings representing the two parents of the progeny. The parents must be included in the VCF file.

### out.file
This is the name of the file where the resulting VCF file should be saved. The function will append the text ".vcf.gz" to indicate that the output file is a gzipped vcf format file, the type of file the vcfR library generates.



## LabyrinthImpute

This is the main function to do an imputation. Its function definition is this:

```r
LabyrinthImpute <- function (vcf, parents, generation, out.file,
                             use.fwd.bkwd=FALSE, use.viterbi=TRUE,
                             calc.posteriors=TRUE,
                             read.err=0.05, geno.err=0.05,
                             recomb.dist=1e6, parallel=TRUE,
                             cores=4)
```

### vcf
This is the input file to be imputed. It must be a either a VCF format file that can be read by the vcfR library or a vcfR object.

### parents
This is a vector containing strings representing the two parents of the progeny. The parents must be included in the VCF file.

### generation
This is the generation of inbreeding. For example, "2" indicates an F2 plant and "5" indicates an F5 plant. LaByRInth must have a corresponding file in the `new-transition-probs` in order to impute the specified generation. Currently, generations 2-6 are supported. In order to impute greater generations, a corresponding file must be generated. See the `README` file in `new-transition-probs` for instructions on how to do this.

### out.file
This is the name of the file where the resulting VCF file should be saved. The function will append the text ".vcf.gz" to indicate that the output file is a gzipped vcf format file, the type of file the vcfR library generates.

### use.fwd.bkwd and use.viterbi
To encourage awareness, one of these parameters must be `TRUE` and the other must be `FALSE`. The Viterbi algorithm and the forward-backward algorithm are two algorithms used for solving hidden Markov models optimally. Each is distinct in what it considers optimal. For a given pair of homologous chromosomes in a given population member, the Viterbi algorithm will find the genotype that is most probable based on the read depths and theoretical probability of each genotype occurring biologically. The forward-backward algorithm will (efficiently) compute the probability of every genotype based on the read depths and theoretical probability of each genotype occurring biologically, then for every SNP/site in the homologous pair it will find which of the three states (homozygous reference, homozygous alternate, or heterozygous) is most probable. For further explanation, <a href="https://stats.stackexchange.com/questions/31746/what-is-the-difference-between-the-forward-backward-and-viterbi-algorithms">this</a> is a good place to start. In practice, we have found that the viterbi algorithm will perform better.

### calc.posteriors
If set to `TRUE`, marginal posterior probabilities will be included in `out.file` in the "GP" field as phred-scaled values according to the <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">VCF specification</a>. The forward-backward algorithm is capable of computing marginal posterior probabilities for each state (homozygous reference, homozygous alternate, or heterozygous) for each SNP/site in a homologous chromosome pair. Regardless of the choice for `use.fwd.bkwd` and `use.viterbi`, these posteriors can be calculated. If `use.fwd.bkwd` is `TRUE` then virtually no extra time will be required to calculate the posterior probabilities. If `use.viterbi` is `TRUE`, the imputation will take approximately twice as long to complete if calculating the posteriors.

### read.err
This is the probability that any given read of an allele is erroneous at the haplotype level.

### genotype.err
This is the probability that any given genotype call is erroneous. This is distinct from `read.err`. The <a href="https://www.ncbi.nlm.nih.gov/pubmed/26715670">LB-Impute paper</a> states that while a read error represents the "likelihood that any one read is erroneous, genotyping errors independent of coverage also may affect a data set. These genotyping errors include misalignment of reads to the reference genome, resulting in the incorrect placement of a genotype, and an unannotated paralogous artifact." How this error rate is used by LaByRInth is different than how it is used in LB-Impute, but the idea is the same.

### recomb.dist
This is the physical distance in base pairs (bp) corresponding to a 50 centimorgan (cM) genetic/map distance.

### parallel
If `TRUE`, the imputation process will run in parallel. If `FALSE`, it will run in serial. Parallelization was implemented with the "parallel" R package available on CRAN (Comprehensive R Archive Network). If you want to run in parallel, you will have to install this package. If you are using a Windows machine, this package will not work[^1] and you may have to set this to `FALSE` and run in serial.

[^1]: The "parallel" package utilizes the Unix fork command. Because Apple operating systems are based on BSD which is based on Unix, and because Linux systems are based on Unix, this package will work for both of these types of operating system. Windows, though, was created in a different manner, and thus it does not implement this functionality. There is little that can be done about this.

### cores
If `parallel` is set to `TRUE`, this indicates how many cores LaByRInth will try to use. If the number is larger than the number of cores on your computer, it will still run fine, but if the value is set excessively large it could cause slow performance.



# Examples



## Example 1

In this example, assume that the parents are named "parent_1" and "parent_2" in the VCF file "example_01.vcf" and that you want the output to be saved as "example_01_result.vcf.gz". Then the following code is all that is needed.

```r
# load the LaByRInth functions
source("labyrinth.R")

# set the parents
parents <- c("parent_1", "parent_2")

# filter the file
LabyrinthFilter("example_01.vcf", parents, "example_01_filtered")

# save the results to "example_01_result.vcf.gz"
LabyrinthImpute("example_01.filtered.vcf.gz", parents, "example_01_result")
```
