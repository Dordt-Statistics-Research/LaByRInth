# What is LaByRInth?

LaByRInth stands for <b>L</b>ow-coverage <b>B</b>iallelic <b>R</b> <b>I</b>mputation. It is a genetic imputaton algorithm based heavily on the <a href="https://github.com/dellaporta-laboratory/LB-Impute">LB-Impute algorithm</a>). It is designed to impute the missing sites in a low-coverage biallelic population where the parental genotypes are known.

The original purpose of the project was to port LB-Impute from Java to R to make it more practical to use in genetics research. Making the algorithm available in R will hopefully increase the ease of its application into larger projects. We also wanted to implement some of the improvements that the authors of LB-Impute mentioned in their <a href="https://www.ncbi.nlm.nih.gov/pubmed/26715670">paper</a>,so we have designed LaByRInth with these in mind. One of these was computational performance, so we have designed LaByRInth to be able to run in parallel on multiple cores to increase the speed of an imputation. A second improvement that they offered was to use the entire genome of a sample when finding the optimal imputation instead of a sliding window accross the genome which we have done by optimizing the implementation of the viterbi algorithm. We intend to make further optimizations to both the theoretical method of the imputation and to its implementation.



# Install

All of the code required to use LaByRInth currently exists in a single file in the project called `functions.R`. After <a href="https://github.com/Dordt-Statistics-Research/LaByRInth/archive/master.zip">downloading</a> and unzipping LaByRInth, simply add `functions.R` into your project directory. Then, in you can load all of the LaByRInth functions into an R script by including this line in your R script:
```r
source("functions.R")
```



# Usage

The following is the list of LaByRInth functions that are designed for you to use. Any of the other functions in `functions.R` are designed for use within the file and should not need to be used.



## LabyrinthImpute

This is the main function to do an imputation. Its function definition is this:

```r
LabyrinthImpute <- function(file, parents, out.file="", use.only.ad=TRUE,
                            leave.all.calls=TRUE, ref.alt.by.parent=TRUE,
                            recomb.double=TRUE, read.err=0.05,
                            genotype.err=0.05, recomb.dist=1e6,
                            write=TRUE, parallel=TRUE, cores=4,
                            quiet=FALSE)
```

## LabyrinthFilter

This function will filter a VCF file. That is, it will remove all of the sites from the file in which the parents are not homozygous within and polymorphic between. This can be used before running LabyrinthImpute, though it is not strictly necessary to do so. Specifically, if `leave.all.calls` is set to false, LabyrinthImpute will produce the same result on the original file and the filtered file (but it will take more time if you choose to cleaning the file first). The difference will be if you have `leave.all.calls` set as true, which may be useful for purposes of comparing the results to the original file. If this is the case, it may make sense to filter the file first. Its function definition is as follows:

```r
LabyrinthFilter <- function(file, parents, out.file="", use.only.ad=TRUE,
                            leave.all.calls=TRUE, ref.alt.by.parent=TRUE,
                            recomb.double=TRUE, read.err=0.05,
                            genotype.err=0.05, recomb.dist=1e6,
                            write=TRUE, parallel=TRUE, cores=4,
                            quiet=FALSE)
```

For both of these functions, the parameters descriptions are below.

### file
This is the input file to be imputed. It must be a VCF format file.

### parents
This is a vector containing strings representing the two parents of the progeny. The parents must be included in the VCF file.

### out.file
This is the name of the file where the resulting VCF file should be saved. If no output file is not specified, LaByRInth will save the results to the same directory as the input file, and it will prepend the text "LaByRInth_" (if LabyrinthImpute is called) or "filtered_" (if LabyrinthFilter is called) to input file name.

### use.only.ad
If true, the genotypes will be inferred from the allelic depth data in the input file rather than directly using the genotypes specified in the input file. The reason for this is that some programs that generate the VCF files my put the most likely genotype instead of the genotype matching the allele reads. There is probably no reason to use false for this parameter.

### leave.all.calls
LaByRInth will only use sites in the input file where the parents are homozygous and polymorphic between (the parents are different genotypes). If this parameter is true, these sites will appear in the output file, but the file will not include calls for any of these sites. It this parameter is false, only the sites that are homozygous within and polymorphic between will be included in the output file.

### ref.alt.by.parent
This parameter is currently not used. In the future it may be possible to specify that in the output file, the reference alleles are swapped so that at all sites that are homozygous within and polymorphic between, parent 1 will be homozygous reference and parent 2 will be homozygous alternate.

### recomb.double
If true, this indicates that probability of a double recombination (homozygous to homozygous) between sites is the square of the probability of a single recombination. If false, double and single recombination events will have equal probability. This may be useful for inbred populations.

### read.err
This is the probability that any given read of an allele is erroneous.

### genotype.err
This is the probability that any given genotype call is erroneous. This is distinct from `read.err`. The <a href="https://www.ncbi.nlm.nih.gov/pubmed/26715670">LB-Impute paper</a> states that while a read error represents the "likelihood that any one read is erroneous, genotyping errors independent of coverage also may affect a data set. These genotyping errors include misalignment of reads to the reference genome, resulting in the incorrect placement of a genotype, and an unannotated paralogous artifact."

### recomb.dist
This is the expected distance of 50cM.

### write
This specifies whether the results should be written to a file. It is probably most useful to leave this set to true. The reason you may want to set it to false is because `LabyrinthImpute` returns an R object that can be used for analysis within an R script; if that is the only object needed, then it may be desirable to save a file.

### parallel
If true, the imputation process will run in parallel. If false, it will run in serial. Parallelization was implemented with the 'parallel' R package available on CRAN (Comprehensive R Archive Network). If you want to run in parallel, you will have to install this package. If you are using a Windows machine, this package will not work[^1] and you may have to set this to false and run in serial.

[^1]: The "parallel" package utilizes the Unix fork command. Because Apple operating systems are based on BSD which is based on Unix, and because Linux systems are based on Unix, this package will work for both of these types of operating system. Windows, though, was created in a vary different manner, and thus it does not implement this functionality. There is little that can be done about this.

### cores
If `parallel` is set to true, this indicates how many cores LaByRInth will try to use. If the number is larger than the number of cores on your computer, it will still run fine, but if the value is set excessively large it could cause slow performance.

### quiet
Currently this option does nothing. Eventually, if it is true, the console output (the information printed to the screen indicating the progress of the imputation) will be limited so that no text appears. False will indicate that the normal messages should appear.



## vcf.to.pgm and result.to.pgm

These functions can be used to generate a <a href="https://en.wikipedia.org/wiki/Netpbm_format">pgm</a> image file from a vcf file or LaByRInth result respectively. This can be useful to easily see a visual representation of the sample genomes without having to use other software. The ouptut will be a pgm image that has a black pixel for any site that is homozygous and matches `parents[1]`; a dark grey pixel for any site that is homozygous and matches `parents[2]`; a light grey pixel for any site that is heterozygous; and a white pixel for any site that is not called. Creating an image from the result of LabyrinthImpute will be much faster than creating an image from a VCF file.

```r
vcf.to.pgm <- function(vcf.file, parents, out.file)

result.to.pgm <- function(labyrinth.result, parents, out.file)
```

### vcf.file
This is a VCF file to make an image of.

### labyrinth.result
This is an object returned as the result of LabyrinthImpute.

### parents
This is a vector containing strings representing the two parents of the progeny. The parents must be included in the VCF file.

### out.file
This is the name of the file where the resulting image will be saved. The file name should end with the extension ".pgm" so that an image viewer on your computer will know how to display it.



## make.diff.ppm

This will produce a red and green <a href="https://en.wikipedia.org/wiki/Netpbm_format">ppm</a> image file showing which pixels of two pgm image files differ. Green indicates the respective pixel is the same, and red indicates that the pixel differs between the images. This is useful for graphically seeing how an imputation looks compared to the original file or comparing imputations with different parameters or even two different programs.

```r
make.diff.ppm <- function(image1, image2, out.image) {
```

### image1
This is a pgm image file produced either by `vcf.to.pgm` or `result.to.pgm`

### image2
This is another pgm image file produced either by `vcf.to.pgm` or `result.to.pgm`

### out.image
This is the name of the file where the resulting image will be saved. The file name should end with the extension ".ppm" so that an image viewer on your computer will know how to display it.



# Examples



## Example 1

In this example, assume that the parents are named "parent_1" and "parent_2" in the VCF file "example_01.vcf" and that you want the output to be saved as "example_01_result.vcf". Then the following code is all that is needed.

```r
# load the LaByRInth functions
source("functions.R")

# set the parents
parents <- c("parent_1", "parent_2")

# save the results to "example_01_result.vcf"
LabyrinthImpute("example_01.vcf", parents, "example_01_result.vcf")
```



## Example 2

Now assume that you first want to filter the file first so that you can compare the results of the imputation with the original VCF file.

```r
source("functions.R")

parents <- c("parent_1", "parent_2")

# save filtered results to "example_01_filtered.vcf"
LabyrinthFilter("example_02.vcf", parents, "example_01_filtered.vcf")

# impute the filtered vcf file and save the results to "example_01_result.vcf"
LabyrinthImpute("example_02_filtered.vcf", parents, "example_01_result.vcf")
```



## Example 3

LaByRInth allows for saving the result as an object in R. Eventually, the exact nature of this structure will be documented, but if you would like to use this feature you can see the details in the source code for `functions.R`.

```r
source("functions.R")

parents <- c("parent_1", "parent_2")

# save the results into the "labyrinth_result_A" variable and prevent the results
# from being saved to a file
labyrinth_result_A <- LabyrinthImpute("example_03.vcf", parents, write=FALSE)

# run the imputation with different parameters
labyrinth_result_B <- LabyrinthImpute("example_03.vcf", parents, write=FALSE, recomb.double=FALSE)

# create images representing the two different imputations
result.to.pgm(labyrinth_result_A, parents, "image_A.pgm")
result.to.pgm(labyrinth_result_B, parents, "image_B.pgm")

# create an image comparison for the results
make.diff.ppm("image_A.pgm", "image_B.pgm", "difference_A_B.ppm")
```



# Known Bugs

1. LaByRInth uses the `sink()` function in R, and if you are using LaByRInth interactively in an R session, and you stop seeing any output at all it is possible that something has gone wrong. Try typing "sink()" and hitting <kbd>Enter</kbd> in the terminal/console/prompt repeatededly until R presents an error message. You may have to do this 3-4 times.

2. Another result of `sink()` is that you will typically see a warning message like the following when LabyrinthImpute finishes. This is "normal" and it should not pose a concern, though it will be fixed.
> Warning message:
> closing unused connection 3 (/tmp/Rtmp08zWwM/file4b7e6a870b6c)
