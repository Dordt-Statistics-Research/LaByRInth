% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labyrinth.R
\name{LabyrinthImputeParents}
\alias{LabyrinthImputeParents}
\title{Parental imputation}
\usage{
LabyrinthImputeParents(vcf, out.file, parents, generation,
  geno.err = 0.01, parent.het = 0.01, parallel = TRUE, cores = 4)
}
\arguments{
\item{vcf}{File path or vcfR object to impute.}

\item{out.file}{File path of output file.}

\item{parents}{Character vector with names of the two parents.}

\item{generation}{Numeric representing generation (e.g. 5 for F5 population).}

\item{geno.err}{Estimate of proportion of calls genotyped incorrectly.}

\item{parent.het}{Estimate of proportion of sites in each parent that are
truly heterozygous.}

\item{parallel}{Logical indicating if imputation should be run in parallel
or serial.}

\item{cores}{Numeric indicating how many sub-processes should be spawned if
running in parallel.}
}
\value{
A vcfR object with parental data replaced by imputation results.
}
\description{
Imputes both of the parents by maximizing likelihood of data across
all progeny at each site.
}
\examples{
input <- system.file(
    "extdata",
    "vcf-files",
    "filtered-lakin-fuller-sample.vcf",
    package = "LaByRInth",
    mustWork = TRUE)
output <- tempfile(fileext="-parental-result.rds")
print(output)
result <- LabyrinthImputeParents(
    vcf = input,
    out.file = output,
    parents = c("LAKIN", "FULLER"),
    generation = 5,
    geno.err = 0.015,
    parent.het = 0.005,
    parallel = FALSE,
    cores = 1)
}
\author{
Jason Vander Woude
}
