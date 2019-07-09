## Copyright 2018 Jason Vander Woude
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


##               __          ____        ____  _____
##              / /         / __ \      / __ \/_  _/
##             / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __
##            / /   / _  \/ _  _/\ \/ / _  _/ / / /  | / /_  __/ /_/ /
##           / /___/ /_/ / /_\ \  \  / / \ \_/ /_/ /||/ / / / / __  /
##          /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/
##
##                  L O W - C O V E R A G E   B I A L L E L I C
##                    R - P A C K A G E   I M P U T A T I O N
##


get.datasets <- function() {c("Lakin-Fuller", "HincII", "RsaI", "IBM-RIL")}

get.all.parents <- function() {
    list("Lakin-Fuller"  = c("LAKIN",      "FULLER"),
         "HincII"        = c("HincII_B73", "HincII_CG"),
         "RsaI"          = c("RsaI_B73",   "RsaI_CG"),
         "IBM-RIL"       = c("B73",        "M17"))
}

## I don't actually know what generation IBM-RIL is, so I'm just saying 7
get.generations <- function() {
    generations <- c(5, 2, 2, 7)
    names(generations) <- get.datasets()
    generations
}

## minimum call depths found with trial and error to get ~1% called sites masked
get.min.mask.depths <- function() {
    min.mask.depths <- c(7, 90, 12, 11)
    names(min.mask.depths) <- get.datasets()
    min.mask.depths
}

## directory associated with extra data files for package
extdata.dir <- function() {system.file("extdata", package = "labyrinth", mustWork = TRUE)}

## directory for imputation result files
result.dir <- function() {paste0(extdata.dir(), "/publication_results/", package = "labyrinth", mustWork = TRUE)}

original.file <- function(dataset.name) {
    file.name <- c("Lakin-Fuller" = "LakinFuller_GBSv2_20170509.vcf.gz",
                   "HincII"       = "HincII_ordered.vcf.gz",
                   "RsaI"         = "RsaI_ordered.vcf.gz",
                   "IBM-RIL"      = "IBM-RIL_ordered.vcf.gz"
      )[dataset.name]
    paste0(extdata.dir(), "datasets/original_files/", dataset.name)
}

## for the masked and filtered datasets
dataset.dir <- function(dataset.name) {
    paste0(extdata.dir(), "datasets/", dataset.name, "/")
}

imputed.dir <- function(dataset.name, algorithm) {
    if (algorithm %in% c("LaByRInth", "LB-Impute")) {
        paste0(result.dir(), "/", dataset.name, "/LaByRInth/")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}

filtered.file <- function(dataset.name) {
    paste0(dataset.dir(dataset.name), "/filtered.vcf.gz")
}

masked.file <- function(dataset.name, mask.ID = 1) {
    paste0(dataset.dir(dataset.name), "/masked_", mask.ID, ".vcf.gz")
}

uncompressed.masked.file <- function(dataset.name, mask.ID = 1) {
    paste0(dataset.dir(dataset.name), "/masked_", mask.ID, ".vcf")
}

imputed.parents.file <- function(dataset.name, algorithm, config = 1) {
    if (algorithm="LaByRInth") {
        paste0(imputed.dir(dataset.name, algorithm), "/parental-imputation-result_", config, ".rds")
    } else if (algorithm="LB-Impute") {
        paste0(imputed.dir(dataset.name, algorithm), "/imputed-parents_", config, ".rds")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}

imputed.progeny.file <- function(dataset.name, parental.config = 1, progeny.config = 1) {
    if (algorithm="LaByRInth") {
        paste0(imputed.dir(dataset.name, algorithm), "/imputed-progeny_", parental.config, "_", progeny.config, ".vcf.gz")
    } else if (algorithm="LB-Impute") {
        paste0(imputed.dir(dataset.name, algorithm), "/imputed-progeny_", parental.config, "_", progeny.config, ".vcf")
    } else {
        stop("Algorithm must be either 'LaByRInth' or 'LB-Impute'")
    }
}
