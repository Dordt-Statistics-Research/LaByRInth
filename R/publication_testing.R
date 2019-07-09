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


datasets <- c("Lakin-Fuller", "HincII", "RsaI", "IBM-RIL")

all.parents <- list("Lakin-Fuller"  = c("LAKIN",      "FULLER"),
                    "HincII"        = c("HincII_B73", "HincII_CG"),
                    "RsaI"          = c("RsaI_B73",   "RsaI_CG"),
                    "IBM-RIL"       = c("B73",        "M17"))
names(all.parents) <- datasets

## I don't actually know what generation IBM-RIL is, so I'm just saying 7
generations <- c(5, 2, 2, 7); names(generations) <- datasets

## minimum call depths found with trial and error to get ~1% called sites masked
min.mask.depths <- c(7, 90, 12, 11); names(min.mask.depths) <- datasets


original.file <- function(dataset.name) {
    c("Lakin-Fuller" = "../data/original_files/LakinFuller_GBSv2_20170509.vcf.gz",
      "HincII"       = "../data/original_files/HincII_ordered.vcf.gz",
      "RsaI"         = "../data/original_files/RsaI_ordered.vcf.gz",
      "IBM-RIL"      = "../data/original_files/IBM-RIL_ordered.vcf.gz"
      )[dataset.name]
}

dataset.dir <- function(dataset.name) {
    paste0("../data/", dataset.name, "/")
}

imputed.dir <- function(dataset.name) {
    paste0(dataset.dir(dataset.name), "/LaByRInth/")
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

parental.file <- function(dataset.name, config = 1) {
    paste0(imputed.dir(dataset.name), "/parental_imputation_result_", config, ".rds")
}

imputed.file <- function(dataset.name, parental.config = 1, progeny.config = 1) {
    paste0(imputed.dir(dataset.name), "/imputed_", parental.config, "_", progeny.config, ".vcf.gz")
}

analysis.file <- function(dataset.name, progeny.config = 1) {
    paste0(imputed.dir(dataset.name), "/analysis_summary_", progeny.config, ".csv")
}

mimicked.file <- function(dataset.name, mimic.ID = 1) {
    paste0(dataset.dir(dataset.name), "/mimicked_", mimic.ID, ".vcf.gz")
}
