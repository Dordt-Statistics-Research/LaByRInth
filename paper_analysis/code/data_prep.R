source("../../labyrinth.R")
source("../../labyrinth_sim.R")
source("../../labyrinth_analysis.R")



################################################################################
#                              DECLARE FILE NAMES                              #
################################################################################

source("./data_definitions.R")


################################################################################
#                        ENSURE THAT DIRECTORIES EXIST                         #
################################################################################

for (dataset in datasets) {
    dir <- dataset.dir(dataset)

    message("Esuring existence of the ", dataset, " directory")

    if (! dir.exists(dir)) {
        dir.create(dir)
    }
}


################################################################################
#                             FILTER THE DATASETS                              #
################################################################################

for (dataset in datasets) {
    f.file <- filtered.file(dataset)

    message("Filtering the ", dataset, " dataset")

    if (file.exists(f.file)) {
        message(f.file, " already exists, so ", dataset,
                " will not be filtered again.\n")
    } else {
        LabyrinthFilter(vcf      = original.file(dataset),
                        parents  = all.parents[[dataset]],
                        out.file = f.file,
                        hom.poly = FALSE)
    }
}


################################################################################
#                              MASK THE DATASETS                               #
################################################################################

for (dataset in datasets) {
    dataset <- "IBM-RIL"
    message("Masking the ", dataset, " dataset")
    set.seed(0)  # for repeatability
    LabyrinthMask(vcf       = filtered.file(dataset),
                  parents   = all.parents[[dataset]],
                  out.file  = masked.file(dataset),
                  depth     = min.mask.depths[dataset],
                  lik.ratio = 100,
                  rerr      = 0.05)
}


## ################################################################################
## ############################## MIMIC POPULATIONS ###############################
## ################################################################################

## LabyrinthMimic(
## set.seed(0)
