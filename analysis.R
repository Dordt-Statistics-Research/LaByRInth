## Basic info
file <- "./analysis/filtered_lakin_fuller.vcf"
parents <- c("LAKIN", "FULLER")
out.file=""

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
prefs$parents           <- parents

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
prefs$out.file          <- out.file


#------------------------------------------------------------------------------#


#lakin.fuller.vcf <- VCF(file, prefs)
#saveRDS(lakin.fuller.vcf, "./analysis/filtered_lakin_fuller.rds")

n.mask <- 100  # how many masked sites per sample
n.samples <- length(lakin.fuller.vcf$GT[, 1, 1])  # length of primary dim
mask <- sample(x=seq(n.samples), size=n.mask)

lakin.fuller.vcf$GT[mask, , ] <- NA
lakin.fuller.vcf$AD[mask, , ] <- 0

lab <- VCF("./analysis/imp_LaByRInth_lf_masked_0001.vcf", prefs)
lb <- VCF("./analysis/imp_LB_lf_masked_0001.vcf", prefs)
#mask <- VCF("./analysis/lf_masked_0001.vcf", prefs)


c1 <- factor(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4), labels=c("correct", "partial", "skipped", "wrong"))
c2 <-    rev(c(1,1,1,1,1,1,2,2,2,3,3,4,1,1,2,2,2,1,2,3,4,3,3,1,1,2,2,2,3,3,4,4,4,4))
d1 <- as.data.frame(c1)
d2 <- as.data.frame(c2)
test <- cbind(d1, d2)
colnames(test) <- c("quality", "depth")



ggplot(diamonds, aes(price, fill=clarity)) + geom_bar(binwidth=1000, position="fill") + scale_fill_brewer(palette = "Blues")
ggplot(diamonds, aes(price, fill=clarity)) + geom_histogram(binwidth=1000, position="fill") + scale_fill_manual(values=c("#000000", "#3333FF","#000000", "#3333FF","#000000", "#3333FF","#000000", "#3333FF"))
ggplot(test, aes(depth, fill=quality)) + geom_histogram(binwidth=1) + scale_fill_manual(values=rev(c("#B70000", "#FFFFFF","#CDFF00","#29ECF3")))
ggplot(test, aes(depth, fill=quality)) + geom_bar(binwidth=1, position="fill") + scale_fill_manual(values=c("#B70000", "#FFFFFF","#CDFF00","#29ECF3"))
ggplot(mpg, aes(fl, fill = drv)) + geom_bar(position = "fill")


df <- AnalyzeImputationsRDS(lab, orig, mask)
df_LB <- AnalyzeImputationsRDS(lb, orig, mask)
## df <- as.data.frame(as.numeric(an))
## colnames(df) <- "depth"

library(scales)
my_tran <- trans_new("l1", function(x) {log10(x+1)}, function(y) {10**y - 1}, domain=c(0, Inf))
ggplot(d, aes(nd)) + geom_histogram(binwidth=1) + scale_y_continuous(trans = my_tran)


ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2) + scale_y_continuous(trans = my_tran, breaks=c(0, 10**c(0:7))) + scale_fill_manual(values=c("#29ECF3", "#CDFF00", "#FFFFFF", "#B70000"), drop=FALSE)

ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2) + scale_fill_manual(values=c("#29ECF3", "#CDFF00", "#FFFFFF", "#B70000"), drop=FALSE)

ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2, position="fill") + scale_fill_manual(values=c("#29ECF3", "#CDFF00", "#FFFFFF", "#B70000"), drop=FALSE)


colors <- c("#29ECF3", "#CDFF00", "#FFFFFF", "#B70000")


## no skipped sites in plot
png("./analysis/imp_LaByRInth_lf_masked_0001_AA.png", width=960, height=540)
ggplot(subset(df, as.numeric(quality)!=3), aes(depth, fill=quality)) + geom_histogram(binwidth=2, position="fill") + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

png("./analysis/imp_LB_lf_masked_0001_AA.png", width=960, height=540)
ggplot(subset(df_LB, as.numeric(quality)!=3), aes(depth, fill=quality)) + geom_histogram(binwidth=2, position="fill") + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

## include skipped sites in plot
png("./analysis/imp_LaByRInth_lf_masked_0001_AB.png", width=960, height=540)
ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2, position="fill") + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

png("./analysis/imp_LB_lf_masked_0001_AB.png", width=960, height=540)
ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2, position="fill") + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

## dont fill plot
png("./analysis/imp_LaByRInth_lf_masked_0001_AC.png", width=960, height=540)
ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2) + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

png("./analysis/imp_LB_lf_masked_0001_AC.png", width=960, height=540)
ggplot(df, aes(depth, fill=quality)) + geom_histogram(binwidth=2) + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

## dont fill plot without skipped sites
png("./analysis/imp_LaByRInth_lf_masked_0001_AD.png", width=960, height=540)
ggplot(subset(df_LB, as.numeric(quality)!=3), aes(depth, fill=quality)) + geom_histogram(binwidth=2) + scale_fill_manual(values=colors, drop=FALSE)
dev.off()

png("./analysis/imp_LB_lf_masked_0001_AD.png", width=960, height=540)
ggplot(subset(df_LB, as.numeric(quality)!=3), aes(depth, fill=quality)) + geom_histogram(binwidth=2) + scale_fill_manual(values=colors, drop=FALSE)
dev.off()


lab <- readRDS("./analysis/imp_LaByRInth_lf_masked_0001.rds")
lb <- readRDS("./analysis/imp_LB_lf_masked_0001.rds")
mask <- readRDS("./analysis/lf_masked_0001.rds")


dir <- "./analysis/filtered_lakin_fuller/mask_0001"
fi <- function(name) {paste0(dir, "/", name)}
f <- fi("imputed_infoLaB.vcf")
v <- VCF(f, prefs)
saveRDS(v, fi("imputed_infoLaB.rds"))


orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_infoLaB.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)

savePlot <- function(analyzed, id, dir, name, width=960, height=540) {
    require(ggplot2)
    colors <- c("#29ECF3", "#CDFF00", "#FFFFFF", "#B70000")
    file <- paste0(dir, "/", id, "_", name, ".png")
    print(file)

    f_aa <- function(){ggplot(analyzed, aes(depth, fill=quality)) + geom_histogram(binwidth=1, position="fill") + scale_fill_manual(values=colors, drop=FALSE)+ guides(fill = guide_legend(reverse = TRUE))}
    f_ab <- function(){ggplot(subset(analyzed, as.numeric(quality)!=3), aes(depth, fill=quality)) + geom_histogram(binwidth=1, position="fill") + scale_fill_manual(values=colors, drop=FALSE)}
    f_ac <- function(){ggplot(analyzed, aes(depth, fill=quality)) + geom_histogram(binwidth=1) + scale_fill_manual(values=colors, drop=FALSE)}
    f_ad <- function(){ggplot(subset(analyzed, as.numeric(quality)!=3), aes(depth, fill=quality)) + geom_histogram(binwidth=1) + scale_fill_manual(values=colors, drop=FALSE)}

    f <- switch(id,
                "AA" = f_aa,
                "AB" = f_ab,
                "AC" = f_ac,
                "AD" = f_ad)

    #png(file, width=960, height=540)
    f()
    ggsave(file)
    #dev.off()
}



#------------------------------------------------------------------------------#
dir <- "./analysis/filtered_lakin_fuller/mask_0001"
fi <- function(name) {paste0(dir, "/", name)}

orig <- readRDS(fi("../orig.rds"))
mask <- readRDS(fi("masked.rds"))

algorithms <- c("LB", "LaByRInth", "infoLaB")
analyses <- lapply(algorithms, function(alg) {
    imp <- readRDS(fi(paste0("imputed_", alg, ".rds")))
    AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
})

for (id in c("AA", "AB", "AC", "AD")) {
    for (i in 1:length(algorithms)) {
        alg <- algorithms[i]
        savePlot(analyses[[i]], id, dir, paste0("imputed_", alg))
    }
}

#savePlot(analyses[[1]], "AD", dir, "junk")

buildMask <- function(vcf, out.file) {
    ## n.samples <- dim(imputed_vcf$GT)[1]  # length of primary dim
    ## n.sites <- dim(imputed_vcf$GT)[2]  # length of primary dim
    ## n.entries <- prod(dim(imputed_vcf$GT))  # length of primary dim
    ## n.mask <- n.samples * 3

    ## print(paste("# of samples:       ", n.samples))
    ## print(paste("# of entries:       ", n.entries))
    ## print(paste("# of masked entries:", n.mask))
    ## entries <- sample(x=seq(n.entries), size=n.mask)

    ## data <- rep(0, n.entries)
    ## data[entries] <- 1
    ## mask <- matrix(data, nrow=n.samples, ncol=n.sites)

    ## mask  # implicit return
    n.mask <- 30  # how many masked sites per sample
    n.samples <- dim(vcf$GT)[1]  # length of primary dim
    mask <- sample(x=seq(n.samples), size=n.mask)

    vcf$GT[mask, , ] <- NA
    vcf$AD[mask, , ] <- 0

    saveRDS(vcf, out.file)
}

null_analysis <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=orig)


orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_LaByRInth.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)


dir <- "./analysis/filtered_lakin_fuller/mask_0002"
buildMask(readRDS(fi("../orig.rds")), fi("masked.rds"))
WriteVCF(readRDS(fi("masked.rds")), fi("masked.vcf"))

infolab <- VCF(fi("imputed_infoLaB.vcf"), prefs)
saveRDS(infolab, fi("imputed_infoLaB.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_infoLaB.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)


#------------------------------------------------------------------------------#
# sandesh mask
#------------------------------------------------------------------------------#
dir <- "./analysis/filtered_lakin_fuller/mask_0003"
sandesh <- VCF(fi("masked.vcf"), prefs)
saveRDS(sandesh, fi("masked.rds"))

temp <- LabyrinthImpute(fi("masked.vcf"), c("LAKIN","FULLER"), fi("imputed_infoLaB.vcf"))

lab <- VCF(fi("imputed_infoLaB.vcf"), prefs)
saveRDS(lab, fi("imputed_infoLaB.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_infoLaB.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
ggplot(analyze, aes(depth, fill=quality)) + geom_histogram(binwidth=1) + scale_fill_manual(values=colors, drop=FALSE)
ggplot(analyze, aes(depth, fill=quality)) + geom_histogram(binwidth=1, position="fill") + scale_fill_manual(values=colors, drop=FALSE)
for (id in c("AA", "AB", "AC", "AD")) {
    savePlot(analyze, id, dir, "imputed_infoLaB")
}


temp <- LabyrinthImpute(fi("masked.vcf"), c("LAKIN","FULLER"), fi("imputed_LaByRInth.vcf"))

lab <- VCF(fi("imputed_LaByRInth.vcf"), prefs)
saveRDS(lab, fi("imputed_LaByRInth.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_LaByRInth.rds"))
analyze2 <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
for (id in c("AA", "AB", "AC", "AD")) {
    savePlot(analyze2, id, dir, "imputed_LaByRInth")
}




lb <- VCF(fi("imputed_LB.vcf"), prefs)
saveRDS(lb, fi("imputed_LB.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_LB.rds"))
analyze3 <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
for (id in c("AA", "AB", "AC", "AD")) {
    savePlot(analyze3, id, dir, "imputed_LB")
}



temp2 <- LabyrinthImpute(fi("masked.vcf"), c("LAKIN","FULLER"), prefs) # run on commit 790fe16
lab2 <- VCF(fi("imputed_LaByRInth_2.vcf"), prefs)
saveRDS(lab, fi("imputed_LaByRInth_2.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_LaByRInth_2.rds"))
analyze4 <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
for (id in c("AA", "AB", "AC", "AD")) {
    savePlot(analyze4, id, dir, "imputed_LaByRInth_2")
}





dir <- "./analysis/filtered_lakin_fuller/mask_0004"
temp <- LabyrinthImpute(fi("masked.vcf"), c("LAKIN","FULLER"), fi("imputed_LaByRInth.vcf")) # commit e22a701 (v0.0.1)

lab <- VCF(fi("imputed_LaByRInth.vcf"), prefs) # commit dev
saveRDS(lab, fi("imputed_LaByRInth.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_LaByRInth.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
for (id in c("AA", "AB", "AC", "AD")) {
    savePlot(analyze, id, dir, "imputed_LaByRInth")
}

temp <- LabyrinthImpute(fi("masked.vcf"), c("LAKIN","FULLER"), fi("imputed_infoLaB.vcf")) # info branch

lab <- VCF(fi("imputed_infoLaB.vcf"), prefs) # commit dev
saveRDS(lab, fi("imputed_infoLaB.rds"))

saveRDS(VCF(fi("imputed_infoLaB.vcf"), prefs), fi("imputed_infoLaB.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_infoLaB.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
for (id in c("AA", "AB", "AC", "AD")) {
    savePlot(analyze, id, dir, "imputed_infoLaB")
}



temp <- LabyrinthImpute(fi("masked.vcf"), c("LAKIN","FULLER"), fi("imputed_infoLaB.vcf"))

lab <- VCF(fi("imputed_infoLaB.vcf"), prefs)
saveRDS(lab, fi("imputed_infoLaB.rds"))

orig = readRDS(fi("../orig.rds"))
mask = readRDS(fi("masked.rds"))
imp  = readRDS(fi("imputed_infoLaB.rds"))
analyze <- AnalyzeImputationsRDS(imp=imp, orig=orig, mask=mask)
