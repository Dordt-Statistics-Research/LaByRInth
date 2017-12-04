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

