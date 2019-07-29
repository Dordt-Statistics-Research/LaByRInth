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


LabyrinthAnalyze <- function(orig, masked, imputed) {

    samples.o <- getSAMPLES(orig)
    samples.m <- getSAMPLES(masked)
    samples.i <- getSAMPLES(imputed)

    site.ids.o <- paste0(vcfR::getCHROM(orig), "\t", vcfR::getPOS(orig))
    site.ids.m <- paste0(vcfR::getCHROM(masked), "\t", vcfR::getPOS(masked))
    site.ids.i <- paste0(vcfR::getCHROM(imputed), "\t", vcfR::getPOS(imputed))

    if (! identical(samples.o, samples.m)) {
        stop("orig contains different progeny than masked")
    }
    if (! identical(samples.o, samples.i)) {
        stop("imputed contains different progeny than orig and masked")
    }


    if (! identical(site.ids.o, site.ids.m)) {
        stop("orig and masked contain different sites")
    }
    if (! identical(site.ids.o, site.ids.i)) {
        warning("orig and imputed contain different sites")
        if (! all(site.ids.i %in% site.ids.o)) {
            stop("progeny in imputed is not a subset of progeny in orig/masked")
        }

        common <- site.ids.o %in% site.ids.i # logical indicating which originals remain
        new.imputed <- orig                  # so all progeny are present
        is.format <- colnames(new.imputed@gt) == "FORMAT"
        new.imputed@gt[ , is.format] <- "GT"
        new.imputed@gt[ , ! is.format] <- "./."
        new.imputed@gt[common, ! is.format] <- imputed@gt[common, ! is.format]
        imputed <- new.imputed
    }


    gt.o <- getGT(orig)
    gt.m <- getGT(masked)
    gt.i <- getGT(imputed)
    gp   <- getGP(imputed)

    ## vector.indices
    masked.sites <-
        (gt.o != "./." & gt.o != ".|.") &
        (gt.m == "./." | gt.m == ".|.")

    gt.to.standard.numeric <- function(gt) {
        switch(gt,
               "0/0"=0,
               "0|0"=0,
               "0/1"=1,
               "0|1"=1,
               "1/0"=1,
               "1|0"=1,
               "1/1"=2,
               "1|1"=2,
                    NA)  #everything else including ./.
    }

    relevant.orig <- gt.o[masked.sites]
    relevant.imputed <- gt.i[masked.sites]
    relevant.posteriors <- gp[masked.sites]

    numeric.relevant.orig <- sapply(relevant.orig, gt.to.standard.numeric)
    numeric.relevant.imputed <- sapply(relevant.imputed, gt.to.standard.numeric)
    numeric.relevant.posteriors <-
        t(sapply(relevant.posteriors, function(gp.string) {
            phred.to.prob(as.numeric(str.split(gp.string, ",")))
        }))

    correct <- numeric.relevant.orig == numeric.relevant.imputed

    df <- data.frame(orig    = numeric.relevant.orig,
                     imputed = numeric.relevant.imputed,
                     correct = correct,
                     p.ref   = numeric.relevant.posteriors[ , 1],
                     p.het   = numeric.relevant.posteriors[ , 2],
                     p.alt   = numeric.relevant.posteriors[ , 3])
    rownames(df) <- NULL

    df  # implicit return
}


get.coverage <- function(df, total.sites) {
    sum(! is.na(df$correct)) / total.sites
}

get.accuracy <- function(df) {
    ## percent correct using boolean addition
    sum(df$correct, na.rm=TRUE) / sum(! is.na(df$correct))
}


get.quality <- function(df, total.sites) {
    sum(df$correct, na.rm=TRUE) / total.sites
}


get.r.sq <- function(df) {
    linear.model <- lm(df$imputed ~ df$orig)
    summary(linear.model)$r.squared
}


LabyrinthPlotPosteriors <- function(labyrinth.df, lb.impute.df, fsfhap.df,
                                    n=21, dataset.name, out.file) {

    ## function to remove low-quality imputation calls from analysis result
    filter.posterior <- function(df, thresh) {
        if (thresh > 1)
            stop("thresh should be in the interval [0, 1]")

        df[mapply(function(r,h,a){max(r,h,a)>=thresh}, df$p.ref, df$p.het, df$p.alt), ]
    }

    ## Because almost no sites will have a posterior probability of 1, make the
    ## last threshold a little less than 1
    thresholds <- seq(0, 1, length.out=n)
    thresholds[length(thresholds)] <- 1 - 1/n/100

    ## filter labyrinth.df
    filtereds <- lapply(thresholds,
                        function(thresh) {filter.posterior(labyrinth.df, thresh)})

    ## accuracy metric
    accuracy <- sapply(filtereds, get.accuracy)

    ## expected accuracy based on averaging all of the posterior probabilities
    expected.accuracy <- sapply(filtereds, function(df) {
        posteriors <-
            apply(df[ , c("p.ref", "p.het", "p.alt")], 1, max)
        mean(posteriors)
    })

    ## quality metric defined as total number of correctly imputed sited divided by the
    n.masked <- nrow(labyrinth.df)
    quality <- sapply(filtereds, function(df) {
        get.quality(df, n.masked)
    })

    r.sq <- sapply(filtereds, get.r.sq)

    my_geom_line <- function(metric, metric.name, linetype) {
        ggplot2::geom_line(
            aes(x = thresholds,
                y = metric,
                color = metric.name),
            linetype = linetype,
            size = 0.25)
    }

    my_geom_point <- function(metric, metric.name) {
        ggplot2::geom_point(
            aes(x = thresholds,
                y = metric,
                color = metric.name,
                shape = metric.name),
            size = 1.2)
    }

    my_trace_factory <- function(metric, metric.name, linetype) {
        list(my_geom_line(metric, metric.name, linetype),
             my_geom_point(metric, metric.name))
    }


    my_theme <- list(
        ggplot2::theme(
            text = element_text(size=10, color="black"),
            axis.text = element_text(color = "black"),
            plot.margin = unit(c(1,1,1,1), "cm"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x=element_text(vjust = -3),
            axis.title.y=element_text(vjust = 3),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.justification=c(0,0), # which corner is positioned
            legend.position=c(0.05,0.05),      # bottom right of graph
            legend.background =
                element_rect(size=0.25, linetype="solid", colour ="black"),
            panel.background =
                element_rect(size=0.5, linetype="solid", colour ="black")
            )
    )


    ## Compute LB-Impute and FSF-Hap accuracy and quality
    if (! is.null(lb.impute.df)) {
        lb.accuracy <- get.accuracy(lb.impute.df)
        lb.quality <- get.quality(lb.impute.df, nrow(fsfhap.df))
        lb.r.sq <- get.r.sq(lb.impute.df)
        lb.gg <- list(
            my_geom_line(rep(lb.accuracy, n), "Accuracy", "longdash"),
            my_geom_line(rep(lb.quality, n), "Quality", "longdash"),
            my_geom_line(rep(lb.r.sq, n), "R²", "longdash")
        )
    } else {
        lb.gg <- list()
    }

    if (! is.null(fsfhap.df)) {
        fsfhap.accuracy <- get.accuracy(fsfhap.df)
        fsfhap.quality <- get.quality(fsfhap.df, nrow(fsfhap.df))
        fsfhap.r.sq <- get.r.sq(fsfhap.df)
        fsfhap.gg <- list(
            my_geom_line(rep(fsfhap.accuracy, n), "Accuracy", "dotdash"),
            my_geom_line(rep(fsfhap.quality, n), "Quality", "dotdash"),
            my_geom_line(rep(fsfhap.r.sq, n), "R²", "dotdash"),

            geom_ribbon(aes(x = thresholds,
                            ymin = 0.72,
                            ymax = 0.77,
                            fill = "Quality"),
                        show.legend = FALSE,
                        alpha=0.4),
            geom_ribbon(aes(x = thresholds,
                            ymin = 0.94,
                            ymax = 0.98,
                            fill = "Accuracy"),
                        show.legend = FALSE,
                        alpha=0.4),
            geom_ribbon(aes(x = thresholds,
                            ymin = 0.7,
                            ymax = 0.8,
                            fill = "R²"),
                        show.legend = FALSE,
                        alpha=0.4)
        )
    }

    p <- ggplot2::ggplot(data.frame(thresholds, accuracy, expected.accuracy, quality)) +
         my_theme +

         ggplot2::ggtitle(paste0("Imputation Metrics for ", dataset.name, "\n")) +
         ggplot2::xlab("Posterior Probability Threshold") +
         ggplot2::ylab("Metric") +

         ggplot2::coord_cartesian(xlim = c(0, 1)) +
         ggplot2::coord_cartesian(ylim = c(0, 1)) +
         ggplot2::scale_x_continuous(breaks=seq(0, 1, 0.2)) +
         ggplot2::scale_y_continuous(breaks=seq(0, 1, 0.1)) +

         ggplot2::scale_shape_discrete(solid=T) +
         ggplot2::labs(linetype="Algorithm") +
         ggplot2::labs(color="Metric") +
         ggplot2::labs(shape="Metric") +

         ggplot2::geom_errorbar(aes(x=thresholds,
                       ymin = mapply(min, accuracy, expected.accuracy),
                       ymax = mapply(max, accuracy, expected.accuracy),
                       color = "Accuracy"),
                       width = 1/n/10) +
         my_trace_factory(accuracy, "Accuracy", "solid") +
         my_trace_factory(quality, "Quality", "solid") +
         my_trace_factory(r.sq, "R²", "solid") +
         lb.gg +
         fsfhap.gg

    ggplot2::ggsave(out.file, plot=p, width=4, height=4, units="in", dpi=500)
}
