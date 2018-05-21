plot.generational.probs <- function(recomb.probs) {
    ps <- recomb.probs

    # F1 probs
    gamete.probs <- f1.gamete.probs(ps)
    n <- length(gamete.probs)
    plot(0:(n-1), gamete.probs, col=1, cex=0.5)

    for (i in 2:5) {
        gamete.probs <- next.gamete.probs(gamete.probs, ps)
        points(0:(n-1), gamete.probs, col=i, cex=0.5)
    }
}


unit.test.2.site <- function(ps) {
    if (length(ps) != 5) stop("length of ps is not 5")
    f1.probs <- initial.6.site(ps)
    plot(0:(n-1), f1.probs, col=1, cex=0.5)
    n <- length(f1.probs)
    probs <- f1.probs
    for (i in 2:5) {
        probs <- recursive.6.site(probs, ps)
        points(0:(n-1), probs, col=i, cex=0.5)
    }
}
