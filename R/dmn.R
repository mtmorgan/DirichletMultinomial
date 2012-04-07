setClass("DMN",
    representation=representation(goodnessOfFit="numeric",
      group="matrix", mixture="list", fit="list"))

.DMN <-
    function(goodnessOfFit, group, mixture, fit, ...)
{
    new("DMN", goodnessOfFit=goodnessOfFit, group=group,
        mixture=mixture, fit=fit, ...)
}

dmn <-
    function(count, k, verbose=FALSE,
             seed=runif(1, 0, .Machine$integer.max))
{
    if (verbose)
        message(sprintf("dmn, k=%d", k))
    mode(count) <- "integer"
    ans <- .Call(.dirichlet_fit, count, as.integer(k),
                 as.logical(verbose), as.integer(seed))
    o <- order(ans$Mixture$Weight, decreasing=TRUE)
    ans <- within(ans, {
        Group <- Group[,o, drop=FALSE]
        Mixture$Weight <- Mixture$Weight[o]
        Fit <- lapply(Fit, function(elt, o) elt[, o, drop=FALSE], o)
    })
    with(ans,
         .DMN(goodnessOfFit=GoodnessOfFit, group=Group,
              mixture=Mixture, fit=Fit))
}

## k-means

mixture <-
    function(object, ..., assign=FALSE)
{
    if (assign) {
        apply(mixture(object), 1, which.max)
    } else {
        object@group
    }
}

## Dirichlet

goodnessOfFit <- function(object, ...) object@goodnessOfFit

laplace <- function(object, ...) goodnessOfFit(object)[["Laplace"]]

.AIC.DMN <- function(object, ...) goodnessOfFit(object)[["AIC"]]

setMethod(AIC, "DMN", .AIC.DMN)

.BIC.DMN <- function(object, ...) goodnessOfFit(object)[["BIC"]]

setMethod(BIC, "DMN", .BIC.DMN)

mixturewt <-
    function(object, ...)
{
    data.frame(pi=object@mixture$Weight, theta=colSums(fitted(object)))
}

.fitted.DMN <- function(object, ..., scale=FALSE)
{
    fit <- object@fit$Estimate
    if (scale)
        fit <- scale(fit, FALSE, mixturewt(object)$theta)
    fit
}

setMethod(fitted, "DMN", .fitted.DMN)

## predict

.neg_log_evidence_i <- 
    function(x, alpha)
{
    .B <- function(x)
        sum(lgamma(x)) - lgamma(sum(x))
    -(.B(x + alpha) - .B(alpha))
}

.predict.DMN <- 
    function(object, newdata, ..., logevidence=FALSE)
{
    if (is.vector(newdata))
        newdata <- matrix(newdata, nrow=1)
    lambda <- fitted(object)

    K <- ncol(lambda)
    alpha <- sapply(seq_len(K), function(k, lamda, x) {
        apply(x, 1, .neg_log_evidence_i, lambda[,k])
    }, lambda, newdata)
    if (is.vector(alpha))
        alpha <- matrix(alpha, nrow=1,
                        dimnames=list(rownames(newdata), NULL))
    if (!logevidence) {
        wt <- mixturewt(object)$pi
        offset <- apply(alpha, 1, min)
        z <- sweep(exp(-(alpha - offset)), 2, wt, "*")
        z / rowSums(z)
    } else {
        alpha
    }
}

setMethod(predict, "DMN", .predict.DMN)

## print / plot

setMethod(show, "DMN",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("k:", ncol(mixture(object)), "\n")
    cat("samples x taxa:", nrow(mixture(object)), "x",
        nrow(fitted(object)), "\n")
    cat("Laplace:", laplace(object), "BIC:", BIC(object),
        "AIC:", AIC(object), "\n")
})

heatmapdmn <- 
    function(count, fit1, fitN, ntaxa=30, ..., transform=sqrt,
             lblwidth=.2 * nrow(count), col=.gradient)
{
    p1 <- fitted(fit1, scale=TRUE)
    pN <- fitted(fitN, scale=TRUE)
    diff <- rowSums(abs(pN - as.vector(p1)))
    taxa <- rev(head(order(diff, decreasing=TRUE), ntaxa))
    pN <- pN[taxa,]

    cl <- mixture(fitN, assign=TRUE)
    ncl <- length(unique(cl))
    nms <- names(cl)
    grp <- factor(cl, levels=as.character(seq(1, ncl)))
    idx <- split(nms, grp)
    ## 2 * ncl + 1 (for labels) panels
    mwd <- .15 * length(cl) / ncl    # 'm's take up 15% of total width
    wd <- c(unlist(Map(c, lapply(idx, length), mwd), use.names=FALSE), lblwidth)
    layout(matrix(seq(1, 2 * ncl + 1), nrow=1), widths=wd)
    op <- par(no.readonly=TRUE)
    on.exit(par(op), add=TRUE)
    par(mar=c(1, 0, 1, 0))
    for (i in seq_along(idx)) {
        image(transform(count[idx[[i]], taxa, drop=FALSE]),
              col=col, xaxt="n", yaxt="n")
        image(t(transform(pN[, i, drop=FALSE])),
              col=col, xaxt="n", yaxt="n")
    }
    xat <- (seq_len(nrow(pN)) - 1) / (nrow(pN) - 1)
    axis(4, xat, labels=rownames(pN), las=1)
}
