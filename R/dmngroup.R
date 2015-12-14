setClass("DMNGroup", contains="SimpleList",
    prototype=prototype(elementType="DMN"))

.DMNGroup <-
    function(...)
{
    new("DMNGroup", listData=list(...))
}

## dmngroup

dmngroup <-
    function(count, group, k, ..., simplify=TRUE, .lapply=parallel::mclapply)
{
    if (length(group) != nrow(count))
        stop("'length(group)' does not equal 'nrow(count)'")
    if (!is.factor(group))
        group <- factor(group)

    counts <- sapply(levels(group), csubset, count, group)
    tasks <- expand.grid(group=names(counts), k=k)
    tid <- seq_len(nrow(tasks))
    ans0 <- .lapply(tid, function(i, tasks, counts, ...) {
        count <- counts[[tasks[i,"group"]]]
        k <- tasks[i,"k"]
        dmn(count, k, ...)
    }, tasks, counts, ...)
    ans <- if (simplify) {
        ans1 <- split(ans0, tasks[,"group"])
        opt <- lapply(ans1, function(ans) {
            which.min(sapply(ans, laplace))
        })
        Map("[[", ans1, opt)
    } else ans0
    do.call(.DMNGroup, ans)
}

## predict

.predict.DMNGroup <-
    function(object, newdata, ..., assign=FALSE)
{
    if (2 < length(object))
        stop("only 2 groups can be used for classification")
    res <- lapply(object, predict, newdata, ..., logevidence=TRUE)
    offset <- apply(do.call(cbind, res), 1, min)
    prClass <- local({
        nClass <- sapply(object, function(x) nrow(mixture(x)))
        nClass / sum(nClass)
    })
    pr <- simplify2array(Map(function(x, alpha, prClass, offset) {
        prMix <- sweep(exp(-(alpha - offset)), 2, mixturewt(x)$pi, "*")
        rowSums(prMix) * prClass
    }, object, res, prClass, MoreArgs=list(offset=offset)))
    if (!is.matrix(pr)) {
        dmnms <- list(rownames(newdata), names(prClass))
        pr <- matrix(pr, nrow=1, dimnames=dmnms)
    }
    
    if (assign)
        names(object)[ifelse((pr[,1] / rowSums(pr)) > .5, 1, 2)]
    else
        pr / rowSums(pr)
}

setMethod(predict, "DMNGroup", .predict.DMNGroup)

## cross-validation

.cv_dmngroup <-
    function(dropidx, count, k, z, ..., verbose=FALSE)
    ## e.g., k = c(Lean=1, Obese=3) --> 1 group for lean, 3 for obese
{
    tryCatch({
        trainz <- z[-dropidx]
        u <- unique(trainz)
        train <- count[-dropidx,,drop=FALSE]
        if (!is.factor(trainz))
            trainz <- factor(trainz, levels=names(k))
        if (any(is.na(trainz)))
            stop("values of 'z' not all in 'names(k)'")
        if (!all(names(k) %in% as.character(trainz)))
            stop("not all names(k) in z subset")

        trains <- sapply(levels(trainz), csubset, train, trainz)
        fits <- Map(dmn, trains, k[levels(trainz)], ...,
                    verbose=verbose)
        fits <- do.call(.DMNGroup, fits)
        predict(fits, count[dropidx,,drop=FALSE], assign=FALSE)
    }, error=function(err) {
        message(".cv_dmngroup error: ", conditionMessage(err))
        matrix(NA_integer_, nrow=length(dropidx), ncol=length(k),
               dimnames=list(rownames(count)[dropidx], names(k)))
    })
}

cvdmngroup <-
    function(ncv, count, k, z, ..., verbose=FALSE, .lapply=parallel::mclapply)
{
    n <- seq_len(nrow(count))
    grp <- split(sample(length(n)), cut(n, ncv))
    names(grp) <- seq_along(grp)
    cvresult <- .lapply(names(grp), function(idx, grp, ..., verbose) {
        if (verbose)
            cat("cross-validation group", names(grp[idx]), "\n")
        .cv_dmngroup(grp[[idx]], ..., verbose=verbose)
    }, grp, count, k, z, ..., verbose=verbose)
    gid <- rep(seq_along(cvresult), sapply(cvresult, nrow))
    cbind(data.frame(group=gid, row.names=NULL),
          do.call(rbind, cvresult))
}

## summary / print / plot

setMethod(summary, "DMNGroup",
    function(object, ...)
{
    k <- data.frame(k=sapply(object, function(elt) ncol(mixture(elt))))
    sxt <- t(sapply(object, function(elt) {
        c(samples=nrow(mixture(elt)), taxa=nrow(fitted(elt)))
    }))
    goodness <- t(sapply(object, goodnessOfFit))
    cbind(k=k, sxt, goodness)
})

setMethod(show, "DMNGroup",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("summary:\n")
    print(summary(object))
})
