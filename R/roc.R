roc <-
    function(exp, obs, ...)
{
    exp0 <- exp[order(obs, decreasing=TRUE)]
    data.frame(TruePostive=cumsum(exp0) / sum(exp0),
               FalsePositive=cumsum(!exp0) / sum(!exp0))
}
