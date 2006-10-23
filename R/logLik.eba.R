## log-likelihood for eba objects
logLik.eba <- function(object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    p <- length(object$estimate) - 1
    val <- object$logL.eba
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
}
