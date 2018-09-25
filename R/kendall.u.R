kendall.u <- function(M, correct = TRUE){
  # Kendall's (1940, 1962) coefficient of agreement
  # Assumptions: equal number of obs. per pair
  #              one obs. per judge and pair
  # Maximum agreement: u = 1, the smaller u the less agreement
  # Chi2 test, H_0: agreement is by chance
  # m observers (judges), n stimuli
  # last mod: Sep/25/2018

  Sigma <- sum( choose(M[upper.tri(M) | lower.tri(M)], 2) )
  u <- 2*Sigma / ( choose(m <- M[1, 2] + M[2, 1], 2) *
                   choose(n <- nrow(M), 2) ) - 1
  min.u <- if(m %% 2) -1/m else -1/(m - 1)
  chi <- 4/(m - 2) * (Sigma - correct -
                      choose(n, 2)/2 * choose(m, 2) * (m - 3)/(m - 2))
  df <- choose(n, 2) * (m*(m - 1))/(m - 2)^2
  out <- list(u=u, min.u=min.u, chi=chi, df=df,
              p.value=pchisq(chi, df, lower.tail = FALSE),
              correct=correct)
  class(out) <- "kendall.u"
  out
}


print.kendall.u <- function(x, digits = max(3L, getOption("digits") - 3L),
                            ...){
  cat("\nKendall's u coefficient of agreement\n\n")
  cat("u = ", format(x$u, digits=digits, ...),
      ", minimum u = ", format(x$min.u, digits=digits, ...),
      "\nchi2 = ", format(x$chi, digits=digits + 1L, ...),
      ", df = ", format(x$df, digits=digits, ...),
      ", p-value = ", format(x$p, digits=digits, ...),
      "\n", sep="")
  cat("alternative hypothesis: between-judges agreement is not by chance\n")
  if(x$correct)
    cat("continuity correction has been applied\n", sep="")
  cat("\n")
  invisible(x)
}

