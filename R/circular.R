# Circular triads (transitivity violations)
# 
# May/10/2010: Bug fix, p-value for zeta is 1 - p-value of chi-square test
#              (reported by Wolfgang Ellermeier)
#
# Sep/21/2018: - exact test for small n
#              - one- and two-sided tests
#              - simulated p-value
#              - list circular triads via strans(), see examples


## User interface
circular <- function(mat, alternative = c("two.sided", "less", "greater"),
                     exact = NULL, correct = TRUE, simulate.p.value = FALSE,
                     nsim = 2000){

  alternative <- match.arg(alternative)
  mat <- as.matrix(mat)
  stopifnot(nrow(mat) > 1,
            nrow(mat) == ncol(mat))
  diag(mat) <- 0                  # remove diagonal
  n <- ncol(mat)

  if(is.null(exact))
    exact <- (n < 11)
  chi2 <- df <- correct.msg <- NULL

  T <- n * (n - 1) * (2*n - 1)/12 - .5 * sum(colSums(mat)^2)
  T.max <- if(n %% 2) n*(n^2 - 1)/24 else n*(n^2 - 4)/24
  T.exp <- choose(n, 3)/4
  zeta <- 1 - T/T.max

  pval <- if(simulate.p.value) {  # MC distribution

    dc <- dcircular(n, T.max = T.max, simulate = TRUE, nsim = nsim)
    switch(alternative,
      "two.sided" = {
        p <- if(T <= T.exp) {
          p1 <- sum(dc[as.character(0:T)])
          p1 + sum(rev(dc)[cumsum(rev(dc)) <= p1])
        } else {
          p1 <- sum(dc[as.character(T:T.max)])
          p1 + sum(dc[cumsum(dc) <= p1])
        }
        min(p, 1)
      },
      "less"      = sum(dc[as.character(0:T)]),
      "greater"   = sum(dc[as.character(T:T.max)])
    )

  } else if(exact && n < 11) {    # exact distribution

    dc <- dcircular(n)
    switch(alternative,
      "two.sided" = {
        p <- if(T <= T.exp) {
          p1 <- sum(dc[as.character(0:T)])
          p1 + sum(rev(dc)[cumsum(rev(dc)) <= p1])
        } else {
          p1 <- sum(dc[as.character(T:T.max)])
          p1 + sum(dc[cumsum(dc) <= p1])
        }
        min(p, 1)
      },
      "less"      = sum(dc[as.character(0:T)]),
      "greater"   = sum(dc[as.character(T:T.max)])
    )

  } else {                        # chi2 approximation

    correct.msg <- correct
    correction <- if(correct) {
      switch(alternative,
        "two.sided" = if(T <= T.exp) -0.5 else 0.5,
        "less"      = -0.5,
        "greater"   =  0.5)
    } else 0

    df   <- if(n < 5) NA else n * (n - 1) * (n - 2) / (n - 4)^2
    chi2 <- if(n < 5) NA else 8/(n - 4) * (T.exp - T + correction) + df

    switch(alternative,
      "two.sided" = 2 * min(pchisq(chi2, df, lower.tail = FALSE),
                            pchisq(chi2, df)),
      "less"      = pchisq(chi2, df, lower.tail = FALSE),
      "greater"   = pchisq(chi2, df)
    )
  }

  if(!simulate.p.value && exact && n > 10)
    warning("cannot compute exact p-value for n > 10")

  if(!simulate.p.value && !exact && n < 11)
    warning("Chi-square approximation may be incorrect")

  z <- list(T=T, T.max=T.max, T.exp=T.exp, zeta=zeta, alternative=alternative,
            chi2 = if(!is.null(chi2)) chi2 else NA,
            df = if(!is.null(df)) df else NA,
            correct = if(!is.null(correct.msg)) correct.msg else FALSE,
            p.value=pval, simulate.p.value=simulate.p.value, nsim=nsim)
  class(z) <- "circular"
  z
}


## Helper function for calculating simulated and exact distributions
dcircular <- function(n, T.max = NULL,
                      simulate = FALSE, nsim = 2000) {
  dc <- if (simulate) {
    M <- matrix(0, n, n)
    dcsim <- setNames(integer(T.max + 1), 0:T.max)
  
    out <- table(replicate(nsim, {
      M[upper.tri(M)] <- sample(0:1, n*(n - 1)/2, replace = TRUE)
      M[lower.tri(M)] <- (1 - t(M))[lower.tri(M)]
      n * (n - 1) * (2*n - 1)/12 - .5 * sum(colSums(M)^2)
    }))
    dcsim[names(out)] <- out
    dcsim
  } else {
    if (n > 10) stop("cannot compute exact density for n > 10")
    list(
      "2" = c(`0` = 2),
      "3" = c(`0` = 6, `1` = 2),
      "4" = c(`0` = 24, `1` = 16, `2` = 24),
      "5" = c(`0` = 120, `1` = 120, `2` = 240, `3` = 240, `4` = 280, `5` = 24),
      "6" = c(`0` = 720, `1` = 960, `2` = 2240, `3` = 2880, `4` = 6240,
              `5` = 3648, `6` = 8640, `7` = 4800, `8` = 2640),
      "7" = c(`0` = 5040, `1` = 8400, `2` = 21840, `3` = 33600, `4` = 75600,
              `5` = 90384, `6` = 179760, `7` = 188160, `8` = 277200,
              `9` = 280560, `10` = 384048, `11` = 244160, `12` = 233520,
              `13` = 72240, `14` = 2640),
      "8" = c(`0` = 40320, `1` = 80640, `2` = 228480, `3` = 403200,
              `4` = 954240, `5` = 1304576, `6` = 3042816, `7` = 3870720,
              `8` = 6926080, `9` = 8332800, `10` = 15821568, `11` = 14755328,
              `12` = 24487680, `13` = 24514560, `14` = 34762240,
              `15` = 29288448, `16` = 37188480, `17` = 24487680,
              `18` = 24312960, `19` = 10402560, `20` = 3230080),
      "9" = c(`0` = 362880, `1` = 846720, `2` = 2580480, `3` = 5093760,
              `4` = 12579840, `5` = 19958400, `6` = 44698752, `7` = 70785792,
              `8` = 130032000, `9` = 190834560, `10` = 361525248,
              `11` = 443931264, `12` = 779950080, `13` = 1043763840,
              `14` = 1529101440, `15` = 1916619264, `16` = 2912257152,
              `17` = 3078407808, `18` = 4506485760, `19` = 4946417280,
              `20` = 6068256768, `21` = 6160876416, `22` = 7730384256,
              `23` = 6292581120, `24` = 6900969600, `25` = 5479802496,
              `26` = 4327787520, `27` = 2399241600, `28` = 1197020160,
              `29` = 163094400, `30` = 3230080),
     "10" = c(`0` = 3628800, `1` = 9676800, `2` = 31449600, `3` = 68275200,
              `4` = 175392000, `5` = 311592960, `6` = 711728640,
              `7` = 1193794560, `8` = 2393475840, `9` = 3784596480,
              `10` = 7444104192, `11` = 10526745600, `12` = 19533696000,
              `13` = 27610168320, `14` = 47107169280, `15` = 64016040960,
              `16` = 107446832640, `17` = 134470425600, `18` = 218941470720,
              `19` = 272302894080, `20` = 417512148480, `21` = 494080834560,
              `22` = 743278970880, `23` = 829743344640, `24` = 1202317401600,
              `25` = 1334577484800, `26` = 1773862272000, `27` = 1878824586240,
              `28` = 2496636103680, `29` = 2406981104640, `30` = 3032021672960,
              `31` = 2841072675840, `32` = 3166378709760, `33` = 2743311191040,
              `34` = 2877794035200, `35` = 2109852702720, `36` = 1840136336640,
              `37` = 1109253196800, `38` = 689719564800, `39` = 230683084800,
              `40` = 48251508480)
    )[[as.character(n)]]
  }
  prop.table(dc)
}


print.circular <- function(x, digits = max(3L, getOption("digits") - 3L),
                           ...){
  cat("\nCircular triads (intransitive cycles)\n\n")
  cat("T = ", x$T, ", max(T) = ", x$T.max, ", E(T) = ", x$T.exp,
      ", zeta = ", format(x$zeta, digits=digits, ...),
      ", p-value = ", format(x$p, digits=digits, ...),
      "\n", sep="")

  cat("alternative hypothesis: T is ")
  cat(switch(x$alternative,
	     two.sided = "smaller or greater",
	     less = "smaller",
	     greater = "greater")
  )
  cat(" than expected by chance\n", sep = "")

  if(x$correct && !x$simulate.p.value)
    cat("continuity correction has been applied\n", sep="")

  if(x$simulate.p.value)
    cat("simulated p-value based on ", x$nsim, " replicates\n", sep="")

  cat("\n")
  invisible(x)
}

