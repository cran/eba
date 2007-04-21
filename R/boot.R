boot <- function(D, R=100, A=1:I, s=rep(1/J,J), constrained=TRUE){
  # performs bootstrapping by resampling the ind. PCMs
  # input: D, a 3d array of ind. PCMs
  # output: bootstrap means, standard errors, and cis
  # author: Florian Wickelmaier (wickelmaier@web.de)
  # last mod: 20/JUL/2004
  #           16/MAR/2007

  # constrained minimization might yield suboptimal estimates
  #  (=> always compare to the eba-estimated values!)
  # try() catches nlm errors, so generally nrow(p) =< R

  I <- ncol(D)  # number of alternatives/stimuli
  J <- max(unlist(A))  # number of eba parameters

  idx1 <- matrix(0, I*(I-1)/2, J)  # index matrices
  idx0 <- idx1
  rdx <- 1
  for(i in 1:(I-1)){
    for(j in (i+1):I){
      idx1[rdx, setdiff(A[[i]], A[[j]])] = 1
      idx0[rdx, setdiff(A[[j]], A[[i]])] = 1
      rdx = rdx + 1
    }
  }

  n <- dim(D)[3]
  p <- numeric()

  if(constrained){  # use eba.boot.constrained
    for(i in 1:R){
      cat(paste("Bootstrap sample ", i, "\n", sep=""), sep="")
      par <- try( eba.boot.constrained(apply(D[,,sample(n, r=TRUE)], 1:2, sum),
                    idx1, idx0, s) )
      if(class(par) != "try-error") p <- rbind(p, par)
    }
  }else{  # use eba.boot
    for(i in 1:R){
      cat(paste("Bootstrap sample ", i, "\n", sep=""), sep="")
      par <- try( eba.boot(apply(D[,,sample(n, r=TRUE)], 1:2, sum),
                    idx1, idx0, s) )
      if(class(par) != "try-error") p <- rbind(p, par)
    }
  }

  rownames(p) <- NULL
  if(nrow(p) < R) cat("Missing bootstrap sample(s) due to nlm errors.\n")
  p.mean <- colMeans(p)
  p.se <- sqrt(apply(p, 2, var))
  p.ci <- apply(p, 2, quantile, c(.025, .975))
  out <- list(p=p, stat=cbind(mean=p.mean, se=p.se, t(p.ci)))
  out
}


eba.boot <- function(M, idx1, idx0, s){
  y1 <- t(M)[lower.tri(t(M))]
  y0 <- M[lower.tri(M)]
  n <- y1 + y0

  out <- nlm(L, s, y1=y1, m=n, i1=idx1, i0=idx0)
  out$est
}


eba.boot.constrained <- function(M, idx1, idx0, s){
  y1 <- t(M)[lower.tri(t(M))]
  y0 <- M[lower.tri(M)]
  n <- y1 + y0

  out <- nlm(L.constrained, s, y1=y1, m=n, i1=idx1, i0=idx0)
  out$est
}
