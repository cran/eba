library(nlme)  # needed for fdHess()
# 2004/MAY/04: removing cov2cor (it's there anyway)
#              changing 'print.coefmat' to 'printCoefmat'
#              'print.matrix' doesn't work anymore (for strans)
# 2004/JUL/20: trying to make it fit for CRAN:
#   removing beta2eba.R (or do you want to explain what it does??)
#   making summary.eba generic: add '...', rename x to object
#   changing all 'T' to 'TRUE' (goddamn why??)


OptiPt = function(M,A=1:I,s=rep(1/J,J)){
  # parameter estimation for BTL/Pretree/EBA models
  # M: paired-comparison matrix
  # A: model specification list(c(1,6),c(2,6),c(3,7),...)
  # s: starting vector (optional)
  # author: Florian Wickelmaier (wickelmaier@web.de)
  # last mod: 18/NOV/2003
  # Reference: Wickelmaier, F. & Schmid, C. (2003). A Matlab function
  #   to estimate choice-model parameters from paired-comparison data.
  #   Behavior Research Methods, Instruments, and Computers. In press.

  I = ncol(M)  # number of alternatives/stimuli
  J = max(unlist(A))  # number of eba parameters

  idx1 = matrix(0,I*(I-1)/2,J)  # index matrices
  idx0 = idx1
  rdx=1
  for(i in 1:(I-1)){                         
    for(j in (i+1):I){
      idx1[rdx,setdiff(A[[i]],A[[j]])] = 1
      idx0[rdx,setdiff(A[[j]],A[[i]])] = 1
      rdx = rdx+1
    }
  }

  y1 = t(M)[lower.tri(t(M))]  # response vectors
  y0 = M[lower.tri(M)]
  n = y1+y0
  logL.sat = sum(dbinom(y1,n,y1/n,log=TRUE))  # likelihood of the sat. model

  out = nlm(L,s,y1=y1,m=n,i1=idx1,i0=idx0)  # Newton-type algorithm
  p = out$est  # optimized parameters
  hes = fdHess(p,L,y1,n,idx1,idx0)$H  # numerical Hessian
  cova = solve(rbind(cbind(hes,1),c(rep(1,J),0)))[1:J,1:J]
  se = sqrt(diag(cova))  # standard error
  ci = qnorm(.975)*se  # 95% confidence interval
  logL.eba = -out$min  # likelihood of the specified model

  fitted = matrix(0,I,I)  # fitted PCM
  fitted[lower.tri(fitted)] = n/(1+idx0%*%p/idx1%*%p)
  fitted = t(fitted)
  fitted[lower.tri(fitted)] = n/(1+idx1%*%p/idx0%*%p)

  chi = 2*(logL.sat-logL.eba)  # goodness-of-fit statistic
  df = I*(I-1)/2 - (J-1)
  pval = 1-pchisq(chi,df)
  gof = c(chi,df,pval)
  names(gof) = c("-2logL","df","pval")
  chi.alt = sum((M-fitted)^2/fitted,na.rm=TRUE)

  u = numeric()  # scale values
  for(i in 1:I) u = c(u,sum(p[A[[i]]]))
  names(u) = colnames(M)

  z = list(estimate=p,se=se,ci95=ci,fitted=fitted,logL.eba=logL.eba,
           logL.sat=logL.sat,goodness.of.fit=gof,u.scale=u,
           hessian=-hes,cov.p=cova,chi.alt=chi.alt,A=A,y1=y1,y0=y0,n=n)
  class(z) = "eba"
  z
}


summary.eba = function(object,...){
  object -> x
  I = length(x$A)
  J = length(x$estimate)
  y = c(x$y1,x$y0)
  chi2 = x$goodness.of.fit[1]
  df = x$goodness.of.fit[2]
  pval = x$goodness.of.fit[3]
  coef = x$estimate
  s.err = x$se
  tvalue = coef/s.err
  pvalue = 2 * pnorm(-abs(tvalue))
  dn = c("Estimate", "Std. Error")
  coef.table = cbind(coef, s.err, tvalue, pvalue)
  dimnames(coef.table) = list(names(coef),c(dn,"z value","Pr(>|z|)"))

  tests = rbind(
    c(1,I*(I-1),sum(dpois(y,mean(y),log=TRUE)),sum(dpois(y,y,log=TRUE)),
      2*(sum(dpois(y,y,log=TRUE))-sum(dpois(y,mean(y),log=TRUE))),
      1-pchisq(2*(sum(dpois(y,y,log=TRUE))-sum(dpois(y,mean(y),log=TRUE))),
        I*(I-1)-1)),
    c(J-1,I*(I-1)/2,x$logL.eba,x$logL.sat,chi2,pval),
    c(0,J-1,sum(dbinom(x$y1,x$n,1/2,log=TRUE)),x$logL.eba,
      2*(x$logL.eba-sum(dbinom(x$y1,x$n,1/2,log=TRUE))),
      1-pchisq(2*(x$logL.eba-sum(dbinom(x$y1,x$n,1/2,log=TRUE))),J-1)),
    c(1,I*(I-1)/2,sum(dpois(x$n,mean(x$n),log=TRUE)),sum(dpois(x$n,x$n,log=TRUE)),
      2*(sum(dpois(x$n,x$n,log=TRUE))-sum(dpois(x$n,mean(x$n),log=TRUE))),
      1-pchisq(2*(sum(dpois(x$n,x$n,log=TRUE))-sum(dpois(x$n,mean(x$n),log=TRUE))),
        I*(I-1)-1))
  )
  rownames(tests) = c('Overall','EBA','Effect','Imbalance')
  colnames(tests) = c('Df1','Df2','logLik1','logLik2','Deviance','Pr(>|Chi|)')

  aic = -2*x$logL.eba + 2*(length(coef)-1)
  ans = list(coefficients=coef.table,chi2=chi2,df=df,pval=pval,aic=aic,
             logL.eba=x$logL.eba,logL.sat=x$logL.sat,tests=tests,
             chi.alt=x$chi.alt)
  class(ans) <- "summary.eba"
  return(ans)
}


print.eba = function(x,digits=max(3,getOption("digits")-3),na.print="",...){
  cat("\nEBA models\n\n")
  cat("Parameter estimates:\n")
  print.default(format(x$estimate, digits = digits), print.gap = 2,
      quote = FALSE)
  chi2 = x$goodness.of.fit[1]
  df = x$goodness.of.fit[2]
  pval = x$goodness.of.fit[3]
  cat("\nGoodness of Fit:\n")
  cat("\tChi2(",df,") = ",format(chi2,digits=digits),", p = ",
      format(pval,digits=digits),"\n",sep="")
  cat("\n")
  invisible(x)
}


print.summary.eba = function(x,digits=max(3,getOption("digits")-3),na.print="",
  symbolic.cor=p>4, signif.stars=getOption("show.signif.stars"),...){
  cat("\nParameter estimates:\n")
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
  cat("\n")
  printCoefmat(x$tests, digits = digits, signif.stars = signif.stars,
    zap.ind = c(1,2), ...)
  cat("\nAIC: ",format(x$aic,digits=max(4,digits+1)),"\n")
  cat("Pearson Chi2:",format(x$chi.alt,digits=digits))
  cat("\n")
  invisible(x)
}


L = function(p,y1,m,i1,i0) -sum(dbinom(y1,m,1/(1+i0%*%p/i1%*%p),log=TRUE))


strans = function(M){
  # Checks stochastic transitivities in a PCM (abs. freq.)
  # last mod: 08/Oct/2003 (works now for unbalanced design)
  # author: Florian Wickelmaier (wickelmaier@web.de)

  I = sqrt(length(as.matrix(M)))  # number of stimuli
  R = as.matrix(M/(M+t(M)))  # pcm rel. freq.
  R[which(is.nan(R),arr.ind=TRUE)] = rep(0,I)
  pre=0
  wst=0; mst=0; sst=0
  wstv = numeric(); mstv = numeric(); sstv = numeric()

  for(ii in 1:(I-2)){ for(jj in (ii+1):(I-1)){ for(kk in (jj+1):I){
    iSST=0; iMST=0; iWST=0
    iwv=1; imv=1; isv=1
    for(i in c(ii,jj,kk)){
      for(j in c(ii,jj,kk)){
        for(k in c(ii,jj,kk)){
          if(i!=j && j!=k && i!=k){
            if(R[i,j]>=.5 && R[j,k]>=.5){
              if(.5-R[i,k]<iwv) iwv = .5-R[i,k]
              if(min(R[i,j],R[j,k])-R[i,k]<imv)
                imv = min(R[i,j],R[j,k])-R[i,k]
              if(max(R[i,j],R[j,k])-R[i,k]<isv)
                isv = max(R[i,j],R[j,k])-R[i,k]
              if(R[i,k]>=.5) iWST = iWST+1
              if(R[i,k]>=min(R[i,j],R[j,k])) iMST = iMST+1
              if(R[i,k]>=max(R[i,j],R[j,k])) iSST = iSST+1
            }
          }
        }
      }
    }
    if(iSST==0){ sst = sst+1; sstv = c(sstv,isv) }
    if(iMST==0){ mst = mst+1; mstv = c(mstv,imv) }
    if(iWST==0){ wst = wst+1; wstv = c(wstv,iwv) }
    pre = pre+1
  } } }
  wv = mv = sv = 0
  if(length(wstv)) wv = wstv
  if(length(mstv)) mv = mstv
  if(length(sstv)) sv = sstv
  z = list(weak=wst,moderate=mst,strong=sst,n.tests=pre,
           wst.violations=wv,mst.violations=mv,sst.violations=sv,pcm=R)
  class(z) = "strans"
  z
}


print.strans = function(x, digits = max(3,getOption("digits")-4), ...){
  cat("\nStochastic Transitivity\n\n")
  tran = c(x$weak,x$moderate,x$strong)
  ntst = x$n.tests
  ttab = cbind(tran/ntst,
               c(mean(x$wst.violations),mean(x$mst.violations),
                 mean(x$sst.violations)),
               c(max(x$wst.violations),max(x$mst.violations),
                 max(x$sst.violations)))
 
  # 2004/MAY/04 new:
  ttran = cbind(tran,ttab)
  rownames(ttran)=c("weak","moderate","strong")
  colnames(ttran)=c("violations","error.ratio","mean.dev","max.dev")
  printCoefmat(ttran, digits = digits, signif.stars = NULL,
    zap.ind = 1, tst.ind = 0, cs.ind = 2:4, ...)

# old:
#  print.matrix(cbind(tran,format(ttab,digits=digits)),
#               rowlab=c("weak","moderate","strong"),
#               collab=c("violations","error.ratio","mean.dev","max.dev"),
#               quote=FALSE,right=TRUE)

  cat("---\nNumber of Tests: ",ntst,"\n")
  cat("\n")
  invisible(x)
}


cov.u = function(eba){
  eba -> x
  A = x$A
  cov.p = x$cov.p
  cov.u = matrix(0,length(A),length(A))
  for(i in 1:length(A)){
    for(j in 1:length(A)){
      cell = 0
      for(k in 1:length(A[[i]]))
        for(l in 1:length(A[[j]]))
          cell = cell + cov.p[A[[i]][k],A[[j]][l]]
      cov.u[i,j] = cell
    }
  }
  cov.u
}