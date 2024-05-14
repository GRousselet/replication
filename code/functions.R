keeporder <- function(x){
  x <- as.character(x)
  x <- factor(x, levels=unique(x))
  x
}

# Formula from Spence & Stanley (2024)
# https://journals.sagepub.com/doi/10.1177/25152459231217932
predint <- function(n1, n2, m1, sd1, alpha = 0.05){
  ci <- vector(mode = "numeric", length = 2)
  se <- sqrt( (sd1^2 / n1) + (sd1^2 / n2))
  df <- n1 - 1
  tval <- qt(1-alpha/2, df)
  ci[1] <- m1 - qt(1-alpha/2,df)*se
  ci[2] <- m1 + qt(1-alpha/2,df)*se
  ci
}

predint.tm <- function(M, sd, orig.n, new.n, tr=0, alpha=.05){
  # Generalization of prediction method in Spence & Stanley 
  # Advances in Methods and Practices in Psychological Science January-March 2024, Vol. 7, No. 1,
  # pp. 1–13
  # by Rand R Wilcox, Feb 2024
  #
  # M = observed mean or can be a trimmed mean
  # tr = amount of trimming
  # sd = Winsorized standard deviation, which is the usual standard deviation when t=0
  # sd should be computed with winsd(x,tr)
  orig.n=orig.n*(1-2*tr)^2
  new.n=new.n*(1-2*tr)^2
  se=sqrt(sd^2/orig.n+sd^2/new.n)
  g=floor(tr*orig.n)
  df=orig.n-2*g-1
  crit=qt(1-alpha/2,df)
  ci=M-crit*se
  ci=c(ci,M+crit*se)
  ci
}

# from Rand R Wilcox's Rallfun-v41.txt
winvar <- function(x, tr=.2, na.rm=FALSE){
  #
  #  Compute the gamma Winsorized variance for the data in the vector x.
  #  tr is the amount of Winsorization which defaults to .2.
  #
  remx=x
  x<-x[!is.na(x)]
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  y<-ifelse(y<=xbot,xbot,y)
  y<-ifelse(y>=xtop,xtop,y)
  wv<-var(y)
  if(!na.rm)if(sum(is.na(remx)>0))wv=NA
  wv
}

winsd<-function(x,tr=.2,na.rm=FALSE){
  val=sqrt(winvar(x,tr=tr,na.rm=na.rm))
  val
}


trimci<-function(x,tr=.2,alpha=.05,null.value=0,nullval=NULL){
  #
  #  Compute a 1-alpha confidence interval for the trimmed mean
  #
  #  The default amount of trimming is tr=.2
  #
  if(!is.null(nullval))null.value=nullval
  x<-x[!is.na(x)]
  se<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
  trimci<-vector(mode="numeric",length=2)
  df<-length(x)-2*floor(tr*length(x))-1
  trimci[1]<-mean(x,tr)-qt(1-alpha/2,df)*se
  trimci[2]<-mean(x,tr)+qt(1-alpha/2,df)*se
  test<-(mean(x,tr)-null.value)/se
  sig<-2*(1-pt(abs(test),df))
  list(estimate=mean(x,tr),ci=trimci,test.stat=test,se=se,p.value=sig,n=length(x))
}


# http://www.pallier.org/pdfs/aprime.pdf
dprime <- function(hit,fa){ 
  qnorm(hit) - qnorm(fa)
}

aprime <-function(hit,fa){ 
  a<-1/2+((hit-fa)*(1+hit-fa) / (4*hit*(1-fa))) 
  b<-1/2-((fa-hit)*(1+fa-hit) / (4*fa*(1-hit)))
  a[fa>hit]<-b[fa>hit]
  a[fa==hit]<-.5
  a }

aprime.vec <-function(hits,falsealarms){ 
  L <- length(hits)
  res <- vector(mode = "numeric", length = L)
  for(E in 1:L){
    hit <- hits[E]
    fa <- falsealarms[E]
    a <- 1/2+((hit-fa)*(1+hit-fa) / (4*hit*(1-fa))) 
    b <- 1/2-((fa-hit)*(1+fa-hit) / (4*fa*(1-hit)))
    a[fa>hit] <- b[fa>hit]
    a[fa==hit] <- .5
    res[E] <- a
  }
  res
}
# --------------------------------------

sim.counter <- function(S, nsim, inc){
  if(S == 1){
    # print(paste(nsim,"iterations:",S))
    cat(nsim,"iterations:",S)
    beep(2)
  }
  if(S %% inc == 0){
    # print(paste("iteration",S,"/",nsim))
    cat(" /",S)
    beep(2)
  }
}

# From Rand Wilcox Rallfun.txt
# https://osf.io/xhe8u/
cnorm<-function(n,epsilon=.1,k=10){
  #
  # generate n observations from a contaminated normal
  # distribution
  # probability 1-epsilon from a standard normal
  # probability epsilon from normal with mean 0 and standard deviation k
  #
  if(epsilon>1)stop("epsilon must be less than or equal to 1")
  if(epsilon<0)stop("epsilon must be greater than or equal to 0")
  if(k<=0)stop("k must be greater than 0")
  val<-rnorm(n)
  uval<-runif(n)
  flag<-(uval<=1-epsilon)
  val[!flag]<-k*val[!flag]
  val
}

trimci<-function(x,tr=.2,alpha=.05,null.value=0,pr=TRUE,nullval=NULL){
  #
  #  Compute a 1-alpha confidence interval for the trimmed mean
  #
  #  The default amount of trimming is tr=.2
  #
  if(pr){
    print("The p-value returned by this function is based on the")
    print("null value specified by the argument null.value, which defaults to 0")
    print('To get a measure of effect size using a Winsorized measure of scale,  use trimciv2')
  }
  if(!is.null(nullval))null.value=nullval
  x<-elimna(x)
  se<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
  trimci<-vector(mode="numeric",length=2)
  df<-length(x)-2*floor(tr*length(x))-1
  trimci[1]<-mean(x,tr)-qt(1-alpha/2,df)*se
  trimci[2]<-mean(x,tr)+qt(1-alpha/2,df)*se
  test<-(mean(x,tr)-null.value)/se
  sig<-2*(1-pt(abs(test),df))
  list(ci=trimci,estimate=mean(x,tr),test.stat=test,se=se,p.value=sig,n=length(x))
}

# MODIFIED TO RETURN ONLY THE P VALUE
trimci.pval <- function(x, tr=.2, null.value){
  #
  #  The default amount of trimming is tr=.2
  #
  x <- elimna(x)
  se <- sqrt(winvar(x,tr)) / ((1-2*tr) * sqrt(length(x)))
  df <- length(x) -2 * floor(tr*length(x)) -1
  test <- (mean(x,tr) - null.value) / se
  pval <- 2*(1-pt(abs(test),df))
  pval
}

elimna<-function(m){
  #
  # remove any rows of data having missing values
  #
  DONE=FALSE
  if(is.list(m) && is.matrix(m)){
    z=pool.a.list(m)
    m=matrix(z,ncol=ncol(m))
    DONE=TRUE
  }
  if(!DONE){
    if(is.list(m) && is.matrix(m[[1]])){
      for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
      for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    #if(!is.list(m)){
    #if(is.null(dim(m)))
    m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    e<-m[ikeep[ikeep>=1],]
    #}
  }
  e
}

winvar<-function(x,tr=.2,na.rm=FALSE,STAND=NULL){
  #
  #  Compute the gamma Winsorized variance for the data in the vector x.
  #  tr is the amount of Winsorization which defaults to .2.
  #
  remx=x
  x<-x[!is.na(x)]
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  y<-ifelse(y<=xbot,xbot,y)
  y<-ifelse(y>=xtop,xtop,y)
  wv<-var(y)
  if(!na.rm)if(sum(is.na(remx)>0))wv=NA
  wv
}

ghdist<-function(n,g=0,h=0){
  #
  # generate n observations from a g-and-h dist.
  #
  x<-rnorm(n)
  if (g>0){
    ghdist<-(exp(g*x)-1)*exp(h*x^2/2)/g
  }
  if(g==0)ghdist<-x*exp(h*x^2/2)
  ghdist
}

ghmean<-function(g,h){
  #
  #Compute the mean and variance of a g-and-h distribution
  #
  val=0
  if(h==0){
    if(g>0){
      val=(exp(g^2/2)-1)/g
      val2=(1-2*exp(g^2/2)+exp(2*g^2))/g^2
      val2=val2-val^2
    }}
  if(g>0 & h!=0){
    if(h<1)
      val=(exp(g^2/(2*(1-h)))-1)/(g*sqrt(1-h))
    val2=NA
    if(h>0){
      if(h<.5)
        val2=(exp(2*g^2/(1-2*h))-2*exp(g^2/(2*(1-2*h)))+1)/(g^2*sqrt(1-2*h))-
          (exp(g^2/(2*(1-h)))-1)^2/(g^2*(1-h))
    }}
  if(g==0){
    val=0
    val2=1/(1-2*h)^1.5   #Headrick et al. (2008)
  }
  list(mean=val,variance=val2)
}

ghtrim<-function(tr=.2,g=0,h=0){
  #
  #  Compute trimmed mean of a g-and-h distribution.
  #
  # 
  if(g==0)val=0
  if(g>0){
    low=qnorm(tr)
    up=-1*low
    val=integrate(ftrim,low,up,tr=tr,g=g,h=h)$value
    val=val/(1-2*tr)
  }
  val
}

ftrim<-function(z,tr,g,h){
  gz=(exp(g*z)-1)*exp(h*z^2/2)/g
  res=dnorm(z)*gz
  res
}




