# Functions to use in reg.Rmd
# From Rallfun-v41.txt (https://osf.io/xhe8u/)
# Only default sub-functions are available here. 
# Source the original Rallfun file to access other options in scor().
scor<-function(x,y=NULL,corfun=pcor,gval=NA,plotit=FALSE,op=TRUE,MM=FALSE,cop=3,xlab='VAR 1',
               ylab='VAR 2',STAND=TRUE,pr=TRUE,SEED=TRUE,MC=FALSE,RAN=FALSE){
  #
  # Compute a skipped correlation coefficient.
  #
  # Eliminate outliers using a projection method
  # That is, compute Donoho-Gasko median, for each point
  # consider the line between it and the median,
  # project all points onto this line, and
  # check for outliers using a boxplot rule.
  # Repeat this for all points. A point is declared
  # an outlier if for any projection it is an outlier
  # using a modification of the usual boxplot rule.
  #
  # For information about the argument cop, see the function
  # outpro.
  #
  # Eliminate any outliers and compute correlation using
  # remaining data.
  #
  #  MC=TRUE, the multicore version of outpro is used
  #
  # corfun=pcor means Pearson's correlation is used.
  # corfun=spear means Spearman's correlation is used.
  # corfun=tau means Kendall tau is used.
  #
  #.  RAN=TRUE uses random projections instead, which results in faster execution time
  #
  if(SEED){
    oldSeed <- .Random.seed
    set.seed(12) # So when using MVE or MCD, get consistent results
  }
  if(identical(corfun,wincor))corfun=winall
  if(is.null(y[1]))m<-x
  if(!is.null(y[1]))m<-cbind(x,y)
  m<-elimna(m)
  if(!RAN){
    if(!MC)temp<-outpro(m,gval=gval,plotit=plotit,op=op,cop=cop,MM=MM,
                        xlab=xlab,ylab=ylab,STAND=STAND,pr=pr)$keep
    if(MC)temp<-outproMC(m,gval=gval,plotit=plotit,op=op,cop=cop,MM=MM,
                         xlab=xlab,ylab=ylab,STAND=STAND,pr=pr)$keep
  }
  if(RAN)temp=outpro.depth(m,MM=MM,plotit=plotit)$keep
  tcor<-corfun(m[temp,])$cor
  if(!is.null(dim((tcor))))tcor<-tcor[1,2]
  test<-abs(tcor*sqrt((nrow(m)-2)/(1-tcor**2)))
  if(ncol(m)!=2)diag(test)<-NA
  crit<-6.947/nrow(m)+2.3197
  if(SEED) {
    assign(x='.Random.seed', value=oldSeed, envir=.GlobalEnv)
  }
  list(cor=tcor,test.stat=test,crit.05=crit)
}

outpro<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
                 xlab="VAR 1",ylab="VAR 2",STAND=TRUE,tr=.2,q=.5,pr=TRUE,...){
  #
  # Detect outliers using a modification of the
  # Stahel-Donoho  projection method.
  #
  # Determine center of data cloud, for each point,
  # connect it with center, project points onto this line
  # and use distances between projected points to detect
  # outliers. A boxplot method is used on the
  # projected distances.
  #
  # plotit=TRUE creates a scatterplot when working with
  # bivariate data.
  #
  # op=T
  # means the .5 depth contour is plotted
  # based on data with outliers removed.
  #
  # op=F
  # means .5 depth contour is plotted without removing outliers.
  #
  #  MM=F  Use interquatile range when checking for outliers
  #  MM=T  uses MAD.
  #
  #  If value for center is not specified,
  #  there are four options for computing the center of the
  #  cloud of points when computing projections:
  #
  #  cop=2 uses MCD center
  #  cop=3 uses median of the marginal distributions.
  #  cop=4 uses MVE center
  #  cop=5 uses TBS
  #  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
  #
  #  args q and tr having are not used by this function. They are included to deal
  #  with situations where smoothers have optional arguments for q and tr
  #
  #  When using cop=2, 3 or 4, default critical value for outliers
  #  is square root of the .975 quantile of a
  #  chi-squared distribution with p degrees
  #  of freedom.
  #
  #  STAND=T means that marginal distributions are standardized before
  #  checking for outliers.
  #
  #  Donoho-Gasko (Tukey) median is marked with a cross, +.
  #
  m<-as.matrix(m)
  if(pr){
    if(!STAND){
      if(ncol(m)>1)print("STAND=FALSE. If measures are on different scales, might want to use STAND=TRUE")
    }}
  library(MASS)
  m=elimna(m)
  m<-as.matrix(m)
  nv=nrow(m)
  if(ncol(m)==1){
    dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
    dis<-sqrt(dis)
    dis[is.na(dis)]=0
    crit<-sqrt(qchisq(.975,1))
    chk<-ifelse(dis>crit,1,0)
    vec<-c(1:nrow(m))
    outid<-vec[chk==1]
    keep<-vec[chk==0]
  }
  if(ncol(m)>1){
    M=m
    if(STAND)m=standm(m,est=median,scat=mad)
    if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
    if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
    if(cop==1 && is.na(center[1])){
      if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
      if(ncol(m)==2){
        tempd<-NA
        for(i in 1:nrow(m))
          tempd[i]<-depth(m[i,1],m[i,2],m)
        mdep<-max(tempd)
        flag<-(tempd==mdep)
        if(sum(flag)==1)center<-m[flag,]
        if(sum(flag)>1)center<-apply(m[flag,],2,mean)
      }}
    if(cop==2 && is.na(center[1])){
      center<-cov.mcd(m)$center
    }
    if(cop==4 && is.na(center[1])){
      center<-cov.mve(m)$center
    }
    if(cop==3 && is.na(center[1])){
      center<-apply(m,2,median)
    }
    if(cop==5 && is.na(center[1])){
      center<-tbs(m)$center
    }
    if(cop==6 && is.na(center[1])){
      center<-rmba(m)$center
    }
    if(cop==7 && is.na(center[1])){
      center<-spat(m)
    }
    flag<-rep(0, nrow(m))
    outid <- NA
    vec <- c(1:nrow(m))
    for (i in 1:nrow(m)){
      B<-m[i,]-center
      dis<-NA
      BB<-B^2
      bot<-sum(BB)
      if(bot!=0){
        for (j in 1:nrow(m)){
          A<-m[j,]-center
          temp<-sum(A*B)*B/bot
          dis[j]<-sqrt(sum(temp^2))
        }
        temp<-idealf(dis)
        if(!MM)cu<-median(dis)+gval*(temp$qu-temp$ql)
        if(MM)cu<-median(dis)+gval*mad(dis)
        outid<-NA
        temp2<-(dis> cu)
        flag[temp2]<-1
      }}
    if(sum(flag) == 0) outid <- NA
    if(sum(flag) > 0)flag<-(flag==1)
    outid <- vec[flag]
    idv<-c(1:nrow(m))
    keep<-idv[!flag]
    if(ncol(m)==2){
      if(plotit){
        m=M # plot data using the original scale.
        plot(m[,1],m[,2],type="n",xlab=xlab,ylab=ylab)
        points(m[keep,1],m[keep,2],pch="*")
        if(length(outid)>0)points(m[outid,1],m[outid,2],pch="o")
        if(op){
          tempd<-NA
          keep<-keep[!is.na(keep)]
          mm<-m[keep,]
          for(i in 1:nrow(mm))tempd[i]<-depth(mm[i,1],mm[i,2],mm)
          mdep<-max(tempd)
          flag<-(tempd==mdep)
          if(sum(flag)==1)center<-mm[flag,]
          if(sum(flag)>1)center<-apply(mm[flag,],2,mean)
          m<-mm
        }
        points(center[1],center[2],pch="+")
        x<-m
        temp<-fdepth(m,plotit=FALSE)
        flag<-(temp>=median(temp))
        xx<-x[flag,]
        xord<-order(xx[,1])
        xx<-xx[xord,]
        temp<-chull(xx)
        lines(xx[temp,])
        lines(xx[c(temp[1],temp[length(temp)]),])
      }}}
  list(n=nv,n.out=length(outid),out.id=outid,keep=keep)
}

pcor<-function(x,y=NA){
  if(!is.na(y[1]))temp<-wincor(x,y,tr=0)
  if(is.na(y[1]))temp<-winall(x,tr=0)
  list(cor=temp$cor,p.value=temp$p.value)
}



wincor<-function(x,y=NULL,tr=.2){
  #   Compute the Winsorized correlation between x and y.
  #
  #   tr is the amount of Winsorization
  #   This function also returns the Winsorized covariance
  #
  #    Pairwise deletion of missing values is performed.
  #
  #   x is a vector, or it can be a matrix with two columns when y=NULL
  #
  
  if(!is.null(y[1])){
    m=cbind(x,y)
  }
  else m=x
  m<-elimna(m)
  nval=nrow(m)
  if(ncol(m)==2){
    a=wincor.sub(m[,1],m[,2],tr=tr)
    wcor=a$cor
    wcov=a$cov
    sig=a$p.value
  }
  if(ncol(m)>2){
    #if(is.data.frame(m))m=as.matrix(m)
    if(!is.matrix(m))stop("The data must be stored in a n by p matrix")
    wcor<-matrix(1,ncol(m),ncol(m))
    wcov<-matrix(0,ncol(m),ncol(m))
    siglevel<-matrix(NA,ncol(m),ncol(m))
    for (i in 1:ncol(m)){
      ip<-i
      for (j in ip:ncol(m)){
        val<-wincor.sub(m[,i],m[,j],tr)
        wcor[i,j]<-val$cor
        wcor[j,i]<-wcor[i,j]
        if(i==j)wcor[i,j]<-1
        wcov[i,j]<-val$cov
        wcov[j,i]<-wcov[i,j]
        if(i!=j){
          siglevel[i,j]<-val$p.value
          siglevel[j,i]<-siglevel[i,j]
        }
      }}
    sig=siglevel
  }
  list(n=nval,cor=wcor,cov=wcov,p.value=sig)
}

wincor.sub<-function(x,y,tr=tr){
  sig<-NA
  g<-floor(tr*length(x))
  xvec<-winval(x,tr)
  yvec<-winval(y,tr)
  wcor<-cor(xvec,yvec)
  wcov<-var(xvec,yvec)
  if(sum(x==y)!=length(x)){
    test<-wcor*sqrt((length(x)-2)/(1.-wcor^2))
    sig<-2*(1-pt(abs(test),length(x)-2*g-2))
  }
  list(cor=wcor,cov=wcov,p.value=sig)
}

standm<-function(x,locfun=lloc,est=mean,scat=var,...){
  # standardize a matrix x
  #
  x=elimna(x)
  x=as.matrix(x)
  m1=lloc(x,est=est)
  v1=apply(x,2,scat)
  p=ncol(x)
  for(j in 1:p)x[,j]=(x[,j]-m1[j])/sqrt(v1[j])
  x
}

lloc<-function(x,est=tmean,...){
  if(is.data.frame(x)){
    x=as.matrix(x)
    x=apply(x,2,as.numeric) # earlier versions of R require this command
  }
  if(!is.list(x))val<-est(x,...)
  if(is.list(x))val=lapply(x,est,...)
  if(is.matrix(x))val<-apply(x,2,est,...)
  val
}

idealf<-function(x,na.rm=FALSE){
  #
  # Compute the ideal fourths for data in x
  #
  if(na.rm)x<-x[!is.na(x)]
  j<-floor(length(x)/4 + 5/12)
  y<-sort(x)
  g<-(length(x)/4)-j+(5/12)
  ql<-(1-g)*y[j]+g*y[j+1]
  k<-length(x)-j+1
  qu<-(1-g)*y[k]+g*y[k-1]
  list(ql=ql,qu=qu)
}

winall<-function(m,tr=.2){
  #
  #    Compute the Winsorized correlation and covariance matrix for the
  #    data in the n by p matrix m.
  #
  #    This function also returns the two-sided significance level
  #
  if(is.data.frame(m))m=as.matrix(m)
  if(!is.matrix(m))stop("The data must be stored in a n by p matrix")
  wcor<-matrix(1,ncol(m),ncol(m))
  wcov<-matrix(0,ncol(m),ncol(m))
  siglevel<-matrix(NA,ncol(m),ncol(m))
  for (i in 1:ncol(m)){
    ip<-i
    for (j in ip:ncol(m)){
      val<-wincor(m[,i],m[,j],tr)
      wcor[i,j]<-val$cor
      wcor[j,i]<-wcor[i,j]
      if(i==j)wcor[i,j]<-1
      wcov[i,j]<-val$cov
      wcov[j,i]<-wcov[i,j]
      if(i!=j){
        siglevel[i,j]<-val$p.value
        siglevel[j,i]<-siglevel[i,j]
      }
    }
  }
  cent=apply(m,2,mean,tr,na.rm=TRUE)
  list(cor=wcor,cov=wcov,center=cent,p.values=siglevel)
}

winval<-function(x,tr=.2){
  #
  #  Winsorize the data in the vector x.
  #  tr is the amount of Winsorization which defaults to .2.
  #
  #  This function is used by several other functions that come with this book.
  #
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  winval<-ifelse(x<=xbot,xbot,x)
  winval<-ifelse(winval>=xtop,xtop,winval)
  winval
}


depth<-function(U,V,m){
  #
  #  Compute the halfspace depth of the point (u,v) for the pairs of points
  #  in the n by 2 matrix m.
  #
  X<-m[,1]
  Y<-m[,2]
  FV<-NA
  NUMS<-0
  NUMH<-0
  SDEP<-0.0
  HDEP<-0.0
  N<-length(X)
  P<-acos(-1)
  P2<-P*2.0
  EPS<-0.000001
  ALPHA<-NA
  NT<-0
  for(i in 1:nrow(m)){
    DV<-sqrt(((X[i]-U)*(X[i]-U)+(Y[i]-V)*(Y[i]-V)))
    if (DV <= EPS){
      NT<-NT+1
    }
    else{
      XU<-(X[i]-U)/DV
      YU<-(Y[i]-V)/DV
      if (abs(XU) > abs(YU)){
        if (X[i] >= U){
          ALPHA[i-NT]<-asin(YU)
          if(ALPHA[i-NT] < 0.0)
            ALPHA[i-NT]<-P2+ALPHA[i-NT]
        }
        else{
          ALPHA[i-NT]<-P-asin(YU)
        }
      }
      else{
        if (Y[i] >= V)
          ALPHA[i-NT]<-acos(XU)
        else
          ALPHA[i-NT]<-P2-acos(XU)
      }
      if (ALPHA[i-NT] >= P2-EPS) ALPHA[i-NT]<-0.0
    }
  }
  NN<-N-NT
  if(NN<=1){
    NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
      depths1(NT,3)
    if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
    NUMH<-NUMH+NT
    HDEP<-(NUMH+0.0)/(N+0.0)
    return(HDEP)
  }
  ALPHA<-sort(ALPHA[1:NN])
  ANGLE<-ALPHA[1]-ALPHA[NN]+P2
  for(i in 2:NN){
    ANGLE<-max(c(ANGLE,ALPHA[i]-ALPHA[i-1]))
  }
  if(ANGLE > (P+EPS)){
    NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
      depths1(NT,3)
    if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
    NUMH<-NUMH+NT
    HDEP<-(NUMH+0.0)/(N+0.0)
    return(HDEP)
  }
  ANGLE<-ALPHA[1]
  NU<-0
  for (i in 1:NN){
    ALPHA[i]<-ALPHA[i]-ANGLE
    if(ALPHA[i]<(P-EPS))NU<-NU+1
  }
  if(NU >= NN){
    NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
      depths1(NT,3)
    if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
    NUMH<-NUMH+NT
    HDEP<-(NUMH+0.0)/(N+0.0)
    return(HDEP)
  }
  #
  #  Mergesort the alpha with their antipodal angles beta,
  #  and at the same time update I, F(I), and NBAD.
  #
  JA<-1
  JB<-1
  ALPHK<-ALPHA[1]
  BETAK<-ALPHA[NU+1]-P
  NN2<-NN*2
  NBAD<-0
  I<-NU
  NF<-NN
  for(J in 1:NN2){
    ADD<-ALPHK+EPS
    if (ADD < BETAK){
      NF<-NF+1
      if(JA < NN){
        JA<-JA+1
        ALPHK<-ALPHA[JA]
      }
      else
        ALPHK<-P2+1.0
    }
    else{
      I<-I+1
      NN1<-NN+1
      if(I==NN1){
        I<-1
        NF<-NF-NN
      }
      FV[I]<-NF
      NFI<-NF-I
      NBAD<-NBAD+depths1(NFI,2)
      if(JB < NN){
        JB<-JB+1
        if(JB+NU <= NN)
          BETAK<-ALPHA[JB+NU]-P
        else
          BETAK<-ALPHA[JB+NU-NN]+P
      }
      else
        BETAK<-P2+1.0
    }
  }
  NUMS<-depths1(NN,3)-NBAD
  #
  #  Computation of NUMH for halfspace depth.
  #
  GI<-0
  JA<-1
  ANGLE<-ALPHA[1]
  dif<-NN-FV[1]
  NUMH<-min(FV[1],dif)
  for(I in 2:NN){
    AEPS<-ANGLE+EPS
    if(ALPHA[I] <= AEPS){
      JA<-JA+1
    }
    else{
      GI<-GI+JA
      JA<-1
      ANGLE<-ALPHA[I]
    }
    KI<-FV[I]-GI
    NNKI<-NN-KI
    NUMH<-min(c(NUMH,min(c(KI,NNKI))))
  }
  NUMS<-NUMS+depths1(NT,1)*depths1(NN,2)+depths1(NT,2)*depths1(NN,1)+
    depths1(NT,3)
  if(N >= 3)SDEP<-(NUMS+0.0)/(depths1(N,3)+0.0)
  NUMH<-NUMH+NT
  HDEP<-(NUMH+0.0)/(N+0.0)
  HDEP
}

depths1<-function(m,j){
  if(m < j)depths1<-0
  else{
    if(j==1)depths1<-m
    if(j==2)depths1<-(m*(m-1))/2
    if(j==3)depths1<-(m*(m-1)*(m-2))/6
  }
  depths1
}


fdepth<-function(m,pts=NA,plotit=TRUE,cop=3,center=NA,xlab="VAR 1",
                 ylab="VAR 2"){
  #
  # Determine depth of points in pts,  relative to
  # points in m. If pts is not specified,
  # depth of all points in m are determined.
  #
  # m and pts can be vectors or matrices with
  # p columns (the number of variables).
  #
  # Determine center, for each point, draw a line
  # connecting it with center, project points onto this line
  # and determine depth of the projected points.
  # The final depth of a point is its minimum depth
  # among all projections.
  #
  # plotit=TRUE creates a scatterplot when working with
  # bivariate data and pts=NA
  #
  #  There are three options for computing the center of the
  #  cloud of points when computing projections, assuming center=NA:
  #
  #  cop=2 uses MCD center
  #  cop=3 uses median of the marginal distributions.
  #  cop=4 uses MVE center
  #
  #  If a value for center is passed to this function,
  #  this value is used to determine depths.
  #
  #  When plotting,
  #  center is marked with a cross, +.
  #
  library(MASS)
  if(cop!=2 && cop!=3 && cop!=4)stop("Only cop=2, 3 or 4 is allowed")
  if(is.list(m))stop("Store data in a matrix; might use function listm")
  m<-as.matrix(m)
  pts<-as.matrix(pts)
  if(!is.na(pts[1]))remm<-m
  nm<-nrow(m)
  nm1<-nm+1
  if(!is.na(pts[1])){
    if(ncol(m)!=ncol(pts))stop("Number of columns of m is not equal to number of columns for pts")
  }
  m<-elimna(m) # Remove missing values
  m<-as.matrix(m)
  if(ncol(m)==1)dep<-unidepth(as.vector(m[,1]),pts=pts)
  if(ncol(m)>1){
    if(is.na(center[1])){
      if(cop==2){
        center<-cov.mcd(m)$center
      }
      if(cop==4){
        center<-cov.mve(m)$center
      }
      if(cop==3){
        center<-apply(m,2,median)
      }}
    if(is.na(pts[1])){
      mdep <- matrix(NA,nrow=nrow(m),ncol=nrow(m))
    }
    if(!is.na(pts[1])){
      mdep <- matrix(NA,nrow=nrow(m),ncol=nrow(pts))
    }
    for (i in 1:nrow(m)){
      B<-m[i,]-center
      dis<-NA
      BB<-B^2
      bot<-sum(BB)
      if(bot!=0){
        if(is.na(pts[1])){
          for (j in 1:nrow(m)){
            A<-m[j,]-center
            temp<-sum(A*B)*B/bot
            dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
          }}
        if(!is.na(pts[1])){
          m<-rbind(remm,pts)
          for (j in 1:nrow(m)){
            A<-m[j,]-center
            temp<-sum(A*B)*B/bot
            dis[j]<-sign(sum(A*B))*sqrt(sum(temp^2))
          }}
        #
        # For ith projection, store depths of
        # points in mdep[i,]
        #
        if(is.na(pts[1]))mdep[i,]<-unidepth(dis)
        if(!is.na(pts[1])){
          mdep[i,]<-unidepth(dis[1:nm],dis[nm1:nrow(m)])
        }}
      if(bot==0)mdep[i,]<-rep(0,ncol(mdep))
    }
    dep<-apply(mdep,2,min)
    if(ncol(m)==2 && is.na(pts[1])){
      flag<-chull(m)
      dep[flag]<-min(dep)
    }
  }
  if(ncol(m)==2){
    if(is.na(pts[1]) && plotit){
      plot(m,xlab=xlab,ylab=ylab)
      points(center[1],center[2],pch="+")
      x<-m
      temp<-dep
      flag<-(temp>=median(temp))
      xx<-x[flag,]
      xord<-order(xx[,1])
      xx<-xx[xord,]
      temp<-chull(xx)
      xord<-order(xx[,1])
      xx<-xx[xord,]
      temp<-chull(xx)
      lines(xx[temp,])
      lines(xx[c(temp[1],temp[length(temp)]),])
    }}
  dep<-round(dep*nrow(m))/nrow(m)
  dep
}

unidepth<-function(x,pts=NA){
  #
  # Determine depth of points in the vector x
  #
  if(!is.vector(x))stop("x should be a vector")
  if(is.na(pts[1]))pts<-x
  pup<-apply(outer(pts,x,FUN="<="),1,sum)/length(x)
  pdown<-apply(outer(pts,x,FUN="<"),1,sum)/length(x)
  pdown<-1-pdown
  m<-matrix(c(pup,pdown),nrow=2,byrow=TRUE)
  dep<-apply(m,2,min)
  dep
}


ols<-function(x,y,xout=FALSE,outfun=outpro,alpha=.05,plotit=FALSE,xlab='X',ylab='Y',zlab='Z',RES=TRUE,...){
  #
  # Performs OLS regression calling built-in R function.
  #
  # xout=T will eliminate any leverage points (outliers among x values)
  # if one predictor,
  # plotit=TRUE will plot the points and the regression line
  #
  m<-elimna(cbind(x,y))
  n=nrow(m)
  n.keep=n
  x<-as.matrix(x)
  p<-ncol(x)
  pp<-p+1
  x<-m[,1:p]
  y<-m[,pp]
  if(xout){
    m<-cbind(x,y)
    flag<-outfun(x,plotit=FALSE,...)$keep
    m<-m[flag,]
    x<-m[,1:p]
    y<-m[,pp]
    n.keep=length(y)
  }
  x<-as.matrix(x)
  temp<-summary(lm(y~x))
  coef<-temp[4]$coefficients
  CI=matrix(NA,nrow(coef),ncol=2)
  CI[,1]=coef[,1]-qt(1-alpha/2,temp[10]$fstatistic[3])*coef[,2]
  CI[,2]=coef[,1]+qt(1-alpha/2,temp[10]$fstatistic[3])*coef[,2]
  dimnames(CI)=list(NULL,c("low.ci","up.ci"))
  coef=cbind(coef,CI)
  if(plotit){
    if(p==1){
      plot(x,y,xlab=xlab,ylab=ylab)
      abline(coef[,1])
    }
    if(p==2){
      regp2plot(x,y,regfun=ols,xlab=xlab,ylab=ylab,zlab=zlab)
    }}
  Ftest<-temp[10]$fstatistic
  Ftest.p.value<-1-pf(Ftest[1],Ftest[2],Ftest[3])
  Rval=Rsq.ols(x,y)
  res=NULL
  if(RES)res=y-x%*%coef[2:pp,1]-coef[1,1]
  list(n=n,n.keep=n.keep,summary=coef,coef=coef[,1],F.test=temp[10]$fstatistic[1],Ftest.p.value=Ftest.p.value,
       F.test.degrees.of.freedom=temp[10]$fstatistic[2:3],R.squared=Rval,residuals=as.vector(res))
}

Rsq.ols<-function(x,y){
  res=lsfit(x,y)$residuals
  yhat=y-res
  rsq=var(yhat)/var(y)
  rsq
}





