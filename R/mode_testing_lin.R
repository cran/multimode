#A function computing the KDE number of modes

nmodes=function(data,bw,lowsup=-Inf,uppsup=Inf,n=2^15){


  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  ndata=length(data)
  if (ndata==0) stop("No observations (at least after removing missing values)")

  if (!is.numeric(bw)) stop("Argument 'bw' must be a positive integer number")
  if ((length(bw)!=1) | (bw<=0)) stop("Argument 'bw' must be a positive number")

  if (!is.numeric(lowsup)){
    warning("Argument 'lowsup' must be numeric. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (!is.numeric(uppsup)){
    warning("Argument 'uppsup' must be numeric. Default value of 'uppsup' was used")
    uppsup=Inf
  }
  if (length(lowsup)==0){
    warning("Argument 'lowsup' must be specified. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (length(uppsup)==0){
    warning("Argument 'uppsup' must be specified. Default value of 'uppsup' was used")
    uppsup=Inf
  }

  if (length(lowsup)>1) warning("Argument 'lowsup' has length > 1 and only the first element will be used")
  if (length(uppsup)>1) warning("Argument 'uppsup' has length > 1 and only the first element will be used")
  lowsup=lowsup[1]
  uppsup=uppsup[1]

  if(lowsup==uppsup){
    warning("Arguments 'lowsup' and 'uppsup' must be different. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(lowsup>uppsup){
    warning("Argument 'uppsup' must be greater than 'lowsup'. They were been interchanged")
    lowsup2=lowsup
    lowsup=uppsup
    uppsup=lowsup2
  }

  if (!is.numeric(n)){
    warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
    n=2^15
  }
  if ((length(n)!=1)|(n%%1!=0) | (n<=0)){
    warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
    n=2^15
  }



  fn=density(data,bw=bw,n=n)

  #Positions of the modes
  z=c(1:(n-1))
  re=z[diff(fn$y)>0]
  z2=c(1:length(re))
  se=z2[diff(re)>1]
  posic=re[se]
  if(re[length(re)]<(n-1)){posic=c(posic,re[length(re)])}

  posic=posic[fn$x[posic]>lowsup]
  posic=posic[fn$x[posic]<uppsup]

  #Number of modes
  num=length(posic)

  return(num)
}



###############################################################
###############################################################



#Critical bandwidth (Silverman by default)

bw.crit=function(data,mod0=1,lowsup=-Inf,uppsup=Inf,n=2^15,tol=10^(-5)){



  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  ndata=length(data)
  if (ndata==0) stop("No observations (at least after removing missing values)")
  if (!is.numeric(mod0)) stop("Argument 'mod0' must be a positive integer number")
  if ((length(mod0)!=1)|(mod0%%1!=0) | (mod0<=0)) stop("Argument 'mod0' must be a positive integer number")

  if (!is.numeric(lowsup)){
    warning("Argument 'lowsup' must be numeric. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (!is.numeric(uppsup)){
    warning("Argument 'uppsup' must be numeric. Default value of 'uppsup' was used")
    uppsup=Inf
  }
  if (length(lowsup)==0){
    warning("Argument 'lowsup' must be specified. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (length(uppsup)==0){
    warning("Argument 'uppsup' must be specified. Default value of 'uppsup' was used")
    uppsup=Inf
  }

  if (length(lowsup)>1) warning("Argument 'lowsup' has length > 1 and only the first element will be used")
  if (length(uppsup)>1) warning("Argument 'uppsup' has length > 1 and only the first element will be used")
  lowsup=lowsup[1]
  uppsup=uppsup[1]

  if(lowsup==uppsup){
    warning("Arguments 'lowsup' and 'uppsup' must be different. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(((lowsup>-Inf)&(uppsup==Inf))|((lowsup==-Inf)&(uppsup<Inf))){
    warning("Both 'lowsup' and 'uppsup' must be finite or infinite. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(lowsup>uppsup){
    warning("Argument 'uppsup' must be greater than 'lowsup'. They were been interchanged")
    lowsup2=lowsup
    lowsup=uppsup
    uppsup=lowsup2
  }

  if (!is.numeric(n)){
    warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
    n=2^15
  }
  if ((length(n)!=1)|(n%%1!=0) | (n<=0)){
    warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
    n=2^15
  }


  if (!is.numeric(tol) ){
    warning("Argument 'tol' must be a positive element. Default value of 'tol' was used")
    tol=10^(-5)
  }
  if ((tol<=0) | (length(tol)!=1)){
    warning("Argument 'tol' must be a positive element. Default value of 'tol' was used")
    tol=10^(-5)
  }



  #Initial bw
  bw0=1
  control0=0

  while(control0==0){
    if(nmodes(data,bw0,lowsup,uppsup,n)>mod0){
      control0=1
    }else{
      bw0=bw0/2
    }
  }


  #Final bw
  bwf=2*bw0
  controlf=0

  while(controlf==0){
    if(nmodes(data,bwf,lowsup,uppsup,n)<=mod0){
      controlf=1
    }else{
      bwf=bwf*2
    }
  }

  #Critical bandwidth
  controli=0

  while(controli==0){
    bwi=bw0+(bwf-bw0)/2
    numi=nmodes(data,bwi,lowsup,uppsup,n)
    if(numi<(mod0+1)){
      bwf=bwi
    }else{
      bw0=bwi
    }
    dif=bwf-bw0
    if(dif<tol){
      cbw=bwf
      controli=1
    }
  }
  return(cbw)
}


###############################################################
###############################################################



#Critical bandwidth test (Silverman, 1981)

cbws=function(data,mod0=1,B=500,n=2^10,tol=10^(-5)){

  #Test statistic
  cbw=bw.crit(data,mod0,n=n,tol=tol)

  #Bootstrap replicas
  cbwB=numeric()
  for(i in 1:B){
    eps=rnorm(length(data),0,cbw)
    samp=sample(data,length(data),replace=T)
    dataB=samp+eps
    cbwB[i]=bw.crit(dataB,mod0,n=n,tol=tol)
  }

  #P-value
  pv=mean(cbw<cbwB)

  return(c(pv,cbw))

}



###############################################################
###############################################################




#Critical bandwidth test (Hall and York, 2001)


cbwhy=function(data,lowsup,uppsup,B=500,methodhy=1,alpha=0.05,n=2^10,tol=10^(-5),nMC=100,BMC=100){

  #Hall and York is for testing H0:j=1
  mod0=1

  #Test statistic
  cbw=bw.crit(data,mod0,lowsup,uppsup,n=n,tol=tol)

  #Bootstrap replicas
  cbwB=numeric()
  for(i in 1:B){
    eps=rnorm(length(data),0,cbw)
    samp=sample(data,length(data),replace=T)
    dataB=samp+eps
    cbwB[i]=bw.crit(dataB,mod0,lowsup,uppsup,n=n,tol=tol)
  }

  #P-value

  #Method 1 of Hall and York (2001)
  if(methodhy==1){
    lambda=(0.94029*alpha^3-1.59914*alpha^2 + 0.17695*alpha+ 0.48971)/(alpha^3-1.77793*alpha^2 + 0.36162*alpha + 0.42423)
    pv=mean((cbw*lambda)<cbwB)
  }

  #Method 2 of Hall and York (2001)
  if(methodhy==2){

    #Computation of the MC p-values

    pvMC=numeric()
    for(i in 1:nMC){
      dataMC=rnorm(length(data))
      cbwMC=bw.crit(dataMC,mod0,-1.5,1.5,n=n,tol=tol)
      cbwBMC=numeric()
      for(j in 1:BMC){
        eps=rnorm(length(dataMC),0,cbwMC)
        samp=sample(dataMC,length(dataMC),replace=T)
        dataBMC=samp+eps
        cbwBMC[j]=bw.crit(dataBMC,mod0,-1.5,1.5,n=n,tol=tol)
      }
      pvMC[i]=mean(cbwMC<cbwBMC)
    }

    pvSil=mean(cbw<cbwB)
    pv=mean(pvMC<pvSil)

  }

  return(c(pv,cbw))

}



###############################################################
###############################################################

#Cramér-von Mises test statistic

cramvm=function(data,bw){

  #Value of the Kernel Distribution Estimation in the data points

  MFg=try(pnorm(outer(data,data,"-"),mean=0,sd=bw))
  if(is.numeric(MFg)){
    Fg=rowSums(MFg)/length(data)
  }else{
    Fg=rep(0,length(data))
    for(k1 in 1:length(data)){
      for (k2 in 1:length(data)){
        Fg[k1]=Fg[k1]+pnorm(data[k1]-data[k2],mean=0,sd=bw)
      }}
    Fg=Fg/length(data)
  }

  U=Fg[order(Fg)]

  #Value of the test statistic
  z3=c(1:length(data))
  sumand=(U-(2*z3-1)/(2*length(data)))^2+1/(12*length(data))
  Tk=sum(sumand)

  return(Tk)
}



###############################################################
###############################################################




#Cramér-von Mises test (Hall and York, 2001)

cbwcvm=function(data,mod0=1,B=500,n=2^10,tol=10^(-5)){

  #Test statistic
  cbw=bw.crit(data,mod0,n=n,tol=tol)
  Tk=cramvm(data,cbw)

  #Bootstrap replicas
  TkB=numeric()
  for(i in 1:B){
    eps=rnorm(length(data),0,cbw)
    samp=sample(data,length(data),replace=T)
    dataB=samp+eps
    cbwB=bw.crit(dataB,mod0,n=n,tol=tol)
    TkB[i]=cramvm(dataB,cbwB)
  }

  #P-value
  pv=mean(Tk<TkB)

  return(c(pv,Tk))

}





###############################################################
###############################################################







#Excess of mass statistic

excessmass=function(data,mod0=1){


  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  ndata=length(data)
  if (ndata==0) stop("No observations (at least after removing missing values)")
  if (!is.numeric(mod0)) stop("Argument 'mod0' must be a positive integer number")
  if ((length(mod0)!=1)|(mod0%%1!=0) | (mod0<=0)) stop("Argument 'mod0' must be a positive integer number")

  #If there are duplicated data a modification is made

  if((sum(duplicated(data))>0)|(sum(duplicated(diff(sort(data)))))){
    distance=dist(data)
    mindist=min(distance[distance>0])
    data=data+runif(ndata,-mindist/2,mindist/2)
    warning("A modification of the data was made in order to compute the excess of mass statistic")
  }


  if(mod0==1){

    em=2*dip(data)
    return(em)

  }else{


    data=sort(data)
    large=10*(data[length(data)]-data[1])
    A=matrix(0,nrow=ndata,ncol=ndata)
    A[upper.tri(A)]=-outer(data,data,"-")[upper.tri(A)]
    distanc=diff(data)

    Aord=matrix(0,nrow=ndata-1,ncol=ndata-1)
    minA=numeric()
    cualminA=numeric()

    for(i in 1:(ndata-1)){
      for(j in 1:(ndata-i)){
        Aord[i,j]=A[j,j+i]
      }
      minA[i]=min(Aord[i,1:(ndata-i)])
      cualminA[i]=which.min(Aord[i,1:(ndata-i)])
    }

    minimos=list()
    iniint=list()
    finint=list()
    iniint2=list()
    finint2=list()

    minimos[[1]]=minA
    iniint[[1]]=cualminA
    finint[[1]]=cualminA+seq(1:(ndata-1))

    for(k in 1:mod0){

      minimos[[k+1]]=minimos[[k]][(1:(length(minimos[[k]])-1))]
      for(l3 in 1:k){
        iniint2[[l3]]=iniint[[l3]][(1:(length(minimos[[k]])-1))]
        finint2[[l3]]=finint[[l3]][(1:(length(minimos[[k]])-1))]
      }
      iniint2[[k+1]]=rep(-1,(length(minimos[[k]])-1))
      finint2[[k+1]]=rep(ndata+1,(length(minimos[[k]])-1))



      for(l in 1:(length(minimos[[k]]))){

        #when we add one interval (outside) ->Bord
        Bord=matrix(large,nrow=(ndata-1),ncol=(ndata-1))
        elemB=matrix(1,nrow=(ndata-1),ncol=(ndata-1))

        #when we quit one interval (inside) ->Cord
        Cord=matrix(-1,nrow=(ndata-1),ncol=(ndata-1))

        for(l2 in 1:k){
          if(l<=(length(minimos[[k]])-2)){
            if(iniint[[l2]][l]>1){
              elemBtemp=((row(Aord)+col(Aord))<iniint[[l2]][l])&(col(Aord)<iniint[[l2]][l])
            }else{
              elemBtemp=matrix(0,nrow=(ndata-1),ncol=(ndata-1))
            }
            if(finint[[l2]][l]<(ndata-1)){
              elemBtemp=elemBtemp|((row(Aord)+col(Aord))<=ndata)&(col(Aord)>finint[[l2]][l])
            }
            if(iniint[[l2]][l]>0){
              elemB=elemB&elemBtemp}
          }
          if(l>1){
            if(iniint[[l2]][l]>0){
              Cord[col(Aord)>=iniint[[l2]][l]&col(Aord)<finint[[l2]][l]&(row(Aord)+col(Aord))<=(finint[[l2]][l])]=
                Aord[col(Aord)>=iniint[[l2]][l]&col(Aord)<finint[[l2]][l]&(row(Aord)+col(Aord))<=(finint[[l2]][l])]}
          }
        }

        if(l<=(length(minimos[[k]])-2)){
          Bord[elemB]=Aord[elemB]

          elemmin=apply(Bord,MARGIN=1,which.min)
          posmin=apply(Bord,MARGIN=1,min)+minimos[[k]][l]

          valmenor=posmin[1:(ndata-l-1-k)]<minimos[[k+1]][(l+1):length(minimos[[k+1]])]
          valmenor=(1:(ndata-l-1-k))[valmenor]
          if(length(valmenor)>0){
            minimos[[k+1]][(valmenor+l)]=posmin[valmenor]
            for(l3 in 1:k){
              iniint2[[l3]][(valmenor+l)]=iniint[[l3]][(l)]
              finint2[[l3]][(valmenor+l)]=finint[[l3]][(l)]
            }
            iniint2[[k+1]][(valmenor+l)]=elemmin[valmenor]
            finint2[[k+1]][(valmenor+l)]=elemmin[valmenor]+valmenor
          }
        }


        if(l>1){
          elemmax=apply(Cord,MARGIN=1,which.max)
          posmax=-apply(Cord,MARGIN=1,max)+minimos[[k]][l]

          valmayor=posmax[(l-1):1]<minimos[[k+1]][1:(l-1)]
          poselemmax=elemmax[(l-1):1]

          valmayor=(1:(l-1))[valmayor]
          if(length(valmayor)>0){

            minimos[[k+1]][(valmayor)]=posmax[l-valmayor]
            for(l4 in 1:length(valmayor)){
              for(l3 in 1:k){
                if((iniint[[l3]][l]>0)&(poselemmax[valmayor[l4]]>=iniint[[l3]][l])&((poselemmax[valmayor[l4]]+l-valmayor[l4])<=finint[[l3]][l])){
                  iniint2[[l3]][valmayor[l4]]=iniint[[l3]][l]
                  finint2[[l3]][valmayor[l4]]=poselemmax[valmayor[l4]]
                  iniint2[[k+1]][valmayor[l4]]=(poselemmax[valmayor[l4]]+l-valmayor[l4])
                  finint2[[k+1]][valmayor[l4]]=finint[[l3]][l]

                  if(iniint2[[l3]][valmayor[l4]]==finint2[[l3]][valmayor[l4]]){
                    iniint2[[l3]][valmayor[l4]]=-1
                    finint2[[l3]][valmayor[l4]]=ndata+1
                  }

                  if(iniint2[[k+1]][valmayor[l4]]==finint2[[k+1]][valmayor[l4]]){
                    iniint2[[k+1]][valmayor[l4]]=-1
                    finint2[[k+1]][valmayor[l4]]=ndata+1
                  }

                }else{
                  iniint2[[l3]][valmayor[l4]]=iniint[[l3]][(l)]
                  finint2[[l3]][valmayor[l4]]=finint[[l3]][(l)]
                }
              }
            }
          }

        }


      }


      iniint=iniint2
      finint=finint2

    }

    #Test statistics

    estadisticos=list()



    for(k2 in 1:(mod0+1)){

      if(k2>1){
        minAfin=minBfin
        pesAfin=pesBfin
        excmastempA=excmastempB
        lambdasA=lambdasB
        minA=minB
      }

      minB=c(0,minimos[[k2]])
      sonlamB=1
      k=1
      lambdasB=numeric()
      while(sonlamB[length(sonlamB)]<(ndata-k2+1)){
        lambtempB=(((ndata-sonlamB[length(sonlamB)]+1)-(ndata-sonlamB[length(sonlamB)]):k2)/ndata)/(minB[(ndata-sonlamB[length(sonlamB)]+2-k2)]-minB[(ndata-k2+1-sonlamB[length(sonlamB)]):1])
        sonlamB=c(sonlamB,which.min(lambtempB)+sonlamB[length(sonlamB)])
        lambdasB[k]=min(lambtempB)
        k=k+1
      }

      minBfin=minB[((ndata-k2+1):1)][sonlamB]
      pesBfin=(ndata:1)[sonlamB]/ndata
      excmastempB=pesBfin[-1]-lambdasB*minBfin[-1]


      if(k2>1){
        lambdastot=unique(sort(c(lambdasA,lambdasB)))

        excmasA=numeric()
        cualenA=lambdastot%in%lambdasA
        queintA=c(0,diff((1:length(lambdastot))[cualenA==0])-1)
        queintA=cumsum(queintA)
        queintA=queintA+((1:length(lambdastot))[(cualenA==0)])[1]

        excmasA[cualenA==0]=pesAfin[queintA]-lambdastot[cualenA==0]*minAfin[queintA]
        excmasA[cualenA]=excmastempA

        excmasB=numeric()
        cualenB=lambdastot%in%lambdasB
        queintB=c(0,diff((1:length(lambdastot))[cualenB==0])-1)
        queintB=cumsum(queintB)
        queintB=queintB+((1:length(lambdastot))[(cualenB==0)])[1]

        excmasB[cualenB==0]=pesBfin[queintB]-lambdastot[cualenB==0]*minBfin[queintB]
        excmasB[cualenB]=excmastempB

        estadisticos[[k2-1]]=max(excmasB-excmasA)

      }


    }

    return(estadisticos[[mod0]])

  }


}



###############################################################
###############################################################



# Excess of mass test (Cheng and Hall, 1998)

emch=function(data,B=500,n=2^15){

  #Cheng and Hall is for testing H0:j=1

  #library(rootSolve)
  #library(diptest)

  #Computation of the d value
  ndata=length(data)
  hest=((4/(3*ndata))^(1/5))*sd(data)
  fest=density(data,bw=hest,n=n)
  xmodest=fest$x[which.max(fest$y)]
  fmodest=max(fest$y)
  hest2=0.94*sd(data)*ndata^(-1/9)
  dataest=(xmodest-data)/hest2
  vecfdmodest=(1/ hest2^3)*((dataest)^2-1)*dnorm(dataest)
  fdmodest=mean(vecfdmodest)

  d=abs(fdmodest)/(fmodest^3)

  #Computation of the test statistic
  em=excessmass(data,1)

  #Computation of the resamples
  emB=numeric()

  if(d<2*pi){
    fun1=function (x) (beta(x,x))^2*2^(4*x-1)*(x-1)-d
    betaest=uniroot.all(fun1, c(1,256.25))
    for (j in 1:B){
      nuevdat=rbeta(ndata,betaest,betaest)
      emB[j]=excessmass(nuevdat,1)
    }
  }


  if(d>2*pi){
    fun2=function (x) 2*(beta(x-0.5,0.5))^2*x-d
    betaest=uniroot.all(fun2, c(0.5,2^50))[1]
    for (j in 1:B){
      nuevdat=rt(ndata,2*betaest-1)
      nuevdat=nuevdat/sqrt(2*betaest-1)
      emB[j]=excessmass(nuevdat,1)
    }
  }

  #P-value
  pv=mean(em<emB)
  return(c(pv,em))

}







###############################################################
###############################################################






locmodes=function(data,mod0=1,lowsup=-Inf,uppsup=Inf,n=2^15,tol=10^(-5),display=F,addplot=NULL,xlab=NULL,ylab=NULL,addLegend=NULL,posLegend=NULL){

  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  ndata=length(data)
  if (ndata==0) stop("No observations (at least after removing missing values)")
  if (!is.numeric(mod0)) stop("Argument 'mod0' must be a positive integer number")
  if ((length(mod0)!=1)|(mod0%%1!=0) | (mod0<=0)) stop("Argument 'mod0' must be a positive integer number")

  if (!is.numeric(lowsup)){
    warning("Argument 'lowsup' must be numeric. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (!is.numeric(uppsup)){
    warning("Argument 'uppsup' must be numeric. Default value of 'uppsup' was used")
    uppsup=Inf
  }
  if (length(lowsup)==0){
    warning("Argument 'lowsup' must be specified. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (length(uppsup)==0){
    warning("Argument 'uppsup' must be specified. Default value of 'uppsup' was used")
    uppsup=Inf
  }

  if (length(lowsup)>1) warning("Argument 'lowsup' has length > 1 and only the first element will be used")
  if (length(uppsup)>1) warning("Argument 'uppsup' has length > 1 and only the first element will be used")
  lowsup=lowsup[1]
  uppsup=uppsup[1]

  if(lowsup==uppsup){
    warning("Arguments 'lowsup' and 'uppsup' must be different. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(((lowsup>-Inf)&(uppsup==Inf))|((lowsup==-Inf)&(uppsup<Inf))){
    warning("Both 'lowsup' and 'uppsup' must be finite or infinite. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(lowsup>uppsup){
    warning("Argument 'uppsup' must be greater than 'lowsup'. They were been interchanged")
    lowsup2=lowsup
    lowsup=uppsup
    uppsup=lowsup2
  }

  if (!is.numeric(n)){
    warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
    n=2^15
  }
  if ((length(n)!=1)|(n%%1!=0) | (n<=0)){
    warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
    n=2^15
  }


  if (!is.numeric(tol) ){
    warning("Argument 'tol' must be a positive element. Default value of 'tol' was used")
    tol=10^(-5)
  }
  if ((tol<=0) | (length(tol)!=1)){
    warning("Argument 'tol' must be a positive element. Default value of 'tol' was used")
    tol=10^(-5)
  }

  if(display!=T&display!=F){
    warning("Argument 'display' must be T or F. Default value of 'display' was used")
    display=F
  }

  if(is.null(addplot)){
    addplot=F
  }else{
    if(addplot!=T&addplot!=F){
      warning("Argument 'addplot' must be T or F. Default value of 'addplot' was used")
      addplot=F
    }else{
      if(display==F){
        warning("Argument 'addplot' is not needed when 'display' is FALSE")
      }
    }
  }

  if(is.null(addLegend)){
    addLegend=T
  }else{
    if(addLegend!=T&addLegend!=F){
      warning("Argument 'addLegend' must be T or F. Default value of 'addLegend' was used")
      addLegend=T
    }else{
      if(display==F){
        warning("Argument 'addLegend' is not needed when 'display' is FALSE")
      }
    }
  }

  if(is.null(posLegend)){
    posLegend="topright"
  }else{
    if(display==F){
      warning("Argument 'posLegend' is not needed when 'display' is FALSE")
    }else{
      possleg= c("bottomright", "bottom", "bottomleft","left", "topleft", "top", "topright", "right", "center")
      if(!is.numeric(posLegend)&!any(posLegend==possleg)){
        warning(paste0("Argument 'posLegend' must be a valid coordinates or one of '",paste(possleg,collapse="','"),"'. Default value of 'posLegend' was used"),sep="")
        posLegend="topright"
      }
      if(is.numeric(posLegend)&(length(posLegend)!=2)){
        warning(paste0("Argument 'posLegend' must be a valid coordinates or one of '",paste(possleg,collapse="','"),"'. Default value of 'posLegend' was used"),sep="")
        posLegend="topright"
      }

    }
  }

  cbw=bw.crit(data,mod0,lowsup,uppsup,n,tol)

  fn=density(data,bw=cbw,n=n)


  #Posiciones en las que están los máximos
  z=c(1:(n-1))
  re=z[diff(fn$y)>0]
  z2=c(1:length(re))
  se=z2[diff(re)>1]
  posic=re[se]
  posicm=re[(se+1)]
  if(re[length(re)]<(n-1)){posic=c(posic,re[length(re)])}
  localization=fn$x[(posic+1)]
  localizationm=fn$x[(posicm)]

  flocalization=fn$y[(posic+1)]
  flocalizationm=fn$y[(posicm)]


  posloc2=which((localization>lowsup)&(localization<uppsup))
  localization=localization[posloc2]
  flocalization=flocalization[posloc2]

  posloc3=which((localizationm>localization[1])&(localizationm<localization[length(localization)]))
  localizationm=localizationm[posloc3]
  flocalizationm=flocalizationm[posloc3]



  localizations=numeric()
  localizations=localization[1]
  if(mod0>1){
    for(l in 1:(mod0-1)){
      localizations[(2*l)]=localizationm[l]
      localizations[(2*l+1)]=localization[(l+1)]
    }
  }

  flocalizations=numeric()
  flocalizations=flocalization[1]
  if(mod0>1){
    for(l in 1:(mod0-1)){
      flocalizations[(2*l)]=flocalizationm[l]
      flocalizations[(2*l+1)]=flocalization[(l+1)]
    }
  }

  if((!is.null(xlab)|!is.null(ylab))&addplot==T){warning("Arguments 'xlab' and 'ylab' are not needed when 'addplot' is TRUE")}

  if(is.null(xlab)){xlab=paste("N = ",ndata," Critical bandwidth = ",formatC(cbw))}
  if(is.null(ylab)){ylab="Density"}

  if(display==T){
    if(addplot==F){
        plot(NA,xlim=c(min(fn$x),max(fn$x)),ylim=c(min(fn$y),max(fn$y)),main="",ylab=ylab,xlab=xlab)
    }
    lines(fn)
    for(i in 1:mod0){
      lines(c(localization[i],localization[i]),c(-1,flocalization[i]),lty=2,lwd=2)
    }
    if(length(localizationm)>0){
      for(i in 1:length(localizationm)){
        lines(c(localizationm[i],localizationm[i]),c(-1,flocalizationm[i]),lty=3,lwd=2)
      }
    }
    lines(c(lowsup,lowsup),c(-1,2*max(fn$y)),lwd=2)
    lines(c(uppsup,uppsup),c(-1,2*max(fn$y)),lwd=2)
    if(addLegend==T){
      if((mod0==1)&(length(localizationm)==0)){
        if((lowsup>min(fn$x))|(uppsup<max(fn$x))){
          legend(posLegend,legend = c("Mode","Support"), lty = c(2,1),lwd= c(2,2))
        }else{
          legend(posLegend,legend = c("Mode"), lty = c(2),lwd= c(2))
        }
      }
      if((mod0==1)&(length(localizationm)>0)){
        if((lowsup>min(fn$x))|(uppsup<max(fn$x))){
          if(length(localizationm)==1){
            legendnames=c("Mode","Antimode","Support")
          }else{
            legendnames=c("Mode","Antimodes","Support")
          }
          legend(posLegend,legend = legendnames, lty = c(2,3,1),lwd= c(2,2,2))
        }else{
          if(length(localizationm)==1){
            legendnames=c("Mode","Antimode")
          }else{
            legendnames=c("Mode","Antimodes")
          }
          legend(posLegend,legend = legendnames, lty = c(2,3),lwd= c(2,2))
        }
      }
      if(mod0>1){
        if((lowsup>min(fn$x))|(uppsup<max(fn$x))){
          if(length(localizationm)==1){
            legendnames=c("Modes","Antimode","Support")
          }else{
            legendnames=c("Modes","Antimodes","Support")
          }
          legend(posLegend,legend = legendnames, lty = c(2,3,1),lwd= c(2,2,2))
        }else{
          if(length(localizationm)==1){
            legendnames=c("Modes","Antimode")
          }else{
            legendnames=c("Modes","Antimodes")
          }
          legend(posLegend,legend = legendnames, lty = c(2,3),lwd= c(2,2))
        }
      }
    }
  }

  loctod=list(localizations,flocalizations,cbw)
  names(loctod)=c("locations","fvalue","cbw")
  message("Estimated location")
  if(mod0>1){
    message("Modes: ",paste(formatC(localization),""))
    if(mod0>2){
      message("Antimodes: ",paste(formatC(localizationm),""))
    }else{
      message("Antimode: ",paste(formatC(localizationm),""))
    }
  }else{
    message("Mode: ",paste(formatC(localization),""))
  }
  message("")
  message("Estimated value of the density")
  if(mod0>1){
    message("Modes: ",paste(formatC(flocalization),""))
    if(mod0>2){
      message("Antimodes: ",paste(formatC(flocalizationm),""))
    }else{
      message("Antimode: ",paste(formatC(flocalizationm),""))
    }
  }else{
    message("Mode: ",paste(formatC(flocalization),""))
  }
  message("")
  message("Critical bandwidth: ", formatC(cbw))
  invisible(loctod)
}






###############################################################
###############################################################


generatorHYbw=function(data,bw,B,lowsup,uppsup,n=2^15,tol2=10^(-5)){

  #library(ks)

  ndata=length(data)
  fn=density(data,bw=bw,n=n)
  fn=density(data,bw=bw,n=n,from=min(c(min(fn$x),lowsup)),to=max(c(max(fn$x),uppsup)))


  #Positions of the modes and antimodes
  z=c(1:(n-1))
  re=z[diff(fn$y)>0]
  z2=c(1:length(re))
  se=z2[diff(re)>1]
  posic=re[se]
  posicm=re[(se+1)]
  if(re[length(re)]<(n-1)){
    posic=c(posic,re[length(re)])
    if(length(posic)>1){
      if(posic[length(posic)]==posic[(length(posic)-1)]){
        posicm=posicm[-length(posicm)]
      }}
  }

  posicm=posicm[is.na(posicm)==0]


  #Are there modes outside the support?
  controltail1=0
  controltail2=0

  localizacion=fn$x[posic]
  localizacionm=fn$x[posicm]
  pos1=which((localizacion>lowsup)&(localizacion<uppsup))
  pos2=which((localizacionm>lowsup)&(localizacionm<uppsup))
  if(pos1[1]>1){controltail1=1}
  if(pos1[length(pos1)]<length(posic)){controltail2=1}
  posic=posic[pos1]
  posicm=posicm[pos2]
  localizacion=localizacion[pos1]
  localizacionm=localizacionm[pos2]

  poscortes=numeric()

  if(length(localizacionm)>0){

    if(localizacionm[1]<localizacion[1]){
      poscortes[1]=posicm[1]+1
    }else{
      poscortes[1]=which(fn$x>(lowsup))[1]
    }

    if(localizacionm[length(localizacionm)]>localizacion[length(localizacion)]){
      poscortes[2]=posicm[length(localizacionm)]-1
    }else{
      poscortes[2]=which(fn$x>=(uppsup))[1]-1
    }

  }else{
    poscortes[1]=which(fn$x>(lowsup))[1]
    poscortes[2]=which(fn$x>=(uppsup))[1]-1
  }

  xcortes=fn$x[poscortes]
  fncortes=fn$y[poscortes]
  intfnmod=kcde(data,h=bw,eval.points=fn$x[poscortes])$estimate
  fpncortes=kdde(data, h=bw, deriv.order=1,eval.points=fn$x[poscortes])$estimate
  if(fpncortes[1]<0){fpncortes[1]=-fpncortes[1]}
  if(fpncortes[(length(fpncortes))]>0){fpncortes[(length(fpncortes))]=-fpncortes[(length(fpncortes))]}

  #Link function to remove the modes in the tails
  functionlink=function(x,x0,hl1,f0,f1,fp0,fp1)(f0+f1)/2+(1+2*((x-x0)/hl1)^3-3*((x-x0)/hl1)^2)*((f0-f1)/2)*exp((x-x0)*fp0/((f0-f1)/2))+((f0-f1)/2)*(2*((x-x0)/hl1)^3-3*((x-x0)/hl1)^2)*exp((hl1-(x-x0))*fp1/((f0-f1)/2))


  intcola=numeric()
  intcola[1]=intfnmod[1]
  intcola[2]=1-intfnmod[(length(intfnmod))]
  epscola=numeric()
  difcola=0

  if(controltail1==1){

    intfuncola1=function(xepscola){
      funcola1=function(x) functionlink(x,xcortes[1]-xepscola,xepscola,0,fncortes[1],0,fpncortes[1])
      return(abs(integrate(funcola1,xcortes[1]-xepscola,xcortes[1])[[1]]-intcola[1]))
    }

    controlepscola1=0
    extepscola1=(uppsup-lowsup)

    while(controlepscola1==0){
      epscola[1]=try(optimize(intfuncola1, interval = c(0, extepscola1), tol = tol2/4)$minimum,silent=TRUE)
      if(is.numeric(epscola[1])==0){
        extepscola1=extepscola1/1.1
        epscola=numeric()
      }else{
        controlepscola1=1
      }
    }

    xepscola=epscola[1]
    funcola1=function(x) functionlink(x,xcortes[1]-xepscola,xepscola,0,fncortes[1],0,fpncortes[1])
    difcola=intcola[1]-integrate(funcola1,xcortes[1]-xepscola,xcortes[1])[[1]]
    intcola[2]=intcola[2]+difcola
  }

  if(controltail2==1){


    intfuncola2=function(xepscola){
      funcola2=function(x) functionlink(x,xcortes[(length(xcortes))],xepscola,fncortes[(length(fncortes))],0,fpncortes[(length(fpncortes))],0)
      return(abs(integrate(funcola2,xcortes[(length(xcortes))],(xcortes[(length(xcortes))]+xepscola))[[1]]-intcola[2]))
    }


    controlepscola2=0
    extepscola2=(uppsup-lowsup)

    while(controlepscola2==0){
      epscola[2]=try(optimize(intfuncola2, interval = c(0, extepscola2), tol = tol2/4)$minimum,silent=TRUE)
      if(is.numeric(epscola[2])==0){
        extepscola2=extepscola2/1.1
        epscola=as.numeric(epscola[1])
      }else{
        controlepscola2=1
      }
    }

    xepscola=epscola[2]
    funcola2=function(x) functionlink(x,xcortes[(length(xcortes))],xepscola,fncortes[(length(fncortes))],0,fpncortes[(length(fpncortes))],0)
    difcola=intcola[2]-integrate(funcola2,xcortes[(length(xcortes))],(xcortes[(length(xcortes))]+xepscola))[[1]]
    intcola[1]=intcola[1]+difcola

  }


  if((difcola>tol2/4)&(controltail2==1)&(controltail1==1)){


    intfuncola1=function(xepscola){
      funcola1=function(x) functionlink(x,xcortes[1]-xepscola,xepscola,0,fncortes[1],0,fpncortes[1])
      return(abs(integrate(funcola1,xcortes[1]-xepscola,xcortes[1])[[1]]-intcola[1]))
    }

    controlepscola1=0
    extepscola1=(uppsup-lowsup)

    while(controlepscola1==0){
      epscola[1]=try(optimize(intfuncola1, interval = c(0, extepscola1), tol = tol2/4)$minimum,silent=TRUE)
      if(is.numeric(epscola[1])==0){
        extepscola1=extepscola1/1.1
        epscola=numeric()
      }else{
        controlepscola1=1
      }
    }

    xepscola=epscola[1]
    funcola1=function(x) functionlink(x,xcortes[1]-xepscola,xepscola,0,fncortes[1],0,fpncortes[1])
    difcola=intcola[1]-integrate(funcola1,xcortes[1]-xepscola,xcortes[1])[[1]]
    intcola[2]=intcola[2]+difcola

  }

  znueran=1
  if(controltail1==1){znueran=c(1,2)}
  if(controltail2==1){znueran=c(znueran,3)}

  qint=numeric()
  qint=1-difcola
  if(controltail1==1){
    qint=c(qint-intcola[1],intcola[1])
  }
  if(controltail2==1){
    qint[1]=qint[1]-intcola[2]
    qint=c(qint,intcola[2])
  }


  ftope=fncortes

  if(controltail1==0){xcortes[1]=-Inf}
  if(controltail2==0){xcortes[2]=Inf}

  fdef=function(x){
    fx=rep(0,length(x))

    fx=fx+replace(((x<=xcortes[1])&(x>=(xcortes[1]-epscola[1])))*functionlink(x,xcortes[1]-epscola[1],epscola[1],0,fncortes[1],0,fpncortes[1]), is.na(((x<=xcortes[1])&(x>=(xcortes[1]-epscola[1])))*functionlink(x,xcortes[1]-epscola[1],epscola[1],0,fncortes[1],0,fpncortes[1])), 0)
    fx=fx+replace(((x>=xcortes[(length(xcortes))])&(x<=(xcortes[(length(xcortes))]+epscola[2])))*functionlink(x,xcortes[(length(xcortes))],epscola[2],fncortes[(length(fncortes))],0,fpncortes[(length(fpncortes))],0), is.na(((x>=xcortes[(length(xcortes))])&(x<=(xcortes[(length(xcortes))]+epscola[2])))*functionlink(x,xcortes[(length(xcortes))],epscola[2],fncortes[(length(fncortes))],0,fpncortes[(length(fpncortes))],0)), 0)
    fx=fx+((x>xcortes[1])&(x<=xcortes[(2)]))*kde(data,h=bw,eval.points=x)$estimate
    return(fx)
  }


  Mnuevdat=matrix(0,nrow=B,ncol=ndata)

  for(j5 in 1:B){
    nuevdat=numeric()
    znue=sample(znueran,size=ndata,replace=T,prob=qint)

    for(j6 in 1:ndata){

      if(znue[j6]==1){

        for(j7 in 1:2^12){
          eps=rnorm(1,0,bw)
          nuev=sample(data,1,replace=T)
          tgen=nuev+eps
          tgenenrango=0
          if((tgen<=xcortes[2])&(tgen>xcortes[1])){tgenenrango=1}

          if(tgenenrango==1){
            nuevdat=c(nuevdat,tgen)
            break
          }
        }


      }else{

        for(j7 in 1:2^12){
          ugen=runif(1)
          vgen=runif(1)

          if(znue[j6]==2){
            tgen=xcortes[1]-epscola[1]+epscola[1]*vgen
          }

          if(znue[j6]==3){
            tgen=xcortes[(length(xcortes))]+epscola[2]*vgen
          }


          if(ftope[(znue[j6]-1)]*ugen <= fdef(tgen)){
            nuevdat=c(nuevdat,tgen)
            break
          }

        }

      }
    }

    Mnuevdat[j5,]=nuevdat
  }

  return(Mnuevdat)




}






###############################################################
###############################################################





#Excess of mass test
#(generating resamples from the KDE using the critical bw)

emcbw=function(data,mod0=1,B=500,lowsup=-Inf,uppsup=Inf,n=2^15,tol=10^(-5),tol2=10^(-5)){
  #library(ks)
  #library(diptest)
  ndata=length(data)

  #Test statistic
  em=excessmass(data,mod0)

  #Generating resamples
  cbw=bw.crit(data,mod0,lowsup,uppsup,n=n,tol=tol)
  emB=numeric()

  if((lowsup>-Inf)|(uppsup<Inf)){
    Mnuevdat=generatorHYbw(data,cbw,B,lowsup,uppsup,n,tol2)
  }

  for (j in 1:B){
    if((lowsup>-Inf)|(uppsup<Inf)){
      nuevdat=Mnuevdat[j,]
    }else{
      nuevdat=numeric()
      eps=rnorm(length(data),0,cbw)
      samp=sample(data,length(data),replace=T)
      nuevdat=samp+eps
    }
    emB[j]=excessmass(nuevdat,mod0)
  }


  #P-value
  pv=mean(em<emB)
  return(c(pv,em))
}



#############################################################
#############################################################

modetest=function(data,mod0=1,method="NP",B=500,full.result=FALSE,lowsup=-Inf,uppsup=Inf,n=NULL,tol=NULL,tol2=NULL,methodhy=NULL,alpha=NULL,nMC=NULL,BMC=NULL){

  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  ndata=length(data)
  if (ndata==0) stop("No observations (at least after removing missing values)")
  if (!any(method==c("SI","HY","FM","HH","CH","NP"))) stop("Value specified for argument 'method' is not valid")
  if (!is.numeric(B)){
    warning("Argument 'B' must be a positive integer number. Default value of 'B' was used")
    B=500
  }
  if ((length(B)!=1)|(B%%1!=0) | (B<=0)){
    warning("Argument 'B' must be a positive integer number. Default value of 'B' was used")
    B=500
  }
  if (!is.numeric(mod0)) stop("Argument 'mod0' must be a positive integer number")
  if ((length(mod0)!=1)|(mod0%%1!=0) | (mod0<=0)) stop("Argument 'mod0' must be a positive integer number")
  if ((mod0!=1)&any(method==c("HY","HH","CH"))) stop("This method is only valid for testing unimodality (mod0=1)")
  if (method=="HY"&(is.infinite(lowsup)|is.infinite(uppsup))){
    warning("Arguments 'lowsup' and 'uppsup' must be finite. Although the compact support should be imposed, a non-compact interval was used")
  }
  if (!is.numeric(lowsup)){
    warning("Argument 'lowsup' must be numeric. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (!is.numeric(uppsup)){
    warning("Argument 'uppsup' must be numeric. Default value of 'uppsup' was used")
    uppsup=Inf
  }
  if (length(lowsup)==0){
    warning("Argument 'lowsup' must be specified. Default value of 'lowsup' was used")
    lowsup=-Inf
  }
  if (length(uppsup)==0){
    warning("Argument 'uppsup' must be specified. Default value of 'uppsup' was used")
    uppsup=Inf
  }

  if (length(lowsup)>1) warning("Argument 'lowsup' has length > 1 and only the first element will be used")
  if (length(uppsup)>1) warning("Argument 'uppsup' has length > 1 and only the first element will be used")
  lowsup=lowsup[1]
  uppsup=uppsup[1]

  if(lowsup==uppsup){
    warning("Arguments 'lowsup' and 'uppsup' must be different. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(((lowsup>-Inf)&(uppsup==Inf))|((lowsup==-Inf)&(uppsup<Inf))){
    warning("Both 'lowsup' and 'uppsup' must be finite or infinite. Default values of 'lowsup' and 'uppsup' were used")
    uppsup=Inf
    lowsup=-Inf
  }

  if(lowsup>uppsup){
    warning("Argument 'uppsup' must be greater than 'lowsup'. They were been interchanged")
    lowsup2=lowsup
    lowsup=uppsup
    uppsup=lowsup2
  }

  if(full.result!=T&full.result!=F){
    warning("Argument 'full.result' must be T or F. Default value of 'full.result' was used")
    full.result=F
  }

  if (is.null(n)){
    if (any(method==c("SI","HY","FM"))){
      n=2^10
    }else{
      n=2^15}
  }else{
    if(method=="HH"){
      warning("Argument 'n' is not needed for this method")
    }else{
      if (!is.numeric(n)){
        warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
        if (any(method==c("SI","HY","FM"))){
          n=2^10
        }else{
          n=2^15}
      }
      if ((length(n)!=1)|(n%%1!=0) | (n<=0)){
        warning("Argument 'n' must be a positive integer number. Default value of 'n' was used")
        if (!any(method==c("SI","HY","FM"))){
          n=2^10
        }else{
          n=2^15}
      }
    }
  }

  if(is.null(tol)){
    tol=10^(-5)
  }else{
    if (any(method==c("HH","CH"))){
      warning("Argument 'tol' is not needed for this method")
    }else{
      if (!is.numeric(tol) ){
        warning("Argument 'tol' must be a positive element. Default value of 'tol' was used")
        tol=10^(-5)
      }
      if ((tol<=0) | (length(tol)!=1)){
        warning("Argument 'tol' must be a positive element. Default value of 'tol' was used")
        tol=10^(-5)
      }
    }
  }

  if(is.null(tol2)){
    tol2=10^(-5)
  }else{
    if (method!="NP"|lowsup==-Inf|uppsup==Inf){
      warning("Argument 'tol2' is not needed for this method")
    }else{
      if (!is.numeric(tol2) ){
        warning("Argument 'tol2' must be a positive element. Default value of 'tol2' was used")
        tol2=10^(-5)
      }
      if ((tol2<=0) | (length(tol2)!=1)){
        warning("Argument 'tol2' must be a positive element. Default value of 'tol2' was used")
        tol2=10^(-5)
      }
    }
  }

  if(is.null(methodhy)){
    methodhy=1
  }else{
    if (method!="HY"){
      warning("Argument 'methodhy' is not needed for this method")
    }else{
      if ((methodhy!=1) & (methodhy!=2)){
        warning("Argument 'methodhy' must be equal to 1 or 2. Default value of 'methodhy' was used")
        methodhy=1
      }
    }
  }

  if(is.null(alpha)){
    alpha=0.05
  }else{
    if (method!="HY"){
      warning("Argument 'alpha' is not needed for this method")
    }else{
      if (!is.numeric(alpha) ){
        warning("Argument 'alpha' must be an element in the interval (0,1). Default value of 'alpha' was used")
        alpha=0.05
      }
      if ((alpha<=0) | (length(alpha)!=1)|(alpha>=1)){
        warning("Argument 'alpha' must be an element in the interval (0,1). Default value of 'alpha' was used")
        alpha=0.05
      }
    }
  }

  if(is.null(nMC)){
    nMC=100
  }else{
    if ((method!="HY")|(methodhy==1)){
      warning("Argument 'nMC' is not needed for this method")
    }else{
      if (!is.numeric(nMC)){
        warning("Argument 'nMC' must be a positive integer number. Default value of 'nMC' was used")
        nMC=100
      }
      if ((length(nMC)!=1)|(nMC%%1!=0) | (nMC<=0)){
        warning("Argument 'nMC' must be a positive integer number. Default value of 'nMC' was used")
        nMC=100
      }
    }
  }

  if(is.null(BMC)){
    BMC=100
  }else{
    if ((method!="HY")|(methodhy==1)){
      warning("Argument 'BMC' is not needed for this method")
    }else{
      if (!is.numeric(BMC) ){
        warning("Argument 'BMC' must be a positive integer number. Default value of 'BMC' was used")
        BMC=100
      }
      if ((length(BMC)!=1)|(BMC%%1!=0) | (BMC<=0)){
        warning("Argument 'BMC' must be a positive integer number. Default value of 'BMC' was used")
        BMC=100
      }
    }
  }

  #In the excess of mass, if there are duplicated data a modification is made

  if (any(method==c("HH","CH","NP"))){
    if(sum(duplicated(data))>0){
      distance=dist(data)
      mindist=min(distance[distance>0])
      data=data+runif(ndata,-mindist/2,mindist/2)
      warning("A modification of the data was made in order to compute the excess of mass or the dip statistic")
    }
  }












  if(method=="SI"){
    pval=cbws(data,mod0,B,n,tol)
    namemethod="Silverman (1981) critical bandwidth test"
    namestatistic="Critical bandwidth"
  }else if (method=="HY"){
    pval=cbwhy(data,lowsup,uppsup,B,methodhy,alpha,n,tol,nMC,BMC)
    namemethod="Hall and York (2001) critical bandwidth test"
    namestatistic="Critical bandwidth"
  }else if (method=="FM"){
    pval=cbwcvm(data,mod0,B,n,tol)
    namemethod="Fisher and Marron (2001) Cramer-von Mises test"
    namestatistic="Cramer-von Mises"
  }else if (method=="HH"){
    pval=as.double(dip.test(data,simulate.p.value=T,B=B)[2:1])
    namemethod="Hartigan and Hartigan (1985) dip test"
    namestatistic="Dip"
  }else if (method=="CH"){
    pval=emch(data,B,n)
    namemethod="Cheng and Hall (1998) excess of mass test"
    namestatistic="Excess of mass"
  }else if (method=="NP"){
    pval=emcbw(data,mod0,B,lowsup,uppsup,n,tol,tol2)
    namemethod="Ameijeiras et al. (2016) excess of mass test"
    namestatistic="Excess of mass"
  }

  message(namemethod)
  message("")
  message(c(namestatistic," statistic: ",formatC(pval[2]),", p-value: ", formatC(pval[1])))
  if(mod0==1){
    message(c("Null hypothesis: unimodality"))
  }else{
    if(method=="NP"){
      message(c("Null hypothesis: ",mod0," modes"))
    }else{
      message(c("Null hypothesis: at most ",mod0," modes"))
    }
  }
  message(c("Alternative hypothesis: at least ",mod0+1," modes"))
  message("")

  names(pval)=c("p.value","statistic")

  if(full.result==T){
    pval=as.list(pval)
  }else{
    pval=pval[1]
  }

  invisible(pval)

}
