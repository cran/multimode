modetree=function(data,bws=NULL,gridsize=NULL,cbw1=NULL,cbw2=NULL,display=TRUE,logbw=FALSE,logbw.regulargrid=NULL,...){

  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  ndata=length(data)
  if (ndata==0) stop("No observations (at least after removing missing values)")

  if(is.null(gridsize)){
    gridsize=c(512,151)
    if(!is.null(bws)){
      if(length(bws)>2){
        gridsize[2]=length(bws)
      }}
  }else{
    if (!is.numeric(gridsize)){
      warning("Argument 'gridsize' must be numeric. Default values of 'gridsize' were used")
      gridsize=c(512,151)
    }else{
      if(length(gridsize)!=2){
        warning("Argument 'gridsize' must be of length two. Default values of 'gridsize' were used")
        gridsize=c(512,151)
      }
    }
  }

  n=gridsize[1]
  n.user=n
  n=max(n, 512)
  if (n > 512){n=2^ceiling(log2(n))}
  range.x=seq.int(min(data),max(data), length.out = n.user)
  range.data=range(data)


  if(logbw!=T&logbw!=F){
    warning("Argument 'logbw' must be T or F. Default value of 'logbw' was used")
    logbw=F
  }


  if(logbw==F&!is.null(logbw.regulargrid)){
    warning("Argument 'logbw.regulargrid' is not employed when 'logbw' is F")
    logbw=T
  }


  if(is.null(logbw.regulargrid)){
    logbw.regulargrid=F
  }

  if(logbw==T&(logbw.regulargrid!=T&logbw.regulargrid!=F)){
    warning("Argument 'logbw.regulargrid' must be T or F. Default value of 'logbw.regulargrid' was used")
    logbw.regulargrid=F
  }



  if(!is.null(bws)&(!is.null(cbw1)|!is.null(cbw2))){warning("Arguments 'cbw1' and 'cbw2' are not needed when 'bws' are given")}

  if(!is.null(cbw1)&!is.null(cbw2)&is.null(bws)){

    if (!is.numeric(cbw1)){
      warning("Argument 'cbw1' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if ((length(cbw1)!=1)|(cbw1%%1!=0) | (cbw1<=0)){
      warning("Argument 'cbw1' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if (!is.numeric(cbw2)){
      warning("Argument 'cbw2' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if ((length(cbw2)!=1)|(cbw2%%1!=0) | (cbw2<=0)){
      warning("Argument 'cbw2' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }

  }

  if((!is.null(cbw1)&is.null(cbw2)&is.null(bws))|(is.null(cbw1)&!is.null(cbw2)&is.null(bws))){
    warning("Both arguments 'cbw1' and 'cbw2' must be specified. Default values of 'bws' were used")
    cbw1=NULL;cbw2=NULL
  }

  if(!is.null(cbw1)&!is.null(cbw2)&is.null(bws)){

    if(cbw1>=cbw2){
      warning("Argument 'cbw2' must be greater than 'cbw1'. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }else{
      hmin=bw.crit(data,cbw2)
      hmax=bw.crit(data,cbw1)
      range.h=seq(hmin,hmax,len=gridsize[2])
      if(logbw.regulargrid==T){
        range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
        range.h=10^range.h
      }
    }
  }

  if(!is.null(bws)){

    if (!is.numeric(bws)){
      warning("Argument 'bws' must be numeric. Default values of 'bws' were used")
      bws=NULL
    }else{
      if(length(bws)<2){
        warning("Argument 'bws' must be of length greater than one. Default values of 'bws' were used")
        bws=NULL
      }
      if(length(bws)==2){
        hmax=bws[2]
        hmin=bws[1]
        warning("A grid of 'bws' between the given ones were used")
        range.h=seq(hmin,hmax,len=gridsize[2])
        if(logbw.regulargrid==T){
          range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
          range.h=10^range.h
        }
      }
      if(length(bws)>2){
        hmax=max(bws)
        hmin=min(bws)
        range.h=bws
        if(length(bws)!=gridsize[2]){
          warning(paste("The 'gridsize' for the 'bws' is equal to",length(bws)))
        }
      }
    }
  }


  if(is.null(bws)&is.null(cbw2)&is.null(cbw1)){
    hmin=2*min(diff(range.x))
    hmax=diff(range.data)
    range.h=seq(hmin,hmax,len=gridsize[2])
    if(logbw.regulargrid==T){
      range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
      range.h=10^range.h
    }
  }

  if(display!=T&display!=F){
    warning("Argument 'display' must be T or F. Default value of 'display' was used")
    display=T
  }


  range.h=sort(range.h,decreasing=T)
  maxmodes=nmodes(data,range.h[length(range.h)],n=n)

  matmodes=matrix(NA,nrow=length(range.h),ncol=maxmodes)
  l=1
  for(bw in range.h){
    fn=density(data,bw=bw,n=n)

    #Positions of the modes
    z=c(1:(n-1))
    re=z[diff(fn$y)>0]
    z2=c(1:length(re))
    se=z2[diff(re)>1]
    posic=re[se]
    if(re[length(re)]<(n-1)){posic=c(posic,re[length(re)])}
    matmodes[l,1:length(posic)]=fn$x[posic]
    l=l+1
  }

  if(logbw==T){
    range.h=log10(range.h)
    hmin=log10(hmin)
    hmax=log10(hmax)
  }

  matmodes2=matmodes
  linemt2=matrix(NA,nrow=length(range.h),ncol=maxmodes)
  l2=1
  for(l in 1:maxmodes){
    controlmt=0
    while(controlmt==0){
      modetemp=matmodes2[l2,]
      if(sum(!is.na(modetemp))>0){
        modetemp=modetemp[!is.na(modetemp)]
        modetemp=modetemp[1]
        matmodes2[l2,which(matmodes2[l2,]==modetemp)]=NA
        if(l>1){
          lh=which.min(abs(modetemp-linemt2[l2,]))
          #lines(c(modetemp,linemt2[l2,lh]),rep(range.h[l2],2),lty=2,col=col.lines[2])
        }
        controlmt=1
      }else{
        l2=l2+1
      }
    }
    linemt=modetemp
    if(l2<length(range.h)){
      for(l3 in (l2+1):length(range.h)){
        whichlmt=which.min(abs(matmodes2[l3,]-linemt[l3-l2]))
        linemt=c(linemt,matmodes2[l3,whichlmt])
        matmodes2[l3,whichlmt]=NA
      }
    }
    #lines(linemt,range.h[l2:length(range.h)],lwd=2,col=col.lines[1])
    linemt2[l2:length(range.h),l]=linemt
  }

  row.names(linemt2)=formatC(range.h)
  colnames(linemt2)=paste("Mode",1:maxmodes)
  callname=match.call()
  modetreeret=list(linemt2,range.data,range.h,logbw,ndata,callname)
  names(modetreeret)=c("locations","range.x","range.bws","logbw","sample.size","call")
  class(modetreeret)="gtmod"
  if(display==T){plot(modetreeret,...)}
  return(modetreeret)
}



###############################################################
###############################################################


print.gtmod<-function(x, digits = getOption("digits"),...){
  stopifnot(is.numeric(sample.size <- x$sample.size),is.numeric(range.x <- x$range.x),is.numeric(range.bws <- x$range.bws),is.logical(logbw<-x$logbw),is.call(callname<-x$call))
  cat("\nCall:\n\t", deparse(callname))
  cat("\n\nn=",sample.size,". Location values range: ",sep="")
  cat(formatC(range(range.x), digits=digits,...))
  cat("\nNumber of employed bandwidths: ",length(x$range.bws),". log10 bandwidths scale: ",logbw,sep="")
  cat("\nBandwidths range:",formatC(range(range.bws), digits=digits,...))
  cat("\n\n")
}



###############################################################
###############################################################

plot.gtmod<-function(x,addplot=FALSE,xlab=NULL,ylab=NULL,col.lines="black",col.sizer=NULL,addlegend=TRUE,poslegend="topright",...){

  stopifnot(is.numeric(ndata <- x$sample.size),is.numeric(range.x <- x$range.x),is.numeric(range.bws <- x$range.bws),is.logical(logbw<-x$logbw))

  if(addplot!=T&addplot!=F){
    warning("Argument 'addplot' must be T or F. Default value of 'addplot' was used")
    addplot=F
  }

  if((!is.null(xlab)|!is.null(ylab))&addplot==T){warning("Arguments 'xlab' and 'ylab' are not needed when 'addplot' is TRUE")}
  if(is.null(xlab)){xlab=paste("N =",ndata)}
  if(is.null(ylab)&logbw==FALSE){ylab="bandwidths"}
  if(is.null(ylab)&logbw==TRUE){ylab=expression(log[10] * "(bandwidths)")}

  if(length(col.lines)!=1&length(col.lines)!=2){
    warning("Argument 'col.lines' must be of length one or two. Default values of 'col.lines' were used")
    col.lines="black"
  }

  if(length(col.lines)==1){
    col.lines=rep(col.lines,2)
  }

  if(!is.null(x$locations)){

    if(addplot==F){
      plot(NA,ylim=range(range.bws),xlim=range(range.x),xlab=xlab,ylab=ylab,...)
    }

    for(j in 1:dim(x$locations)[2]){
      lines(x$locations[,j],x$range.bws,lwd=2,col=col.lines[1])
    }
    if(dim(x$locations)[2]>1){
      for(j in 2:dim(x$locations)[2]){
        whichfirst=which(!is.na(x$locations[,j]))[1]
        whichclose=which.min(abs(x$locations[whichfirst,j]-x$locations[whichfirst,1:(j-1)]))
        lines(c(x$locations[whichfirst,j],x$locations[whichfirst,whichclose]),rep(x$range.bws[whichfirst],2),lty=2,col=col.lines[2])
      }
    }
  }



  if(!is.null(x$modeforest)){
    pixdif=length(unique(as.double(x$modeforest)))-1
    col2=grey((pixdif:0)/pixdif)
    image(range.x,range.bws,x$modeforest,col=col2,xlab=xlab,ylab=ylab,...)

    box()

  }



  if(!is.null(x$sizer)){


    if(!is.null(col.sizer)){
      if(length(col.sizer)!=4 & x$method>1){
        warning("Argument 'col.sizer' must be of length four. Default values of 'col.sizer' were used")
        col.sizer=NULL
      }
      if(length(col.sizer)!=3 & x$method==1){
        if(length(col.sizer)==4){
          warning("The locations where the data are not dense are not plotted for this method")
        }else{
          warning("Argument 'col.sizer' must be of length three. Default values of 'col.sizer' were used")
          col.sizer=NULL
        }
      }
    }

    if(is.null(col.sizer)){
      col.sizer=c("red","orchid","blue","grey")
    }


    if(addlegend!=T&addlegend!=F){
      warning("Argument 'addlegend' must be T or F. Default value of 'addlegend' was used")
      addlegend=T
    }


    col2=col.sizer[min(which(tabulate(x$sizer)>0)):max(which(tabulate(x$sizer)>0))]

    image(range.x,range.bws,x$sizer,col=col2,xlab=xlab,ylab=ylab,...)
    box()

    legendnames=numeric()
    legendnames[3]=expression("Sign. incr." * ""%up%"")
    legendnames[1]=expression("Sign. dec." * ""%down%"")
    legendnames[2]=expression("Not sign."!=0)
    legendnames[4]="Sparse data"
    elemnonzero=c(3,2,1,4)[c(3,2,1,4)%in%which(tabulate(x$sizer,4)>0)]

    colleg=col.sizer[elemnonzero]
    legendnames2=legendnames[elemnonzero]

    if(addlegend==T){
      legend(poslegend,legend = legendnames2, pch=rep(22,length(colleg)),col=rep("white",length(colleg)),pt.bg=colleg,pt.cex=1.5,pt.lwd=1.5,bty="n")
    }
  }


}



###############################################################
###############################################################


summary.gtmod<-function(object, bandwidths=TRUE, digits = getOption("digits"), width=1, levelmf=NULL, ...){
  stopifnot(is.numeric(object$range.x),is.numeric(object$range.bws),is.logical(object$logbw))

  totnmodes=0

  if(!is.null(object$locations)){

    if(!is.null(levelmf)){
      warning("Argument 'levelmf' is not employed for this method")
    }

   totnmodes=dim(object$locations)[2]
   if(totnmodes>0){
   modesval=matrix(0,nrow=totnmodes,ncol=2)
   bandval=matrix(0,nrow=totnmodes,ncol=2)
   for(i in 1:totnmodes){
     modesval[i,]=range(object$locations[,i],na.rm=T)
     bandval[i,]=range(object$range.bws[!is.na(object$locations[,i])])
   }
   }
  }

  if(!is.null(object$modeforest)){

    if(is.null(levelmf)){
      levelmf=0.5
    }

    estmodes=matrix(0,dim(object$modeforest)[1],dim(object$modeforest)[2])
    modecand=rowSums(object$modeforest>levelmf)>0
    wmodecand=which(modecand)
    totnmodes=0
    if(length(wmodecand)>0){
    j=1
    labelmode=1
    i=wmodecand[1]
    while(j<length(wmodecand)){
     i=wmodecand[j]
     estmodes[i,which(object$modeforest[i,]>levelmf)]=labelmode
     j=j+1
     dim2mf=dim(object$modeforest)[2]
     cond1=(object$modeforest[i,]>levelmf)&(object$modeforest[i+1,]>levelmf)
     cond2=(object$modeforest[i,1:(dim2mf-1)]>levelmf)&(object$modeforest[i+1,2:dim2mf]>levelmf)
     cond3=(object$modeforest[i,2:dim2mf]>levelmf)&(object$modeforest[i+1,1:(dim2mf-1)]>levelmf)

     while((sum(cond1)+sum(cond2)+sum(cond3))>0){
       i=wmodecand[j]
       estmodes[i,which(object$modeforest[i,]>levelmf)]=labelmode
       j=j+1
       dim2mf=dim(object$modeforest)[2]
       cond1=(object$modeforest[i,]>levelmf)&(object$modeforest[i+1,]>levelmf)
       cond2=(object$modeforest[i,1:(dim2mf-1)]>levelmf)&(object$modeforest[i+1,2:dim2mf]>levelmf)
       cond3=(object$modeforest[i,2:dim2mf]>levelmf)&(object$modeforest[i+1,1:(dim2mf-1)]>levelmf)
     }
     labelmode=labelmode+1
    }
    totnmodes=max(estmodes)
    modesval=matrix(0,nrow=totnmodes,ncol=2)
    bandval=matrix(0,nrow=totnmodes,ncol=2)
    diffrx=diff(object$range.x)[1]/2
    for(i in 1:totnmodes){
      modesval[i,]=range(object$range.x[which(rowSums(estmodes==i)>0)])
      modesval[i,1]=modesval[i,1]-diffrx
      modesval[i,2]=modesval[i,2]+diffrx
      bandval[i,]=range(object$range.bws[which(colSums(estmodes==i)>0)])
    }
    orderestmodes=order(bandval[,2]-bandval[,1],decreasing=T)
    modesval=modesval[orderestmodes,]
    bandval=bandval[orderestmodes,]
    }
  }


  if(!is.null(object$sizer)){

    if(!is.null(levelmf)){
      warning("Argument 'levelmf' is not employed for this method")
    }

    estmodes=matrix(0,dim(object$sizer)[1],dim(object$sizer)[2])
    flatsizer=matrix(0,dim(object$sizer)[1],dim(object$sizer)[2])
    for(i in 1:(dim(object$sizer)[2])){
      wincreasing=which(object$sizer[,i]==3)
      wdecreasing=which(object$sizer[,i]==1)
      lastwinc=c(wincreasing[which(diff(wincreasing)>1)],wincreasing[length(wincreasing)])
      wddec=which(diff(wdecreasing)>1)+1
      wddec=c(1,wddec)
      firstwdec=wdecreasing[wddec]
      if(sum(is.na(lastwinc))){lastwinc=numeric()}
      if(sum(is.na(firstwdec))){firstwdec=numeric()}
      if((length(lastwinc)>0)&(length(firstwdec)>0)){
        j=1
        while(j<=length(firstwdec)){
          if(sum(lastwinc<firstwdec[j])>0){
            wlf=which(lastwinc<firstwdec[j])
            startmode=lastwinc[wlf[length(wlf)]]
            flatsizer[(startmode+1):(firstwdec[j]-1),i]=1
            flw=firstwdec>lastwinc[wlf[length(wlf)]+1]
            if(sum(flw,na.rm=T)==0){
              j=length(firstwdec)+1
            }else{
              j=which(flw)[1]
            }
          }else{
            j=j+1
          }
        }

      }
    }

    modecand=rowSums(flatsizer)>0
    wmodecand=which(modecand)
    totnmodes=0
    if(length(wmodecand)>0){
      j=1
      labelmode=1
      i=wmodecand[1]
      dim2mf=dim(object$sizer)[2]
      while(j<length(wmodecand)){
        i=wmodecand[j]
        estmodes[i,which(flatsizer[i,]==1)]=labelmode
        j=j+1
        cond1=(flatsizer[i,]==1)&(flatsizer[i+1,]==1)
        cond2=(flatsizer[i,1:(dim2mf-1)]==1)&(flatsizer[i+1,2:dim2mf]==1)
        cond3=(flatsizer[i,2:dim2mf]==1)&(flatsizer[i+1,1:(dim2mf-1)]==1)

        while((sum(cond1)+sum(cond2)+sum(cond3))>0){
          i=wmodecand[j]
          estmodes[i,which(flatsizer[i,]==1)]=labelmode
          j=j+1
          cond1=(flatsizer[i,]==1)&(flatsizer[i+1,]==1)
          cond2=(flatsizer[i,1:(dim2mf-1)]==1)&(flatsizer[i+1,2:dim2mf]==1)
          cond3=(flatsizer[i,2:dim2mf]==1)&(flatsizer[i+1,1:(dim2mf-1)]==1)

        }
        labelmode=labelmode+1
      }
      totnmodes=max(estmodes)
      modesval=matrix(0,nrow=totnmodes,ncol=2)
      bandval=matrix(0,nrow=totnmodes,ncol=2)
      diffrx=diff(object$range.x)[1]/2
      for(i in 1:totnmodes){
        modesval[i,]=range(object$range.x[which(rowSums(estmodes==i)>0)])
        modesval[i,1]=modesval[i,1]-diffrx
        modesval[i,2]=modesval[i,2]+diffrx
        bandval[i,]=range(object$range.bws[which(colSums(estmodes==i)>0)])
      }
      orderestmodes=order(bandval[,2]-bandval[,1],decreasing=T)
      modesval=modesval[orderestmodes,]
      bandval=bandval[orderestmodes,]
    }
  }

  if(totnmodes>0){
  cat("\nThe estimated modes are located between the following values")
  if(bandwidths==TRUE){
  if(object$logbw==T){
  cat("\n(for the log10 bandwidths range indicated in parenthesis)\n")
  }
  if(object$logbw==F){
    cat("\n(for the bandwidths range indicated in parenthesis)\n")
  }
  }
    if(!is.matrix(modesval)){modesval=matrix(modesval,ncol=2)}
    if(!is.matrix(bandval)){bandval=matrix(bandval,ncol=2)}
  for(i in 1:totnmodes){
    cat(paste("\nMode ", i,": ",sep=""),
    paste(formatC(modesval[i,1], digits=digits,width=width,...),"|",formatC(modesval[i,2], digits=digits,width=width,...)))
    if(bandwidths==TRUE){
    cat(" (",formatC(bandval[i,1], digits=digits,width=width,...)," ",formatC(bandval[i,2], digits=digits,width=width,...),")",sep="")
    }
  }

  cat("\n\n")
  }else{
    cat("\nNo modes were found according to the given criteria")
  }
}






###############################################################
###############################################################





sizer=function(data,method=2,bws=NULL,gridsize=NULL,alpha=0.05,
               B=NULL,n0=NULL,cbw1=NULL,cbw2=NULL,display=TRUE,
               logbw=TRUE,from=NULL,to=NULL,logbw.regulargrid=NULL,...){

  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  N=as.integer(length(data))
  nx=N
  if (nx==0) stop("No observations (at least after removing missing values)")

  if(display!=T&display!=F){
    warning("Argument 'display' must be T or F. Default value of 'display' was used")
    display=T
  }

  if(logbw!=T&logbw!=F){
    warning("Argument 'logbw' must be T or F. Default value of 'logbw' was used")
    logbw=T
  }

  if(logbw==F&!is.null(logbw.regulargrid)){
    warning("Argument 'logbw.regulargrid' is not employed when 'logbw' is F")
    logbw=T
  }


  if(is.null(logbw.regulargrid)){
    logbw.regulargrid=F
  }

  if(logbw==T&(logbw.regulargrid!=T&logbw.regulargrid!=F)){
    warning("Argument 'logbw.regulargrid' must be T or F. Default value of 'logbw.regulargrid' was used")
    logbw.regulargrid=F
  }


  if(is.null(gridsize)){
    gridsize=c(512,151)
    if(!is.null(bws)){
      if(length(bws)>2){
        gridsize[2]=length(bws)
      }}
  }else{
    if (!is.numeric(gridsize)){
      warning("Argument 'gridsize' must be numeric. Default values of 'gridsize' were used")
      gridsize=c(512,151)
    }else{
      if(length(gridsize)!=2){
        warning("Argument 'gridsize' must be of length two. Default values of 'bws' were used")
        gridsize=c(512,151)
      }
    }
  }

  if(is.null(from)){from=min(data)}
  if(is.null(to)){to=max(data)}

  if (!is.numeric(from)){
    warning("Argument 'from' must be numeric. Default value of 'from' was used")
    from=min(data)
  }
  if (!is.numeric(to)){
    warning("Argument 'to' must be numeric. Default value of 'to' was used")
    to=max(data)
  }

  if (length(from)>1) warning("Argument 'from' has length > 1 and only the first element will be used")
  if (length(to)>1) warning("Argument 'to' has length > 1 and only the first element will be used")
  from=from[1]
  to=to[1]

  if(from==-Inf&to==Inf){
    warning("Both 'from' and 'to' must be finite. Default values of 'from' and 'to' were used")
    to=max(data)
    from=min(data)
  }

  if(from==to){
    warning("Arguments 'from' and 'to' must be different. Default values of 'from' and 'to' were used")
    to=max(data)
    from=min(data)
  }


  if(from>to){
    warning("Argument 'to' must be greater than 'from'. They were been interchanged")
    from2=from
    from=to
    to=from2
  }


  n=gridsize[1]
  n.user=n
  n=max(n, 512)
  if (n > 512){n=2^ceiling(log2(n))}
  range.x=seq.int(from, to, length.out = n.user)
  range.data=range(data)

  if(!is.null(bws)&(!is.null(cbw1)|!is.null(cbw2))){warning("Arguments 'cbw1' and 'cbw2' are not needed when 'bws' are given")}

  if(!is.null(cbw1)&!is.null(cbw2)&is.null(bws)){

    if (!is.numeric(cbw1)){
      warning("Argument 'cbw1' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if ((length(cbw1)!=1)|(cbw1%%1!=0) | (cbw1<=0)){
      warning("Argument 'cbw1' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if (!is.numeric(cbw2)){
      warning("Argument 'cbw2' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if ((length(cbw2)!=1)|(cbw2%%1!=0) | (cbw2<=0)){
      warning("Argument 'cbw2' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }

  }

  if((!is.null(cbw1)&is.null(cbw2)&is.null(bws))|(is.null(cbw1)&!is.null(cbw2)&is.null(bws))){
    warning("Both arguments 'cbw1' and 'cbw2' must be specified. Default values of 'bws' were used")
    cbw1=NULL;cbw2=NULL
  }

  if(!is.null(cbw1)&!is.null(cbw2)&is.null(bws)){

    if(cbw1>=cbw2){
      warning("Argument 'cbw2' must be greater than 'cbw1'. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }else{
      hmin=bw.crit(data,cbw2)
      hmax=bw.crit(data,cbw1)
      range.h=seq(hmin,hmax,len=gridsize[2])
      if(logbw.regulargrid==T){
        range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
        range.h=10^range.h
      }
    }
  }


  if(!is.null(bws)){

    if (!is.numeric(bws)){
      warning("Argument 'bws' must be numeric. Default values of 'bws' were used")
      bws=NULL
    }else{
      if(length(bws)<2){
        warning("Argument 'bws' must be of length greater than one. Default values of 'bws' were used")
        bws=NULL
      }
      if(length(bws)==2){
        hmax=bws[2]
        hmin=bws[1]
        warning("A grid of 'bws' between the given ones were used")
        range.h=seq(hmin,hmax,len=gridsize[2])
        if(logbw.regulargrid==T){
          range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
          range.h=10^range.h
        }
      }
      if(length(bws)>2){
        hmax=max(bws)
        hmin=min(bws)
        range.h=bws
        if(length(bws)!=gridsize[2]){
          warning(paste("The 'gridsize' for the 'bws' is equal to",length(bws)))
        }
      }
    }
  }


  if(is.null(bws)&is.null(cbw2)&is.null(cbw1)){
    hmin=2*min(diff(range.x))
    hmax=diff(range.data)
    range.h=seq(hmin,hmax,len=gridsize[2])
    if(logbw.regulargrid==T){
      range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
      range.h=10^range.h
    }
  }


  if ((method!=1) & (method!=2) & (method!=3) & (method!=4)){
    warning("Argument 'method' must be equal to 1, 2, 3 or 4. Default value of 'method' was used")
    method=2
  }


  if(!is.null(B)&method<3){
    warning("Argument 'B' is not needed for this method")
  }

  if(is.null(B)){B=100}
  if (!is.numeric(B)){
    warning("Argument 'B' must be a positive integer number. Default value of 'B' was used")
    B=100
  }
  if ((length(B)!=1)|(B%%1!=0) | (B<=0)){
    warning("Argument 'B' must be a positive integer number. Default value of 'B' was used")
    B=100
  }


  if(!is.null(n0)&method==1){
    warning("Argument 'n0' is not needed for this method")
  }

  if(is.null(n0)){n0=5}
  if (!is.numeric(n0)){
    warning("Argument 'n0' must be a positive number. Default value of 'n0' was used")
    n0=5
  }
  if ((length(n0)!=1)|(n0<0)){
    warning("Argument 'n0' must be a positive number. Default value of 'n0' was used")
    n0=5
  }





  kernel="gaussian"
  weights=rep.int(1/nx, nx)
  totMass=nx/N

  fhatp=matrix(0,gridsize[1],gridsize[2])
  cfkords2=list()

  for(ibw in 1:length(range.h)){
    bw=range.h[ibw]
    lo=from - 4 * bw
    up=to + 4 * bw
    y=.Call(BinDist2, data, weights, lo, up, n) * totMass
    kords2=seq.int(0, 2 * (up - lo), length.out = 2L * n)
    kords2[(n + 2):(2 * n)]=-kords2[n:2]
    kords2=switch(kernel, gaussian = -kords2*dnorm(kords2, sd = bw))
    cfkords2[[ibw]]=Conj(fft(kords2))
    kords=fft(fft(y) * cfkords2[[ibw]], inverse = TRUE)
    kords=Re(kords)[1L:n]/length(y)
    xords=seq.int(lo, up, length.out = n)
    y = approx(xords, kords, range.x)$y
    fhatp[,ibw] = -y/bw^2
  }

  fhatp2=matrix(0,gridsize[1],gridsize[2])

  for(ibw in 1:length(range.h)){
    bw=range.h[ibw]
    lo=from - 4 * bw
    up=to + 4 * bw
    y=.Call(BinDist2, data, weights, lo, up, n) * totMass
    kords2=seq.int(0, 2 * (up - lo), length.out = 2L * n)
    kords2[(n + 2):(2 * n)]=-kords2[n:2]
    kords=switch(kernel, gaussian = (kords2*dnorm(kords2, sd = bw))^2)
    kords=fft(fft(y) * Conj(fft(kords)), inverse = TRUE)
    kords=Re(kords)[1L:n]/length(y)
    xords=seq.int(lo, up, length.out = n)
    y = approx(xords, kords, range.x)$y
    fhatp2[,ibw] = y/bw^4
  }

  varfhatp=t(t((fhatp2-(fhatp)^2))/range.h)/nx
  eps <- 10*.Machine$double.eps
  varfhatp[varfhatp<=0]=eps
  dtfhatp=sqrt(varfhatp)


  #pointwise Gaussian quantiles
  if(method==1){
    qalpha=qnorm(1-alpha/2)
    upquan=fhatp+qalpha*dtfhatp
    downquan=fhatp-qalpha*dtfhatp
  }


  if(method>1){
    #ESS
    ESS=matrix(0,gridsize[1],gridsize[2])
    for(ibw in 1:length(range.h)){
      bw=range.h[ibw]
      ESS[,ibw]=density(data,bw=bw,from=from,to=to,n=gridsize[1])$y
    }
    ESS=t(t(ESS)*range.h)*nx*sqrt(2*pi)

  }


  #Approximate simultaneous over x Gaussian quantiles
  if(method==2){
    mh=nx/colMeans(ESS)
    qalpha=(1+(1-alpha)^(1/mh))/2
    upquan=fhatp+t(qalpha*t(dtfhatp))
    downquan=fhatp-t(qalpha*t(dtfhatp))
  }

  #bootstrap simultaneous over x (and h)

  if(method==3|method==4){

    if(method==3){
      qalphamatrix=matrix(0,B,gridsize[2])
    }
    if(method==4){
      qalphavector=rep(0,B)
    }

    for(boot in 1:B){

      databoot=sample(data,replace=T)

      fhatpboot=matrix(0,gridsize[1],gridsize[2])

      for(ibw in 1:length(range.h)){
        bw=range.h[ibw]
        lo=from - 4 * bw
        up=to + 4 * bw
        y=.Call(BinDist2, databoot, weights, lo, up, n) * totMass
        kords=fft(fft(y) * cfkords2[[ibw]], inverse = TRUE)
        kords=Re(kords)[1L:n]/length(y)
        xords=seq.int(lo, up, length.out = n)
        y = approx(xords, kords, range.x)$y
        fhatpboot[,ibw] = -y/bw^2
      }

      Zxh=abs((fhatp-fhatpboot)/dtfhatp)
      Zxh[ESS<n0]=0

      if(method==3){
        qalphamatrix[boot,]=apply(Zxh,2,"max")
      }

      if(method==4){
        qalphavector[boot]=max(Zxh)
      }

    }

  }

  if(method==3){
    quanalpha=function(x) quantile(x,alpha/2)
    qalpha=apply(qalphamatrix,2,quanalpha)
    upquan=fhatp+t(qalpha*t(dtfhatp))
    downquan=fhatp-t(qalpha*t(dtfhatp))
  }

  if(method==4){
    qalpha=quantile(qalphavector,alpha/2)
    upquan=fhatp+qalpha*dtfhatp
    downquan=fhatp-qalpha*dtfhatp
  }

  matrixSiZer=sign(upquan)+sign(downquan)

  #Plot in grey ESS<n0 (the set of x locations where the data are dense)
  if(method>1){
    matrixSiZer[ESS<n0]=4
  }

  matrixSiZer=(matrixSiZer/2)+2


  if(method==1){
    ESS="ESS is not computed for this method"
  }

  if(logbw==T){
    range.h=log10(range.h)
  }

  colnames(matrixSiZer)=formatC(range.h)
  colnames(downquan)=formatC(range.h)
  colnames(fhatp)=formatC(range.h)
  colnames(upquan)=formatC(range.h)
  if(method>1){colnames(ESS)=formatC(range.h)}

  rownames(matrixSiZer)=formatC(range.x)
  rownames(downquan)=formatC(range.x)
  rownames(fhatp)=formatC(range.x)
  rownames(upquan)=formatC(range.x)
  if(method>1){rownames(ESS)=formatC(range.x)}

  callname=match.call()
  sizerret=list(matrixSiZer,method,downquan,fhatp,upquan,ESS,range.x,range.h,logbw,nx,callname)
  names(sizerret)=c("sizer","method","lower.CI","estimate","upper.CI","ESS","range.x","range.bws","logbw","sample.size","call")

  class(sizerret)="gtmod"
  if(display==T){plot(sizerret,...)}
  return(sizerret)

}







modeforest=function(data,bws=NULL,gridsize=NULL,B=99,n=512,cbw1=NULL,cbw2=NULL,display=TRUE,logbw=FALSE,from=NULL,to=NULL,logbw.regulargrid=NULL,...){

  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data))>0) warning("Missing values were removed")
  data=data[!is.na(data)]
  N=as.integer(length(data))
  nx=N
  if (nx==0) stop("No observations (at least after removing missing values)")

  if(display!=T&display!=F){
    warning("Argument 'display' must be T or F. Default value of 'display' was used")
    display=T
  }

  if(logbw!=T&logbw!=F){
    warning("Argument 'logbw' must be T or F. Default value of 'logbw' was used")
    logbw=T
  }


  if(logbw==F&!is.null(logbw.regulargrid)){
    warning("Argument 'logbw.regulargrid' is not employed when 'logbw' is F")
    logbw=T
  }


  if(is.null(logbw.regulargrid)){
    logbw.regulargrid=F
  }

  if(logbw==T&(logbw.regulargrid!=T&logbw.regulargrid!=F)){
    warning("Argument 'logbw.regulargrid' must be T or F. Default value of 'logbw.regulargrid' was used")
    logbw.regulargrid=F
  }

  if(is.null(gridsize)){
    gridsize=c(100,151)
    if(!is.null(bws)){
      if(length(bws)>2){
        gridsize[2]=length(bws)
      }}
  }else{
    if (!is.numeric(gridsize)){
      warning("Argument 'gridsize' must be numeric. Default values of 'gridsize' were used")
      gridsize=c(100,151)
    }else{
      if(length(gridsize)!=2){
        warning("Argument 'gridsize' must be of length two. Default values of 'bws' were used")
        gridsize=c(100,151)
      }
    }
  }

  nmf=gridsize[1]
  n.user=n
  n=max(n, 512)
  if (n > 512){n=2^ceiling(log2(n))}
  range.data=range(data)

  slack=2*diff(range.data)/n
  if(is.null(from)){from=min(data)-slack}
  if(is.null(to)){to=max(data)+slack}

  if (!is.numeric(from)){
    warning("Argument 'from' must be numeric. Default value of 'from' was used")
    from=min(data)
  }
  if (!is.numeric(to)){
    warning("Argument 'to' must be numeric. Default value of 'to' was used")
    to=max(data)
  }

  if (length(from)>1) warning("Argument 'from' has length > 1 and only the first element will be used")
  if (length(to)>1) warning("Argument 'to' has length > 1 and only the first element will be used")
  from=from[1]
  to=to[1]

  if(from==-Inf&to==Inf){
    warning("Both 'from' and 'to' must be finite. Default values of 'from' and 'to' were used")
    to=max(data)
    from=min(data)
  }

  if(from==to){
    warning("Arguments 'from' and 'to' must be different. Default values of 'from' and 'to' were used")
    to=max(data)
    from=min(data)
  }


  if(from>to){
    warning("Argument 'to' must be greater than 'from'. They were been interchanged")
    from2=from
    from=to
    to=from2
  }


  range.x=seq.int(from, to, length.out = n.user)
  range.mf=seq.int(from, to, length.out = nmf)
  range.mf2=c(range.mf[1]-diff(range.mf)[1]/2,range.mf+diff(range.mf)[1]/2)

  if(!is.null(bws)&(!is.null(cbw1)|!is.null(cbw2))){warning("Arguments 'cbw1' and 'cbw2' are not needed when 'bws' are given")}

  if(!is.null(cbw1)&!is.null(cbw2)&is.null(bws)){

    if (!is.numeric(cbw1)){
      warning("Argument 'cbw1' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if ((length(cbw1)!=1)|(cbw1%%1!=0) | (cbw1<=0)){
      warning("Argument 'cbw1' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if (!is.numeric(cbw2)){
      warning("Argument 'cbw2' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }
    if ((length(cbw2)!=1)|(cbw2%%1!=0) | (cbw2<=0)){
      warning("Argument 'cbw2' must be a positive integer number. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }

  }

  if((!is.null(cbw1)&is.null(cbw2)&is.null(bws))|(is.null(cbw1)&!is.null(cbw2)&is.null(bws))){
    warning("Both arguments 'cbw1' and 'cbw2' must be specified. Default values of 'bws' were used")
    cbw1=NULL;cbw2=NULL
  }

  if(!is.null(cbw1)&!is.null(cbw2)&is.null(bws)){

    if(cbw1>=cbw2){
      warning("Argument 'cbw2' must be greater than 'cbw1'. Default values of 'bws' were used")
      cbw1=NULL;cbw2=NULL
    }else{
      hmin=bw.crit(data,cbw2)
      hmax=bw.crit(data,cbw1)
      range.h=seq(hmin,hmax,len=gridsize[2])
      if(logbw.regulargrid==T){
        range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
        range.h=10^range.h
      }
    }
  }


  if(!is.null(bws)){

    if (!is.numeric(bws)){
      warning("Argument 'bws' must be numeric. Default values of 'bws' were used")
      bws=NULL
    }else{
      if(length(bws)<2){
        warning("Argument 'bws' must be of length greater than one. Default values of 'bws' were used")
        bws=NULL
      }
      if(length(bws)==2){
        hmax=bws[2]
        hmin=bws[1]
        warning("A grid of 'bws' between the given ones were used")
        range.h=seq(hmin,hmax,len=gridsize[2])
        if(logbw.regulargrid==T){
          range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
          range.h=10^range.h
        }
      }
      if(length(bws)>2){
        hmax=max(bws)
        hmin=min(bws)
        range.h=bws
        if(length(bws)!=gridsize[2]){
          warning(paste("The 'gridsize' for the 'bws' is equal to",length(bws)))
        }
      }
    }
  }


  if(is.null(bws)&is.null(cbw2)&is.null(cbw1)){
    hmin=2*min(diff(range.x))
    hmax=diff(range.data)
    range.h=seq(hmin,hmax,len=gridsize[2])
    if(logbw.regulargrid==T){
      range.h=seq(log10(hmin),log10(hmax),len=gridsize[2])
      range.h=10^range.h
    }
  }

  if (!is.numeric(B)){
    warning("Argument 'B' must be a positive integer number. Default value of 'B' was used")
    B=99
  }
  if ((length(B)!=1)|(B%%1!=0) | (B<=0)){
    warning("Argument 'B' must be a positive integer number. Default value of 'B' was used")
    B=99
  }

  kernel="gaussian"
  weights=rep.int(1/nx, nx)
  totMass=nx/N

  modehat=matrix(0,gridsize[1],gridsize[2])
  cfkords2=list()

  for(ibw in 1:length(range.h)){
    bw=range.h[ibw]
    lo=from - 4 * bw
    up=to + 4 * bw
    y=.Call(BinDist2, data, weights, lo, up, n) * totMass
    kords2=seq.int(0, 2 * (up - lo), length.out = 2L * n)
    kords2[(n + 2):(2 * n)]=-kords2[n:2]
    kords2=switch(kernel, gaussian = dnorm(kords2, sd = bw))
    cfkords2[[ibw]]=Conj(fft(kords2))
    kords=fft(fft(y) * cfkords2[[ibw]], inverse = TRUE)
    kords=Re(kords)[1L:n]/length(y)
    xords=seq.int(lo, up, length.out = n)
    y = approx(xords, kords, range.x)$y
    z=c(1:(n-1))
    re=z[diff(y)>0]
    z2=c(1:length(re))
    se=z2[diff(re)>1]
    posic=re[se]
    if(re[length(re)]<(n-1)){posic=c(posic,re[length(re)])}
    temp=range.x[posic]
    for(k in 1:length(posic)){
      if(k==1){
        modehat[(which(temp[k]<range.mf2)[1]-1),ibw]=1+modehat[(which(temp[k]<range.mf2)[1]-1),ibw]
      }else{
        if(which(temp[k]<range.mf2)[1]!=which(temp[k-1]<range.mf2)[1]){
          modehat[(which(temp[k]<range.mf2)[1]-1),ibw]=1+modehat[(which(temp[k]<range.mf2)[1]-1),ibw]
        }
      }
    }
  }

  #bootstrap

    for(boot in 1:B){

      databoot=sample(data,replace=T)

      fhatboot=matrix(0,gridsize[1],gridsize[2])

      for(ibw in 1:length(range.h)){
        bw=range.h[ibw]
        lo=from - 4 * bw
        up=to + 4 * bw
        y=.Call(BinDist2, databoot, weights, lo, up, n) * totMass
        kords=fft(fft(y) * cfkords2[[ibw]], inverse = TRUE)
        kords=Re(kords)[1L:n]/length(y)
        xords=seq.int(lo, up, length.out = n)
        y = approx(xords, kords, range.x)$y
        re=z[diff(y)>0]
        z2=c(1:length(re))
        se=z2[diff(re)>1]
        posic=re[se]
        if(re[length(re)]<(n-1)){posic=c(posic,re[length(re)])}
        temp=range.x[posic]
        for(k in 1:length(posic)){
          if(k==1){
            modehat[(which(temp[k]<range.mf2)[1]-1),ibw]=1+modehat[(which(temp[k]<range.mf2)[1]-1),ibw]
          }else{
            if(which(temp[k]<range.mf2)[1]!=which(temp[k-1]<range.mf2)[1]){
              modehat[(which(temp[k]<range.mf2)[1]-1),ibw]=1+modehat[(which(temp[k]<range.mf2)[1]-1),ibw]
            }
          }
        }
      }

    }



  modehat=modehat/(B+1)



  if(logbw==T){
    range.h=log10(range.h)
  }
  colnames(modehat)=formatC(range.h)
  row.names(modehat)=formatC(range.mf)
  callname=match.call()
  modehatret=list(modehat,range.mf,range.h,logbw,nx,callname)
  names(modehatret)=c("modeforest","range.x","range.bws","logbw","sample.size","call")
  class(modehatret)="gtmod"
  if(display==T){plot(modehatret,...)}
  return(modehatret)



}








