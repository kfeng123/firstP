# one sample
library(xtable)
library(ggplot2)
library(tikzDevice)
library(highD2pop)
myPvalue=function(X,B){
    myStat=function(X){
        myMean=colMeans(X)
        myEigen=eigen(var(X),symmetric=TRUE)
        myRhat=which.max(myEigen$values[1:(n-2)]/myEigen$values[2:(n-1)])
        # restrict myRhat
        myRhat=min(myRhat,10)
        myTildeV=myEigen$vectors[,-(1:myRhat)]
        mytemp=X%*%myTildeV
        chenStat(mytemp)
    }
    chenStat=function(X){
        myMean=colMeans(X)
        (n^2*sum(myMean^2)-sum(apply(X,1,function(x){sum(x^2)})))/(n*(n-1))
    }
    fastStat=function(X){
        myMean=colMeans(X)
        myEigen=eigen(var(X),symmetric=TRUE)
        myRhat=which.max(myEigen$values[1:(n-2)]/myEigen$values[2:(n-1)])
        # restrict myRhat
        myRhat=min(myRhat,10)
        n*sum(myMean^2)-sum(myEigen$values)-sum(myEigen$values[1:myRhat])
    }
    sriStat=function(X){
        myMean=colMeans(X)
        n*sum(myMean^2/diag(var(X)))
    }
    theNewTemp=rep(0,B)
    theChenTemp=rep(0,B)
    theFastTemp=rep(0,B)
    theSriTemp=rep(0,B)
    
    theNewStat=myStat(X)
    theChenStat=chenStat(X)
    theFastStat=fastStat(X)
    theSriStat=sriStat(X)
    for(i in 1:B){
        tempSample=as.logical(rbinom(n,1,0.5))
        temp=rbind(X[tempSample,],-X[!tempSample,])
        theNewTemp[i]=myStat(temp)
        theChenTemp[i]=chenStat(temp)
        theFastTemp[i]=fastStat(temp)
        theSriTemp[i]=sriStat(temp)
    }
    list(new=sum(theNewStat<theNewTemp)/B,
         chen=sum(theChenStat<theChenTemp)/B,
         fast=sum(theFastStat<theFastTemp)/B,
         sri=sum(theSriStat<theSriTemp)/B)
}

mySimulation=function(n,p,r,alpha){
    c=1
    sigmaSquare=1
    if(r==0){
        r=1
        D=diag(0,nrow=1,ncol=1)
    }else{
        D=diag(sqrt(rep(c*p^alpha,r)+runif(r,0,1)),nrow=r,ncol=r)
    }
    V=svd(matrix(rnorm(r*p,0,1),c(r,p)))$v
    
    mu0=rnorm(p,0,1)
    mu0=mu0/sqrt(sum(mu0^2))
    mu0=mu0/n*sqrt(p)
    
    myOut=NULL
    chenOut=NULL
    fastOut=NULL
    sriOut=NULL
    for(cc in seq(0,6,1)){
        mu=mu0*cc
        hh=100
        myResult=rep(0,hh)
        chenResult=rep(0,hh)
        fastResult=rep(0,hh)
        sriResult=rep(0,hh)
        for(i in 1:hh){
            U=rnorm(n*r,0,1)
            dim(U)=c(n,r)
            Z=rnorm(n*p,0,sigmaSquare)
            dim(Z)=c(n,p)
            X=outer(rep(1,n),mu)+U%*%D%*%t(V)+Z
            jj=myPvalue(X,B=100)
            myResult[i]=jj$new
            chenResult[i]=jj$chen
            fastResult[i]=jj$fast
            sriResult[i]=jj$sri
        }
        myOut=c(myOut,sum(myResult<0.05)/hh)
        chenOut=c(chenOut,sum(chenResult<0.05)/hh)
        fastOut=c(fastOut,sum(fastResult<0.05)/hh)
        sriOut=c(sriOut,sum(sriResult<0.05)/hh)
    }
    list(my=myOut,chen=chenOut,fast=fastOut,sri=sriOut)
}


n=50
p=100
r=1
alpha=1
jjj1=mySimulation(n,p,r,alpha)


n=50
p=100
r=0
alpha=1
jjj2=mySimulation(n,p,r,alpha)

n=50
p=100
r=1
alpha=0.5
jjj3=mySimulation(n,p,r,alpha)


n=50
p=100
r=5
alpha=2
jjj4=mySimulation(n,p,r,alpha)

n=100
p=150
r=1
alpha=1
jjj5=mySimulation(n,p,r,alpha)

n=100
p=150
r=1
alpha=0.5
jjj6=mySimulation(n,p,r,alpha)

myPlot=function(uio){
    cc=seq(0,6,1)
    temp1=data.frame(h=cc,Power=uio$my,Method="New")
    temp2=data.frame(h=cc,Power=uio$chen,Method="CQ")
    temp3=data.frame(h=cc,Power=uio$fast,Method="FAST")
    temp4=data.frame(h=cc,Power=uio$sri,Method="S")
    myD=rbind(temp1,temp2,temp3,temp4)
    ggplot(data=myD,aes(x=h,y=Power,color=Method,linetype=Method))+
        geom_line()+xlab(expression( paste("||",mu,"||"^2)))
}
jpeg("fig1.jpeg",width=300,height=300)
myPlot(jjj1)+ggtitle(expression(paste("n=50, p=100, r=1, ",alpha,"=1")))+theme_linedraw()
dev.off()
jpeg("fig2.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste("n=50, p=100, r=0, ",alpha,"=1")))+theme_linedraw()
dev.off()
jpeg("fig3.jpeg",width=300,height=300)
myPlot(jjj3)+ggtitle(expression(paste("n=50, p=100, r=1, ",alpha,"=0.5")))+theme_linedraw()
dev.off()
jpeg("fig4.jpeg",width=300,height=300)
myPlot(jjj4)+ggtitle(expression(paste("n=50, p=100, r=5, ",alpha,"=2")))+theme_linedraw()
dev.off()
jpeg("fig5.jpeg",width=300,height=300)
myPlot(jjj5)+ggtitle(expression(paste("n=100, p=150, r=1, ",alpha,"=1")))+theme_linedraw()
dev.off()
jpeg("fig6.jpeg",width=300,height=300)
myPlot(jjj6)+ggtitle(expression(paste("n=100, p=150, r=1, ",alpha,"=0.5")))+theme_linedraw()
dev.off()
