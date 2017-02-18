source("./stat.R")
library(ggplot2)
# beta=1
# difference simulation
#pSeq=c(100,200,300,400)
pSeq=seq(20,1000,by=40)
myDif=NULL
for(p in pSeq){
    n=p
    n1=n
    n2=n
    r=3
    beta=1
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:1){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
        myEigen=eigen(S,symmetric = TRUE)
        myRhat=which.max(myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)])
        myRhat=min(myRhat,rmax)
        myTildeV=myEigen$vectors[,(1:myRhat)]
        myCT[i]=(chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2)-chenStat(X1%*%(temp[[1]]),X2%*%(temp[[1]]),n1,n2))*
            n1*n2/(n1+n2)/sqrt(2*p)
    }
    # empirical variance
    myDif=c(myDif,median(abs(myCT)))
}

jpeg("difference1.jpeg",width=300,height=300)
ggplot(data=data.frame(p=pSeq,dif=myDif),aes(x=p,y=dif))+
    geom_line()+
    xlab("p")+
    #ylab(expression(paste("|",T[dif],"|")))+
    ylab("Absolute Error")+
    ylim(0,NA)+
    ggtitle(expression(paste("p=n, r=3, ",beta,"=1")))+
    theme_linedraw()
dev.off()


nSeq=seq(10,40,by=1)
#pSeq=c(100,225,400,625,900)
pSeq=nSeq^2
myDif=NULL
for(p in pSeq){
    n=sqrt(p)
    n1=n
    n2=n
    r=3
    beta=1
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:1){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
        myEigen=eigen(S,symmetric = TRUE)
        myRhat=which.max(myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)])
        myRhat=min(myRhat,rmax)
        myTildeV=myEigen$vectors[,(1:myRhat)]
        myCT[i]=(chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2)-chenStat(X1%*%(temp[[1]]),X2%*%(temp[[1]]),n1,n2))*
            n1*n2/(n1+n2)/sqrt(2*p)
    }
    # empirical variance
    myDif=c(myDif,mean(abs(myCT)))
}

jpeg("difference2.jpeg",width=300,height=300)#,res=100)
ggplot(data=data.frame(p=pSeq,dif=myDif),aes(x=p,y=dif))+
    geom_line()+
    xlab("p")+
    #ylab(expression(paste("|",T[dif],"|")))+
    ylab("Absolute Error")+
    xlim(0,NA)+ylim(0,NA)+
    ggtitle(expression(paste("p=",n^2,", r=3, ",beta,"=1")))+
    theme_linedraw()
dev.off()



# beta=0.5
# difference simulation
#pSeq=c(100,200,300,400)
pSeq=seq(20,1000,by=40)
myDif=NULL
for(p in pSeq){
    n=p
    n1=n
    n2=n
    r=3
    beta=0.5
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:1){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
        myEigen=eigen(S,symmetric = TRUE)
        myRhat=which.max(myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)])
        myRhat=min(myRhat,rmax)
        myTildeV=myEigen$vectors[,(1:myRhat)]
        myCT[i]=(chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2)-chenStat(X1%*%(temp[[1]]),X2%*%(temp[[1]]),n1,n2))*
            n1*n2/(n1+n2)/sqrt(2*p)
    }
    # empirical variance
    myDif=c(myDif,median(abs(myCT)))
}

jpeg("difference3.jpeg",width=300,height=300)
ggplot(data=data.frame(p=pSeq,dif=myDif),aes(x=p,y=dif))+
    geom_line()+
    xlab("p")+
    #ylab(expression(paste("|",T[dif],"|")))+
    ylab("Absolute Error")+
    ylim(0,NA)+
    ggtitle(expression(paste("p=n, r=3, ",beta,"=0.5")))+
    theme_linedraw()
dev.off()


nSeq=seq(10,40,by=1)
#pSeq=c(100,225,400,625,900)
pSeq=nSeq^2
myDif=NULL
for(p in pSeq){
    n=sqrt(p)
    n1=n
    n2=n
    r=3
    beta=0.5
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:1){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
        myEigen=eigen(S,symmetric = TRUE)
        myRhat=which.max(myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)])
        myRhat=min(myRhat,rmax)
        myTildeV=myEigen$vectors[,(1:myRhat)]
        myCT[i]=(chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2)-chenStat(X1%*%(temp[[1]]),X2%*%(temp[[1]]),n1,n2))*
            n1*n2/(n1+n2)/sqrt(2*p)
    }
    # empirical variance
    myDif=c(myDif,mean(abs(myCT)))
}

jpeg("difference4.jpeg",width=300,height=300)#,res=100)
ggplot(data=data.frame(p=pSeq,dif=myDif),aes(x=p,y=dif))+
    geom_line()+
    xlab("p")+
    #ylab(expression(paste("|",T[dif],"|")))+
    ylab("Absolute Error")+
    xlim(0,NA)+ylim(0,NA)+
    ggtitle(expression(paste("p=",n^2,", r=3, ",beta,"=0.5")))+
    theme_linedraw()
dev.off()







# beta=2
# difference simulation
#pSeq=c(100,200,300,400)
pSeq=seq(20,1000,by=40)
myDif=NULL
for(p in pSeq){
    n=p
    n1=n
    n2=n
    r=3
    beta=2
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:1){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
        myEigen=eigen(S,symmetric = TRUE)
        myRhat=which.max(myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)])
        myRhat=min(myRhat,rmax)
        myTildeV=myEigen$vectors[,(1:myRhat)]
        myCT[i]=(chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2)-chenStat(X1%*%(temp[[1]]),X2%*%(temp[[1]]),n1,n2))*
            n1*n2/(n1+n2)/sqrt(2*p)
    }
    # empirical variance
    myDif=c(myDif,median(abs(myCT)))
}

jpeg("difference5.jpeg",width=300,height=300)
ggplot(data=data.frame(p=pSeq,dif=myDif),aes(x=p,y=dif))+
    geom_line()+
    xlab("p")+
    #ylab(expression(paste("|",T[dif],"|")))+
    ylab("Absolute Error")+
    ylim(0,NA)+
    ggtitle(expression(paste("p=n, r=3, ",beta,"=2")))+
    theme_linedraw()
dev.off()


nSeq=seq(10,40,by=1)
#pSeq=c(100,225,400,625,900)
pSeq=nSeq^2
myDif=NULL
for(p in pSeq){
    n=sqrt(p)
    n1=n
    n2=n
    r=3
    beta=2
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:1){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
        myEigen=eigen(S,symmetric = TRUE)
        myRhat=which.max(myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)])
        myRhat=min(myRhat,rmax)
        myTildeV=myEigen$vectors[,(1:myRhat)]
        myCT[i]=(chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2)-chenStat(X1%*%(temp[[1]]),X2%*%(temp[[1]]),n1,n2))*
            n1*n2/(n1+n2)/sqrt(2*p)
    }
    # empirical variance
    myDif=c(myDif,mean(abs(myCT)))
}

jpeg("difference6.jpeg",width=300,height=300)#,res=100)
ggplot(data=data.frame(p=pSeq,dif=myDif),aes(x=p,y=dif))+
    geom_line()+
    xlab("p")+
    #ylab(expression(paste("|",T[dif],"|")))+
    ylab("Absolute Error")+
    xlim(0,NA)+ylim(0,NA)+
    ggtitle(expression(paste("p=",n^2,", r=3, ",beta,"=2")))+
    theme_linedraw()
dev.off()