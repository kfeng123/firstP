library(ggplot2)
# do permutation
myPermutation = function(X1, X2, n1, n2, rmax = 10, B) {
    chen = chenStat(X1, X2, n1, n2)
    my = myStat(X1, X2, n1, n2, rmax)$stat
    chenTemp = NULL
    myTemp = NULL
    for (i in 1:B) {
        myOrder = sample(n1 + n2, size = n1 + n2, replace = FALSE)
        tempX1 = rbind(X1,X2)[myOrder[1:n1], ]
        tempX2 = rbind(X1,X2)[myOrder[(n1+1):(n1+n2)], ]
        chenTemp[i] = chenStat(tempX1, tempX2, n1, n2)
        myTemp[i] = myStat(tempX1, tempX2, n1, n2)$stat
    }
    return(list(
        my = sum(myTemp > my) / length(myTemp),
        chen = sum(chenTemp > chen) / length(chenTemp)
    ))
}

# main simulation
myMainSimulation=function(p,n1,n2,r3,beta){
    modelSimulator=modelGenerator(p=p,r=r,beta=beta)$modelSimulator
    temp1=rnorm(p,0,1)
    temp2=rnorm(p,0,1)
    tau=(n1+n2)/n1/n2
    theoryPower=function(mu1,mu2){
        sum((mu1-mu2)^2)/sqrt(2*tau^2*p)
    }
    
    outMy=NULL
    outChen=NULL
    for(hh in c(0,1,2,3,4,5)){
        myC=sqrt(hh*sqrt(2*tau^2*p)/sum((temp1-temp2)^2))
        mu1=temp1*myC
        mu2=temp2*myC
        tempMy=NULL
        tempChen=NULL
        for(i in 1:100){
            X1=modelSimulator(n=n1)+outer(rep(1,n1),mu1)
            X2=modelSimulator(n=n2)+outer(rep(1,n2),mu2)
            temp=myPermutation(X1,X2,n1,n2,B=50)
            tempMy[i]=temp$my
            tempChen[i]=temp$chen
        }
        outMy=c(outMy,sum(tempMy<0.05)/length(tempMy))
        outChen=c(outChen,sum(tempChen<0.05)/length(tempChen))
    }
    list(my=outMy,chen=outChen)
}

myPlot=function(uio){
    cc=seq(0,5,1)
    temp1=data.frame(h=cc,Power=uio$my,Method="New")
    temp2=data.frame(h=cc,Power=uio$chen,Method="CQ")
    #temp3=data.frame(h=cc,Power=uio$fast,Method="FAST")
    #temp4=data.frame(h=cc,Power=uio$sri,Method="S")
    myD=rbind(temp1,temp2)#,temp3,temp4)
    ggplot(data=myD,aes(x=h,y=Power,color=Method,linetype=Method))+
        geom_line()+
        #xlab(expression( paste("||",mu[1]-mu[2],"||"^2)))
        xlab("d")
}


p=50
n1=30
n2=35
r=3
beta=0.5
jjj1=myMainSimulation(p=100,n1=30,n2=30,r=3,beta=0.5)
jpeg("fig1.jpeg",width=300,height=300)
myPlot(jjj1)+ggtitle(expression(paste(n[1],"=30, ",n[2],"=35, p=50, r=3, ",beta,"=0.5")))+theme_linedraw()
dev.off()

p=50
n1=30
n2=35
r=3
beta=1
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("fig2.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=30, ",n[2],"=35, p=50, r=3, ",beta,"=1")))+theme_linedraw()
dev.off()

p=50
n1=30
n2=35
r=3
beta=2
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("fig3.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=30, ",n[2],"=35, p=50, r=3, ",beta,"=2")))+theme_linedraw()
dev.off()

p=50
n1=30
n2=35
r=0
beta=2
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("fig4.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=30, ",n[2],"=35, p=50, r=0")))+theme_linedraw()
dev.off()

p=150
n1=30
n2=35
r=3
beta=1
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("fig5.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=30, ",n[2],"=35, p=150, r=3, ",beta,"=1")))+theme_linedraw()
dev.off()

p=150
n1=30
n2=35
r=3
beta=2
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("fig6.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=30, ",n[2],"=35, p=150, r=3, ",beta,"=2")))+theme_linedraw()
dev.off()

# variance simulation
n=50
myRatio=NULL
pSeq=c(50,100,150,200,250,300)
for( p in pSeq){
    n1=n
    n2=n
    r=3
    beta=1
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp[[2]]
    for(i in 1:1000){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        myCT[i]=myStat(X1,X2,n1,n2)
    }
    # empirical variance
    myEmp=var(myCT)
    # asymptotic variance
    tau=(n1+n2)/n1/n2
    myAsy=2*tau^2*p
    myRatio=c(myRatio,myEmp/myAsy)
}
plot(pSeq,myRatio)

jpeg("varianceRatio.jpeg",width=300,height=300)
ggplot(data=data.frame(p=pSeq,R=myRatio),aes(x=p,y=R))+
    geom_line()+
    xlab("p")+ylab("R")+
    ggtitle(expression(paste(n,"=50, r=3, ",beta,"=1")))+
    theme_linedraw()
dev.off()




