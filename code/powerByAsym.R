# two sample
source("./stat.R")
# main simulation
myMainSimulation=function(p,n1,n2,r=3,beta){
    
    theEig <- rep(1,p)
    theEig[1:r] <- rep(p^beta,r)
    temp <- newModelGenerator(theEig)
    normalModelSimulator <- temp$normalModelSimulator
    V <- temp$V
    
    tau=(n1+n2)/n1/n2
    theoryPower=function(mu1,mu2){
        sum((mu1-mu2)^2)/sqrt(2*tau^2*p)
    }
    
    
    temp1=rnorm(p,0,1)
    temp2=rnorm(p,0,1)
    outMy=NULL
    for(hh in c(0,1,2,3,4,5)){
        myC=sqrt(hh*sqrt(2*tau^2*p)/sum((temp1-temp2)^2))
        mu1=temp1*myC
        mu2=temp2*myC
        tempMy=NULL
        for(i in 1:300){
            X1 <- normalModelSimulator(n=n1)+outer(rep(1,n1),mu1)
            X2 <- normalModelSimulator(n=n2)+outer(rep(1,n2),mu2)
            temp=doTest(X1,X2,n1,n2,p)
            tempMy[i]=temp
        }
        outMy=c(outMy,sum(tempMy<0.05)/length(tempMy))
    }
    list(my=outMy)
}

myPlot=function(uio){
    cc=seq(0,5,1)
    temp1=data.frame(h=cc,Power=uio$my,Method="New")
    #temp2=data.frame(h=cc,Power=uio$chen,Method="CQ")
    #temp3=data.frame(h=cc,Power=uio$fast,Method="FAST")
    #temp4=data.frame(h=cc,Power=uio$sri,Method="S")
    #myD=rbind(temp1,temp2)#,temp3,temp4)
    myD=temp1
    ggplot(data=myD,aes(x=h,y=Power,color=Method,linetype=Method))+
        geom_line()+ylim(0,1)+
        #xlab(expression( paste("||",mu[1]-mu[2],"||"^2)))
        xlab("d")
}


p=400;n1=300;n2=310;r=3;beta=0.5
jjj1=myMainSimulation(p,n1,n2,r,beta)
jpeg("newfig1.jpeg",width=300,height=300)
myPlot(jjj1)+ggtitle(expression(paste(n[1],"=300, ",n[2],"=310, p=400, r=3, ",beta,"=0.5")))+theme_linedraw()
dev.off()

p=400;n1=300;n2=310;r=3;beta=1
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("newfig2.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=300, ",n[2],"=310, p=400, r=3, ",beta,"=1")))+theme_linedraw()
dev.off()

p=400;n1=300;n2=310;r=3;beta=2
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("newfig3.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=300, ",n[2],"=310, p=400, r=3, ",beta,"=2")))+theme_linedraw()
dev.off()

p=400;n1=300;n2=310;r=0;beta=2
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("newfig4.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=300, ",n[2],"=310, p=400, r=0")))+theme_linedraw()
dev.off()

p=1000;n1=500;n2=510;r=3;beta=1
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("newfig5.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=500, ",n[2],"=510, p=1000, r=3, ",beta,"=1")))+theme_linedraw()
dev.off()

p=1000;n1=500;n2=510;r=3;beta=2
jjj2=myMainSimulation(p=p,n1=n1,n2=n2,r=r,beta=beta)
jpeg("newfig6.jpeg",width=300,height=300)
myPlot(jjj2)+ggtitle(expression(paste(n[1],"=500, ",n[2],"=510, p=1000, r=3, ",beta,"=2")))+theme_linedraw()
dev.off()
