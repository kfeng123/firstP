source("./stat.R")
# variance simulation
n=50
myRatio=NULL
pSeq=c(50,100,150)#,200,250,300)
for( p in pSeq){
    n1=n
    n2=n
    r=3
    beta=1
    myCT=NULL
    rmax=10
    temp=modelGenerator(p=p,r=r,beta=beta)
    modelSimulator=temp$modelSimulator
    for(i in 1:500){
        X1=modelSimulator(n=n1)
        X2=modelSimulator(n=n2)
        myCT[i]=myStat(X1,X2,n1,n2,rmax)$stat
    }
    # empirical variance
    myEmp=var(myCT)
    # asymptotic variance
    tau=(n1+n2)/n1/n2
    myAsy=2*tau^2*p
    myRatio=c(myRatio,myEmp/myAsy)
}
plot(pSeq,myRatio)

sum(pnorm(myCT/sqrt(myEmp),0,1,lower.tail = FALSE)<0.05)/length(myCT)


jpeg("varianceRatio.jpeg",width=300,height=300)
ggplot(data=data.frame(p=pSeq,R=myRatio),aes(x=p,y=R))+
    geom_line()+
    xlab("p")+ylab("R")+
    ggtitle(expression(paste(n,"=50, r=3, ",beta,"=1")))+
    theme_linedraw()
dev.off()
