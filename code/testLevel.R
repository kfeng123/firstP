source("./stat.R")
library(ggplot2)
library(xtable)

simulateLevel=function(n1,n2,p,r,beta,B=1000,rmax=10){
    temp=modelGenerator(p=p,r=r,beta=beta) 
    V=temp$V
    D=temp$D
    modelSimulator=temp$modelSimulator
    
    myPvalue=NULL
    oraclePvalue=NULL
    for(i in 1:B){
        X1=modelSimulator(n1)
        X2=modelSimulator(n2)
        myPvalue[i]=doTest(X1,X2,n1,n2,p,rmax=10)
        # oracle: V known
        temp=chenStat(X1%*%Null(V),X2%*%Null(V),n1,n2)
        tempOracleStat=n1*n2*temp/(sqrt(2*p)*(n1+n2)*1)
        oraclePvalue[i]=pnorm(tempOracleStat,0,1,lower.tail = FALSE)
    }
    myLevel=sum(myPvalue<0.05)/B
    oracleLevel=sum(oraclePvalue<0.05)/B
    list(myLevel=myLevel,oracleLevel=oracleLevel)
} 

Out=data.frame(n1=0,n2=0,p=0,r=0,beta=0,myLevel=0,oracleLevel=0)

for(beta in c(0.5,1,2))
for(r in 2)
for(p in c(200,400,600,800))
for(n in c(300,600)){
    n1=n
    n2=n
    level=simulateLevel(n1,n2,p,r,beta)
    temp=data.frame(n1=n1,n2=n2,p=p,r=r,beta=beta,myLevel=level$myLevel,oracleLevel=level$oracleLevel)
    Out=rbind(Out,temp)
}
Out=Out[-1,]
row.names(Out)=NULL
write.csv(Out,"level.csv",row.names = FALSE)


Out=Out[,-2]
Temp1=Out[Out$beta==0.5,]
Temp2=Out[Out$beta==1,]
Temp3=Out[Out$beta==2,]

TTT1=merge(Temp1,Temp2,by=c("n1","p","r"))
TTT2=merge(TTT1,Temp3,by=c("n1","p","r"))
TTT2=TTT2[,-c(3,4,7,10)]
names(TTT2)=c("n","p","myLevelbeta0.5","oracleLevelbeta0.5","myLevelbeta1","oracleLevelbeta1","myLevelbeta2","oracleLevelbeta2")

TTT2=TTT2[order(TTT2[,1],TTT2[,2]),]


myTable1=xtable(TTT2,digits=c(0,0,0,3,3,3,3,3,3),caption="Test level simulation",label="biaoge1")
align(myTable1) <- "rrccccccc"
print(myTable1,file="level.tex",include.rownames=FALSE)


# myTable1=xtable(Out[1:15,],digits=c(0,0,0,0,0,0,3),caption="Test level simulation",label="biaoge1")
# align(myTable1) <- "|r|rrrrr|r|"
# print(myTable1,file="level1.tex",include.rownames=FALSE)
