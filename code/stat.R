library(ggplot2)
library(MASS)


newModelGenerator <- function(eigValue){
    # generate V
    V = svd(matrix(rnorm(p * p, 0, 1), c(p, p)))$v
    normalModelSimulator <- function(n){
        t(
            V%*%diag(eigValue)%*%matrix(rnorm(p*n),c(p,n))
        )   
    }
    return(
        list(
            V=V,
            normalModelSimulator=normalModelSimulator
        )
    )
}


chenStat=function(X1,X2,n1,n2){
    T1=sum(colMeans(X1)^2)*n1/(n1-1)-sum(X1^2)/n1/(n1-1)
    T2=sum(colMeans(X2)^2)*n2/(n2-1)-sum(X2^2)/n2/(n2-1)
    T3=sum(colMeans(X1)*colMeans(X2))
    return(T1+T2-2*T3)
}
sriStat <- function(X1,X2,n1,n2){
    S1 <- var(X1)
    S2 <- var(X2)
    S <- ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
    sum((colMeans(X1)-colMeans(X2))^2/diag(S))
}

myStat=function(X1,X2,n1,n2,rmax=10){
    S1 <- var(X1)
    S2 <- var(X2)
    S=((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
    myEigen=eigen(S,symmetric = TRUE)
    theTemp=myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)]
    myRhat=which.max(theTemp[1:rmax])
    myTildeV=myEigen$vectors[,-(1:myRhat)]
    
    trace1=sum(eigen(S1,symmetric = TRUE)$values[(myRhat+1):p])
    trace2=sum(eigen(S2,symmetric = TRUE)$values[(myRhat+1):p])
     
    # variance estimator
    sigmaSqEst <- mean(myEigen$values[(myRhat+1):p])
    
                     
    stat <- sum((t(myTildeV)%*%(colMeans(X1)- colMeans(X2)))^2)-trace1/n1-trace2/n2
    
    tau <- 1/n1+1/n2
    studentStat <- stat/sigmaSqEst/sqrt(2*tau^2*p)
    
    return(
        list(
            stat=stat,
            studentStat=studentStat
        )
    )
}
# do test
doTest = function(X1, X2, n1, n2, p, rmax = 10) {
    my = myStat(X1, X2, n1, n2, rmax)
    return(pnorm(my$studentStat,0,1,lower.tail = FALSE))
}

doUneqTest <- function(X1, X2, n1, n2, p, rmax = 10) {
    S1 <- var(X1)
    S2 <- var(X2)
    myEigen=eigen(S1,symmetric = TRUE)
    theTemp=myEigen$values[1:(n1-2)]/myEigen$values[2:(n1-1)]
    myRhat1=which.max(theTemp[1:rmax])
    myV1=myEigen$vectors[,(1:myRhat1)]
    
    myEigen=eigen(S2,symmetric = TRUE)
    theTemp=myEigen$values[1:(n2-2)]/myEigen$values[2:(n2-1)]
    myRhat2=which.max(theTemp[1:rmax])
    myV2 <- myEigen$vectors[,(1:myRhat2)]
    
    myV <- svd(cbind(myV1,myV2))$u
    myTildeV <- Null(myV)
    
    
    trace1=sum(eigen(S1,symmetric = TRUE)$values[(myRhat1+1):p])
    trace2=sum(eigen(S2,symmetric = TRUE)$values[(myRhat2+1):p])
     
    # variance estimator
    sigmaSqEst <- 2*(p-myRhat1-myRhat2)*(
        (trace1/(p-myRhat1))^2/n1/(n1-1)+
            (trace2/(p-myRhat2))^2/n2/(n2-1)+
            2*(trace1/(p-myRhat1))*(trace2/(p-myRhat2))/n1/n2
    )
    
                     
    stat <- sum((t(myTildeV)%*%(colMeans(X1)- colMeans(X2)))^2)-trace1/n1-trace2/n2
    
    studentStat <- stat/sqrt(sigmaSqEst)
    
    return(pnorm(studentStat,0,1,lower.tail = FALSE))
}

# n1 <- 30
# n2 <- 30
# p <- 100
# r <- 3
# beta <- 1
# theEig <- rep(1,p)
# theEig[1:r] <- rep(p^beta,r)
# normalModelSimulator <- newModelGenerator(theEig)
# X1 <- normalModelSimulator(n1)
# X2 <- normalModelSimulator(n2)
# myStat(X1,X2,n1,n2)