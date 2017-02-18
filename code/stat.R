library(ggplot2)
library(MASS)

# model generator
# that is, generate V,D
# return a function which can generate sample
modelGenerator = function(c = 1,
                          sigmaSquare = 1,
                          p,
                          r,
                          beta = 1)
{
    # generate D
    if (r == 0) {
        r = 1
        D = diag(0, nrow = 1, ncol = 1)
    } else{
        D = diag(sqrt(rep(c * p ^ beta, r) + runif(r, 0, 1)), nrow = r, ncol = r)
    }
    # generate V
    V = svd(matrix(rnorm(r * p, 0, 1), c(r, p)))$v
    # modelSimulator 
    list(V=V,D=D,
    modelSimulator=function(n){
        t(
        #outer(mu,rep(1,n))+
            V%*%D%*%matrix(rnorm(r*n),c(r,n))+
            matrix(rnorm(p*n),c(p,n))*sigmaSquare
        )
    })
}

# single sample
chenStat=function(X1,X2,n1,n2){
    T1=sum(colMeans(X1)^2)*n1/(n1-1)-sum(X1^2)/n1/(n1-1)
    T2=sum(colMeans(X2)^2)*n2/(n2-1)-sum(X2^2)/n2/(n2-1)
    T3=sum(colMeans(X1)*colMeans(X2))
    return(T1+T2-2*T3)
}
myStat=function(X1,X2,n1,n2,rmax=10){
    S=((n1-1)*var(X1)+(n2-1)*var(X2))/(n1+n2-2)
    myEigen=eigen(S,symmetric = TRUE)
    theTemp=myEigen$values[1:(n1+n2-3)]/myEigen$values[2:(n1+n2-2)]
    myRhat=which.max(theTemp[1:rmax])
    myTildeV=myEigen$vectors[,-(1:myRhat)]
    # variance estimator
    varEstimator1=mean(myEigen$values[(myRhat+1):p])
    varEstimator2=mean(myEigen$values[(2*myRhat+1):(p-2*myRhat)])
                       
    return(list(stat=chenStat(X1%*%myTildeV,X2%*%myTildeV,n1,n2),
                varEstimator1=varEstimator1,varEstimator2=varEstimator2))
}
# do test
doTest = function(X1, X2, n1, n2, p, rmax = 10) {
    my = myStat(X1, X2, n1, n2, rmax)
    myT1=n1*n2*my$stat/(sqrt(2*p)*(n1+n2)*my$varEstimator1)
    return(pnorm(myT1,0,1,lower.tail = FALSE))
}
