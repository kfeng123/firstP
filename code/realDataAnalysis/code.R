source("../stat.R")
data=read.table("data")
label=read.table("label")
data=t(data)
data=log(data)
X1=data[label>0,]
X2=data[label<0,]
row.names(X1)=NULL
row.names(X2)=NULL

# do permutation
myPermutation = function(X1, X2, n1, n2, rmax = 10, B) {
    chen = chenStat(X1, X2, n1, n2)
    my = myStat(X1, X2, n1, n2, rmax)
    chenTemp = NULL
    myTemp = NULL
    for (i in 1:B) {
        myOrder = sample(n1 + n2, size = n1 + n2, replace = FALSE)
        tempX1 = rbind(X1,X2)[myOrder[1:n1], ]
        tempX2 = rbind(X1,X2)[myOrder[(n1+1):(n1+n2)], ]
        chenTemp[i] = chenStat(tempX1, tempX2, n1, n2)
        myTemp[i] = myStat(tempX1, tempX2, n1, n2)
    }
    return(list(
        my = sum(myTemp > my) / length(myTemp),
        chen = sum(chenTemp > chen) / length(chenTemp)
    ))
}
