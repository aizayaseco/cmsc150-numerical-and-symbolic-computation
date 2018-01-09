"
Author: Aizaya L. Seco
Date: October 17, 2017
Description: solves unkown independent variables within the given range in a problem using Neville and Lagrange. 
"

options(scipen=999)
Neville <-function(mat,verbose=TRUE,x)
{
  #default of true for verbose

  matCopy <- mat
  matCopy<- (matCopy[order(abs(x-matCopy[,1])),]) #sort via nearest to x
  
  r <- nrow(matCopy)
  c <- ncol(matCopy)
  iterations = r-1 #maximum P i,n
  
  # creating of new matrix given maximum iterations
  mat2 <- matrix(, nrow = r, ncol = c+iterations+1)
  colnames(mat2) <- c("i","x",1:((c+iterations)-1))
  rownames(mat2) <- 1:r
  mat2[ ,1]<- 1:r
  mat2[1:r, 2:(c+1)]<- matCopy
  
  #Neville mechanism application
  n=0
  for(i in 4:ncol(mat2)){
    for(j in 1:(r-(n+1))){
      mat2[j,i] <- ((x-mat2[j,'x'])*mat2[j+1,i-1] +(mat2[(j+1+n),'x']-x)*mat2[j,i-1])/(mat2[(j+1+n),'x']-mat2[j,'x'])
    }
    n=n+1
  }
  
  #table printing
  if(verbose==TRUE){
    print(mat2)  
  }
  result <- list(fx=mat2[1,ncol(mat2)])
  return(result)
}

Lagrange<- function(mat){
  #sort via x
  matCopy <- mat
  matCopy <- matCopy[order(mat[,1]),]
  r <- nrow(matCopy)
  c <- ncol(matCopy)
  
  #function concatenation
  func= "function(x)"
  for(i in 1:r){
    func= paste(func,"(", matCopy[i,2], sep=" ")
    for(j in 1:r){
      if(i!=j){
        func= paste(func, " * (x-",matCopy[j,1],")/(",matCopy[i,1],"-",matCopy[j,1],")", sep="")
      }
    }
    func=paste(func, ")", sep=" ")
    if(i!=r)
      func=paste(func, "+", sep=" ")
  }
  
  result <- list(f=eval(parse(text=func)))
  return(result)
}

"Commented Section for plotting"
#mat = matrix( 
#  c(1990, 61229622, 1995, 68349452, 2000, 75505061, 2005, 82079348, 2010, 87940171, 2015,93440274, 2020, 99417214),  " the data elements " 
#  ncol=2,              " number of columns " 
#  byrow = TRUE
#)

#C<- Lagrange(mat)
#C
#F1<-C$f
#plot(2004,80824322,xlim=range(c(1990,2020)),ylim=range(61229622,99417214), col="red",type="p", xlab="Year", ylab="Population Count", main="Population Count from Year 1990 to 2020",pch=20)
#par(new=TRUE)
#plot(mat,xlim=range(c(1990,2020)),ylim=range(61229622,99417214),col="blue",log="x", pch=3,xlab="", ylab="")
#curve(F1, from=1990, to=2020, add=TRUE)