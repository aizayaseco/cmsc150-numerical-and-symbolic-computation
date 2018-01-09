"
Author: Aizaya L. Seco
Date: October 14, 2017
Description: program for Newton's Divided Difference using augmented coefficient matrix
"

NDD <-function(mat)
{
  #sort via x
  matCopy <- mat
  matCopy <- matCopy[order(mat[,1]),]
  r <- nrow(matCopy)
  c <- ncol(matCopy)
  iterations = r-1 #maximum iterations
  
  # creating of new matrix given maximum iterations
  mat2 <- matrix(, nrow = r, ncol = c+iterations)
  colnames(mat2) <- c("x",1:((c+iterations)-1))
  rownames(mat2) <- 1:r
  mat2[1:r, 1:c]<- matCopy
  
  #NDD mechanism application
  n=0
  for(i in 3:(r+1)){
    for(j in 1:(r-(n+1))){
      mat2[j,i] <- (mat2[j+1,i-1]-mat2[j,i-1])/(mat2[(j+1+n),'x']-mat2[j,'x'])
    }
    n=n+1
  }
  
  #getting of as
  as= as.vector(mat2[1,(1:r+1)])
  
  l= length(as)
  
  #getting f
  func= ""
  for(i in 1:l){
    #update(func,as[i])
    func= paste(func, as[i], sep=" ")
    for(j in 1:i){
     if(j==i) break;
     #temp= ""
     #temp=paste(func, "(x-", mat2[j,'x'],")", sep="")
     func= paste(func, " * (x-", mat2[j,'x'],")", sep="");
     #update(func, (+ x-mat2[j,'x']))
    }
    if(i==l) break;
    func=paste(func, "+ ", sep=" ")
    
  }
  result <- list(table=mat2, as=as, f=force((parse(text=func))))
  return(result)
}


"Commented out code are for graphing and for comparison"

#F1 <- function(x) (768.832372457) -(27.696213846)*x  + (0.335236551)*x^2 -(0.001333299)*x^3;
#F2 <- function(x) 15.4 +  -2.33333333333336 * (x-69.4) +  2.87037037037043 * (x-69.4) * (x-69.7) +  -1.49875952507535 * (x-69.4) * (x-69.7) * (x-70.6) +  0.268798349363337 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) +  -0.040618012859624 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) +  0.00417479695662574 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) +  -0.000397535849995408 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) +  3.12268684530735e-05 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) +  -1.19741434405002e-07 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) * (x-81.6) +  -1.05126139772842e-06 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) * (x-81.6) * (x-82.6) +  7.02987725909163e-07 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) * (x-81.6) * (x-82.6) * (x-83.3) +  -3.58298924262608e-07 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) * (x-81.6) * (x-82.6) * (x-83.3) * (x-83.5) +  5.85135308625051e-08 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) * (x-81.6) * (x-82.6) * (x-83.3) * (x-83.5) * (x-84.3) +  -5.44581322031705e-09 * (x-69.4) * (x-69.7) * (x-70.6) * (x-71.6) * (x-75.2) * (x-76.3) * (x-79.6) * (x-80.6) * (x-81.6) * (x-82.6) * (x-83.3) * (x-83.5) * (x-84.3) * (x-88.6);
#plot(mat, col="red", xlab="Temperature (℉)",ylab="Chirps/Second")
#curve(F2, from=0, to=100, add=TRUE)
#curve(F1, from=0, to=100, add=TRUE, col="blue")