"
Author: Aizaya L. Seco
Date: September 19, 2017
Description: Using augmented coefficient matrix, it solves for variables using Gauss Jordan or Gaussian Elimination.
"

AugCoeffMatrix<-function(system)
{
  result <- list()
  
  #gets the list of variables
  fxn <- deparse(system[[1]])[1]
  fxn1 <- unlist(strsplit(fxn, split=" "))
  len <- length(fxn1)
  fxn2 <- (fxn1)[2:length(fxn1)]
  variables <- gsub(pattern = "[[:punct:]]", replacement = "", fxn2)
  
  #gets the values based on the gathered variables and list of systems
  n <- length(system)
  m <- length(variables)
  
  #matrix initialization
  mat <- matrix(, nrow = n, ncol = m+1)
  colnames(mat) <- c(variables,"RHS")
  rownames(mat) <- 1:n
  
  
  for(i in 1:n){
    equation <- deparse(system[[i]])[2:length(deparse(system[[i]]))]
    values1 <- unlist(strsplit(equation, split="\\+"))
    values2 <- unlist(strsplit(values1, split="\\*"))
    values2 <- gsub(pattern = " ", replacement = "", values2)
    values2 <- values2[values2!=""]
    o <- length(values2)
    print(values2)
    for(j in 1:m){
      for(k in 1:o){
        if(variables[[j]]==values2[[k]]){
          mat[i,variables[[j]]] <- as.numeric(values2[[k-1]])
          values2 <- values2[-c(k,(k-1))]
          break;
        }
      }
    }
    
    #checker for unaccepted acm
    if(length(values2)==1){
      if(is.numeric(as.numeric(values2[[1]])) && !is.na(values2[[1]]) && !is.na(as.numeric(values2[[1]]))){
        mat[i,"RHS"] <- -1 * as.numeric(values2[[1]])
      }
      else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }
  result <- list(variables=variables, augcoeffmatrix=mat)
  return(result)
}

GaussJordan<- function(augcoeffmatrix)
{
  acm <- augcoeffmatrix
  rows <- nrow(acm) #3
  cols <- ncol(acm) #4
  #Pivot Element acm[n,n]
  #Pivot Row acm[n, ]
  #Pivot Column acm[ ,n ]
  n <- 1
  n1 <-1
  
  for(i in 1:rows){
    #pivot checker
    m<- abs(acm[n,n])
    for(k in n:rows){
      if(abs(acm[k,n]) > m){
        n1 <- k
        m <- abs(acm[k,n])
      }
    }
    if(n1!=n){
      temp <- acm[n1, ]
      acm[n1, ] <- acm[n, ]
      acm[n, ] <- temp    
      n1<- n
    }
    #normalization
    acm[n, ] <- (acm[n, ]/acm[n,n])
    for(j in 1:rows){
      if(j!=n){
        #m=acm[j,n]
        #m*PR
        temp <- (acm[j,n]*acm[n, ])
        acm[j, ] <- acm[j, ] - temp
      } 
    }
    n <- n+1
    n1 <- n
  }
  
  xs <- numeric(cols-1)
  for(i in 1:(cols-1)){
    xs[i] <- acm[i,"RHS"]
  }
  result1 <- list()
  result1[["xs"]]<- xs
  result1[["forelem"]]<- acm
  return(result1)
}

Gaussian<- function(augcoeffmatrix)
{
  acm <- augcoeffmatrix
  rows <- nrow(acm) #3
  cols <- ncol(acm) #4
  #Pivot Element acm[n,n]
  #Pivot Row acm[n, ]
  n <- 1
  n1 <-1
  for(i in 1:rows){
    m<- abs(acm[n,n])
    #pivot checker
    for(k in n:rows){
      if (abs(acm[k,n]) > m){
        n1 <- k
        m <- abs(acm[k,n])
      }
    }
    if(n1!=n){
      temp <- acm[n1, ]
      acm[n1, ] <- acm[n, ]
      acm[n, ] <- temp    
      n1<- n
      #print("performed pivot")
      #print(acm)
    }
    
    
    for(j in 1:rows){
      if(j>n){
        #multiplier
        #multiplier*Pivot Row acm[n, ]
        multiplier <- (acm[j,n]/acm[n,n])
        temp <- (multiplier*acm[n, ])
        acm[j, ] <- acm[j, ] - temp
      } 
    }
    n <- n+1
    n1 <- n
  }
  
  xs <- numeric(cols-1)
  #backward substition
  for(i in rows:1){
    xs[i]<- acm[i,cols]
    for(j in (cols-1):1){
      if(j==i){
        if(acm[i,j]!=0){
          xs[i]<-(xs[i]/acm[i,j])
        }
        else{
          return(NA)
        }
      }
      else{
        if(j>i){
          xs[i]<- (xs[i]-(acm[i,j]*xs[j]))
        }
      }
    }
  }
  
  result1 <- list()
  result1[["xs"]]<- xs
  result1[["forelem"]]<- acm
  return(result1)
}

E1 <- function (x0, x1, x2) 15 * x0 + 1190.2 * x1 + 95128.42 * x2 + -248.5;
E2 <- function (x0, x1, x2) 1190.2 * x0 + 95128.42 * x1 + 7658515.36 * x2 + -19857.72;
E3 <- function (x0, x1, x2) 95128.42 * x0 + 7658515.36 * x1 + 620995159.4 * x2 + -1598754.258;
system <- list(E1, E2, E3);
result1 = AugCoeffMatrix(system);
result1
GaussJordan(result1$augcoeffmatrix)