p = 30

condition1 = function(A,B,C){
  f = 4*n*(B^2)-2*(p-1)*(A^2)+n*p*(C^2)+(-1-2*n)*2*(p-1)*A*B+(2-n*p)*(p-1)*A*C-n*p*(1+2*n)*B*C
  return(f)
} 

condition2 = function(A,B,C){
  f = n*p*C-2*(p-1)*A
  return(f)
}

condition3 = function(A,B,C){
  f = n*B-(((n+1)/2)-(2/(2*p+n*p-2)))*n*p*C-(n+2)*(p-1)*A
  return(f)
}

condition4 = function(A,B,C){
  f = (n*B-(p-1)*A)*(2*p-2+n*p)+n*p*C-2*(p-1)*A
  return(f)
}

condition5 = function(A,B,C){
  f = ((p-1)^2)*(A^2)+n*(1-(6/(p*(n+2))))*(B^2)+(n/(n+2))*(C^2)+(p-1)*(n+1+((2*n)/(p*(n+2))))*A*B-((2*(p-1)*(n+1))/(n+2))*A*C-((n*(p-2-n*p))/(p*(n+2)))*B*C
  return(f)
}

condition6 = function(A,B,C){
  f = 2*n*B+(n+1)*n*p*C-2*(p-1)*(n+2)*A
  return(f)
}


for (n in 30:200){
  for (A in 0.00001:100){
    for (B in 0.0001:100){
      for (C in 0.0001:100){
        # if ((2*B>n*p*C+2*(p-1)*A) & (n*p*C>2*(p-1)*A)){
        if ((B-C+(p-1)*A>0) & ((p-1)*(A+C)>B)){
          if ((condition1(A,B,C)>0 & condition2(A,B,C)>0) & (condition3(A,B,C)<0 & condition4(A,B,C)<0)){
          # if ((condition1(A,B,C)<0 & condition2(A,B,C)>0) & (condition4(A,B,C)>0 & condition5(A,B,C)>0)){
          # if ((condition1(A,B,C)<0 & condition2(A,B,C)<0) & (condition3(A,B,C)<0 & condition6(A,B,C)>0)){
          # if ((condition1(A,B,C)>0 & condition2(A,B,C)<0) & (condition6(A,B,C)>0 & condition5(A,B,C)>0)){
          # if (condition1(A,B,C)>0 & condition2(A,B,C)<0){
              print(paste0('n:',n))
              print(paste0('A:',A))
              print(paste0('B:',B))
              print(paste0('C:',C)) 
          }
        }
      }
    }
  }
}


formula1 = function(A,B,C){
  f = 2*(p-1)*(2-p)*A^2+2*n*(n*p-2)*B^2+(p-1)*2*n*(n*p+2)*A*B
  return(f)
} 
for (n in 30:200){
  for (A in 0.00001:10){
    for (B in 0.0001:10){
      for (C in 0.0001:10){
        # if ((2*B>n*p*C+2*(p-1)*A) & (n*p*C>2*(p-1)*A)){
          if (formula1(A,B,C)<0){
            print(paste0('n:',n))
            print(paste0('A:',A))
            print(paste0('B:',B))
            print(paste0('C:',C)) 
        }
      }
    }
  }
}
