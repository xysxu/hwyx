#' solves a linear system using Gauss-Seidel or Jacobi(parallel) method
#' @param A matrix
#' @param b vector
#' @param tol tolarence
#' @param p number of cores to use for parallel
#' @return estimation of parameter x in Ax=b
#' @export

solve_ols = function(A, b, tol, methods, p){
    D = diag(A)
    D = diag(D)
    L = A
    L[lower.tri(L,diag=TRUE)] = 0
    U = A
    U[upper.tri(U,diag=TRUE)] = 0
    if(methods == "GS"){
        C_g=norm(solve(L + D)%*%U, type = "2")
        if(C_g < 1){
            GS(A, b, tol)
        }else{
            cat("Gauss-Seidel doesn't converge")
            stop()
        }
    }else{
        if(methods == "JS"){
            C_j=norm(-solve(D)%*%(L + U),type = "2")
            if(C_j < 1){
                JS(A, b, tol)
            }else{
                cat("Jacob Sequential doesn't converge")
                stop()
            }
        }else{
            C_j=norm(-solve(D)%*%(L + U),type = "2")
            if(C_j < 1){
                JP(A, b, tol, p)
            }else{
                cat("Jacob Paralell doesn't converge")
                stop()
            }
        }
    }
}

# #TEST
# a1=2
# n1=10
# D1= matrix(0, n1, n1)
# L1= matrix(0, n1, n1)
# U1= matrix(0, n1, n1)
# for(i in 1:n1-1){
#     D1[i,i]=a1;
#     L1[i,i+1]=-1;
#     U1[i+1,i]=-1
# }
# D1[n1,n1]=a1
# A=D1+L1+U1
# m1=n1/2
# v=rep(c(1,0),times=m1)
# b=A%*%v
# p=2
# tol = 1e-10
# solve_ols(A, b, tol = 1e-10, methods="JS", p=2)



#Gauss-Seidel
GS = function(A, b, tol){
    N=100000
    n=length(b)
    x=matrix(0, nrow=n, ncol=N)
    x[,1]=rep(0,n)
    k=1
    while(k < N+1){
        for(i in 1:n){
            if(i>1){
                if (i == n) {
                    x[i,k] = (-sum(A[i,1:(i-1)] * x[1:(i-1),k])+ b[i])/ A[i,i]
                } else {
                    x[i,k] = (-sum(A[i,1:(i-1)] * x[1:(i-1),k]) -
                                  sum(A[i,(i+1):n] * x[(i+1):n,k-1]) + b[i]) / A[i,i]
                }
            }else{
                x[i,k] = ( -sum(A[i,(i+1):n] * x[(i+1):n,k-1]) + b[i]) / A[i,i]
            }

        }
        if(k > 1 & sum((x[,k]-x[,k-1])^2)<tol){
            cat("stop at", k, "th iteration")
            return(x[,k])
        }else{
            k=k+1
        }
    }
}


#jacobi sequential "JS"
JS = function(A, b, tol){
n=length(b)
N=100000
x=matrix(0, nrow=n, ncol=N)
x[,1]=rep(0,n)
k=1
while(k < N+1){
    for(i in 1:n){
        if(i>1){
            if (i == n) {
                x[i,k] = (-sum(A[i,1:(i-1)] * x[1:(i-1),k-1])+ b[i])/ A[i,i]
            } else {
                x[i,k] = (-sum(A[i,1:(i-1)] * x[1:(i-1),k-1]) -
                              sum(A[i,(i+1):n] * x[(i+1):n,k-1]) + b[i]) / A[i,i]
            }
        }else{
            x[i,k] = ( -sum(A[i,(i+1):n] * x[(i+1):n,k-1]) + b[i]) / A[i,i]
        }

    }
    if(k > 1 & sum((x[,k]-x[,k-1])^2)<tol){
        cat("stop at", k, "th iteration")
        return(x[,k])
    }else{
        k=k+1
    }
}

}


#Jacobi parallel
JP =  function(A, b, tol, p){
l=detectCores(all.tests = FALSE, logical = TRUE)
if(p > l){
    cat("number of parallel overflow")
    stop()
}else{
cl = makeCluster(p)
registerDoParallel(cl)
N=100000
n=length(b)
x=matrix(0, nrow=n, ncol=N)
x[,1]=rep(0,n)
k=1
while(k < N+1){
    x[,k] = foreach(i = 1:n, .combine = "c") %dopar% {
        if(i>1){
            if (i == n) {
               (-sum(A[i,1:(i-1)] * x[1:(i-1),k-1])+ b[i])/ A[i,i]
            } else {
                (-sum(A[i,1:(i-1)] * x[1:(i-1),k-1]) -
                              sum(A[i,(i+1):n] * x[(i+1):n,k-1]) + b[i]) / A[i,i]
            }
        }else{
            (-sum(A[i,(i+1):n] * x[(i+1):n,k-1]) + b[i]) / A[i,i]
        }
    }
    if(k > 1 & sum((x[,k]-x[,k-1])^2)<tol){
        cat("stop at", k, "th iteration")
        return(x[,k])
    }else{
        k=k+1
    }
}

}
}










