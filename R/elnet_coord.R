#' fits elastic net to data using coordinate descent algorithm
#' @param x independent variable
#' @param y dependent variable
#' @param lambda LASSO tuning
#' @param alpha proportion of usage in l1 penalty comparing with l2
#' @return estimation of parameter beta in y=beta*x
#' @export

elnet_coord = function(x, y, lambda, alpha, tol){
    n = length(y)
    p = ncol(x)
    N =100000
    cda(x, y, n, p, lambda, alpha, N, tol)
}

#TEST
# rm(list = ls())
# p = 20
# n = 100
# x=matrix(0,nrow = n, ncol = p)
# beta=c(2,0,-2,0,1,0,-1,0,rep(0,12))
# x[,1:2] = mvrnorm(n, mu = c(0,0), Sigma = matrix(c(1,0.8,0.8,1), ncol = 2))
# x[,5:6] = mvrnorm(n, mu = c(0,0), Sigma = matrix(c(1,0.8,0.8,1), ncol = 2))
# x[,3:4] = mvrnorm(n, mu = c(0,0), Sigma = matrix(c(1,0,0,1), ncol = 2))
# for(i in 7:p){
#     x[,i] = rnorm(n, mean = 0, sd = 1)
# }
# e = rnorm(n, mean = 0, sd = 1)
# y = x%*%beta+e
# elnet_coord(x, y, lambda=1, alpha=0.5, tol=1e-6)





#coordinate descent algorithm
cda = function(x, y, n, p, lambda, alpha, N, tol = 1e-6){
    beta=matrix(0, p, N)
    for(k in 1:(N-1)){
        for(j in 1:p){
            if(j > 1){
                if(j == p){
                    r = y - x[,1:(j-1)] %*% as.matrix(beta[1:(j-1),k+1])
                }else{
                    r = y - x[,1:(j-1)] %*% as.matrix(beta[1:(j-1),k+1]) - x[,(j+1):p] %*% as.matrix(beta[(j+1):p, k])
                }
            }else{
                r = y - x[,(j+1):p] %*% as.matrix(beta[(j+1):p,k])
            }
            t = sum(x[,j] * r)
            pos_part = ifelse(abs(t) - alpha * lambda > 0, abs(t) - alpha * lambda, 0)
            beta[j,k+1] = sign(t) * pos_part / (sum(x[,j] * x[,j]) + 2 * lambda * (1-alpha))
        }

        if (sum((beta[,k+1] - beta[,k]) ** 2) < tol)
            return(beta[,k+1])
    # cat("stop at",k, "\n")
    }
    return(beta[,k+1])
}



