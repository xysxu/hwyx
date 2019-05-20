#' fits elastic net to data using coordinate descent algorithm
#' @param x independent variable
#' @param y dependent variable
#' @param lambda LASSO tuning
#' @param alpha proportion of usage in l1 penalty comparing with l2
#' @param tol tolerance
#' @return estimation of parameter beta in y=beta*x
#' @export
#' @examples
#' y = x%*%beta+e
#' elnet_coord(x, y, lambda=1, alpha=0.5, tol=1e-6)
#'

elnet_coord = function(x, y, lambda, alpha, tol){
    n = length(y)
    p = ncol(x)
    N =100000
    cda(x, y, n, p, lambda, alpha, N, tol)
}









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



