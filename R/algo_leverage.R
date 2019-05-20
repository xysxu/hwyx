#' algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows
#' @param x independent variable
#' @param y dependent variable
#' @param r subsample size
#' @return estimation of parameter beta in y=beta*x
#' @export
#' @examples e = rnorm(500,0,1)
#' x = rt(500, df = 6)
#' y = -x+e
#' algo_leverage(x, y, r=10, methods="UNIF")

algo_leverage = function(x, y, r, methods){
    nu =c(1:length(x))
    if(methods == "UNIF"){
        sim_UNIF(x, y, nu, r)
    }else{
        sim_BLEV(x, y, nu, r)
    }
}

#UNIF
sim_UNIF = function(x, y, nu, r) {
    s = sample(nu,r,replace = T)
    y_r = y[s]
    x_r = x[s]
    phi = rep(1/r, r)
    lmodel = lm(y_r ~ 0 + x_r, weights = phi)
    coef(lmodel)

}

#BLEV
sim_BLEV = function(x, y, nu, r) {
    s = sample(nu,r,replace = T,prob=w)
    y_r = y[s]
    x_r = x[s]
    phi = w[s]
    P = sqrt(1/phi)
    lmodel = lm(y_r ~ 0 + x_r, weights = P)
    coef(lmodel)

}




