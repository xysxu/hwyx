% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_ols.R
\name{solve_ols}
\alias{solve_ols}
\title{solves a linear system using Gauss-Seidel or Jacobi method}
\usage{
solve_ols(A, b, tol, methods, p)
}
\arguments{
\item{A}{matrix}

\item{b}{vector}

\item{tol}{tolarence}

\item{methods}{"GS", "JS", "JP" for Gauss-Seide, Jacobi(sequential), Jacobi(parallel)}

\item{p}{number of cores to use for parallel}
}
\value{
estimation of parameter x in Ax=b
}
\description{
solves a linear system using Gauss-Seidel or Jacobi method
}
\examples{
a1=2
n1=10
D1= matrix(0, n1, n1)
L1= matrix(0, n1, n1)
U1= matrix(0, n1, n1)
for(i in 1:n1-1){
    D1[i,i]=a1;
    L1[i,i+1]=-1;
    U1[i+1,i]=-1
}
D1[n1,n1]=a1
A=D1+L1+U1
m1=n1/2
v=rep(c(1,0),times=m1)
b=A\%*\%v
p=2
solve_ols(A, b, tol = 1e-10, methods="JP", p=2)
}
