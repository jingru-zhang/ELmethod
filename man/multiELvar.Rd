\name{multiELvar}
\alias{multiELvar}
\title{Empirical Likelihood Inference of Variance Components at multiple time points}
\description{This function provides an empirical likelihood method for the inference of variance components at multiple time points in linear mixed-effects models.}
\usage{
    multiELvar(X,Y.all,Philist,theta0=0,beta.all=NA,other=FALSE)
}
\arguments{
    \item{X}{design matrix for all observations, in which each row represents a p-dimentional covariates.}
    \item{Y.all}{response matrix, in which each column is the response vector at time t.}
    \item{Philist}{list of design matrices of variance components. Its i-th element is an ni by d*ni matrix that combines design matrices of variance components by columns for the i-th subject, where ni is the number of repeated measures for the i-th subject and d is the number of variance components.}
    \item{theta0}{value of the first variance component under the null. Its default value is 0.}
    \item{beta.all}{fixed effects. Each column is the fixed effects at time t. Its default value is NA (unknown fixed effects).}
    \item{other}{logical; if TRUE, the function gives auxiliary terms. Its default value is FALSE.}
}
\value{
    \item{stat.all}{vector of test statistics at multiple time points.}
    \item{pvalue.all}{vector of approximated p-value at multiple time points based on asymptotic theory.}
    \item{Z.all,D.all,M.all,nv1sq.all}{auxiliary terms if other=TRUE.}
}
\references{
    Zhang J., Guo W., Carpenter J.S., Leroux A., Merikangas K.R., Martin N.G., Hickie I.B., Shou H., and Li H. (2022). Empirical likelihood tests for variance components in linear
mixed-effects models. 
}
\examples{

# Datasets "exampleNE0" and "exampleNE1" contain normal distributed longitudinal data.
# Datasets "exampleTE0" and "exampleTE1" contain t distributed longitudinal data.
# The fist variance components in the datasets "exampleNE0" and "exampleTE0" are zero.
# The fist variance components in the datasets "exampleNE1" and "exampleTE1" are nonzero at the 24, 25, 26, 27 time points.

# X is an N by p matrix with N being the number of all observations and p being the dimension of covariates.
# Y.all is an N by T matrix with T being the number of time points.
# Philist is an n list of design matrices of variance components with n being the number of subjects. Its $i$th element Philist[[i]] is an $n_i$ by $n_id$ matrix that combines design matrices of variance components by columns for the $i$th subject, where $n_i$ is the number of repeated measures for the $i$th subject and $d$ is the number of variance components.
# beta.all is a p by T matrix. Each column is the fixed effects at time t.
# thetastar is a d by T matrix. Each column is the variance components at time t.

    data(exampleNE0)
    re = multiELvar(X,Y.all,Philist,theta0=0)
}

\seealso{
    \code{\link{GELvar}}
}













