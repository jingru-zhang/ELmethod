# ELmethodVar
Empirical Likelihood Inference of Variance Components in Linear Mixed-Effects Models

## Description
This package provides empirical likelihood-based methods for the inference of  variance components in linear mixed-effects models.

## Maintainer
Jingru Zhang (jingru.zhang@pennmedicine.upenn.edu)

## Examples
```{r}
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

# Empirical likelihood method for the inference of a local variance component
t = 1 # consider the local problem at time t
re = ELvar(X,Y.all[,t],Philist,theta0=0) # with unknown fixed effects
re = ELvar(X,Y.all[,t],Philist,theta0=0,beta=beta.all[,t]) # with known fixed effects
    
# Empirical likelihood method for the inference of variance components over an interval
re = GELvar(X,Y.all,Philist,theta0=0)
```

## Depends
R (>= 3.5.0)

## References
Zhang J., Guo W., Carpenter J.S., Leroux A., Merikangas K.R., Martin N.G., Hickie I.B., Shou H., and Li H. (2022). Empirical likelihood tests for variance components in linear mixed-effects models. 
