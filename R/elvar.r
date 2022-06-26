ELvar = function(X,Y,Philist,theta0=0,beta=NA,other=FALSE){ # X: obs*p; Y: obs*1; Philist: n-list; Philist[[1]]: ni*(ni*d)

    n = length(Philist)
    if(is.na(beta[1])){
        beta = solve(t(X)%*%X,t(X)%*%Y)
    }
    r = Y - X%*%beta

    d = ncol(Philist[[1]])/nrow(Philist[[1]])
    matE = matrix(0,d,d)
    vecLam = rep(0,d)
    vecni = unlist(lapply(Philist,nrow))
    cumni = c(0,cumsum(vecni))
    R = vector("list",n)

    for(i in 1:n){
        ni = vecni[i]
        Phi = matrix(Philist[[i]],ni^2) # ni^2*d matrix
        matE = matE + t(Phi)%*%Phi
        st = cumni[i] + 1
        ed = cumni[i+1]
        R[[i]] = as.vector(r[st:ed]%*%t(r[st:ed])) # ni^2 vector
        vecLam = vecLam + t(Phi)%*%R[[i]] # d*1
    }
    
    
    Ft = solve(matE[-1,-1],matE[-1,1])
    alpha = 1 - sum(matE[1,-1]*Ft)/matE[1,1]
    theta = solve(matE,vecLam)
    Zi = Di = Mi= rep(0,n)
    nv1sq = 0
        
    theta1 = theta[-1]
    
    for (i in 1:n) { 
        ni = vecni[i]
        Phi = matrix(Philist[[i]],ni^2)
        thetanow = c(theta0,theta1)
        Hnow = Phi%*%thetanow # ni^2 vector
        Rnow = R[[i]]
        deltai = Rnow-Hnow
        Zi[i] = sum(Phi[,1]*deltai)

        temp1 = Phi[,1]-Phi[,-1]%*%Ft
        Di[i] = sum(temp1*(Rnow-theta0*Phi[,1]))/alpha
        Mi[i] = sum(temp1*deltai)/alpha
        nv1sq = nv1sq + Mi[i]^2
    }
    nv1sq = nv1sq/n
         
    if(theta0==0){
        stat = (max(sum(Zi),0))^2/n/nv1sq
        if(stat==0){
            pvalue = 0.5
        }else{
            pvalue = pnorm(sqrt(stat),lower.tail=FALSE)
        }
    }else{
        stat = (sum(Zi))^2/n/nv1sq
        pvalue = pnorm(sqrt(stat),lower.tail=FALSE)*2
    }

    if(other){
        return(list(stat=stat,pvalue=pvalue,Zi=Zi,Di=Di,Mi=Mi,nv1sq=nv1sq))
    }else{
        return(list(stat=stat,pvalue=pvalue))
    }
}

multiELvar = function(X,Y.all,Philist,theta0=0,beta.all=NA,other=FALSE){ # X: obs*p; Y.all: obs*lenT; Philist: n-list; beta.all: d*lenT
    lenT = ncol(Y.all)
    n = length(Philist)
    if(other){
        Z.all = D.all = M.all = matrix(0,n,lenT)
        nv1sq.all = rep(0,lenT)
    }
    stat.all = pvalue.all = rep(0,lenT)
    for(t in 1:lenT){
        Y = Y.all[,t]
        if(is.na(beta.all[1])){
            beta = NA
        }else{
            beta = beta.all[,t]   
        }
        re = ELvar(X,Y,Philist,theta0,beta,other)
        if(other){
            Z.all[,t] = re$Zi
            D.all[,t] = re$Di
            M.all[,t] = re$Mi
            nv1sq.all[t] = re$nv1sq
        }
        stat.all[t] = re$stat
        pvalue.all[t] = re$pvalue
        # print(t)
    }
    if(other){
        return(list(stat.all=stat.all,pvalue.all=pvalue.all,Z.all=Z.all,D.all=D.all,M.all=M.all,nv1sq.all=nv1sq.all))
    }else{
        return(list(stat.all=stat.all,pvalue.all=pvalue.all))
    }
}

GELvar = function(X,Y.all,Philist,theta0=0,beta.all=NA,permnum=1e3){ 

    re = multiELvar(X,Y.all,Philist,theta0,beta.all,other=TRUE)
    M.all = re$M.all
    stat.all = re$stat.all
    nv1sq.all = re$nv1sq.all
    n = nrow(M.all)
    lenT = ncol(M.all)
    permmat = matrix(rnorm(permnum*n),permnum,n)
    statperm = matrix(0,permnum,lenT)

    for (perm in 1:permnum) {
        xi = permmat[perm,]
        Mperm = M.all*matrix(xi,n,lenT)        
        statperm[perm,] = (pmax(colSums(Mperm),0))^2/n/nv1sq.all     
    }

    stat.global = max(stat.all)
    pvalue.global = length(which(apply(statperm,1,max)>max(stat.all)))/permnum
   
    return(list(stat.global=stat.global,pvalue.global=pvalue.global)) 
}



