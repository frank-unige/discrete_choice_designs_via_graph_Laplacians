# This file contains the algorithm and help functions

optimal_Gamma = function(p,k,algorithm="NLOPT_LD_SLSQP",tol=1e-15, maxeval=500,print_level=0){
  m=length(p)
  p=p/min(p)
  Smatrix = S_norm(p,k)
  negCM = function(Gvec){
    G = matrix(0,nrow = m,ncol = m)
    G[lower.tri(G, diag=FALSE)] <- Gvec 
    Gam = G+t(G)
    CM= log(det(rbind(c(0,rep(-1,m)),cbind(rep(1,m),-Gam/2))))
    return(-CM)
  }
  grad = function(Gvec){
    grad = -Theta2Qvec(Gamma2Theta(Gvec2Gamma(Gvec)))
    return(grad)
  }
  dirder = function(Gvec){
    ineq = Smatrix %*% Gvec - rep(m-1,choose(m,k))
    return(ineq)
  }
  dirder_jacob = function(Gvec){
    jacob= Smatrix
    return(jacob)
  }
  G = nloptr(x0 = rep(2,choose(m,2)),lb=rep(0,choose(m,2)), eval_f = negCM, eval_grad_f = grad, eval_g_ineq = dirder,eval_jac_g_ineq = dirder_jacob ,opts = list("algorithm"=algorithm,maxeval=maxeval, xtol_rel=tol,print_level=print_level))
  G_new = Gvec2Gamma(G$solution)
  return(G_new)
}

design_from_Gamma = function(Gamma,p,k,algorithm="NLOPT_LD_SLSQP",tol=1e-15){
  m=length(p)
  p=p/min(p)
  G_new = Gamma
  Q_new= Theta2Qvec(Gamma2Theta(G_new))
  A= diag(L(p=p))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p,k=k)))
  fn = function(w){
    return(sum((A%*%w-Q_new)^2))
  }
  gr = function(w){
    2*t(A)%*%(A%*%w-Q_new)
  }
  heq=function(w){
    return(sum(w)-1)
  }
  heqjac=function(w){
    return(rep(1,choose(m,k)))
  }
  result=slsqp(x0=rep(1/choose(m,k),choose(m,k)),fn=fn,gr=gr, lower = rep(0,choose(m,k)),heq = heq, heqjac = heqjac,control = list(xtol_rel = tol))
  w_new=result$par
  check_dirder_new=S_norm(p,k)%*%Gamma2vec(Theta2Gamma(Qvec2Theta(A%*%w_new)))-rep(m-1,choose(m,k))
  return(list(w = w_new,dirder_max = max(check_dirder_new), sol_distance = result$value))
}


discrete_choice_design = function(p,k,algorithm="NLOPT_LD_SLSQP",tol=1e-15, maxeval=500,print_level=0){
  G_new = optimal_Gamma(p=p,k=k,algorithm = algorithm, tol = tol, maxeval = maxeval, print_level = print_level)
  result=design_from_Gamma(p=p,k=k,Gamma=G_new,algorithm = algorithm, tol = tol)
  w_new=result$w
  return(list(w = w_new,Gamma_opt = G_new, dirder_max = max(result$dirder_max), sol_distance = result$sol_distance))
}




hyper_incidence = function(m,k){
  a= combn(m,2,simplify = TRUE)
  b= combn(m,k,simplify = TRUE)
  Z=matrix(nrow = choose(m,2) ,ncol = choose(m,k))
  for (i in 1:choose(m,2)) {
    for (j in 1:choose(m,k)) {
      Z[i,j]=as.integer(all(a[,i] %in% b[,j]))
    }
  }
  return(Z)
}

R= function(p,k){
  m=length(p)
  R=vector(length = choose(m,k))
  b= combn(m,k,simplify = TRUE)
  for (i in 1:choose(m,k)) {
    R[i] = sum(p[b[,i]])^2
  }
  return(R)
}

L=function(p){
  m=length(p)
  L=vector(length = choose(m,2))
  a= combn(m,2,simplify = TRUE)
  for (i in 1:choose(m,2)) {
    L[i] = p[a[1,i]]*p[a[2,i]]
  }
  return(L)
}

G_bar <- function(p,k){
  m = length(p)
  Y= solve(diag(L(p=p)))%*% corpcor::pseudoinverse(t(hyper_incidence(m,k))) %*% diag(R(p=p,k=k)) %*%rep(m-1,choose(m,k))
  Gb= matrix(0, m, m)
  Gb[lower.tri(Gb, diag=FALSE)] <- Y
  return(t(Gb)+Gb)
}

CM = function(Gvec){
  n=length(Gvec)
  m=sqrt(2*n+1/4)+1/2
  G = matrix(0,nrow = m,ncol = m)
  G[lower.tri(G, diag=FALSE)] <- Gvec 
  Gam = G+t(G)
  CM = det(rbind(c(0,rep(-1,m)),cbind(rep(1,m),-Gam/2)))
  return(CM)}

Theta2Qvec <-function(Theta){
  m = ncol(Theta)
  Q <- vector(length = choose(m,2))
  Q<- -Theta[lower.tri(Theta, diag=FALSE)]  
  return(Q)
}

Gamma2vec <- function(Gamma){
  m = ncol(Gamma)
  G <- vector(length = choose(m,2))
  G<- Gamma[lower.tri(Gamma, diag=FALSE)]  
  return(G)
}

Gvec2Gamma <-function(Gvec){
  n=length(Gvec)
  m=sqrt(2*n+1/4)+1/2
  G = matrix(0,nrow = m,ncol = m)
  G[lower.tri(G, diag=FALSE)] <- Gvec
  return(G+t(G))
}

Qvec2Theta <- function(Qvec){
  n=length(Qvec)
  m=sqrt(2*n+1/4)+1/2
  Theta = matrix(0,nrow = m,ncol = m)
  Theta[lower.tri(Theta, diag=FALSE)]=-Qvec
  Theta = Theta +t(Theta)
  diag(Theta)=-colSums(Theta)
  return(Theta)
}

Qvec2Thetawithm <- function(Qvec,m){
  Theta = matrix(0,nrow = m,ncol = m)
  Theta[lower.tri(Theta, diag=FALSE)]=-Qvec
  Theta = Theta +t(Theta)
  diag(Theta)=-colSums(Theta)
  return(Theta)
}

S_norm = function(p,k){
  m=length(p)
  return(solve(diag(R(p=p,k=k)))%*%t(hyper_incidence(m,k))%*%diag(L(p=p)))
}