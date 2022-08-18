library(nloptr)
library(graphicalExtremes)
library(latex2exp)
source("algorithm_functions.R")

#########################################
#### Efficiency comparison for BIBDs ####
#########################################



### Set parameters

m=6
k=3

### Define BIBD designs

w_bibd = c(0,1,0,1,0,1,1,1,0,0,1,1,0,0,0,1,0,1,0,1)/10
w_bibd2 = rep(1/10,20)-w_bibd

### Compute efficiencies

eff=matrix(0,nrow=100,ncol = 3)

for (ell in (0:99)) {
  pi1=exp(ell/10)
  p=c(pi1,pi1^(1/2),pi1^(5/4),pi1^(7/4),pi1^(3/4),1)
  
  Q_complete= diag(L(p=p))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p,k=k)))%*%rep(1/choose(m,k),choose(m,k))
  Theta_complete_m = Qvec2Theta(Q_complete)[-m,-m]
  complete_det = det(Theta_complete_m)
  
  Q_bibd= diag(L(p=p))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p,k=k)))%*%w_bibd
  Theta_bibd_m = Qvec2Theta(Q_bibd)[-m,-m]
  BIBD_det = det(Theta_bibd_m)
  
  Q_bibd2= diag(L(p=p))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p,k=k)))%*%w_bibd2
  Theta_bibd2_m = Qvec2Theta(Q_bibd2)[-m,-m]
  BIBD2_det = det(Theta_bibd2_m)
  
  G_opt = discrete_choice_design(p=p,k=k)$Gamma_opt
  OD_det = 1/CM(Gamma2vec(G_opt))
  
  eff[ell+1,]=c((complete_det/OD_det)^(1/(m-1)),(BIBD_det/OD_det)^(1/(m-1)),(BIBD2_det/OD_det)^(1/(m-1)))
}

### Plot results

plot((0:99)/10,eff[,1],xlab = TeX("$ \\log(\\pi_1) $"),type="l",ylab = "D-eff",ylim = c(0.2,1),lwd=2)
lines((0:99)/10,eff[,2],lty=2,lwd=2)
lines((0:99)/10,eff[,3],lty=3,lwd=2)

### Limiting design

ell=100
pi1=exp(ell/10)
p=c(pi1,pi1^(1/2),pi1^(5/4),pi1^(7/4),pi1^(3/4),1)
w_lim = discrete_choice_design(p=p,k=k)$w
round(w_lim,digits = 5)

