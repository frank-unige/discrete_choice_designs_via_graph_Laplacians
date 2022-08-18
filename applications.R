library(nloptr)
library(graphicalExtremes)
library(hyper2)
library(xtable)
source("algorithm_functions.R")


#####################################
#### Climate change perception ######
#####################################

#load icons data and set parameters

data(icons)
p_est=icons_maxp
m=6
k=4

### Study design

w_hankin=c(0,0,0,11/124,15/124,9/124,18/124,18/124,8/124,0,11/124,16/124,18/124,0,0)

### Compute BIBD information determinant

Q_bibd= diag(L(p=p_est))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p_est,k=k)))%*%rep(1/choose(m,k),choose(m,k))
Theta_m_bibd = Qvec2Theta(Q_bibd)[-m,-m]
BIBD_det = det(Theta_m_bibd)

### Compute study design information determinant

Q_hankin= diag(L(p=p_est))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p_est,k=k)))%*%w_hankin
Theta_m_hankin = Qvec2Theta(Q_hankin)[-m,-m]
hankin_det = det(Theta_m_hankin)

### Compute D-optimal design

result_opt = discrete_choice_design(p=p_est,k=k)
G_opt=result_opt$Gamma_opt
w_opt=result_opt$w

### Compute rounded design from optimal design

(round(w_opt*124)-w_opt*124)
w_opt_rounded=round(w_opt*124)
w_opt_rounded[8]=11
w_opt_rounded=w_opt_rounded/124
sum(w_opt_rounded)

### Compute rounded design information determinant

Q_opt_rounded= diag(L(p=p_est))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p_est,k=k)))%*%w_opt_rounded
Theta_m_rounded = Qvec2Theta(Q_opt_rounded)[-m,-m]
rounded_det = det(Theta_m_rounded)

### Compute efficiencies

OD_det = 1/CM(Gamma2vec(G_opt))

efficiencies = matrix(0,nrow = 1,ncol = 3)
colnames(efficiencies) = c("study","complete","rounded")
rownames(efficiencies) = c("D-efficiency")

efficiencies[1,]=c((hankin_det/OD_det)^(1/(m-1)),
                   (BIBD_det/OD_det)^(1/(m-1)), 
                   (rounded_det/OD_det)^(1/(m-1)))

xtable(efficiencies,digits = 5)



######################
#### Cricket data ####
######################

### Load data and parameters

data(T20)
p_est=T20_maxp
m=length(p_est)
k=2

### Read study design from data

nam=names(T20_maxp)
pair = combn(nam,2,simplify = TRUE)
w_cricket = rep(0,ncol(pair))

for (i in 1:ncol(pair)) {
  for (j in 1:nrow(T20_table)) {
    if((pair[1,i]==T20_table[j,1]&&pair[2,i]==T20_table[j,2] )|(pair[1,i]==T20_table[j,2] &&pair[2,i]==T20_table[j,1])) w_cricket[i]= w_cricket[i]+1
  }
}

w_cricket=w_cricket/nrow(T20_table)
sum(w_cricket)

### Compute Study design information determinant

Q_cricket= diag(L(p=p_est))%*%solve(diag(R(p=p_est,k=k)))%*%w_cricket
Theta_m_cricket = Qvec2Theta(Q_cricket)[-m,-m]
cricket_det = det(Theta_m_cricket)

### Compute constant design information determinant

Q_bal= diag(L(p=p_est))%*%hyper_incidence(m,k)%*%solve(diag(R(p=p_est,k=k)))%*%rep(1/choose(m,k),choose(m,k))
Theta_m_bal = Qvec2Theta(Q_bal)[-m,-m]
bal_det = det(Theta_m_bal)

### Compute D-optimal design

result_opt = discrete_choice_design(p=p_est,k=k)
G_opt=result_opt$Gamma_opt
w_opt=result_opt$w

### Compute efficiencies

OD_det = 1/CM(Gamma2vec(G_opt))

efficiencies = matrix(0,nrow = 1,ncol = 2)
colnames(efficiencies) = c("study","complete")
rownames(efficiencies) = c("D-efficiency")

efficiencies[1,]=c((cricket_det/OD_det)^(1/(m-1)),
                   (bal_det/OD_det)^(1/(m-1)))

xtable(efficiencies,digits = 5)

