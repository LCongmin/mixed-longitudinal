##########################################################################
######  带方差CP，自动窗宽，分开估+同时估，2015+2019cao,h_first     ######
##########################################################################
rm(list = ls (all = TRUE))
library(MASS)
n<-400
K<-1000
t<-0.3

#分开估beta的窗宽
h.can.beta<-seq(n^(-0.8),n^(-0.6),length=10)
#分开估gamma的窗宽，同时估betagamma的窗宽
h.can.gamma<-seq(n^(-0.5),n^(-0.4),length=10)

beta.hat_p_2019<-array(0,c(length(h.can.beta),4,K))
beta.hat1_p_2019<-array(0,c(length(h.can.beta),4,K))
beta.hat2_p_2019<-array(0,c(length(h.can.beta),4,K))

gamma.hat_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat1_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat2_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hatd_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat1d_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat2d_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))

mbeta_p_2019<-matrix(0:0,nrow = K,ncol = length(h.can.beta))
mbeta1_p_2019<-matrix(0:0,nrow = K,ncol = length(h.can.beta))
mbeta2_p_2019<-matrix(0:0,nrow = K,ncol = length(h.can.beta))

mgamma_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
mgamma1_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
mgamma2_p_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))

htvar_beta_p_2019<-matrix(0:0,ncol = K,nrow = length(h.can.beta))
htvar_gamma_p_2019<-array(0,c(length(h.can.beta),length(h.can.gamma),K))

beta.hat_2019<-array(0,c(length(h.can.beta),2,K))
beta.hat1_2019<-array(0,c(length(h.can.beta),2,K))
beta.hat2_2019<-array(0,c(length(h.can.beta),2,K))

gamma.hat_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat1_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat2_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hatd_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat1d_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
gamma.hat2d_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))

mbeta_2019<-matrix(0:0,nrow = K,ncol = length(h.can.beta))
mbeta1_2019<-matrix(0:0,nrow = K,ncol = length(h.can.beta))
mbeta2_2019<-matrix(0:0,nrow = K,ncol = length(h.can.beta))

mgamma_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
mgamma1_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))
mgamma2_2019<-array(0,dim = c(length(h.can.beta),length(h.can.gamma),K))

htvar_beta_2019<-matrix(0:0,ncol = K,nrow = length(h.can.beta))
htvar_gamma_2019<-array(0,c(length(h.can.beta),length(h.can.gamma),K))

betagamma.hat_2019<-array(0,c(length(h.can.gamma),4,K))
betagamma.hat1_2019<-array(0,c(length(h.can.gamma),4,K))
betagamma.hat2_2019<-array(0,c(length(h.can.gamma),4,K))

mbetagamma_2019<-matrix(0:0,nrow = K,ncol = length(h.can.gamma))
mbetagamma1_2019<-matrix(0:0,nrow = K,ncol = length(h.can.gamma))
mbetagamma2_2019<-matrix(0:0,nrow = K,ncol = length(h.can.gamma))

hvarb_2019<-matrix(0:0,ncol = K,nrow = length(h.can.gamma))
hvarg_2019<-matrix(0:0,ncol = K,nrow = length(h.can.gamma))


#fun.beta<-function(t){
#  Beta<-0.4*t+0.5
#  return(Beta)
#}
fun.beta<-function(t){
  Beta<-3*((t-0.4)^2)
  return(Beta)
}

#fun.gamma<-function(t){
#  Gamma<-sqrt(t)
#  return(Gamma)
#}
fun.gamma<-function(t){
  Gamma<-sin(2*pi*t)
  return(Gamma)
}
epanechnikov <- function(t){
  tst <- (-1.0 <= t) & (t <= 1.0 )
  kt <- 0.75*( 1.0 - t*t )
  kt[ !tst ] <- 0.0
  return(kt)
}
#核函数#
local_kernel <- function(t, h){
  kt <- epanechnikov(t/h)
  return(kt/h)
}

beta_est_omitgamma<-function(x,z,y,T.x,S.z,n,nx,nz,t,h){
  kw<-list()
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  
  #omitted估计x的系数beta
  Q<-matrix(0:0,ncol=2,nrow=2)#求解系数的2×2矩阵
  q<-matrix(0:0,ncol = 1,nrow = 2)#求解系数的2×1矩阵
  for (i in 1:n) {
    q<-q+matrix(c(sum(kw[[i]]*x[[i]]*y[[i]])/n,
                  sum(kw[[i]]*x[[i]]*y[[i]]*(T.x[[i]]-t))/n),ncol = 1,nrow = 2)
    Q<-Q+matrix(c(sum(kw[[i]]*x[[i]]*x[[i]])/n,
                  sum(kw[[i]]*x[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*x[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n),
                ncol = 2,nrow = 2)
  }
  beta_omitgamma<-solve(Q,q)#l代表第几次模拟
  
  sig_beta<-matrix(0:0,ncol = 2,nrow = 2)
  
  for (i in 1:n) {
    sig_beta[1,1]<-sig_beta[1,1]+sum((kw[[i]]*x[[i]]*(y[[i]]-x[[i]]*beta_omitgamma[1,1]-x[[i]]*beta_omitgamma[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(y[[i]]-x[[i]]*beta_omitgamma[1,1]-x[[i]]*beta_omitgamma[2,1]*(T.x[[i]]-t))))/n
    sig_beta[1,2]<-sig_beta[1,2]+sum((kw[[i]]*x[[i]]*(y[[i]]-x[[i]]*beta_omitgamma[1,1]-x[[i]]*beta_omitgamma[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-x[[i]]*beta_omitgamma[1,1]-x[[i]]*beta_omitgamma[2,1]*(T.x[[i]]-t))))/n
    sig_beta[2,2]<-sig_beta[2,2]+sum((kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-x[[i]]*beta_omitgamma[1,1]-x[[i]]*beta_omitgamma[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-x[[i]]*beta_omitgamma[1,1]-x[[i]]*beta_omitgamma[2,1]*(T.x[[i]]-t))))/n
  }
  sig_beta[2,1]<-sig_beta[1,2]
  hvar_beta_omitgamma<-ginv(Q)%*%sig_beta%*%ginv(Q)
  var_beta_omitgamma<-hvar_beta_omitgamma[1,1]
  return(list("beta_omitgamma"=beta_omitgamma,"var_beta_omitgamma"=var_beta_omitgamma))
}
#omit gamma beta均方预测误差#
m.beta_omitgamma<-function(x,z,y,T.x,S.z,n,nx,nz,t,h,beta1,beta2){
  kw<-list()
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  #先将所有核权重求和
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  betam_omitgamma<-0
  for (i in 1:n) {
    betam_omitgamma<-betam_omitgamma+sum(kw[[i]]*((y[[i]]-x[[i]]*(beta1+beta2*(T.x[[i]]-t)))^2))/kw.all
  }
  return(betam_omitgamma)
}

#beta部分线性估计函数
beta.est.p<- function(x,z,y,T.x,S.z,n,nx,nz,t,h) {
  #部分线性方法估计x的系数beta
  kw<-list()
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  #先将所有核权重求和
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  Q<-matrix(0:0,ncol=4,nrow=4)#求解系数的2×2矩阵
  q<-matrix(0:0,ncol = 1,nrow = 4)#求解系数的2×1矩阵
  for (i in 1:n) {
    q<-q+matrix(c(sum(kw[[i]]*y[[i]])/n,
                  sum(kw[[i]]*(T.x[[i]]-t)*y[[i]])/n,
                  sum(kw[[i]]*x[[i]]*y[[i]])/n,
                  sum(kw[[i]]*x[[i]]*y[[i]]*(T.x[[i]]-t))/n),ncol = 1,nrow = 4)
    Q<-Q+matrix(c(sum(kw[[i]])/n,
                  sum(kw[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]])/n,
                  sum(kw[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]])/n,
                  sum(kw[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*x[[i]])/n,
                  sum(kw[[i]]*x[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*x[[i]]*(T.x[[i]]-t))/n,
                  sum(kw[[i]]*x[[i]]*x[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n),
                ncol = 4,nrow = 4,byrow = TRUE)
  }
  beta.hat_p<-ginv(Q)%*%q#l代表第几次模拟
  
  sig_beta_p<-matrix(0:0,ncol = 4,nrow = 4)
  for (i in 1:n) {
    sig_beta_p[1,1]<-sig_beta_p[1,1]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    sig_beta_p[1,2]<-sig_beta_p[1,2]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    sig_beta_p[1,3]<-sig_beta_p[1,3]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    sig_beta_p[1,4]<-sig_beta_p[1,4]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    
    sig_beta_p[2,2]<-sig_beta_p[2,2]+sum((kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    sig_beta_p[2,3]<-sig_beta_p[2,3]+sum((kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    sig_beta_p[2,4]<-sig_beta_p[2,4]+sum((kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t)*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    
    sig_beta_p[3,3]<-sig_beta_p[3,3]+sum((kw[[i]]*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    sig_beta_p[3,4]<-sig_beta_p[3,4]+sum((kw[[i]]*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
    
    sig_beta_p[4,4]<-sig_beta_p[4,4]+sum((kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
  }
  sig_beta_p[2,1]<-sig_beta_p[1,2]
  sig_beta_p[3,1]<-sig_beta_p[1,3]
  sig_beta_p[4,1]<-sig_beta_p[1,4]
  sig_beta_p[3,2]<-sig_beta_p[2,3]
  sig_beta_p[4,2]<-sig_beta_p[2,4]
  sig_beta_p[4,3]<-sig_beta_p[3,4]
  
  hbe_p<-ginv(Q)%*%sig_beta_p%*%ginv(Q)
  hb_p<-hbe_p[3,3]
  return(list("beta_p"=beta.hat_p,"hb_p"=hb_p))
  
}

#gamma部分线性估计函数与中心化一样

#beta部分线性预测误差
m.beta.p<-function(x,z,y,T.x,S.z,n,nx,nz,t,h,beta1,beta2,beta3,beta4){
  kw<-list()
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  #先将所有核权重求和
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  
  betam<-0
  for (i in 1:n) {
    betam<-betam+sum(kw[[i]]*((y[[i]]-(beta1+beta2*(T.x[[i]]-t))-x[[i]]*(beta3+beta4*(T.x[[i]]-t)))^2))/kw.all
  }
  return(betam)
}

#gamma部分线性预测误差与中心化一样




#beta中心化估计函数#
beta.est<-function(x,z,y,T.x,S.z,n,nx,nz,t,h){
  kw<-list()
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  #先将所有核权重求和
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  km.x<-0#代表所有数据x的核估计
  for (i in 1:n) {
    km.x<-km.x+sum(x[[i]]*kw[[i]]/kw.all)#每个个体的均值核估计
  }
  #中心化之后的x的数据，记为x.tilde
  for (i in 1:n) {
    x.tilde[[i]]<-x[[i]]-km.x
  }
  #得到数据y的核估计，这里应该是所有个体放在一起的核估计吧？
  km.y<-0#代表所有数据y的核估计
  for (i in 1:n) {
    km.y<-km.y+sum(y[[i]]*kw[[i]])/kw.all#每个个体的均值核估计
  }
  #中心化之后的x的数据，记为x.tilde
  for (i in 1:n) {
    y.tilde[[i]]<-y[[i]]-km.y
  }
  #估计x的系数beta
  Q.beta<-matrix(0:0,ncol=2,nrow=2)#求解系数的2×2矩阵
  q.beta<-matrix(0:0,ncol = 1,nrow = 2)#求解系数的2×1矩阵
  for (i in 1:n) {
    q.beta<-q.beta+matrix(c(sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]])/n,
                            sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]]*(T.x[[i]]-t))/n),ncol = 1,nrow = 2)
    Q.beta<-Q.beta+matrix(c(sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]])/n,
                            sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                            sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                            sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n),
                          ncol = 2,nrow = 2)
  }
  beta<-ginv(Q.beta)%*%q.beta#l代表第几次模拟
  sig_beta<-matrix(0:0,ncol = 2,nrow = 2)
  for (i in 1:n) {
    sig_beta[1,1]<-sig_beta[1,1]+sum((kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta[1,1]-x.tilde[[i]]*beta[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta[1,1]-x.tilde[[i]]*beta[2,1]*(T.x[[i]]-t))))/n
    sig_beta[1,2]<-sig_beta[1,2]+sum((kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta[1,1]-x.tilde[[i]]*beta[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(y.tilde[[i]]-x.tilde[[i]]*beta[1,1]-x.tilde[[i]]*beta[2,1]*(T.x[[i]]-t))))/n
    sig_beta[2,2]<-sig_beta[2,2]+sum((kw[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(y.tilde[[i]]-x.tilde[[i]]*beta[1,1]-x.tilde[[i]]*beta[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(y.tilde[[i]]-x.tilde[[i]]*beta[1,1]-x.tilde[[i]]*beta[2,1]*(T.x[[i]]-t))))/n
  }
  sig_beta[2,1]<-sig_beta[1,2]
  hbe<-ginv(Q.beta)%*%sig_beta%*%ginv(Q.beta)
  hb<-hbe[1,1]
  return(list("beta"=beta,"hb"=hb))
}
#gamma中心化估计函数#
gamma.est<-function(x,z,y,T.x,S.z,n,nx,nz,t,beta1,beta2,h){
  ker<-list(matrix)
  ker_xy<-list(matrix)
  for (i in 1:n) {
    ker[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_xy[[i]]<-matrix(rep((y[[i]]-x[[i]]*(beta1+beta2*(T.x[[i]]-t)))*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
  }
  Q.gamma<-matrix(0:0,ncol = 2,nrow = 2)
  q.gamma<-matrix(0:0,ncol = 1,nrow = 2)
  for (i in 1:n) {
    Q.gamma[1,1]<-Q.gamma[1,1]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*z[[i]])))/n
    Q.gamma[1,2]<-Q.gamma[1,2]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
    Q.gamma[2,2]<-Q.gamma[2,2]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
    q.gamma[1,]<-q.gamma[1,]+colSums(ker_xy[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]])))/n
    q.gamma[2,]<-q.gamma[2,]+colSums(ker_xy[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*(S.z[[i]]-t))))/n
  }
  Q.gamma[2,1]<-Q.gamma[1,2]
  gamma<-ginv(Q.gamma)%*%q.gamma
  p<-matrix(0:0,ncol = 2,nrow = 2)
  for (i in 1:n) {
    p1<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p1<-t(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta1+beta2*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma[1,1]+gamma[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p2<-matrix(0, ncol=nx[i], nrow=nz[i])
    p3<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p3<-t(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)*(S.z[[i]]-t)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta1+beta2*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma[1,1]+gamma[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p4<-matrix(0, ncol=nx[i], nrow=nz[i])
    for (j in 1:nx[i]) {
      for (k in 1:nz[i]) {
        if(is.matrix(p1[,,j,k]==TRUE)) {p2[k,j]<-p1[,,j,k][k,j]}
        else {p2[k,j]<-p1[,,j,k][k]}
        if(is.matrix(p3[,,j,k]==TRUE)) {p4[k,j]<-p3[,,j,k][k,j]}
        else {p4[k,j]<-p3[,,j,k][k]}
      }
    }
    p[1,1]<-p[1,1]+sum(p2%o%p2)/n
    p[1,2]<-p[1,2]+sum(p2%o%p4)/n
    p[2,2]<-p[2,2]+sum(p4%o%p4)/n
  }
  p[2,1]<-p[1,2]
  hgam<-ginv(Q.gamma)%*%p%*%ginv(Q.gamma)
  hg<-hgam[1,1]
  return(list("gamma"=gamma,"hg"=hg))
}
#beta中心化均方预测误差#
m.beta<-function(x,z,y,T.x,S.z,n,nx,nz,t,h,beta1,beta2){
  kw<-list()
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  #先将所有核权重求和
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  km.x<-0#代表所有数据x的核估计
  for (i in 1:n) {
    km.x<-km.x+sum(x[[i]]*kw[[i]]/kw.all)#每个个体的均值核估计
  }
  #中心化之后的x的数据，记为x.tilde
  for (i in 1:n) {
    x.tilde[[i]]<-x[[i]]-km.x
  }
  #得到数据y的核估计，这里应该是所有个体放在一起的核估计吧？
  km.y<-0#代表所有数据y的核估计
  for (i in 1:n) {
    km.y<-km.y+sum(y[[i]]*kw[[i]])/kw.all#每个个体的均值核估计
  }
  #中心化之后的x的数据，记为x.tilde
  for (i in 1:n) {
    y.tilde[[i]]<-y[[i]]-km.y
  }
  betam<-0
  for (i in 1:n) {
    betam<-betam+sum(kw[[i]]*((y.tilde[[i]]-x.tilde[[i]]*(beta1+beta2*(T.x[[i]]-t)))^2))/kw.all
  }
  return(betam)
}
#gamma中心化均方预测误差#
m.gamma<-function(x,z,y,T.x,S.z,n,nx,nz,t,h,beta1,beta2,gamma1,gamma2){
  gammam1<-0
  gammam2<-0
  gammam<-0
  for (i in 1:n) {
    gammam1<-gammam1+sum(diag(t(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)))%*%((matrix(rep(y[[i]]-x[[i]]*(beta1+beta2*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma1+gamma2*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))^2)))
    gammam2<-gammam2+sum(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)))
  }
  gammam<-gammam1/gammam2
  return(gammam)
}

#betagamma同时估计
fun.betagamma<-function(x,z,y,T.x,S.z,n,nx,nz,t,h){
  ker<-list(matrix)
  ker_x1<-list(matrix)
  ker_x2<-list(matrix)
  ker_x1t1<-list(matrix)
  ker_x2t1<-list(matrix)
  ker_x2t2<-list(matrix)
  ker_y<-list(matrix)
  ker_xy<-list(matrix)
  ker_xyt<-list(matrix)
  
  for (i in 1:n) {
    ker[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_x1[[i]]<-matrix(rep(x[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_x2[[i]]<-matrix(rep(x[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_x1t1[[i]]<-matrix(rep((T.x[[i]]-t)*x[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_x2t1[[i]]<-matrix(rep((T.x[[i]]-t)*x[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_x2t2[[i]]<-matrix(rep((T.x[[i]]-t)*(T.x[[i]]-t)*x[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    
    ker_y[[i]]<-matrix(rep(y[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_xy[[i]]<-matrix(rep(y[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_xyt[[i]]<-matrix(rep(y[[i]]*x[[i]]*(T.x[[i]]-t)*local_kernel(T.x[[i]]-t,h),nz[i]),ncol = nz[i],nrow = nx[i])
  }
  
  Q<-matrix(0:0,ncol = 4,nrow = 4)
  q<-matrix(0:0,ncol = 1,nrow = 4)
  
  for (i in 1:n) {
    Q[1,1]<-Q[1,1]+colSums(ker_x2[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h))))/n
    Q[1,2]<-Q[1,2]+colSums(ker_x2t1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h))))/n
    Q[1,3]<-Q[1,3]+colSums(ker_x1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]])))/n
    Q[1,4]<-Q[1,4]+colSums(ker_x1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*(S.z[[i]]-t))))/n
    Q[2,2]<-Q[2,2]+colSums(ker_x2t2[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h))))/n
    Q[2,3]<-Q[2,3]+colSums(ker_x1t1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]])))/n
    Q[2,4]<-Q[2,4]+colSums(ker_x1t1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*(S.z[[i]]-t))))/n
    Q[3,3]<-Q[3,3]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*z[[i]])))/n
    Q[3,4]<-Q[3,4]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
    Q[4,4]<-Q[4,4]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
    
    q[1,]<-q[1,]+colSums(ker_xy[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h))))/n
    q[2,]<-q[2,]+colSums(ker_xyt[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h))))/n
    q[3,]<-q[3,]+colSums(ker_y[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]])))/n
    q[4,]<-q[4,]+colSums(ker_y[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h)*z[[i]]*(S.z[[i]]-t))))/n
    
  }
  Q[2,1]<-Q[1,2]
  Q[3,1]<-Q[1,3]
  Q[4,1]<-Q[1,4]
  Q[3,2]<-Q[2,3]
  Q[4,2]<-Q[2,4]
  Q[4,3]<-Q[3,4]
  betagamma<-ginv(Q)%*%q#l代表第几次模拟
  
  p<-matrix(0:0,ncol = 4,nrow = 4)
  for (i in 1:n) {
    p11<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p11<-t(as.vector(local_kernel(T.x[[i]]-t,h)*x[[i]])%o%as.vector(local_kernel(S.z[[i]]-t,h)))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma[1,1]+betagamma[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma[3,1]+betagamma[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p1<-matrix(0, ncol=nx[i], nrow=nz[i])
    p22<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p22<-t(as.vector(local_kernel(T.x[[i]]-t,h)*(T.x[[i]]-t)*x[[i]])%o%as.vector(local_kernel(S.z[[i]]-t,h)))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma[1,1]+betagamma[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma[3,1]+betagamma[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p2<-matrix(0, ncol=nx[i], nrow=nz[i])
    p33<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p33<-t(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma[1,1]+betagamma[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma[3,1]+betagamma[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p3<-matrix(0, ncol=nx[i], nrow=nz[i])
    p44<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p44<-t(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)*(S.z[[i]]-t)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma[1,1]+betagamma[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma[3,1]+betagamma[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p4<-matrix(0, ncol=nx[i], nrow=nz[i])
    for (j in 1:nx[i]) {
      for (k in 1:nz[i]) {
        if(is.matrix(p11[,,j,k]==TRUE)) {p1[k,j]<-p11[,,j,k][k,j]}
        else {p1[k,j]<-p11[,,j,k][k]}
        if(is.matrix(p22[,,j,k]==TRUE)) {p2[k,j]<-p22[,,j,k][k,j]}
        else {p2[k,j]<-p22[,,j,k][k]}
        if(is.matrix(p33[,,j,k]==TRUE)) {p3[k,j]<-p33[,,j,k][k,j]}
        else {p3[k,j]<-p33[,,j,k][k]}
        if(is.matrix(p44[,,j,k]==TRUE)) {p4[k,j]<-p44[,,j,k][k,j]}
        else {p4[k,j]<-p44[,,j,k][k]}
      }
    }
    p[1,1]<-p[1,1]+sum(p1%o%p1)/n
    p[2,2]<-p[2,2]+sum(p2%o%p2)/n
    p[3,3]<-p[3,3]+sum(p3%o%p3)/n
    p[4,4]<-p[4,4]+sum(p4%o%p4)/n
    p[1,2]<-p[1,2]+sum(p1%o%p2)/n
    p[1,3]<-p[1,3]+sum(p1%o%p3)/n
    p[1,4]<-p[1,4]+sum(p1%o%p4)/n
    p[2,3]<-p[2,3]+sum(p2%o%p3)/n
    p[2,4]<-p[2,4]+sum(p2%o%p4)/n
    p[3,4]<-p[3,4]+sum(p3%o%p4)/n
  }
  p[2,1]<-p[1,2]
  p[3,1]<-p[1,3]
  p[4,1]<-p[1,4]
  p[3,2]<-p[2,3]
  p[4,2]<-p[2,4]
  p[4,3]<-p[3,4]
  hatvar<-ginv(Q)%*%p%*%ginv(Q)
  return(list("betagamma"=betagamma,"hatvar"=hatvar))
}

#betagamma均方预测误差#
m.betagamma<-function(x,z,y,T.x,S.z,n,nx,nz,t,h,beta1,beta2,gamma1,gamma2){
  m1<-0
  m2<-0
  m<-0
  for (i in 1:n) {
    m1<-m1+sum(diag(t(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)))%*%((matrix(rep(y[[i]]-x[[i]]*(beta1+beta2*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma1+gamma2*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))^2)))
    m2<-m2+sum(as.vector(local_kernel(T.x[[i]]-t,h))%o%as.vector(local_kernel(S.z[[i]]-t,h)))
  }
  m<-m1/m2
  return(m)
}


for (h.ind.gamma in 1:length(h.can.gamma)) {
  print(paste(h.ind.gamma))
  for (l in 1:K) {
    T.x<-list()#代表X和Y的观测时间，代表每个subject的不等长度观测时间向量的集合，可以为任何数据类型向量矩阵数组，不同长度的列表
    S.z<-list()#代表实际Z的n个不等长度观测时间向量的集合,但是生成Y时，需要一起使用XY的观测时间去生成
    kw<-list()
    kw1<-list()
    kw2<-list()
    mu_x<-list()#x多元正态均值，n个不等长度向量的合集
    sigma_x<-list(matrix)#x多元正态sigma，n个矩阵的合集
    mu_z.x<-list()#z.x多元正态均值，n个不等长度向量的合集
    sigma_z.x<-list(matrix)#z.x多元正态sigma，n个矩阵的合集
    mu_z<-list()#z多元正态均值，n个不等长度向量的合集
    sigma_z<-list(matrix)#z多元正态sigma，n个矩阵的合集
    mu_epsilon<-list()
    sigma_epsilon<-list(matrix)
    x<-list()
    x.tilde<-list()#代表中心化之后的数据x
    z.x<-list()#代表按x观测时间生成的数据z，用来生成数据y使用
    z<-list()#代表实际按z观测时间生成的数据z，用来实际模拟使用
    y<-list()
    y.tilde<-list()#代表中心化之后的数据y
    epsilon<-list()
    nx<-rpois(n,5)+1#the number of each subject,由每个观测对象的观测次数形成的一个向量
    nz<-rpois(n,5)+1
    nx1<-rep(0,n/2)
    nz1<-rep(0,n/2)
    nx2<-rep(0,n/2)
    nz2<-rep(0,n/2)
    x1<-list()
    z1<-list()#代表实际按z观测时间生成的数据z，用来实际模拟使用
    y1<-list()
    x2<-list()
    z2<-list()#代表实际按z观测时间生成的数据z，用来实际模拟使用
    y2<-list()
    T.x1<-list()
    S.z1<-list()
    T.x2<-list()
    S.z2<-list()
    x.tilde1<-list()#代表中心化之后的数据x
    x.tilde2<-list()#代表中心化之后的数据x
    y.tilde1<-list()#代表中心化之后的数据y
    y.tilde2<-list()#代表中心化之后的数据y
    #定义协变量x和z的观测时间
    for (i in 1:n) {
      T.x[[i]]<-as.matrix(runif(nx[i]))
      S.z[[i]]<-as.matrix(runif(nz[i]))
      t.temp=rbind(T.x[[i]],S.z[[i]])
      n.temp=nx[i]+nz[i]
      corr=exp(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
      corr.e=2^(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
      #corr.x[[i]]=exp(-abs(rep(1,nx[i])%*%t(tx[[i]])-tx[[i]]%*%t(rep(1,nx[i]))))
      corr.v=exp(-abs(rep(1,n.temp)%*%t(t.temp)-t.temp%*%t(rep(1,n.temp))))
      v=mvrnorm(1,rep(0,n.temp),corr.v)
      
      EX=rep(0,n.temp)
      EZ=2*(t.temp-0.5)^2
      #rep(2,n.temp),0.5+t.temp,2*sin(2*pi*t.temp)
      #按不独立，但不相关生成X，Z，epsilon
      #x.temp=rnorm(1,0,1)*v+EX
      #z.temp=v+EZ
      #epsilon.temp=as.matrix(rnorm(1,0,1)*v)
      x.temp=mvrnorm(1,EX,corr)
      z.temp=mvrnorm(1,EZ,corr)
      epsilon.temp=as.matrix(mvrnorm(1,rep(0,n.temp),corr.e))
      
      x[[i]]=as.matrix(x.temp[1:nx[i]])
      z[[i]]=as.matrix(z.temp[-(1:nx[i])])
      y.temp=x.temp*fun.beta(t.temp)+z.temp*fun.gamma(t.temp)+epsilon.temp
      y[[i]]=as.matrix(y.temp[1:nx[i]])
    }
    x1<-list()
    z1<-list()
    y1<-list()
    nx1<-list()
    T.x1<-list()
    nz1<-list()
    S.z1<-list()
    V<-2#分成V组#
    gr=sample(rep(seq(1,V,1), length =n ))
    for (v in 1:V) {
      x1[[v]]<-x[gr==v]
      z1[[v]]<-z[gr==v]
      y1[[v]]<-y[gr==v]
      nx1[[v]]<-nx[gr==v]
      T.x1[[v]]<-T.x[gr==v]
      nz1[[v]]<-nz[gr==v]
      S.z1[[v]]<-S.z[gr==v]
    }
    
    #同时估betagamma2019
    h<-h.can.gamma[h.ind.gamma]
    bg.s_2019<-fun.betagamma(x,z,y,T.x,S.z,n,nx,nz,t,h)
    betagamma.hat_2019[h.ind.gamma,,l]<-t(bg.s_2019$betagamma)
    hvarb_2019[h.ind.gamma,l]<-(bg.s_2019$hatvar)[1,1]
    hvarg_2019[h.ind.gamma,l]<-(bg.s_2019$hatvar)[3,3]
    bg.s1_2019<-fun.betagamma(x1[[1]],z1[[1]],y1[[1]],T.x1[[1]],S.z1[[1]],n/2,nx1[[1]],nz1[[1]],t,h)
    bg.s2_2019<-fun.betagamma(x1[[2]],z1[[2]],y1[[2]],T.x1[[2]],S.z1[[2]],n/2,nx1[[2]],nz1[[2]],t,h)
    betagamma.hat1_2019[h.ind.gamma,,l]<-t(bg.s1_2019$betagamma)
    betagamma.hat2_2019[h.ind.gamma,,l]<-t(bg.s2_2019$betagamma)
    mbetagamma_2019[l,h.ind.gamma]<-m.betagamma(x,z,y,T.x,S.z,n,nx,nz,t,h,betagamma.hat_2019[h.ind.gamma,1,l],betagamma.hat_2019[h.ind.gamma,2,l],betagamma.hat_2019[h.ind.gamma,3,l],betagamma.hat_2019[h.ind.gamma,4,l])
    #去掉某一份
    mbetagamma1_2019[l,h.ind.gamma]<-m.betagamma(x1[[1]],z1[[1]],y1[[1]],T.x1[[1]],S.z1[[1]],n/2,nx1[[1]],nz1[[1]],t,h,betagamma.hat2_2019[h.ind.gamma,1,l],betagamma.hat2_2019[h.ind.gamma,2,l],betagamma.hat2_2019[h.ind.gamma,3,l],betagamma.hat2_2019[h.ind.gamma,4,l])
    mbetagamma2_2019[l,h.ind.gamma]<-m.betagamma(x1[[2]],z1[[2]],y1[[2]],T.x1[[2]],S.z1[[2]],n/2,nx1[[2]],nz1[[2]],t,h,betagamma.hat1_2019[h.ind.gamma,1,l],betagamma.hat1_2019[h.ind.gamma,2,l],betagamma.hat1_2019[h.ind.gamma,3,l],betagamma.hat1_2019[h.ind.gamma,4,l])
    
  }
}


#同时估2019
rm_2019<-rep(0,length(h.can.gamma))
for (l in 1:K) {
  rm_2019<-rm_2019+(mbetagamma1_2019[l,]/K+mbetagamma2_2019[l,]/K)/2
}
h.best_2019<-which(rm_2019==min(rm_2019),arr.ind = T)

bemsame_2019<-matrix(0:0,nrow = K,ncol = length(h.can.gamma))
gamsame_2019<-matrix(0:0,nrow = K,ncol = length(h.can.gamma))
for (l in 1:K) {
  bemsame_2019[l,]<-betagamma.hat_2019[,1,l]
  gamsame_2019[l,]<-betagamma.hat_2019[,3,l]
}
bias.betasame_2019<-colMeans(bemsame_2019)-fun.beta(t)
bias.gammasame_2019<-colMeans(gamsame_2019)-fun.gamma(t)
mse.betasame_2019<-colMeans(bemsame_2019-fun.beta(t))^2
mse.gammasame_2019<-colMeans(gamsame_2019-fun.gamma(t))^2
SD.betasame_2019<-rep(0,length(h.can.gamma))
SD.gammasame_2019<-rep(0,length(h.can.gamma))
for (h.ind in 1:length(h.can.gamma)) {
  SD.betasame_2019[h.ind]<-sqrt(sum((bemsame_2019[,h.ind]-mean(bemsame_2019[,h.ind]))^2)/(K-1))
  SD.gammasame_2019[h.ind]<-sqrt(sum((gamsame_2019[,h.ind]-mean(gamsame_2019[,h.ind]))^2)/(K-1))
}
SE.betasame_2019<-sqrt(rowMeans(hvarb_2019))/sqrt(n)
SE.gammasame_2019<-sqrt(rowMeans(hvarg_2019))/sqrt(n)

Lbsame_2019<-matrix(0,ncol=length(h.can.gamma),nrow = K)
Rbsame_2019<-matrix(0,ncol=length(h.can.gamma),nrow = K)
Lgsame_2019<-matrix(0,ncol=length(h.can.gamma),nrow = K)
Rgsame_2019<-matrix(0,ncol=length(h.can.gamma),nrow = K)
cumbsame_2019<-rep(0,length(h.can.gamma))
cumgsame_2019<-rep(0,length(h.can.gamma))

for (l in 1:K) {
  for(k in 1:length(h.can.gamma)){
    Lbsame_2019[l,k]<-betagamma.hat_2019[k,1,l]-qnorm(0.975)*sqrt(hvarb_2019[k,l])/sqrt(n)
    Rbsame_2019[l,k]<-betagamma.hat_2019[k,1,l]+qnorm(0.975)*sqrt(hvarb_2019[k,l])/sqrt(n)
    Lgsame_2019[l,k]<-betagamma.hat_2019[k,3,l]-qnorm(0.975)*sqrt(hvarg_2019[k,l])/sqrt(n)
    Rgsame_2019[l,k]<-betagamma.hat_2019[k,3,l]+qnorm(0.975)*sqrt(hvarg_2019[k,l])/sqrt(n)
    if(fun.beta(t)<=Rbsame_2019[l,k] & fun.beta(t)>=Lbsame_2019[l,k]) 
    {cumbsame_2019[k]<-cumbsame_2019[k]+1}
    if(fun.gamma(t)<=Rgsame_2019[l,k] & fun.gamma(t)>=Lgsame_2019[l,k]) 
    {cumgsame_2019[k]<-cumgsame_2019[k]+1}
  }
}

fun.beta
n;t
h.best_2019
bias.betasame_2019;mse.betasame_2019;SD.betasame_2019;SE.betasame_2019;cumbsame_2019/K;colMeans(Rbsame_2019);colMeans(Lbsame_2019)
bias.gammasame_2019;mse.gammasame_2019;SD.gammasame_2019;SE.gammasame_2019;cumgsame_2019/K;colMeans(Rgsame_2019);colMeans(Lgsame_2019)
rm_2019
#Rbsame_2019;Lbsame_2019
#Rgsame_2019;Lgsame_2019
#bemsame_2019
#gamsame_2019
#hvarb_2019;hvarg_2019



result<-rbind(Rbsame_2019,Lbsame_2019,bemsame_2019,Rgsame_2019,Lgsame_2019,gamsame_2019,t(hvarb_2019),t(hvarg_2019))
write.csv(result,file="E:/papercao/2020.12.14_datafirst_1K2t/auto0.3kw_sin_n400_ez2t_2019.csv")

#result<-rbind(bem_p_2019,bem_p_2015,bem_2019,bem_2015,bemsame_2019,bemsame_2015,
#gam_p_2019,gam_p_2015,gam_2019,gam_2015,gamsame_2019,gamsame_2015,
#htvar_beta_p_2019,htvar_beta_p_2015,htvar_gamma_p_2019[locationbeta_p_2019,,],htvar_gamma_p_2015[locationbeta_p_2015,,],
#htvar_beta_2019,htvar_beta_2015,htvar_gamma_2019[locationbeta_2019,,],htvar_gamma_2015[locationbeta_2015,,],
#hvarb_2019,hvarb_2015,hvarg_2019,hvarg_2015)
#write.csv(result,file="E:/papercao/2020.8.17跑程序/auto0.25line_independent_n200.csv")



#plot(h.can.beta,rmbeta,type="o",xlim=c(0.00,0.35),xaxt="n",xlab="Bandwidth",ylab="Squred prediction error")
#axis(1,at=seq(0.00,0.35,length=8),labels = seq(0.00,0.35,length=8),cex.axis=0.8)

#plot(h.can.gamma,rmgamma,type="o",xlim=c(0.05,0.40),xaxt="n",xlab="Bandwidth",ylab="Squred prediction error")
#axis(1,at=seq(0.05,0.40,length=8),labels = seq(0.05,0.40,length=8),cex.axis=0.8)
