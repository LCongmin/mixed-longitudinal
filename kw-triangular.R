###############################################################################################
###     带方差CP,分开同时都有估计，固定窗宽h=n^(-0.8),h1=h2=n^(-0.5),t=(0.05,0.95)          ###
###############################################################################################
rm(list = ls (all = TRUE))
yxt<-proc.time()
library(MASS)
p<-seq(0.3,0.9,0.3)


betat.hat_kw<-rep(0,length(p))
bias.beta_kw<-rep(0,length(p))
mse.beta_kw<-rep(0,length(p))
SD.beta_kw<-rep(0,length(p))
SE.beta_kw<-rep(0,length(p))
cum_beta_kw1<-rep(0,length(p))
hv_beta_kw<-rep(0,length(p))
L_beta_kw<-rep(0,length(p))
R_beta_kw<-rep(0,length(p))


gammat.hat_kw<-rep(0,length(p))
bias.gamma_kw<-rep(0,length(p))
mse.gamma_kw<-rep(0,length(p))
SD.gamma_kw<-rep(0,length(p))
SE.gamma_kw<-rep(0,length(p))
cum_gamma_kw1<-rep(0,length(p))
hv_gamma_kw<-rep(0,length(p))
L_gamma_kw<-rep(0,length(p))
R_gamma_kw<-rep(0,length(p))


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

epa <- function(t){
  
  tst <- (-1.0 <= t) & (t <= 1.0 )
  
  kt <- (0.75*( 1.0 - t*t ))^2
  kt[ !tst ] <- 0.0
  
  return(kt)
}
#integrate(epa,-Inf,Inf)

uniform <- function(t){
  
  tst <- (-1.0 <= t) & (t <= 1.0 )
  kt <- t
  kt[  tst ] <- 0.5
  kt[ !tst ] <- 0.0
  
  return(kt)
}


triangular=function(t){
  tst <- (-1.0 <= t) & (t <= 1.0 )
  kt <- 1-abs(t)
  kt[ !tst ] <- 0.0
  return(kt)
}
local_kernel <- function(t, h){
  kt <- triangular(t/h)
  return(kt/h)
}

n<-400
K<-1000#循环次数
h1_kw<-n^(-0.45)
h2_kw<-n^(-0.45)
#t<-0.25
#t<-0.5
#t<-0.75



#存储每个时间点的每一次模拟结果91*K的矩阵

#核加权方法下beta的估计
hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_kw<-matrix(0,ncol = K,nrow = length(p))#beta方差
hatbetadot_var_kw<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta95%上界

#核加权方法下gamma的估计
hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma估计值
hatgamma_var_kw<-matrix(0,ncol = K,nrow = length(p))#gamma方差
hatgammadot_var_kw<-matrix(0,ncol = K,nrow = length(p))#gamma方差
L_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界

for (l in 1:K) {
  
  print(paste(l))
  
  
  T.x<-list()#代表X和Y的观测时间，代表每个subject的不等长度观测时间向量的集合，可以为任何数据类型向量矩阵数组，不同长度的列表
  S.z<-list()#代表实际Z的n个不等长度观测时间向量的集合,但是生成Y时，需要一起使用XY的观测时间去生成
  kw<-list()
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
  v<-list()
  nx<-rpois(n,5)+1#the number of each subject,由每个观测对象的观测次数形成的一个向量
  nz<-rpois(n,5)+1
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
  
  
  for (t_ind in 1:length(p)) {
    t<-0.3+0.3*(t_ind-1)
    
    #核加权同时估计系数beta,gamma
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
      ker[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_x1[[i]]<-matrix(rep(x[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_x2[[i]]<-matrix(rep(x[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_x1t1[[i]]<-matrix(rep((T.x[[i]]-t)*x[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_x2t1[[i]]<-matrix(rep((T.x[[i]]-t)*x[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_x2t2[[i]]<-matrix(rep((T.x[[i]]-t)*(T.x[[i]]-t)*x[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      
      ker_y[[i]]<-matrix(rep(y[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_xy[[i]]<-matrix(rep(y[[i]]*x[[i]]*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_xyt[[i]]<-matrix(rep(y[[i]]*x[[i]]*(T.x[[i]]-t)*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
    }
    
    Q<-matrix(0:0,ncol = 4,nrow = 4)
    q<-matrix(0:0,ncol = 1,nrow = 4)
    
    for (i in 1:n) {
      Q[1,1]<-Q[1,1]+colSums(ker_x2[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw))))/n
      Q[1,2]<-Q[1,2]+colSums(ker_x2t1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw))))/n
      Q[1,3]<-Q[1,3]+colSums(ker_x1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]])))/n
      Q[1,4]<-Q[1,4]+colSums(ker_x1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*(S.z[[i]]-t))))/n
      Q[2,2]<-Q[2,2]+colSums(ker_x2t2[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw))))/n
      Q[2,3]<-Q[2,3]+colSums(ker_x1t1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]])))/n
      Q[2,4]<-Q[2,4]+colSums(ker_x1t1[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*(S.z[[i]]-t))))/n
      Q[3,3]<-Q[3,3]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]])))/n
      Q[3,4]<-Q[3,4]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
      Q[4,4]<-Q[4,4]+colSums(ker[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
      
      q[1,]<-q[1,]+colSums(ker_xy[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw))))/n
      q[2,]<-q[2,]+colSums(ker_xyt[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw))))/n
      q[3,]<-q[3,]+colSums(ker_y[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]])))/n
      q[4,]<-q[4,]+colSums(ker_y[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*(S.z[[i]]-t))))/n
      
    }
    Q[2,1]<-Q[1,2]
    Q[3,1]<-Q[1,3]
    Q[4,1]<-Q[1,4]
    Q[3,2]<-Q[2,3]
    Q[4,2]<-Q[2,4]
    Q[4,3]<-Q[3,4]
    betagamma.hat<-ginv(Q)%*%q#l代表第几次模拟
    hatbeta_kw[t_ind,l]<-betagamma.hat[1,1]
    hatgamma_kw[t_ind,l]<-betagamma.hat[3,1]
    
    
    
    #sig<-matrix(0:0,ncol = 4,nrow = 4)
    P<-matrix(0:0,ncol = 4,nrow = 4)
    for (i in 1:n) {
      p11<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p11<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw)*x[[i]])%o%as.vector(local_kernel(S.z[[i]]-t,h2_kw)))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma.hat[1,1]+betagamma.hat[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma.hat[3,1]+betagamma.hat[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p1<-matrix(0, ncol=nx[i], nrow=nz[i])
      p22<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p22<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw)*(T.x[[i]]-t)*x[[i]])%o%as.vector(local_kernel(S.z[[i]]-t,h2_kw)))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma.hat[1,1]+betagamma.hat[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma.hat[3,1]+betagamma.hat[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p2<-matrix(0, ncol=nx[i], nrow=nz[i])
      p33<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p33<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw))%o%as.vector(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma.hat[1,1]+betagamma.hat[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma.hat[3,1]+betagamma.hat[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p3<-matrix(0, ncol=nx[i], nrow=nz[i])
      p44<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p44<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw))%o%as.vector(local_kernel(S.z[[i]]-t,h2_kw)*(S.z[[i]]-t)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(betagamma.hat[1,1]+betagamma.hat[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(betagamma.hat[3,1]+betagamma.hat[4,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
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
      P[1,1]<-P[1,1]+sum(p1%o%p1)/n
      P[2,2]<-P[2,2]+sum(p2%o%p2)/n
      P[3,3]<-P[3,3]+sum(p3%o%p3)/n
      P[4,4]<-P[4,4]+sum(p4%o%p4)/n
      P[1,2]<-P[1,2]+sum(p1%o%p2)/n
      P[1,3]<-P[1,3]+sum(p1%o%p3)/n
      P[1,4]<-P[1,4]+sum(p1%o%p4)/n
      P[2,3]<-P[2,3]+sum(p2%o%p3)/n
      P[2,4]<-P[2,4]+sum(p2%o%p4)/n
      P[3,4]<-P[3,4]+sum(p3%o%p4)/n
    }
    P[2,1]<-P[1,2]
    P[3,1]<-P[1,3]
    P[4,1]<-P[1,4]
    P[3,2]<-P[2,3]
    P[4,2]<-P[2,4]
    P[4,3]<-P[3,4]
    hatvar<-ginv(Q)%*%P%*%ginv(Q)
    hatbeta_var_kw[t_ind,l]<-hatvar[1,1]
    hatgamma_var_kw[t_ind,l]<-hatvar[3,3]
    hatbetadot_var_kw[t_ind,l]<-hatvar[2,2]
    hatgammadot_var_kw[t_ind,l]<-hatvar[4,4]
    
    
    #正确模型下核加权估计
    L_hatbeta_kw[t_ind,l]<-hatbeta_kw[t_ind,l]-qnorm(0.975)*sqrt(hatbeta_var_kw[t_ind,l])/sqrt(n)
    R_hatbeta_kw[t_ind,l]<-hatbeta_kw[t_ind,l]+qnorm(0.975)*sqrt(hatbeta_var_kw[t_ind,l])/sqrt(n)
    L_hatgamma_kw[t_ind,l]<-hatgamma_kw[t_ind,l]-qnorm(0.975)*sqrt(hatgamma_var_kw[t_ind,l])/sqrt(n)
    R_hatgamma_kw[t_ind,l]<-hatgamma_kw[t_ind,l]+qnorm(0.975)*sqrt(hatgamma_var_kw[t_ind,l])/sqrt(n)
    
    if(fun.beta(t)<=R_hatbeta_kw[t_ind,l] & fun.beta(t)>=L_hatbeta_kw[t_ind,l]) 
    {cum_beta_kw1[t_ind]<-cum_beta_kw1[t_ind]+1}
    
    
    if(fun.gamma(t)<=R_hatgamma_kw[t_ind,l] & fun.gamma(t)>=L_hatgamma_kw[t_ind,l]) 
    {cum_gamma_kw1[t_ind]<-cum_gamma_kw1[t_ind]+1}
    
  } 
}


#正确模型下核加权估计
hv_beta_kw<-rowMeans(hatbeta_var_kw)
hv_gamma_kw<-rowMeans(hatgamma_var_kw)

betat.hat_kw<-rowMeans(hatbeta_kw)#p代表第几个t
gammat.hat_kw<-rowMeans(hatgamma_kw)

SE.beta_kw<-sqrt(hv_beta_kw)/sqrt(n)
SE.gamma_kw<-sqrt(hv_gamma_kw)/sqrt(n)

for (t_ind in 1:length(p)) {
  t<-0.3+0.3*(t_ind-1)
  bias.beta_kw[t_ind]<-mean(hatbeta_kw[t_ind,]-fun.beta(t))
  mse.beta_kw[t_ind]<-mean((hatbeta_kw[t_ind,]-fun.beta(t))^2)
  SD.beta_kw[t_ind]<-sqrt(sum((hatbeta_kw[t_ind,]-betat.hat_kw[t_ind])^2)/(K-1))
  
  bias.gamma_kw[t_ind]<-mean(hatgamma_kw[t_ind,]-fun.gamma(t))
  mse.gamma_kw[t_ind]<-mean((hatgamma_kw[t_ind,]-fun.gamma(t))^2)
  SD.gamma_kw[t_ind]<-sqrt(sum((hatgamma_kw[t_ind,]-gammat.hat_kw[t_ind])^2)/(K-1))
  
}

L_beta_kw_point<-rowMeans(L_hatbeta_kw)
R_beta_kw_point<-rowMeans(R_hatbeta_kw)

L_gamma_kw_point<-rowMeans(L_hatgamma_kw)
R_gamma_kw_point<-rowMeans(R_hatgamma_kw)

fun.beta
n^{-0.5};n^{-0.5}
n;h1_kw;h2_kw
cum_beta_kw1/K
cum_gamma_kw1/K
betat.hat_kw
bias.beta_kw
mse.beta_kw
SD.beta_kw
SE.beta_kw
gammat.hat_kw
bias.gamma_kw
mse.gamma_kw
SD.gamma_kw
SE.gamma_kw
L_beta_kw_point
R_beta_kw_point
L_gamma_kw_point
R_gamma_kw_point

proc.time()-yxt

