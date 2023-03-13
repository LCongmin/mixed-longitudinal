###############################################################################################
###     带方差CP,分开同时都有估计，固定窗宽h=n^(-0.8),h1=h2=n^(-0.5),t=(0.05,0.95)          ###
###############################################################################################
rm(list = ls (all = TRUE))
yxt<-proc.time()
library(MASS)
p<-seq(0.3,0.9,0.3)
#部分线性下beta的估计
betat.hat_p<-rep(0,length(p))
bias.beta_p<-rep(0,length(p))
mse.beta_p<-rep(0,length(p))
SD.beta_p<-rep(0,length(p))
SE.beta_p<-rep(0,length(p))
cum_beta_p1<-rep(0,length(p))
hv_beta_p<-rep(0,length(p))
L_beta_p<-rep(0,length(p))
R_beta_p<-rep(0,length(p))

#部分线性gamma
gammat.hat_p<-rep(0,length(p))
bias.gamma_p<-rep(0,length(p))
mse.gamma_p<-rep(0,length(p))
SD.gamma_p<-rep(0,length(p))
SE.gamma_p<-rep(0,length(p))
cum_gamma_p1<-rep(0,length(p))
hv_gamma_p<-rep(0,length(p))
L_gamma_p<-rep(0,length(p))
R_gamma_p<-rep(0,length(p))

#中心化beta
betat.hat_center<-rep(0,length(p))
bias.beta_center<-rep(0,length(p))
mse.beta_center<-rep(0,length(p))
SD.beta_center<-rep(0,length(p))
SE.beta_center<-rep(0,length(p))
cum_beta_center1<-rep(0,length(p))
hv_beta_center<-rep(0,length(p))
L_beta_center<-rep(0,length(p))
R_beta_center<-rep(0,length(p))


#中心化gamma
gammat.hat_center<-rep(0,length(p))
bias.gamma_center<-rep(0,length(p))
mse.gamma_center<-rep(0,length(p))
SD.gamma_center<-rep(0,length(p))
SE.gamma_center<-rep(0,length(p))
cum_gamma_center1<-rep(0,length(p))
hv_gamma_center<-rep(0,length(p))
L_gamma_center<-rep(0,length(p))
R_gamma_center<-rep(0,length(p))

gammat.hat_lvcf<-rep(0,length(p))
bias.gamma_lvcf<-rep(0,length(p))
mse.gamma_lvcf<-rep(0,length(p))
SD.gamma_lvcf<-rep(0,length(p))
SE.gamma_lvcf<-rep(0,length(p))
cum_gamma_lvcf1<-rep(0,length(p))
L_gamma_lvcf<-rep(0,length(p))
R_gamma_lvcf<-rep(0,length(p))


#kw betagamma
betagammat.hat<-matrix(0:0,nrow=length(p),ncol=2)
betagammatdot.hat<-matrix(0:0,nrow=length(p),ncol=2)
bias.betagamma<-matrix(0:0,nrow=length(p),ncol=2)
mse.betagamma<-matrix(0:0,nrow=length(p),ncol=2)
SD.betagamma<-matrix(0:0,nrow=length(p),ncol=2)
SE.betagamma<-matrix(0:0,nrow=length(p),ncol=2)
hv<-matrix(0:0,ncol = 4,nrow = length(p))
L<-matrix(0:0,ncol = 2,nrow = length(p))
R<-matrix(0:0,ncol = 2,nrow = length(p))
cum_beta_kw1<-rep(0,length(p))
cum_gamma_kw1<-rep(0,length(p))


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
h<-n^{-0.6}
#hv=n^{-0.5}
h1_kw<-n^{-0.5}
h2_kw<-n^{-0.5}
#t<-0.25
#t<-0.5
#t<-0.75


#存储每个时间点的每一次模拟结果90*K的矩阵
#部分线性下beta
hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_p<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
#部分线性下gamma
hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatgamma_var_p<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%上界

#中心化方法下beta的估计
hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_center<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta95%上界

#中心化方法下gamma的估计
hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma估计值
hatgamma_var_center<-matrix(0,ncol = K,nrow = length(p))#gamma方差
L_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界

hatgamma_lvcf<-matrix(0,ncol = K,nrow = length(p))#gamma估计值
hatgamma_SE_lvcf<-matrix(0,ncol = K,nrow = length(p))#gamma方差
L_hatgamma_lvcf<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_lvcf<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界

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
    #kw代表核权重
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
    
    
    
    #中心化方法估计x的系数beta
    Q_center<-matrix(0:0,ncol=2,nrow=2)#求解系数的2×2矩阵
    q_center<-matrix(0:0,ncol = 1,nrow = 2)#求解系数的2×1矩阵
    for (i in 1:n) {
      q_center<-q_center+matrix(c(sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]])/n,
                                  sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]]*(T.x[[i]]-t))/n),ncol = 1,nrow = 2)
      Q_center<-Q_center+matrix(c(sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]])/n,
                                  sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                                  sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                                  sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n),
                                ncol = 2,nrow = 2)
    }
    beta.hat_center<-ginv(Q_center)%*%q_center#l代表第几次模拟
    hatbeta_center[t_ind,l]<-beta.hat_center[1,1]
    
    sig_beta_center<-matrix(0:0,ncol = 2,nrow = 2)
    for (i in 1:n) {
      sig_beta_center[1,1]<-sig_beta_center[1,1]+sum((kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta.hat_center[1,1]-x.tilde[[i]]*beta.hat_center[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta.hat_center[1,1]-x.tilde[[i]]*beta.hat_center[2,1]*(T.x[[i]]-t))))/n
      sig_beta_center[1,2]<-sig_beta_center[1,2]+sum((kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta.hat_center[1,1]-x.tilde[[i]]*beta.hat_center[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(y.tilde[[i]]-x.tilde[[i]]*beta.hat_center[1,1]-x.tilde[[i]]*beta.hat_center[2,1]*(T.x[[i]]-t))))/n
      sig_beta_center[2,2]<-sig_beta_center[2,2]+sum((kw[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(y.tilde[[i]]-x.tilde[[i]]*beta.hat_center[1,1]-x.tilde[[i]]*beta.hat_center[2,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(y.tilde[[i]]-x.tilde[[i]]*beta.hat_center[1,1]-x.tilde[[i]]*beta.hat_center[2,1]*(T.x[[i]]-t))))/n
    }
    sig_beta_center[2,1]<-sig_beta_center[1,2]
    hatvar_beta_center<-ginv(Q_center)%*%sig_beta_center%*%ginv(Q_center)
    hatbeta_var_center[t_ind,l]<-hatvar_beta_center[1,1]
    
    ker_center<-list(matrix)
    ker_xy_center<-list(matrix)
    for (i in 1:n) {
      ker_center[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_xy_center[[i]]<-matrix(rep((y[[i]]-x[[i]]*(beta.hat_center[1,1]+beta.hat_center[2,1]*(T.x[[i]]-t)))*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      #ker_xy[[i]]<-matrix(rep((y[[i]]-x[[i]]*fun.beta(T.x[[i]]))*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nz[i])
    }
    Q_center<-matrix(0:0,ncol = 2,nrow = 2)
    q_center<-matrix(0:0,ncol = 1,nrow = 2)
    for (i in 1:n) {
      Q_center[1,1]<-Q_center[1,1]+colSums(ker_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]])))/n
      Q_center[1,2]<-Q_center[1,2]+colSums(ker_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
      Q_center[2,2]<-Q_center[2,2]+colSums(ker_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
      q_center[1,]<-q_center[1,]+colSums(ker_xy_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]])))/n
      q_center[2,]<-q_center[2,]+colSums(ker_xy_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]*(S.z[[i]]-t))))/n
    }
    Q_center[2,1]<-Q_center[1,2]
    gamma.hat_center<-ginv(Q_center)%*%q_center
    hatgamma_center[t_ind,l]<-gamma.hat_center[1,1]
    
    p_center<-matrix(0:0,ncol = 2,nrow = 2)
    for (i in 1:n) {
      p1_center<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p1_center<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw))%o%as.vector(local_kernel(S.z[[i]]-t,h2_kw)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta.hat_center[1,1]+beta.hat_center[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma.hat_center[1,1]+gamma.hat_center[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p2_center<-matrix(0, ncol=nx[i], nrow=nz[i])
      p3_center<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p3_center<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw))%o%as.vector(local_kernel(S.z[[i]]-t,h2_kw)*(S.z[[i]]-t)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta.hat_center[1,1]+beta.hat_center[2,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma.hat_center[1,1]+gamma.hat_center[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p4_center<-matrix(0, ncol=nx[i], nrow=nz[i])
      for (j in 1:nx[i]) {
        for (k in 1:nz[i]) {
          if(is.matrix(p1_center[,,j,k]==TRUE)) {p2_center[k,j]<-p1_center[,,j,k][k,j]}
          else {p2_center[k,j]<-p1_center[,,j,k][k]}
          if(is.matrix(p3_center[,,j,k]==TRUE)) {p4_center[k,j]<-p3_center[,,j,k][k,j]}
          else {p4_center[k,j]<-p3_center[,,j,k][k]}
        }
      }
      p_center[1,1]<-p_center[1,1]+sum(p2_center%o%p2_center)/n
      p_center[1,2]<-p_center[1,2]+sum(p2_center%o%p4_center)/n
      p_center[2,2]<-p_center[2,2]+sum(p4_center%o%p4_center)/n
    }
    p_center[2,1]<-p_center[1,2]
    hatvar_gamma_center<-ginv(Q_center)%*%p_center%*%ginv(Q_center)
    hatgamma_var_center[t_ind,l]<-hatvar_gamma_center[1,1]
    
    #正确模型下中心化beta估计,用每次模拟的方差求覆盖率
    L_hatbeta_center[t_ind,l]<-hatbeta_center[t_ind,l]-qnorm(0.975)*sqrt(hatbeta_var_center[t_ind,l])/sqrt(n)
    R_hatbeta_center[t_ind,l]<-hatbeta_center[t_ind,l]+qnorm(0.975)*sqrt(hatbeta_var_center[t_ind,l])/sqrt(n)
    
    if(fun.beta(t)<=R_hatbeta_center[t_ind,l] & fun.beta(t)>=L_hatbeta_center[t_ind,l]) 
    {cum_beta_center1[t_ind]<-cum_beta_center1[t_ind]+1}
    
    #正确模型下中心化gamma估计,用每次模拟的方差求覆盖率
    L_hatgamma_center[t_ind,l]<-hatgamma_center[t_ind,l]-qnorm(0.975)*sqrt(hatgamma_var_center[t_ind,l])/sqrt(n)
    R_hatgamma_center[t_ind,l]<-hatgamma_center[t_ind,l]+qnorm(0.975)*sqrt(hatgamma_var_center[t_ind,l])/sqrt(n)
    
    if(fun.gamma(t)<=R_hatgamma_center[t_ind,l] & fun.gamma(t)>=L_hatgamma_center[t_ind,l]) 
    {cum_gamma_center1[t_ind]<-cum_gamma_center1[t_ind]+1}
    
  }
}



hv_beta_center<-rowMeans(hatbeta_var_center)
hv_gamma_center<-rowMeans(hatgamma_var_center)
betat.hat_center<-rowMeans(hatbeta_center)
gammat.hat_center<-rowMeans(hatgamma_center)
SE.beta_center<-sqrt(hv_beta_center)/sqrt(n)
SE.gamma_center<-sqrt(hv_gamma_center)/sqrt(n)


for (t_ind in 1:length(p)) {
  t<-0.3+0.3*(t_ind-1)
  #正确模型下中心化估计,用平均的方差求覆盖率
  bias.beta_center[t_ind]<-mean(hatbeta_center[t_ind,]-fun.beta(t))
  mse.beta_center[t_ind]<-mean((hatbeta_center[t_ind,]-fun.beta(t))^2)
  SD.beta_center[t_ind]<-sqrt(sum((hatbeta_center[t_ind,]-betat.hat_center[t_ind])^2)/(K-1))
  bias.gamma_center[t_ind]<-mean(hatgamma_center[t_ind,]-fun.gamma(t))
  mse.gamma_center[t_ind]<-mean((hatgamma_center[t_ind,]-fun.gamma(t))^2)
  SD.gamma_center[t_ind]<-sqrt(sum((hatgamma_center[t_ind,]-gammat.hat_center[t_ind])^2)/(K-1))
}

L_beta_center_point<-rowMeans(L_hatbeta_center)
R_beta_center_point<-rowMeans(R_hatbeta_center)

L_gamma_center_point<-rowMeans(L_hatgamma_center)
R_gamma_center_point<-rowMeans(R_hatgamma_center)

fun.beta
local_kernel
n^{-0.6};n^{-0.5}
n;h;h1_kw;h2_kw
cum_beta_center1/K
cum_gamma_center1/K
betat.hat_center
bias.beta_center
mse.beta_center
SD.beta_center
SE.beta_center
gammat.hat_center
bias.gamma_center
mse.gamma_center
SD.gamma_center
SE.gamma_center
L_beta_center_point
R_beta_center_point
L_gamma_center_point
R_gamma_center_point



