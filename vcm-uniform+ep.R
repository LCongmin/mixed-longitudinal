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
gammat.hat_lvcf<-rep(0,length(p))
bias.gamma_lvcf<-rep(0,length(p))
mse.gamma_lvcf<-rep(0,length(p))
SD.gamma_lvcf<-rep(0,length(p))
SE.gamma_lvcf<-rep(0,length(p))
cum_gamma_lvcf1<-rep(0,length(p))



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

gauss <- function(t){
  
  kt <- exp(-t*t*0.5)/sqrt(2.0*pi)
  
  return(kt)
  
}

local_kernel <- function(t, h){
  kt <- uniform(t/h)
  return(kt/h)
}
local_kernel2 <- function(t, h){
  kt <- epanechnikov(t/h)
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


#存储每个时间点的每一次模拟结果91*K的矩阵
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
hatgamma_lvcf<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatgamma_SE_lvcf<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatgamma_lvcf<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_lvcf<-matrix(0,ncol = K,nrow = length(p))#beta95%上界


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
    
    #得到数据x的核估计，这里应该是所有个体放在一起的核估计吧？
    #kw代表核权重
    for (i in 1:n) {
      kw[[i]]<-rep(0,nx[i])
    }
    
    for (i in 1:n) {
      for (j in 1:nx[i]) {
        kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
      }
    } 
    
    
    
    #部分线性方法估计x的系数beta
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
    hatbeta_p[t_ind,l]<-beta.hat_p[3,1]
    
    sig_beta_p<-matrix(0:0,ncol = 4,nrow = 4)
    for (i in 1:n) {
      sig_beta_p[1,1]<-sig_beta_p[1,1]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      sig_beta_p[1,2]<-sig_beta_p[1,2]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      sig_beta_p[1,3]<-sig_beta_p[1,3]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      sig_beta_p[1,4]<-sig_beta_p[1,4]+sum((kw[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      
      sig_beta_p[2,2]<-sig_beta_p[2,2]+sum((kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      sig_beta_p[2,3]<-sig_beta_p[2,3]+sum((kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      sig_beta_p[2,4]<-sig_beta_p[2,4]+sum((kw[[i]]*(T.x[[i]]-t)*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t)))%o%(kw[[i]]*(T.x[[i]]-t)*x[[i]]*(y[[i]]-beta.hat_p[1,1]-beta.hat_p[2,1]*(T.x[[i]]-t)-x[[i]]*beta.hat_p[3,1]-x[[i]]*beta.hat_p[4,1]*(T.x[[i]]-t))))/n
      
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
    hatvar_beta_p<-ginv(Q)%*%sig_beta_p%*%ginv(Q)
    hatbeta_var_p[t_ind,l]<-hatvar_beta_p[3,3]
    
    
    ker_p<-list(matrix)
    ker_xy_p<-list(matrix)
    for (i in 1:n) {
      ker_p[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_xy_p[[i]]<-matrix(rep((y[[i]]-x[[i]]*(beta.hat_p[3,1]+beta.hat_p[4,1]*(T.x[[i]]-t)))*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nx[i])
      #ker_xy[[i]]<-matrix(rep((y[[i]]-x[[i]]*fun.beta(T.x[[i]]))*local_kernel(T.x[[i]]-t,h1_kw),nz[i]),ncol = nz[i],nrow = nz[i])
    }
    Q_p<-matrix(0:0,ncol = 2,nrow = 2)
    q_p<-matrix(0:0,ncol = 1,nrow = 2)
    for (i in 1:n) {
      Q_p[1,1]<-Q_p[1,1]+colSums(ker_p[[i]]%*%t(t(local_kernel2(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]])))/n
      Q_p[1,2]<-Q_p[1,2]+colSums(ker_p[[i]]%*%t(t(local_kernel2(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
      Q_p[2,2]<-Q_p[2,2]+colSums(ker_p[[i]]%*%t(t(local_kernel2(S.z[[i]]-t,h2_kw)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
      q_p[1,]<-q_p[1,]+colSums(ker_xy_p[[i]]%*%t(t(local_kernel2(S.z[[i]]-t,h2_kw)*z[[i]])))/n
      q_p[2,]<-q_p[2,]+colSums(ker_xy_p[[i]]%*%t(t(local_kernel2(S.z[[i]]-t,h2_kw)*z[[i]]*(S.z[[i]]-t))))/n
    }
    Q_p[2,1]<-Q_p[1,2]
    gamma.hat_p<-ginv(Q_p)%*%q_p
    hatgamma_p[t_ind,l]<-gamma.hat_p[1,1]
    
    
    p_p<-matrix(0:0,ncol = 2,nrow = 2)
    for (i in 1:n) {
      p1_p<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p1_p<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw))%o%as.vector(local_kernel2(S.z[[i]]-t,h2_kw)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta.hat_p[3,1]+beta.hat_p[4,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma.hat_p[1,1]+gamma.hat_p[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p2_p<-matrix(0, ncol=nx[i], nrow=nz[i])
      p3_p<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p3_p<-t(as.vector(local_kernel(T.x[[i]]-t,h1_kw))%o%as.vector(local_kernel2(S.z[[i]]-t,h2_kw)*(S.z[[i]]-t)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta.hat_p[3,1]+beta.hat_p[4,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma.hat_p[1,1]+gamma.hat_p[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p4_p<-matrix(0, ncol=nx[i], nrow=nz[i])
      for (j in 1:nx[i]) {
        for (k in 1:nz[i]) {
          if(is.matrix(p1_p[,,j,k]==TRUE)) {p2_p[k,j]<-p1_p[,,j,k][k,j]}
          else {p2_p[k,j]<-p1_p[,,j,k][k]}
          if(is.matrix(p3_p[,,j,k]==TRUE)) {p4_p[k,j]<-p3_p[,,j,k][k,j]}
          else {p4_p[k,j]<-p3_p[,,j,k][k]}
        }
      }
      p_p[1,1]<-p_p[1,1]+sum(p2_p%o%p2_p)/n
      p_p[1,2]<-p_p[1,2]+sum(p2_p%o%p4_p)/n
      p_p[2,2]<-p_p[2,2]+sum(p4_p%o%p4_p)/n
    }
    p_p[2,1]<-p_p[1,2]
    hatvar_gamma_p<-ginv(Q_p)%*%p_p%*%ginv(Q_p)
    hatgamma_var_p[t_ind,l]<-hatvar_gamma_p[1,1]
    
    #正确模型下部分线性beta估计,用每次模拟的方差求覆盖率
    L_hatbeta_p[t_ind,l]<-hatbeta_p[t_ind,l]-qnorm(0.975)*sqrt(hatbeta_var_p[t_ind,l])/sqrt(n)
    R_hatbeta_p[t_ind,l]<-hatbeta_p[t_ind,l]+qnorm(0.975)*sqrt(hatbeta_var_p[t_ind,l])/sqrt(n)
    
    if(fun.beta(t)<=R_hatbeta_p[t_ind,l] & fun.beta(t)>=L_hatbeta_p[t_ind,l]) 
    {cum_beta_p1[t_ind]<-cum_beta_p1[t_ind]+1}
    
    #正确模型下部分线性gamma估计,用每次模拟的方差求覆盖率
    L_hatgamma_p[t_ind,l]<-hatgamma_p[t_ind,l]-qnorm(0.975)*sqrt(hatgamma_var_p[t_ind,l])/sqrt(n)
    R_hatgamma_p[t_ind,l]<-hatgamma_p[t_ind,l]+qnorm(0.975)*sqrt(hatgamma_var_p[t_ind,l])/sqrt(n)
    
    if(fun.gamma(t)<=R_hatgamma_p[t_ind,l] & fun.gamma(t)>=L_hatgamma_p[t_ind,l]) 
    {cum_gamma_p1[t_ind]<-cum_gamma_p1[t_ind]+1}
    
  } 
}


#正确模型下部分线性估计,用平均的方差求覆盖率
hv_beta_p<-rowMeans(hatbeta_var_p)
hv_gamma_p<-rowMeans(hatgamma_var_p)

betat.hat_p<-rowMeans(hatbeta_p)
SE.beta_p<-sqrt(hv_beta_p)/sqrt(n)
gammat.hat_p<-rowMeans(hatgamma_p)
SE.gamma_p<-sqrt(hv_gamma_p)/sqrt(n)

for (t_ind in 1:length(p)) {
  t<-0.3+0.3*(t_ind-1)
  bias.beta_p[t_ind]<-mean(hatbeta_p[t_ind,]-fun.beta(t))
  mse.beta_p[t_ind]<-mean((hatbeta_p[t_ind,]-fun.beta(t))^2)
  SD.beta_p[t_ind]<-sqrt(sum((hatbeta_p[t_ind,]-betat.hat_p[t_ind])^2)/(K-1))
  
  bias.gamma_p[t_ind]<-mean(hatgamma_p[t_ind,]-fun.gamma(t))
  mse.gamma_p[t_ind]<-mean((hatgamma_p[t_ind,]-fun.gamma(t))^2)
  SD.gamma_p[t_ind]<-sqrt(sum((hatgamma_p[t_ind,]-gammat.hat_p[t_ind])^2)/(K-1))
}


L_beta_p_point<-rowMeans(L_hatbeta_p)
R_beta_p_point<-rowMeans(R_hatbeta_p)

L_gamma_p_point<-rowMeans(L_hatgamma_p)
R_gamma_p_point<-rowMeans(R_hatgamma_p)


local_kernel
local_kernel2
fun.beta
n^{-0.6};n^{-0.5}
n;h;h1_kw;h2_kw
cum_beta_p1/K
cum_gamma_p1/K
betat.hat_p
bias.beta_p
mse.beta_p
SD.beta_p
SE.beta_p
L_beta_p_point
R_beta_p_point
gammat.hat_p
bias.gamma_p
mse.gamma_p
SD.gamma_p
SE.gamma_p
L_gamma_p_point
R_gamma_p_point


proc.time()-yxt

