###############################################################################################
###     带方差CP,分开同时都有估计，固定窗宽h=n^(-0.8),h1=h2=n^(-0.5),t=(0.05,0.95)          ###
###############################################################################################
rm(list = ls (all = TRUE))
yxt<-proc.time()
library(MASS)
p<-seq(0.05,0.95,0.005)
#部分线性下beta的估计
betat.hat_p<-rep(0,length(p))
bias.beta_p<-rep(0,length(p))
mse.beta_p<-rep(0,length(p))
SD.beta_p<-rep(0,length(p))
SE.beta_p<-rep(0,length(p))
cum_beta_p1<-rep(0,length(p))
cum_beta_p_bonf1<-rep(0,length(p))
cum_beta_p_scb1<-rep(0,length(p))
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
cum_gamma_p_bonf1<-rep(0,length(p))
cum_gamma_p_scb1<-rep(0,length(p))
hv_gamma_p<-rep(0,length(p))
L_gamma_p<-rep(0,length(p))
R_gamma_p<-rep(0,length(p))



fun.beta<-function(t){
  Beta<-0.4*t+0.5
  return(Beta)
}

#fun.beta<-function(t){
#  Beta<-3*((t-0.4)^2)
#  return(Beta)
#}

fun.gamma<-function(t){
  Gamma<-sqrt(t)
  return(Gamma)
}

#fun.gamma<-function(t){
#  Gamma<-sin(2*pi*t)
#  return(Gamma)
#}
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

local_kernel <- function(t, h){
  kt <- epanechnikov(t/h)
  return(kt/h)
}

n<-400
K<-1000#循环次数
h<-n^{-0.6}
h1<-n^{-0.45}
h2<-n^{-0.45}
#t<-0.25
#t<-0.5
#t<-0.75

#评估估计曲线
bi_beta_p<-matrix(0, nrow = length(p), ncol = K)
bi_beta_center<-matrix(0, nrow = length(p), ncol = K)
bi_beta_kw<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_p<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_center<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_kw<-matrix(0, nrow = length(p), ncol = K)

#寻找最大值分布，求同时置信带
sim_beta_p<-matrix(0, nrow = length(p), ncol = K)
sim_beta_center<-matrix(0, nrow = length(p), ncol = K)
sim_beta_kw<-matrix(0, nrow = length(p), ncol = K)
sim_gamma_p<-matrix(0, nrow = length(p), ncol = K)
sim_gamma_center<-matrix(0, nrow = length(p), ncol = K)
sim_gamma_kw<-matrix(0, nrow = length(p), ncol = K)

#提取最大值
max_beta_p<-rep(0,K)
max_beta_center<-rep(0,K)
max_beta_kw<-rep(0,K)
max_gamma_p<-rep(0,K)
max_gamma_center<-rep(0,K)
max_gamma_kw<-rep(0,K)


#存储每个时间点的每一次模拟结果91*K的矩阵
#部分线性下beta
hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_p<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatbeta_p_bonf1<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_p_bonf1<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
abs_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatbeta-beta|/sqrt(var)
L_scb_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))
#部分线性下gamma
hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatgamma_var_p<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatgamma_p_bonf1<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_p_bonf1<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
abs_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatbeta-beta|/sqrt(var)
L_scb_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))


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
    
    t<-0.05+0.005*(t_ind-1)
    
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
      ker_p[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nx[i])
      ker_xy_p[[i]]<-matrix(rep((y[[i]]-x[[i]]*(beta.hat_p[3,1]+beta.hat_p[4,1]*(T.x[[i]]-t)))*local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nx[i])
      #ker_xy[[i]]<-matrix(rep((y[[i]]-x[[i]]*fun.beta(T.x[[i]]))*local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nz[i])
    }
    Q_p<-matrix(0:0,ncol = 2,nrow = 2)
    q_p<-matrix(0:0,ncol = 1,nrow = 2)
    for (i in 1:n) {
      Q_p[1,1]<-Q_p[1,1]+colSums(ker_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]])))/n
      Q_p[1,2]<-Q_p[1,2]+colSums(ker_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
      Q_p[2,2]<-Q_p[2,2]+colSums(ker_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
      q_p[1,]<-q_p[1,]+colSums(ker_xy_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]])))/n
      q_p[2,]<-q_p[2,]+colSums(ker_xy_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*(S.z[[i]]-t))))/n
    }
    Q_p[2,1]<-Q_p[1,2]
    gamma.hat_p<-ginv(Q_p)%*%q_p
    hatgamma_p[t_ind,l]<-gamma.hat_p[1,1]
    
    
    p_p<-matrix(0:0,ncol = 2,nrow = 2)
    for (i in 1:n) {
      p1_p<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p1_p<-t(as.vector(local_kernel(T.x[[i]]-t,h1))%o%as.vector(local_kernel(S.z[[i]]-t,h2)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta.hat_p[3,1]+beta.hat_p[4,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma.hat_p[1,1]+gamma.hat_p[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
      p2_p<-matrix(0, ncol=nx[i], nrow=nz[i])
      p3_p<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
      p3_p<-t(as.vector(local_kernel(T.x[[i]]-t,h1))%o%as.vector(local_kernel(S.z[[i]]-t,h2)*(S.z[[i]]-t)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta.hat_p[3,1]+beta.hat_p[4,1]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma.hat_p[1,1]+gamma.hat_p[2,1]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
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
    
    
    bi_beta_p[t_ind,l]<-(hatbeta_p[t_ind,l]-fun.beta(t))^2
    bi_gamma_p[t_ind,l]<-(hatgamma_p[t_ind,l]-fun.gamma(t))^2
    
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
    
    
    L_hatbeta_p_bonf1[t_ind,l]<-hatbeta_p[t_ind,l]-qnorm(1-0.05/(2*length(p)))*sqrt(hatbeta_var_p[t_ind,l])/sqrt(n)
    R_hatbeta_p_bonf1[t_ind,l]<-hatbeta_p[t_ind,l]+qnorm(1-0.05/(2*length(p)))*sqrt(hatbeta_var_p[t_ind,l])/sqrt(n)
    
    if(fun.beta(t)<=R_hatbeta_p_bonf1[t_ind,l] & fun.beta(t)>=L_hatbeta_p_bonf1[t_ind,l]) 
    {cum_beta_p_bonf1[t_ind]<-cum_beta_p_bonf1[t_ind]+1}
    
    #正确模型下部分线性gamma估计,用每次模拟的方差求覆盖率
    L_hatgamma_p_bonf1[t_ind,l]<-hatgamma_p[t_ind,l]-qnorm(1-0.05/(2*length(p)))*sqrt(hatgamma_var_p[t_ind,l])/sqrt(n)
    R_hatgamma_p_bonf1[t_ind,l]<-hatgamma_p[t_ind,l]+qnorm(1-0.05/(2*length(p)))*sqrt(hatgamma_var_p[t_ind,l])/sqrt(n)
    
    if(fun.gamma(t)<=R_hatgamma_p_bonf1[t_ind,l] & fun.gamma(t)>=L_hatgamma_p_bonf1[t_ind,l]) 
    {cum_gamma_p_bonf1[t_ind]<-cum_gamma_p_bonf1[t_ind]+1}
    
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
  t<-0.05+0.005*(t_ind-1)
  bias.beta_p[t_ind]<-mean(hatbeta_p[t_ind,]-fun.beta(t))
  mse.beta_p[t_ind]<-mean((hatbeta_p[t_ind,]-fun.beta(t))^2)
  SD.beta_p[t_ind]<-sqrt(sum((hatbeta_p[t_ind,]-betat.hat_p[t_ind])^2)/(K-1))
  
  bias.gamma_p[t_ind]<-mean(hatgamma_p[t_ind,]-fun.gamma(t))
  mse.gamma_p[t_ind]<-mean((hatgamma_p[t_ind,]-fun.gamma(t))^2)
  SD.gamma_p[t_ind]<-sqrt(sum((hatgamma_p[t_ind,]-gammat.hat_p[t_ind])^2)/(K-1))
}


#MB
K_simu<-10000
abs_G_beta_p<-rep(0,K_simu)
abs_G_gamma_p<-rep(0,K_simu)


beta_p<-array(0,c(4,1,length(p)))
gamma_p<-array(0,c(2,1,length(p)))
G_p_beta<-matrix(0,ncol = length(p),nrow = K_simu)
G_p_gamma<-matrix(0,ncol = length(p),nrow = K_simu)


for (t_ind in 1:length(p)) {
  print(paste(t_ind))
  t<-0.05+0.005*(t_ind-1)
  
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
  beta_p[,,t_ind]<-ginv(Q)%*%q#l代表第几次模拟
  
  
  sig_beta_p<-rep(0,n)
  for (i in 1:n) {
    sig_beta_p[i]<-sum(kw[[i]]*x[[i]]*(y[[i]]-beta_p[1,1,t_ind]-beta_p[2,1,t_ind]*(T.x[[i]]-t)-x[[i]]*beta_p[3,1,t_ind]-x[[i]]*beta_p[4,1,t_ind]*(T.x[[i]]-t)))
  }
  
  ker_p<-list(matrix)
  ker_xy_p<-list(matrix)
  for (i in 1:n) {
    ker_p[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_xy_p[[i]]<-matrix(rep((y[[i]]-x[[i]]*(beta_p[3,1,t_ind]+beta_p[4,1,t_ind]*(T.x[[i]]-t)))*local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nx[i])
    #ker_xy[[i]]<-matrix(rep((y[[i]]-x[[i]]*fun.beta(T.x[[i]]))*local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nz[i])
  }
  Q_p<-matrix(0:0,ncol = 2,nrow = 2)
  q_p<-matrix(0:0,ncol = 1,nrow = 2)
  for (i in 1:n) {
    Q_p[1,1]<-Q_p[1,1]+colSums(ker_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]])))/n
    Q_p[1,2]<-Q_p[1,2]+colSums(ker_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
    Q_p[2,2]<-Q_p[2,2]+colSums(ker_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
    q_p[1,]<-q_p[1,]+colSums(ker_xy_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]])))/n
    q_p[2,]<-q_p[2,]+colSums(ker_xy_p[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*(S.z[[i]]-t))))/n
  }
  Q_p[2,1]<-Q_p[1,2]
  gamma_p[,,t_ind]<-ginv(Q_p)%*%q_p
  
  
  p_p<-rep(0,n)
  for (i in 1:n) {
    p1_p<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p1_p<-t(as.vector(local_kernel(T.x[[i]]-t,h1))%o%as.vector(local_kernel(S.z[[i]]-t,h2)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta_p[3,1,t_ind]+beta_p[4,1,t_ind]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma_p[1,1,t_ind]+gamma_p[2,1,t_ind]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p2_p<-matrix(0, ncol=nx[i], nrow=nz[i])
    for (j in 1:nx[i]) {
      for (k in 1:nz[i]) {
        if(is.matrix(p1_p[,,j,k]==TRUE)) {p2_p[k,j]<-p1_p[,,j,k][k,j]}
        else {p2_p[k,j]<-p1_p[,,j,k][k]}
      }
    }
    p_p[i]<-sum(p2_p)
  }
  
  for (k_simu in 1:K_simu) {
    #for (i in 1:n) {
    tau<-rbinom(n,1,0.5)
    tau[which(tau==0)]=-1
    G_p_beta[k_simu,t_ind]<-ginv(Q[3,3])*sum(tau*sig_beta_p)/n
    G_p_gamma[k_simu,t_ind]<-ginv(Q_p[1,1])*sum(tau*p_p)/n
      
    #}
  }
}


for (k_simu in 1:K_simu) {
  abs_G_beta_p[k_simu]<-max(abs(G_p_beta[k_simu,]))
  abs_G_gamma_p[k_simu]<-max(abs(G_p_gamma[k_simu,]))
}

quantile(abs_G_beta_p,0.95)
quantile(abs_G_gamma_p,0.95)

realbeta<-rep(0,length(p))
for (t_ind in 1:length(p)) {
  realbeta[t_ind]<-fun.beta(0.05+(t_ind-1)*0.005)
}

realgamma<-rep(0,length(p))
for (t_ind in 1:length(p)) {
  realgamma[t_ind]<-fun.gamma(0.05+(t_ind-1)*0.005)
}

#beta_p所有的逐点覆盖率
cum_beta_p2<-0
for (l in 1:K) {
  if(all(realbeta<=R_hatbeta_p[,l] & realbeta>=L_hatbeta_p[,l]))
    cum_beta_p2<-cum_beta_p2+1
}

#gamma_p所有的逐点覆盖率
cum_gamma_p2<-0
for (l in 1:K) {
  if(all(realgamma<=R_hatgamma_p[,l] & realgamma>=L_hatgamma_p[,l]))
    cum_gamma_p2<-cum_gamma_p2+1
}


#beta_p,gamma_p所有的同时覆盖率
for (t_ind in 1:length(p)) {
  for (l in 1:K) {
    L_scb_hatbeta_p[t_ind,l]<-hatbeta_p[t_ind,l]-quantile(abs_G_beta_p,0.95)
    R_scb_hatbeta_p[t_ind,l]<-hatbeta_p[t_ind,l]+quantile(abs_G_beta_p,0.95)
    L_scb_hatgamma_p[t_ind,l]<-hatgamma_p[t_ind,l]-quantile(abs_G_gamma_p,0.95)
    R_scb_hatgamma_p[t_ind,l]<-hatgamma_p[t_ind,l]+quantile(abs_G_gamma_p,0.95)
  }
}



cum_beta_p_scb1<-rep(0,length(p))
cum_gamma_p_scb1<-rep(0,length(p))

for (t_ind in 1:length(p)) {
  t<-0.05+0.005*(t_ind-1)
  for (l in 1:K) {
    
    if(fun.beta(t)<=R_scb_hatbeta_p[t_ind,l] & fun.beta(t)>=L_scb_hatbeta_p[t_ind,l]) 
    {cum_beta_p_scb1[t_ind]<-cum_beta_p_scb1[t_ind]+1}
    
    if(fun.gamma(t)<=R_scb_hatgamma_p[t_ind,l] & fun.gamma(t)>=L_scb_hatgamma_p[t_ind,l]) 
    {cum_gamma_p_scb1[t_ind]<-cum_gamma_p_scb1[t_ind]+1}
    
  }
}

cum_beta_p_scb2<-0
for (l in 1:K) {
  if(all(realbeta<=R_scb_hatbeta_p[,l] & realbeta>=L_scb_hatbeta_p[,l]))
    cum_beta_p_scb2<-cum_beta_p_scb2+1
}

cum_gamma_p_scb2<-0
for (l in 1:K) {
  if(all(realgamma<=R_scb_hatgamma_p[,l] & realgamma>=L_scb_hatgamma_p[,l]))
    cum_gamma_p_scb2<-cum_gamma_p_scb2+1
}

cum_beta_p_bonf2<-0
for (l in 1:K) {
  if(all(realbeta<=R_hatbeta_p_bonf1[,l] & realbeta>=L_hatbeta_p_bonf1[,l]))
    cum_beta_p_bonf2<-cum_beta_p_bonf2+1
}
cum_gamma_p_bonf2<-0
for (l in 1:K) {
  if(all(realgamma<=R_hatgamma_p_bonf1[,l] & realgamma>=L_hatgamma_p_bonf1[,l]))
    cum_gamma_p_bonf2<-cum_gamma_p_bonf2+1
}

#RASE,SD评估整条曲线
RASE_beta_p<-sqrt(colMeans(bi_beta_p))
MEAN_beta_p<-mean(RASE_beta_p)
SD_beta_p<-sqrt(sum((RASE_beta_p-MEAN_beta_p)^2)/(K-1))


RASE_gamma_p<-sqrt(colMeans(bi_gamma_p))
MEAN_gamma_p<-mean(RASE_gamma_p)
SD_gamma_p<-sqrt(sum((RASE_gamma_p-MEAN_gamma_p)^2)/(K-1))


beta_p_quantile<-quantile(abs_G_beta_p,0.95)
gamma_p_quantile<-quantile(abs_G_gamma_p,0.95)

L_beta_p_point<-rowMeans(L_hatbeta_p)
R_beta_p_point<-rowMeans(R_hatbeta_p)
L_beta_p_scb<-rowMeans(L_scb_hatbeta_p)
R_beta_p_scb<-rowMeans(R_scb_hatbeta_p)
L_beta_p_bonf<-rowMeans(L_hatbeta_p_bonf1)
R_beta_p_bonf<-rowMeans(R_hatbeta_p_bonf1)

L_gamma_p_point<-rowMeans(L_hatgamma_p)
R_gamma_p_point<-rowMeans(R_hatgamma_p)
L_gamma_p_scb<-rowMeans(L_scb_hatgamma_p)
R_gamma_p_scb<-rowMeans(R_scb_hatgamma_p)
L_gamma_p_bonf<-rowMeans(L_hatbeta_p_bonf1)
R_gamma_p_bonf<-rowMeans(R_hatbeta_p_bonf1)

beta_p_quantile
cum_beta_p1/K;cum_beta_p2/K
cum_beta_p_bonf1/K;cum_beta_p_bonf2/K
cum_beta_p_scb1/K;cum_beta_p_scb2/K

gamma_p_quantile
cum_gamma_p1/K;cum_gamma_p2/K
cum_gamma_p_bonf1/K;cum_gamma_p_bonf2/K
cum_gamma_p_scb1/K;cum_gamma_p_scb2/K



betat.hat_p
bias.beta_p
mse.beta_p
SD.beta_p
SE.beta_p
L_beta_p_point
R_beta_p_point
L_beta_p_scb
R_beta_p_scb
L_beta_p_bonf
R_beta_p_bonf

gammat.hat_p
bias.gamma_p
mse.gamma_p
SD.gamma_p
SE.gamma_p
L_gamma_p_point
R_gamma_p_point
L_gamma_p_scb
R_gamma_p_scb
L_gamma_p_bonf
R_gamma_p_bonf



MEAN_beta_p;SD_beta_p
#RASE_beta_p
MEAN_gamma_p;SD_gamma_p
#RASE_gamma_p

result=rbind(beta_p_quantile,
             cum_beta_p1/K,cum_beta_p2/K,
             cum_beta_p_bonf1/K,cum_beta_p_bonf2/K,
             cum_beta_p_scb1/K,cum_beta_p_scb2/K,
             gamma_p_quantile,
             cum_gamma_p1/K,cum_gamma_p2/K,
             cum_gamma_p_bonf1/K,cum_gamma_p_bonf2/K,
             cum_gamma_p_scb1/K,cum_gamma_p_scb2/K,
             betat.hat_p,
             bias.beta_p,
             mse.beta_p,
             SD.beta_p,
             SE.beta_p,
             L_beta_p_point,
             R_beta_p_point,
             L_beta_p_scb,
             R_beta_p_scb,
             L_beta_p_bonf,
             R_beta_p_bonf,
             gammat.hat_p,
             bias.gamma_p,
             mse.gamma_p,
             SD.gamma_p,
             SE.gamma_p,
             L_gamma_p_point,
             R_gamma_p_point,
             L_gamma_p_scb,
             R_gamma_p_scb,
             L_gamma_p_bonf,
             R_gamma_p_bonf,
             MEAN_beta_p,SD_beta_p,
             MEAN_gamma_p,SD_gamma_p)
write.csv(result,file = "E:/papercao/2020.12.14_datafirst_1K2t/p_data_line_n400h06045_ez2t.csv")

a<-quantile(abs_G_beta_p,0.95)
b<-quantile(abs_G_gamma_p,0.95)
result<-rbind(a,b)


write.csv(result,file = "E:/papercao/2020.12.14_datafirst_1K2t/p_line_n400h06045_scb_ez2t.csv")

v<-rbind(abs_G_beta_p,abs_G_gamma_p)
write.csv(v,file = "E:/papercao/2020.12.14_datafirst_1K2t/G_p_line_n400h06045_ez2t.csv")

plot_result<-rbind(betat.hat_p,gammat.hat_p,
                   L_beta_p_point,R_beta_p_point,
                   L_beta_p_scb,R_beta_p_scb,
                   L_beta_p_bonf,R_beta_p_bonf,
                   L_gamma_p_point,R_gamma_p_point,
                   L_gamma_p_scb,R_gamma_p_scb,
                   L_gamma_p_bonf,R_gamma_p_bonf)
write.csv(plot_result,row.names = c("beta","gamma",
                                    "L_beta_p_point","R_beta_p_point",
                                    "L_beta_p_scb","R_beta_p_scb",
                                    "L_beta_p_bonf","R_beta_p_bonf",
                                    "L_gamma_p_point","R_gamma_p_point",
                                    "L_gamma_p_scb","R_gamma_p_scb",
                                    "L_gamma_p_bonf","R_gamma_p_bonf"),
          file="E:/papercao/2020.12.14_datafirst_1K2t/plot_p_line_n400_h06045_181point_ez2t.csv")



result<-rbind(hatbeta_p,hatbeta_var_p,
              hatgamma_p,hatgamma_var_p,
              bi_beta_p,bi_gamma_p,
              R_hatbeta_p,L_hatbeta_p,
              R_hatgamma_p,L_hatgamma_p,
              R_scb_hatbeta_p,L_scb_hatbeta_p,
              R_scb_hatgamma_p,L_scb_hatgamma_p)
write.csv(result,row.names = c(rep("hatbeta_p",181),rep("hatbeta_var_p",181),
                               rep("hatgamma_p",181),rep("hatgamma_var_p",181),
                               rep("bi_beta_p",181),rep("bi_gamma_p",181),
                               rep("R_hatbeta_p",181),rep("L_hatbeta_p",181),
                               rep("R_hatgamma_p",181),rep("L_hatgamma_p",181),
                               rep("R_scb_hatbeta_p",181),rep("L_scb_hatbeta_p",181),
                               rep("R_scb_hatgamma_p",181),rep("L_scb_hatgamma_p",181)),
          file="E:/papercao/2020.12.14_datafirst_1K2t/p_line_n400_h06045_181point_ez2t.csv")

#图片大小4.5×4
#a<-seq(0.05,0.95,0.01)
#两步估计，变系数模型+核加权
#plot(a,fun.beta(a),xaxt = "n", yaxt = "n",xlim=c(0,1),ylim=c(-1.5,2.0),xlab="t",ylab=expression(paste(beta, "(t)")),type="l",col="black",lwd=2)
#axis(1,at=c("0.0","0.2","0.4","0.6","0.8","1.0"),label=c("0.0","0.2","0.4","0.6","0.8","1.0")) 
#axis(2,at=c("-2.0","-1.0","0.0","1.0","2.0"),label=c("-2.0","-1.0","0.0","1.0","2.0")) 
#lines(a,betat.hat_p,col="red")
#lines(a,rowMeans(L_hatbeta_p),col="blue",lty=2)
#lines(a,rowMeans(R_hatbeta_p),col="blue",lty=2)
#lines(a,rowMeans(L_scb_hatbeta_p),lty=2,col="magenta")
#lines(a,rowMeans(R_scb_hatbeta_p),lty=2,col="magenta")
#legend(-0.04,0.0, lty=c(1,1,2,2), col=c("black", "red","blue","magenta"),
#       legend=c("true function", "estimated function","95% point-wise CI","95% scb"),
#       merge = TRUE,cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6,
#       text.width=0.5)
#legend(-0.12,2.5, text.col="black",legend="Method: Two-step (VCM)",cex = 0.75,
#       bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6)
#


#plot(a,fun.gamma(a),xaxt = "n", yaxt = "n",xlim=c(0,1),ylim=c(-1.5,2.0),xlab="t",ylab=expression(paste(gamma, "(t)")),type="l",col="black",lwd=2)
#axis(1,at=c("0.0","0.2","0.4","0.6","0.8","1.0"),label=c("0.0","0.2","0.4","0.6","0.8","1.0")) 
#axis(2,at=c("-2.0","-1.0","0.0","1.0","2.0"),label=c("-2.0","-1.0","0.0","1.0","2.0")) 
#lines(a,gammat.hat_p,col="red")
#lines(a,rowMeans(L_hatgamma_p),col="blue",lty=2)
#lines(a,rowMeans(R_hatgamma_p),col="blue",lty=2)
#lines(a,rowMeans(L_scb_hatgamma_p),lty=2,col="magenta")
#lines(a,rowMeans(R_scb_hatgamma_p),lty=2,col="magenta")
#legend(-0.04,0.0, lty=c(1,1,2,2), col=c("black", "red","blue","magenta"),
#       legend=c("true function", "estimated function","95% point-wise CI","95% scb"),
#       merge = TRUE,cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6,
#       text.width=0.5)
#legend(-0.12,2.5, text.col="black",legend="Method: Two-step (VCM+KW)",cex = 0.75,
#       bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6)
#


#两步估计，中心化+核加权
#plot(a,fun.beta(a),xaxt = "n", yaxt = "n",xlim=c(0,1),ylim=c(-1.5,2.0),xlab="t",ylab=expression(paste(beta, "(t)")),type="l",col="black",lwd=2)
#axis(1,at=c("0.0","0.2","0.4","0.6","0.8","1.0"),label=c("0.0","0.2","0.4","0.6","0.8","1.0")) 
#axis(2,at=c("-2.0","-1.0","0.0","1.0","2.0"),label=c("-2.0","-1.0","0.0","1.0","2.0")) 
#lines(a,betat.hat_center,col="red")
#lines(a,rowMeans(L_hatbeta_center),col="blue",lty=2)
#lines(a,rowMeans(R_hatbeta_center),col="blue",lty=2)
#lines(a,rowMeans(L_scb_hatbeta_center),lty=2,col="magenta")
#lines(a,rowMeans(R_scb_hatbeta_center),lty=2,col="magenta")
#legend(-0.04,0.0, lty=c(1,1,2,2), col=c("black", "red","blue","magenta"),
#       legend=c("true function", "estimated function","95% point-wise CI","95% scb"),
#       merge = TRUE,cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6,
#       text.width=0.5)
#legend(-0.12,2.5, text.col="black",legend="Method: Two-step (Center)",cex = 0.75,
#       bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6)
#


#plot(a,fun.gamma(a),xaxt = "n", yaxt = "n",xlim=c(0,1),ylim=c(-1.5,2.0),xlab="t",ylab=expression(paste(gamma, "(t)")),type="l",col="black",lwd=2)
#axis(1,at=c("0.0","0.2","0.4","0.6","0.8","1.0"),label=c("0.0","0.2","0.4","0.6","0.8","1.0")) 
#axis(2,at=c("-2.0","-1.0","0.0","1.0","2.0"),label=c("-2.0","-1.0","0.0","1.0","2.0")) 
#lines(a,gammat.hat_center,col="red")
#lines(a,rowMeans(L_hatgamma_center),col="blue",lty=2)
#lines(a,rowMeans(R_hatgamma_center),col="blue",lty=2)
#lines(a,rowMeans(L_scb_hatgamma_center),lty=2,col="magenta")
#lines(a,rowMeans(R_scb_hatgamma_center),lty=2,col="magenta")
#legend(-0.04,0.0, lty=c(1,1,2,2), col=c("black", "red","blue","magenta"),
#       legend=c("true function", "estimated function","95% point-wise CI","95% scb"),
#       merge = TRUE,cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6,
#       text.width=0.5)
#legend(-0.12,2.5, text.col="black",legend="Method: Two-step (Center+KW)",cex = 0.75,
#       bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6)
#




#一步估计，核加权
#plot(a,fun.beta(a),xaxt = "n", yaxt = "n",xlim=c(0,1),ylim=c(-1.5,2.0),xlab="t",ylab=expression(paste(beta, "(t)")),type="l",col="black",lwd=2)
#axis(1,at=c("0.0","0.2","0.4","0.6","0.8","1.0"),label=c("0.0","0.2","0.4","0.6","0.8","1.0")) 
#axis(2,at=c("-2.0","-1.0","0.0","1.0","2.0"),label=c("-2.0","-1.0","0.0","1.0","2.0")) 
#lines(a,betagammat.hat[,1],col="red")
#lines(a,rowMeans(L_hatbeta_kw),col="blue",lty=2)
#lines(a,rowMeans(R_hatbeta_kw),col="blue",lty=2)
#lines(a,rowMeans(L_scb_hatbeta_kw),lty=2,col="magenta")
#lines(a,rowMeans(R_scb_hatbeta_kw),lty=2,col="magenta")
#legend(-0.04,0.0, lty=c(1,1,2,2), col=c("black", "red","blue","magenta"),
#       legend=c("true function", "estimated function","95% point-wise CI","95% scb"),
#       merge = TRUE,cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6,
#       text.width=0.5)
#legend(-0.12,2.5, text.col="black",legend="Method: One-step (KW)",cex = 0.75,
#       bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6)
#


#plot(a,fun.gamma(a),xaxt = "n", yaxt = "n",xlim=c(0,1),ylim=c(-1.5,2.0),xlab="t",ylab=expression(paste(gamma, "(t)")),type="l",col="black",lwd=2)
#axis(1,at=c("0.0","0.2","0.4","0.6","0.8","1.0"),label=c("0.0","0.2","0.4","0.6","0.8","1.0")) 
#axis(2,at=c("-2.0","-1.0","0.0","1.0","2.0"),label=c("-2.0","-1.0","0.0","1.0","2.0")) 
#lines(a,betagammat.hat[,2],col="red")
#lines(a,rowMeans(L_hatgamma_kw),col="blue",lty=2)
#lines(a,rowMeans(R_hatgamma_kw),col="blue",lty=2)
#lines(a,rowMeans(L_scb_hatgamma_kw),lty=2,col="magenta")
#lines(a,rowMeans(R_scb_hatgamma_kw),lty=2,col="magenta")
#legend(-0.04,0.0, lty=c(1,1,2,2), col=c("black", "red","blue","magenta"),
#       legend=c("true function", "estimated function","95% point-wise CI","95% scb"),
#       merge = TRUE,cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6,
#       text.width=0.5)
#legend(-0.12,2.5, text.col="black",legend="Method: One-step (KW)",cex = 0.75,
#       bty = "n",x.intersp=0.2,y.intersp=0.3,seg.len=0.6)
#



proc.time()-yxt

