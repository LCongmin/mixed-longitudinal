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
cum_beta_center_bonf1<-rep(0,length(p))
L_beta_center_bonf<-rep(0,length(p))
R_beta_center_bonf<-rep(0,length(p))
cum_beta_center_sche1<-rep(0,length(p))
L_beta_center_sche<-rep(0,length(p))
R_beta_center_sche<-rep(0,length(p))



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
cum_gamma_center_bonf1<-rep(0,length(p))
L_gamma_center_bonf<-rep(0,length(p))
R_gamma_center_bonf<-rep(0,length(p))
cum_gamma_center_sche1<-rep(0,length(p))
L_gamma_center_sche<-rep(0,length(p))
R_gamma_center_sche<-rep(0,length(p))



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
cum_beta_kw_bonf1<-rep(0,length(p))
cum_gamma_kw_bonf1<-rep(0,length(p))
cum_beta_kw_sche1<-rep(0,length(p))
cum_gamma_kw_sche1<-rep(0,length(p))




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

local_kernel <- function(t, h){
  kt <- epanechnikov(t/h)
  return(kt/h)
}

n<-400
K<-1000#循环次数
h<-n^{-0.6}
h1<-n^{-0.5}
h2<-n^{-0.5}

#评估估计曲线
bi_beta_p<-matrix(0, nrow = length(p), ncol = K)
bi_beta_center<-matrix(0, nrow = length(p), ncol = K)
bi_beta_kw<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_p<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_center<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_kw<-matrix(0, nrow = length(p), ncol = K)


#存储每个时间点的每一次模拟结果90*K的矩阵
#部分线性下beta
hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_p<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatbeta_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatbeta_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
abs_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatbeta-beta|/sqrt(var)
L_scb_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))
#部分线性下gamma
hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatgamma_var_p<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatgamma_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatgamma_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatgamma_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
abs_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatbeta-beta|/sqrt(var)
L_scb_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))

#中心化方法下beta的估计
hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_center<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_scb_hatbeta_center<-matrix(0,ncol=K,nrow=length(p))
R_scb_hatbeta_center<-matrix(0,ncol=K,nrow=length(p))
L_hatbeta_center_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_center_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatbeta_center_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_center_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
abs_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatbeta-beta|/sqrt(var)

#中心化方法下gamma的估计
hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma估计值
hatgamma_var_center<-matrix(0,ncol = K,nrow = length(p))#gamma方差
L_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界
L_hatgamma_center_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_center_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界
L_hatgamma_center_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_center_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界
abs_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatgamma-gamma|/sqrt(var)
L_scb_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))

#核加权方法下beta的估计
hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta估计值
hatbeta_var_kw<-matrix(0,ncol = K,nrow = length(p))#beta方差
hatbetadot_var_kw<-matrix(0,ncol = K,nrow = length(p))#beta方差
L_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatbeta_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
L_hatbeta_kw_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%下界
R_hatbeta_kw_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%上界
abs_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatbeta-beta|/sqrt(var)
L_scb_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))

#核加权方法下gamma的估计
hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma估计值
hatgamma_var_kw<-matrix(0,ncol = K,nrow = length(p))#gamma方差
hatgammadot_var_kw<-matrix(0,ncol = K,nrow = length(p))#gamma方差
L_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界
L_hatgamma_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界
L_hatgamma_kw_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%下界
R_hatgamma_kw_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%上界
abs_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#绝对值|hatgamma-gamma|/sqrt(var)
L_scb_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))

#MB
K_simu<-10000
abs_G_beta_center_N01<-rep(0,K_simu)
abs_G_gamma_center_N01<-rep(0,K_simu)
abs_G_beta_center_E1<-rep(0,K_simu)
abs_G_gamma_center_E1<-rep(0,K_simu)
abs_G_beta_center_P<-rep(0,K_simu)
abs_G_gamma_center_P<-rep(0,K_simu)
abs_G_beta_center_R<-rep(0,K_simu)
abs_G_gamma_center_R<-rep(0,K_simu)
abs_G_beta_center_NN<-rep(0,K_simu)
abs_G_gamma_center_NN<-rep(0,K_simu)
abs_G_beta_center_22<-rep(0,K_simu)
abs_G_gamma_center_22<-rep(0,K_simu)


beta_center<-array(0,c(2,1,length(p)))
gamma_center<-array(0,c(2,1,length(p)))
G_center_beta_N01<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma_N01<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_beta_E1<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma_E1<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_beta_P<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma_P<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_beta_R<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma_R<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_beta_NN<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma_NN<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_beta_22<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma_22<-matrix(0,ncol = length(p),nrow = K_simu)


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
  Q_beta_center<-matrix(0:0,ncol=2,nrow=2)#求解系数的2×2矩阵
  q_beta_center<-matrix(0:0,ncol = 1,nrow = 2)#求解系数的2×1矩阵
  for (i in 1:n) {
    q_beta_center<-q_beta_center+matrix(c(sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]])/n,
                                          sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]]*(T.x[[i]]-t))/n),ncol = 1,nrow = 2)
    Q_beta_center<-Q_beta_center+matrix(c(sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]])/n,
                                          sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                                          sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                                          sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n),
                                        ncol = 2,nrow = 2)
  }
  beta_center[,,t_ind]<-ginv(Q_beta_center)%*%q_beta_center#l代表第几次模拟
  
  sig_beta<-rep(0,n)
  for (i in 1:n) {
    sig_beta[i]<-sum(kw[[i]]*x.tilde[[i]]*(y.tilde[[i]]-x.tilde[[i]]*beta_center[1,1,t_ind]-x.tilde[[i]]*beta_center[2,1,t_ind]*(T.x[[i]]-t)))
  }
  
  ker_center<-list(matrix)
  ker_xy_center<-list(matrix)
  for (i in 1:n) {
    ker_center[[i]]<-matrix(rep(local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nx[i])
    ker_xy_center[[i]]<-matrix(rep((y[[i]]-x[[i]]*(beta_center[1,1,t_ind]+beta_center[2,1,t_ind]*(T.x[[i]]-t)))*local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nx[i])
    #ker_xy[[i]]<-matrix(rep((y[[i]]-x[[i]]*fun.beta(T.x[[i]]))*local_kernel(T.x[[i]]-t,h1),nz[i]),ncol = nz[i],nrow = nz[i])
  }
  Q_center<-matrix(0:0,ncol = 2,nrow = 2)
  q_center<-matrix(0:0,ncol = 1,nrow = 2)
  for (i in 1:n) {
    Q_center[1,1]<-Q_center[1,1]+colSums(ker_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]])))/n
    Q_center[1,2]<-Q_center[1,2]+colSums(ker_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]]*(S.z[[i]]-t))))/n
    Q_center[2,2]<-Q_center[2,2]+colSums(ker_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*z[[i]]*(S.z[[i]]-t)*(S.z[[i]]-t))))/n
    q_center[1,]<-q_center[1,]+colSums(ker_xy_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]])))/n
    q_center[2,]<-q_center[2,]+colSums(ker_xy_center[[i]]%*%t(t(local_kernel(S.z[[i]]-t,h2)*z[[i]]*(S.z[[i]]-t))))/n
  }
  Q_center[2,1]<-Q_center[1,2]
  gamma_center[,,t_ind]<-ginv(Q_center)%*%q_center
  
  p1_1<-rep(0,length(n))
  for (i in 1:n) {
    p1_center<-array(0,c(nz[i],nx[i],nx[i],nz[i]))
    p1_center<-t(as.vector(local_kernel(T.x[[i]]-t,h1))%o%as.vector(local_kernel(S.z[[i]]-t,h2)*z[[i]]))%o%(matrix(rep(y[[i]]-x[[i]]*(beta_center[1,1,t_ind]+beta_center[2,1,t_ind]*(T.x[[i]]-t)),nz[i]),ncol = nz[i],nrow = nx[i])-matrix(rep(z[[i]]*(gamma_center[1,1,t_ind]+gamma_center[2,1,t_ind]*(S.z[[i]]-t)),nx[i]),ncol = nz[i],nrow = nx[i],byrow = T))
    p2_center<-matrix(0, ncol=nx[i], nrow=nz[i])
    for (j in 1:nx[i]) {
      for (k in 1:nz[i]) {
        if(is.matrix(p1_center[,,j,k]==TRUE)) {p2_center[k,j]<-p1_center[,,j,k][k,j]}
        else {p2_center[k,j]<-p1_center[,,j,k][k]}
      }
    }
    p1_1[i]<-sum(p2_center)
  }
  
  
  for (k_simu in 1:K_simu) {
    #for (i in 1:n) {
    tau_N01<-rnorm(n,0,1)
    G_center_beta_N01[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_N01*sig_beta)/n
    G_center_gamma_N01[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_N01*p1_1)/n
    
    #tau_N11<-rnorm(n,1,1)
    #G_center_beta_N11[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_N11*sig_beta)/n
    #G_center_gamma_N11[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_N11*p1_1)/n
    
    #tau_N02<-rnorm(n,0,2)
    #G_center_beta_N02[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_N02*sig_beta)/n
    #G_center_gamma_N02[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_N02*p1_1)/n
    
    #tau_N12<-rnorm(n,1,2)
    #G_center_beta_N12[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_N12*sig_beta)/n
    #G_center_gamma_N12[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_N12*p1_1)/n
    
    tau_E1<-rexp(n,1)-1
    G_center_beta_E1[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_E1*sig_beta)/n
    G_center_gamma_E1[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_E1*p1_1)/n
    
    tau_P<-rpois(n,1)-1
    G_center_beta_P[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_P*sig_beta)/n
    G_center_gamma_P[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_P*p1_1)/n
    
    #tau_B<-rbeta(n,sqrt(2)-1,1)
    #G_center_beta_B[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_B*sig_beta)/n
    #G_center_gamma_B[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_B*p1_1)/n
    
    tau_B<-rbinom(n,1,0.5)
    tau_B[which(tau_B==0)]=-1
    G_center_beta_R[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_B*sig_beta)/n
    G_center_gamma_R[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_B*p1_1)/n
    
    tau1=rnorm(n,(sqrt(17/6)+sqrt(1/6))/2,sqrt(1/2))
    tau2=rnorm(n,(sqrt(17/6)-sqrt(1/6))/2,sqrt(1/2))
    tau_NN<-tau1*tau2-((sqrt(17/6)+sqrt(1/6))/2)*((sqrt(17/6)-sqrt(1/6))/2)
    G_center_beta_NN[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau_NN*sig_beta)/n
    G_center_gamma_NN[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau_NN*p1_1)/n
    
    tau22=rbinom(n,1,(sqrt(5)+1)/(2*sqrt(5)))
    tau22[which(tau22==1)]=-(sqrt(5)-1)/2
    tau22[which(tau22==0)]=(sqrt(5)+1)/2
    G_center_beta_22[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau22*sig_beta)/n
    G_center_gamma_22[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau22*p1_1)/n
    
    #}
  }
}


for (k_simu in 1:K_simu) {
  abs_G_beta_center_N01[k_simu]<-max(abs(G_center_beta_N01[k_simu,]))
  abs_G_gamma_center_N01[k_simu]<-max(abs(G_center_gamma_N01[k_simu,]))
  #abs_G_beta_center_N02[k_simu]<-max(abs(G_center_beta_N02[k_simu,]))
  #abs_G_gamma_center_N02[k_simu]<-max(abs(G_center_gamma_N02[k_simu,]))
  #abs_G_beta_center_N12[k_simu]<-max(abs(G_center_beta_N12[k_simu,]))
  #abs_G_gamma_center_N12[k_simu]<-max(abs(G_center_gamma_N12[k_simu,]))
  #abs_G_beta_center_N11[k_simu]<-max(abs(G_center_beta_N11[k_simu,]))
  #abs_G_gamma_center_N11[k_simu]<-max(abs(G_center_gamma_N11[k_simu,]))
  abs_G_beta_center_E1[k_simu]<-max(abs(G_center_beta_E1[k_simu,]))
  abs_G_gamma_center_E1[k_simu]<-max(abs(G_center_gamma_E1[k_simu,]))
  abs_G_beta_center_R[k_simu]<-max(abs(G_center_beta_R[k_simu,]))
  abs_G_gamma_center_R[k_simu]<-max(abs(G_center_gamma_R[k_simu,]))
  abs_G_beta_center_P[k_simu]<-max(abs(G_center_beta_P[k_simu,]))
  abs_G_gamma_center_P[k_simu]<-max(abs(G_center_gamma_P[k_simu,]))
  #abs_G_beta_center_1[k_simu]<-max(abs(G_center_beta_1[k_simu,]))
  #abs_G_gamma_center_1[k_simu]<-max(abs(G_center_gamma_1[k_simu,]))
  abs_G_beta_center_NN[k_simu]<-max(abs(G_center_beta_NN[k_simu,]))
  abs_G_gamma_center_NN[k_simu]<-max(abs(G_center_gamma_NN[k_simu,]))
  abs_G_beta_center_22[k_simu]<-max(abs(G_center_beta_22[k_simu,]))
  abs_G_gamma_center_22[k_simu]<-max(abs(G_center_gamma_22[k_simu,]))
}

for (k_simu in 1:K_simu) {
  abs_G_beta_center_N01[k_simu]<-max(abs(G_center_beta_N01[k_simu,]))
  abs_G_gamma_center_N01[k_simu]<-max(abs(G_center_gamma_N01[k_simu,]))
  abs_G_beta_center_E1[k_simu]<-max(abs(G_center_beta_E1[k_simu,]))
  abs_G_gamma_center_E1[k_simu]<-max(abs(G_center_gamma_E1[k_simu,]))
  abs_G_beta_center_R[k_simu]<-max(abs(G_center_beta_R[k_simu,]))
  abs_G_gamma_center_R[k_simu]<-max(abs(G_center_gamma_R[k_simu,]))
  abs_G_beta_center_P[k_simu]<-max(abs(G_center_beta_P[k_simu,]))
  abs_G_gamma_center_P[k_simu]<-max(abs(G_center_gamma_P[k_simu,]))
  abs_G_beta_center_NN[k_simu]<-max(abs(G_center_beta_NN[k_simu,]))
  abs_G_gamma_center_NN[k_simu]<-max(abs(G_center_gamma_NN[k_simu,]))
  abs_G_beta_center_22[k_simu]<-max(abs(G_center_beta_22[k_simu,]))
  abs_G_gamma_center_22[k_simu]<-max(abs(G_center_gamma_22[k_simu,]))
}
quantile(abs_G_beta_center_R,0.95);quantile(abs_G_gamma_center_R,0.95)
quantile(abs_G_beta_center_22,0.95);quantile(abs_G_gamma_center_22,0.95)
quantile(abs_G_beta_center_N01,0.95);quantile(abs_G_gamma_center_N01,0.95)
quantile(abs_G_beta_center_P,0.95);quantile(abs_G_gamma_center_P,0.95)
quantile(abs_G_beta_center_NN,0.95);quantile(abs_G_gamma_center_NN,0.95)
quantile(abs_G_beta_center_E1,0.95);quantile(abs_G_gamma_center_E1,0.95)


par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(1,2))
plot(density(abs_G_beta_center_N01),xlim=c(0,2.5),ylim=c(0,6), main = "",xlab = "Values")
#lines(density(abs_G_beta_center_N11),col="red")
lines(density(abs_G_beta_center_R),col="green",lty=1)
lines(density(abs_G_beta_center_22),col="gray",lty=1)
lines(density(abs_G_beta_center_P),col="red",lty=1)
lines(density(abs_G_beta_center_NN),col="purple",lty=1)
lines(density(abs_G_beta_center_E1),col="blue")
#lines(density(abs_G_beta_center_N12),col="green")
#lines(density(abs_G_beta_center_N02),col="red",lty=2)
#lines(density(abs_G_beta_center_1),col="black",lty=2)
legend("topright", lty=c(1,1,1,1,1,1), col=c("green","gray","black","red","purple","blue"),
       legend=c("Rademacher","Two-point","N(0,1)","Poisson(1)-1","NN","Exp(1)-1"))

plot(density(abs_G_gamma_center_N01),xlim=c(0,2.5),ylim=c(0,6),main = "",xlab = "Values")
#lines(density(abs_G_gamma_center_N11),col="red")
lines(density(abs_G_gamma_center_R),col="green",lty=1)
lines(density(abs_G_gamma_center_22),col="gray",lty=1)
lines(density(abs_G_gamma_center_P),col="red",lty=1)
lines(density(abs_G_gamma_center_NN),col="purple",lty=1)
lines(density(abs_G_gamma_center_E1),col="blue")
#lines(density(abs_G_gamma_center_N12),col="green")
#lines(density(abs_G_gamma_center_N02),col="red",lty=2)
#lines(density(abs_G_gamma_center_1),col="black",lty=2)
legend("topright", lty=c(1,1,1,1,1,1), col=c("green","gray","black","red","purple","blue"),
       legend=c("Rademacher","Two-point","N(0,1)","Poisson(1)-1","NN","Exp(1)-1"))





#t=0.5
par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(1,2))
plot(density(G_center_beta_N01[,91]),xlim=c(-1,1),ylim=c(0,5))
lines(density(G_center_beta_N11[,91]),col="red")
lines(density(G_center_beta_E1[,91]),col="blue")
lines(density(G_center_beta_N02[,91]),col="green")
lines(density(G_center_beta_N12[,91]),col="red",lty=2)
lines(density(G_center_beta_G[,91]),col="blue",lty=2)
lines(density(G_center_beta_B[,91]),col="green",lty=2)
lines(density(G_center_beta_1[,91]),col="black",lty=2)
lines(density(G_center_beta_NN[,91]),col="gray",lty=2)
lines(density(G_center_beta_22[,91]),col="gray",lty=1)
legend("topright", lty=c(1,1,1,1,2,2,2,2,2,1), col=c("black","red","blue","green","red","blue","green","black","gray","gray"),
       legend=c("N(0,1)","N(1,1)","EXP(1)","N(0,2)","N(1,2)","Ga(0.25,0.5)","Be(sqrt(2)-1,1)","Rademacher","NN","two-point"))
legend("topleft",legend = c("beta","t=0.5"))

plot(density(G_center_gamma_N01[,91]),xlim=c(-1,1),ylim=c(0,7))
lines(density(G_center_gamma_N11[,91]),col="red")
lines(density(G_center_gamma_E1[,91]),col="blue")
lines(density(G_center_gamma_N02[,91]),col="green")
lines(density(G_center_gamma_N12[,91]),col="red",lty=2)
lines(density(G_center_gamma_G[,91]),col="blue",lty=2)
lines(density(G_center_gamma_B[,91]),col="green",lty=2)
lines(density(G_center_gamma_1[,91]),col="black",lty=2)
lines(density(G_center_gamma_NN[,91]),col="gray",lty=2)
lines(density(G_center_gamma_22[,91]),col="gray",lty=1)
legend("topright", lty=c(1,1,1,1,2,2,2,2,2,1), col=c("black","red","blue","green","red","blue","green","black","gray","gray"),
       legend=c("N(0,1)","N(1,1)","EXP(1)","N(0,2)","N(1,2)","Ga(0.25,0.5)","Be(sqrt(2)-1,1)","Rademacher","NN","two-point"))
legend("topleft",legend = c("gamma","t=0.5"))







#t=0.4
par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(1,2))
plot(density(G_center_beta_N01[,71]),xlim=c(-1,1),ylim=c(0,5))
lines(density(G_center_beta_N11[,71]),col="red")
lines(density(G_center_beta_E1[,71]),col="blue")
lines(density(G_center_beta_N02[,71]),col="green")
lines(density(G_center_beta_N12[,71]),col="red",lty=2)
lines(density(G_center_beta_G[,71]),col="blue",lty=2)
lines(density(G_center_beta_B[,71]),col="green",lty=2)
lines(density(G_center_beta_1[,71]),col="black",lty=2)
lines(density(G_center_beta_NN[,71]),col="gray",lty=2)
lines(density(G_center_beta_22[,71]),col="gray",lty=1)
legend("topright", lty=c(1,1,1,1,2,2,2,2,2,1), col=c("black","red","blue","green","red","blue","green","black","gray","gray"),
       legend=c("N(0,1)","N(1,1)","EXP(1)","N(0,2)","N(1,2)","Ga(0.25,0.5)","Be(sqrt(2)-1,1)","Rademacher","NN","two-point"))
legend("topleft",legend = c("beta","t=0.4"))

plot(density(G_center_gamma_N01[,71]),xlim=c(-1,1),ylim=c(0,8))
lines(density(G_center_gamma_N11[,71]),col="red")
lines(density(G_center_gamma_E1[,71]),col="blue")
lines(density(G_center_gamma_N02[,71]),col="green")
lines(density(G_center_gamma_N12[,71]),col="red",lty=2)
lines(density(G_center_gamma_G[,71]),col="blue",lty=2)
lines(density(G_center_gamma_B[,71]),col="green",lty=2)
lines(density(G_center_gamma_1[,71]),col="black",lty=2)
lines(density(G_center_gamma_NN[,71]),col="gray",lty=2)
lines(density(G_center_gamma_22[,71]),col="gray",lty=1)
legend("topright", lty=c(1,1,1,1,2,2,2,2,2,1), col=c("black","red","blue","green","red","blue","green","black","gray","gray"),
       legend=c("N(0,1)","N(1,1)","EXP(1)","N(0,2)","N(1,2)","Ga(0.25,0.5)","Be(sqrt(2)-1,1)","Rademacher","NN","two-point"))
legend("topleft",legend = c("gamma","t=0.4"))




#t=0.6
par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(1,2))
plot(density(G_center_beta_N01[,111]),xlim=c(-1,1),ylim=c(0,6))
lines(density(G_center_beta_N11[,111]),col="red")
lines(density(G_center_beta_E1[,111]),col="blue")
lines(density(G_center_beta_N02[,111]),col="green")
lines(density(G_center_beta_N12[,111]),col="red",lty=2)
lines(density(G_center_beta_G[,111]),col="blue",lty=2)
lines(density(G_center_beta_B[,111]),col="green",lty=2)
lines(density(G_center_beta_1[,111]),col="black",lty=2)
lines(density(G_center_beta_NN[,111]),col="gray",lty=2)
lines(density(G_center_beta_22[,111]),col="gray",lty=1)
legend("topright", lty=c(1,1,1,1,2,2,2,2,2,1), col=c("black","red","blue","green","red","blue","green","black","gray","gray"),
       legend=c("N(0,1)","N(1,1)","EXP(1)","N(0,2)","N(1,2)","Ga(0.25,0.5)","Be(sqrt(2)-1,1)","Rademacher","NN","two-point"))
legend("topleft",legend = c("beta","t=0.6"))

plot(density(G_center_gamma_N01[,111]),xlim=c(-1,1),ylim=c(0,8))
lines(density(G_center_gamma_N11[,111]),col="red")
lines(density(G_center_gamma_E1[,111]),col="blue")
lines(density(G_center_gamma_N02[,111]),col="green")
lines(density(G_center_gamma_N12[,111]),col="red",lty=2)
lines(density(G_center_gamma_G[,111]),col="blue",lty=2)
lines(density(G_center_gamma_B[,111]),col="green",lty=2)
lines(density(G_center_gamma_1[,111]),col="black",lty=2)
lines(density(G_center_gamma_NN[,111]),col="gray",lty=2)
lines(density(G_center_gamma_22[,111]),col="gray",lty=1)
legend("topright", lty=c(1,1,1,1,2,2,2,2,2,1), col=c("black","red","blue","green","red","blue","green","black","gray","gray"),
       legend=c("N(0,1)","N(1,1)","EXP(1)","N(0,2)","N(1,2)","Ga(0.25,0.5)","Be(sqrt(2)-1,1)","Rademacher","NN","two-point"))
legend("topleft",legend = c("gamma","t=0.6"))

