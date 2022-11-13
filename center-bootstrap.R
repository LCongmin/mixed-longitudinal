###############################################################################################
###     ������CP,�ֿ�ͬʱ���й��ƣ��̶�����h=n^(-0.8),h1=h2=n^(-0.5),t=(0.05,0.95)          ###
###############################################################################################
rm(list = ls (all = TRUE))
yxt<-proc.time()
library(MASS)
p<-seq(0.05,0.95,0.005)
#����������beta�Ĺ���
betat.hat_p<-rep(0,length(p))
bias.beta_p<-rep(0,length(p))
mse.beta_p<-rep(0,length(p))
SD.beta_p<-rep(0,length(p))
SE.beta_p<-rep(0,length(p))
cum_beta_p1<-rep(0,length(p))
hv_beta_p<-rep(0,length(p))
L_beta_p<-rep(0,length(p))
R_beta_p<-rep(0,length(p))

#��������gamma
gammat.hat_p<-rep(0,length(p))
bias.gamma_p<-rep(0,length(p))
mse.gamma_p<-rep(0,length(p))
SD.gamma_p<-rep(0,length(p))
SE.gamma_p<-rep(0,length(p))
cum_gamma_p1<-rep(0,length(p))
hv_gamma_p<-rep(0,length(p))
L_gamma_p<-rep(0,length(p))
R_gamma_p<-rep(0,length(p))

#���Ļ�beta
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



#���Ļ�gamma
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

n<-900
K<-1000#ѭ������
h<-n^{-0.5}
h1<-n^{-0.5}
h2<-n^{-0.5}

#������������
bi_beta_p<-matrix(0, nrow = length(p), ncol = K)
bi_beta_center<-matrix(0, nrow = length(p), ncol = K)
bi_beta_kw<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_p<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_center<-matrix(0, nrow = length(p), ncol = K)
bi_gamma_kw<-matrix(0, nrow = length(p), ncol = K)


#�洢ÿ��ʱ����ÿһ��ģ����90*K�ľ���
#����������beta
hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta����ֵ
hatbeta_var_p<-matrix(0,ncol = K,nrow = length(p))#beta����
L_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatbeta_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatbeta_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
abs_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))#����ֵ|hatbeta-beta|/sqrt(var)
L_scb_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatbeta_p<-matrix(0,ncol = K,nrow = length(p))
#����������gamma
hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta����ֵ
hatgamma_var_p<-matrix(0,ncol = K,nrow = length(p))#beta����
L_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatgamma_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatgamma_p_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatgamma_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatgamma_p_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
abs_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))#����ֵ|hatbeta-beta|/sqrt(var)
L_scb_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_p<-matrix(0,ncol = K,nrow = length(p))

#���Ļ�������beta�Ĺ���
hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta����ֵ
hatbeta_var_center<-matrix(0,ncol = K,nrow = length(p))#beta����
L_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_scb_hatbeta_center<-matrix(0,ncol=K,nrow=length(p))
R_scb_hatbeta_center<-matrix(0,ncol=K,nrow=length(p))
L_hatbeta_center_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_center_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatbeta_center_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_center_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
abs_hatbeta_center<-matrix(0,ncol = K,nrow = length(p))#����ֵ|hatbeta-beta|/sqrt(var)

#���Ļ�������gamma�Ĺ���
hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma����ֵ
hatgamma_var_center<-matrix(0,ncol = K,nrow = length(p))#gamma����
L_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma95%�½�
R_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#gamma95%�Ͻ�
L_hatgamma_center_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%�½�
R_hatgamma_center_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%�Ͻ�
L_hatgamma_center_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%�½�
R_hatgamma_center_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%�Ͻ�
abs_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))#����ֵ|hatgamma-gamma|/sqrt(var)
L_scb_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_center<-matrix(0,ncol = K,nrow = length(p))

#�˼�Ȩ������beta�Ĺ���
hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta����ֵ
hatbeta_var_kw<-matrix(0,ncol = K,nrow = length(p))#beta����
hatbetadot_var_kw<-matrix(0,ncol = K,nrow = length(p))#beta����
L_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatbeta_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
L_hatbeta_kw_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�½�
R_hatbeta_kw_sche<-matrix(0,ncol = K,nrow = length(p))#beta95%�Ͻ�
abs_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))#����ֵ|hatbeta-beta|/sqrt(var)
L_scb_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatbeta_kw<-matrix(0,ncol = K,nrow = length(p))

#�˼�Ȩ������gamma�Ĺ���
hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma����ֵ
hatgamma_var_kw<-matrix(0,ncol = K,nrow = length(p))#gamma����
hatgammadot_var_kw<-matrix(0,ncol = K,nrow = length(p))#gamma����
L_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma95%�½�
R_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#gamma95%�Ͻ�
L_hatgamma_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%�½�
R_hatgamma_kw_bonf<-matrix(0,ncol = K,nrow = length(p))#gamma95%�Ͻ�
L_hatgamma_kw_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%�½�
R_hatgamma_kw_sche<-matrix(0,ncol = K,nrow = length(p))#gamma95%�Ͻ�
abs_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))#����ֵ|hatgamma-gamma|/sqrt(var)
L_scb_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))
R_scb_hatgamma_kw<-matrix(0,ncol = K,nrow = length(p))

#MB
K_simu<-10000
abs_G_beta_center<-rep(0,K_simu)
abs_G_gamma_center<-rep(0,K_simu)


beta_center<-array(0,c(2,1,length(p)))
gamma_center<-array(0,c(2,1,length(p)))
G_center_beta<-matrix(0,ncol = length(p),nrow = K_simu)
G_center_gamma<-matrix(0,ncol = length(p),nrow = K_simu)


for (t_ind in 1:length(p)) {
  #print(paste(t_ind))
  t<-0.05+0.005*(t_ind-1)
  
  T.x<-list()#����X��Y�Ĺ۲�ʱ�䣬����ÿ��subject�Ĳ��ȳ��ȹ۲�ʱ�������ļ��ϣ�����Ϊ�κ��������������������飬��ͬ���ȵ��б�
  S.z<-list()#����ʵ��Z��n�����ȳ��ȹ۲�ʱ�������ļ���,��������Yʱ����Ҫһ��ʹ��XY�Ĺ۲�ʱ��ȥ����
  kw<-list()
  mu_x<-list()#x��Ԫ��̬��ֵ��n�����ȳ��������ĺϼ�
  sigma_x<-list(matrix)#x��Ԫ��̬sigma��n������ĺϼ�
  mu_z.x<-list()#z.x��Ԫ��̬��ֵ��n�����ȳ��������ĺϼ�
  sigma_z.x<-list(matrix)#z.x��Ԫ��̬sigma��n������ĺϼ�
  mu_z<-list()#z��Ԫ��̬��ֵ��n�����ȳ��������ĺϼ�
  sigma_z<-list(matrix)#z��Ԫ��̬sigma��n������ĺϼ�
  mu_epsilon<-list()
  sigma_epsilon<-list(matrix)
  x<-list()
  x.tilde<-list()#�������Ļ�֮�������x
  z.x<-list()#������x�۲�ʱ�����ɵ�����z��������������yʹ��
  z<-list()#����ʵ�ʰ�z�۲�ʱ�����ɵ�����z������ʵ��ģ��ʹ��
  y<-list()
  y.tilde<-list()#�������Ļ�֮�������y
  epsilon<-list()
  v<-list()
  nx<-rpois(n,5)+1#the number of each subject,��ÿ���۲����Ĺ۲�����γɵ�һ������
  nz<-rpois(n,5)+1
  #����Э����x��z�Ĺ۲�ʱ��
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
    #���������������������X��Z��epsilon
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
  
  #�õ�����x�ĺ˹��ƣ�����Ӧ�������и������һ��ĺ˹��ưɣ�
  #kw������Ȩ��
  for (i in 1:n) {
    kw[[i]]<-rep(0,nx[i])
  }
  
  for (i in 1:n) {
    for (j in 1:nx[i]) {
      kw[[i]][j]<-local_kernel(T.x[[i]][j]-t,h)
    }
  } 
  
  #�Ƚ����к�Ȩ�����
  kw.all<-0
  for (i in 1:n) {
    kw.all<-kw.all+sum(kw[[i]])
  }
  
  km.x<-0#������������x�ĺ˹���
  for (i in 1:n) {
    km.x<-km.x+sum(x[[i]]*kw[[i]]/kw.all)#ÿ������ľ�ֵ�˹���
  }
  
  #���Ļ�֮���x�����ݣ���Ϊx.tilde
  for (i in 1:n) {
    x.tilde[[i]]<-x[[i]]-km.x
  }
  
  #�õ�����y�ĺ˹��ƣ�����Ӧ�������и������һ��ĺ˹��ưɣ�
  km.y<-0#������������y�ĺ˹���
  for (i in 1:n) {
    km.y<-km.y+sum(y[[i]]*kw[[i]])/kw.all#ÿ������ľ�ֵ�˹���
  }
  
  
  
  
  #���Ļ�֮���x�����ݣ���Ϊx.tilde
  for (i in 1:n) {
    y.tilde[[i]]<-y[[i]]-km.y
  }
  
  
  
  #���Ļ���������x��ϵ��beta
  Q_beta_center<-matrix(0:0,ncol=2,nrow=2)#���ϵ����2��2����
  q_beta_center<-matrix(0:0,ncol = 1,nrow = 2)#���ϵ����2��1����
  for (i in 1:n) {
    q_beta_center<-q_beta_center+matrix(c(sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]])/n,
                                          sum(kw[[i]]*x.tilde[[i]]*y.tilde[[i]]*(T.x[[i]]-t))/n),ncol = 1,nrow = 2)
    Q_beta_center<-Q_beta_center+matrix(c(sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]])/n,
                                          sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                                          sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t))/n,
                                          sum(kw[[i]]*x.tilde[[i]]*x.tilde[[i]]*(T.x[[i]]-t)*(T.x[[i]]-t))/n),
                                        ncol = 2,nrow = 2)
  }
  beta_center[,,t_ind]<-ginv(Q_beta_center)%*%q_beta_center#l�����ڼ���ģ��
  
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
    tau<-rbinom(n,1,0.5)
    tau[which(tau==0)]=-1
    G_center_beta[k_simu,t_ind]<-ginv(Q_beta_center[1,1])*sum(tau*sig_beta)/n
    G_center_gamma[k_simu,t_ind]<-ginv(Q_center[1,1])*sum(tau*p1_1)/n
      
    #}
  }
}


for (k_simu in 1:K_simu) {
  abs_G_beta_center[k_simu]<-max(abs(G_center_beta[k_simu,]))
  abs_G_gamma_center[k_simu]<-max(abs(G_center_gamma[k_simu,]))
}

quantile(abs_G_beta_center,0.95)
quantile(abs_G_gamma_center,0.95)

#a<-quantile(abs_G_beta_center,0.95)
#b<-quantile(abs_G_gamma_center,0.95)
#result<-rbind(a,b)

#write.csv(result,file = "E:/papercao/2021.2.17center3000/center_line_h0.0150.09_1.csv")

