#model: mmse=alpha(t)+gamma(t)*fa+beta1(t)*age+beta2(t)*ad+beta3(t)*edu+beta4(t)*apeo12
#带有年龄变量
rm(list = ls (all = TRUE))
library(MASS)

dataset2<-read.csv(choose.files())#FA_hazard_WB
fa<-dataset2[,38]#u=0.65

dataset1<-read.csv(choose.files())#z_covariate_apoe12
#作为同步的X
z_id<-dataset1[,1]
z_time<-dataset1[,2]
mci<-dataset1[,4]
ad<-dataset1[,5]
edu<-dataset1[,7]
apoe12<-dataset1[,8]+dataset1[,9]


dataset<-read.csv(choose.files())#mmse,读入文档y.csv
#MMSE作为非同步y
y_id<-dataset[,1]
unique(y_id)
length(unique(y_id))
y_time<-dataset[,2]
y_value<-dataset[,3]
y_age<-dataset[,4]

zz_id<-unique(z_id)
length(zz_id)

yy_id<-unique(y_id)
length(yy_id)

n<-256

count_zz<-rep(0,n)
count_yy<-rep(0,n)

for (i in 1:n) {
  count_zz[i]<-length(which(z_id==zz_id[i]))
  count_yy[i]<-length(which(y_id==yy_id[i]))
}


count_zz
count_yy

z_fa<-list()
z_fa[[1]]<-fa[1:count_zz[1]]

x_age<-list()
x_age[[1]]<-y_age[1:count_yy[1]]

x_mci<-list()
x_mci[[1]]<-rep(mci[1],count_yy[1])

x_ad<-list()
x_ad[[1]]<-rep(ad[1],count_yy[1])

x_edu<-list()
x_edu[[1]]<-rep(edu[1],count_yy[1])

x_apoe12<-list()
x_apoe12[[1]]<-rep(apoe12[1],count_yy[1])


y_mmse<-list()
y_mmse[[1]]<-y_value[1:count_yy[1]]

T.y<-list()
T.y[[1]]<-y_time[1:count_yy[1]]

T.z<-list()
T.z[[1]]<-z_time[1:count_zz[1]]

for (i in 2:n) {
  z_fa[[i]]<-fa[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))]
  
  x_age[[i]]<-y_age[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  x_mci[[i]]<-rep(mci[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_ad[[i]]<-rep(ad[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_edu[[i]]<-rep(edu[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_apoe12[[i]]<-rep(apoe12[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  
  y_mmse[[i]]<-y_value[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  
  T.y[[i]]<-y_time[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  T.z[[i]]<-z_time[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))]
}




#首先将所有非01变量标准化，包括age,edu,fa065,y_mmse
epanechnikov <- function(t){
  
  tst <- (-1.0 <= t) & (t <= 1.0 )
  
  kt <- 0.75*( 1.0 - t*t )
  kt[ !tst ] <- 0.0
  
  return(kt)
}

local_kernel <- function(t, h){
  kt <- epanechnikov(t/h)
  return(kt/h)
}


#kw变系数模型
a=quantile(c(unlist(T.y),unlist(T.z)))
Q3=a[[4]]
Q1=a[[2]]
#2*(Q3-Q1)*n^{-1/6}
h1_kw<-0.144
h2_kw<-0.144

o<-seq(0.15,0.9,0.01)
theta_kw<-matrix(0, nrow = length(o), ncol = 10)
hatvar_kw<-array(0,dim = c(10,10,length(o)))
diag_hatvar_theta<-matrix(0, nrow = length(o), ncol = 10)

K_simu<-3000

L_alpha_kw<-rep(0,length(o))
R_alpha_kw<-rep(0,length(o))

L_age_kw<-rep(0,length(o))
R_age_kw<-rep(0,length(o))

L_edu_kw<-rep(0,length(o))
R_edu_kw<-rep(0,length(o))

L_apoe12_kw<-rep(0,length(o))
R_apoe12_kw<-rep(0,length(o))

L_fa_kw<-rep(0,length(o))
R_fa_kw<-rep(0,length(o))


se_kw<-matrix(0,ncol=10,nrow=length(o))

G_kw_alpha<-matrix(0,ncol=length(o),nrow=K_simu)
G_kw_age<-matrix(0,ncol=length(o),nrow=K_simu)
G_kw_edu<-matrix(0,ncol=length(o),nrow=K_simu)
G_kw_apoe12<-matrix(0,ncol=length(o),nrow=K_simu)
G_kw_fa<-matrix(0,ncol=length(o),nrow=K_simu)



z_fa<-list()
z_fa[[1]]<-fa[1:count_zz[1]]

x_age<-list()
x_age[[1]]<-y_age[1:count_yy[1]]

x_mci<-list()
x_mci[[1]]<-rep(mci[1],count_yy[1])

x_ad<-list()
x_ad[[1]]<-rep(ad[1],count_yy[1])

x_edu<-list()
x_edu[[1]]<-rep(edu[1],count_yy[1])

x_apoe12<-list()
x_apoe12[[1]]<-rep(apoe12[1],count_yy[1])


y_mmse<-list()
y_mmse[[1]]<-y_value[1:count_yy[1]]

T.y<-list()
T.y[[1]]<-y_time[1:count_yy[1]]

T.z<-list()
T.z[[1]]<-z_time[1:count_zz[1]]

for (i in 2:n) {
  z_fa[[i]]<-fa[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))]
  
  x_age[[i]]<-y_age[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  x_mci[[i]]<-rep(mci[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_ad[[i]]<-rep(ad[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_edu[[i]]<-rep(edu[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_apoe12[[i]]<-rep(apoe12[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  y_mmse[[i]]<-y_value[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  
  T.y[[i]]<-y_time[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  T.z[[i]]<-z_time[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))]
}


#2*(Q3-Q1)*n^{-0.21}
h=0.108
x.tilde_age<-list()
x.tilde_edu<-list()
z.tilde_fa<-list()
y.tilde<-list()

Eage<-list()
Efa<-list()
Ey<-list()
Eageage<-list()
Efafa<-list()
Eyy<-list()
Eedu=0
Eeduedu=0
for (i in 1:n) {
  Eedu=Eedu+x_edu[[i]][1]/n
  Eeduedu=Eeduedu+(x_edu[[i]][1]*x_edu[[i]][1])/n
}

#标准化
for (i in 1:n) {
  x.tilde_edu[[i]]<-(x_edu[[i]]-Eedu)/sqrt(Eeduedu-Eedu^2)
  
  Eage[[i]]<-rep(0,count_yy[i])
  Ey[[i]]<-rep(0,count_yy[i])
  
  Eageage[[i]]<-rep(0,count_yy[i])
  Eyy[[i]]<-rep(0,count_yy[i])
  
  x.tilde_age[[i]]<-rep(0,count_yy[i])
  y.tilde[[i]]<-rep(0,count_yy[i])
  for (j in 1:count_yy[i]) {
    t=T.y[[i]][j]
    Eage[[i]][j]=sum(unlist(x_age)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    Eageage[[i]][j]=sum(unlist(x_age)*unlist(x_age)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    Ey[[i]][j]=sum(unlist(y_mmse)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    Eyy[[i]][j]=sum(unlist(y_mmse)*unlist(y_mmse)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    
    x.tilde_age[[i]][j]<-(x_age[[i]][j]-Eage[[i]][j])/sqrt(Eageage[[i]][j]-Eage[[i]][j]^2)
    y.tilde[[i]][j]<-(y_mmse[[i]][j]-Ey[[i]][j])/sqrt(Eyy[[i]][j]-Ey[[i]][j]^2)
    
  }
} 


#t2<-unlist(T.z)
for (i in 1:n) {
  Efa[[i]]<-rep(0,count_zz[i])
  Efafa[[i]]<-rep(0,count_zz[i])
  z.tilde_fa[[i]]<-rep(0,count_zz[i])
  for (j in 1:count_zz[i]) {
    t=T.z[[i]][j]
    Efa[[i]][j]=sum(unlist(z_fa)*local_kernel(unlist(T.z)-t,h))/sum(local_kernel(unlist(T.z)-t,h))
    Efafa[[i]][j]=sum(unlist(z_fa)*unlist(z_fa)*local_kernel(unlist(T.z)-t,h))/sum(local_kernel(unlist(T.z)-t,h))
    z.tilde_fa[[i]][j]<-(z_fa[[i]][j]-Efa[[i]][j])/sqrt(Efafa[[i]][j]-Efa[[i]][j]^2)
  }
} 
#为了计算方便，对标准化的变量重新赋值
x_age<-list()
x_edu<-list()
z_fa<-list()
y_mmse<-list()

for (i in 1:n) {
  x_age[[i]]<-x.tilde_age[[i]]
  x_edu[[i]]<-x.tilde_edu[[i]]
  z_fa[[i]]<-z.tilde_fa[[i]]
  y_mmse[[i]]<-y.tilde[[i]]
}


for (t_ind in 1:length(o)) {
  
  t<-0.15+(t_ind-1)/100
  print(paste(t))
  
  #核加权方法估计
  ker<-list(matrix)
  ker_t1<-list(matrix)
  
  ker_x1_age<-list(matrix)
  ker_x2_age<-list(matrix)
  ker_x1t1_age<-list(matrix)
  ker_x2t1_age<-list(matrix)
  ker_x2t2_age<-list(matrix)
  
  ker_x1_edu<-list(matrix)
  ker_x2_edu<-list(matrix)
  ker_x1t1_edu<-list(matrix)
  ker_x2t1_edu<-list(matrix)
  ker_x2t2_edu<-list(matrix)
  
  ker_x1_apoe12<-list(matrix)
  ker_x2_apoe12<-list(matrix)
  ker_x1t1_apoe12<-list(matrix)
  ker_x2t1_apoe12<-list(matrix)
  ker_x2t2_apoe12<-list(matrix)
  
  ker_y<-list(matrix)
  
  ker_xy_age<-list(matrix)
  ker_xyt_age<-list(matrix)
  
  ker_xy_edu<-list(matrix)
  ker_xyt_edu<-list(matrix)
  
  ker_xy_apoe12<-list(matrix)
  ker_xyt_apoe12<-list(matrix)
  
  for (i in 1:n) {
    ker[[i]]<-matrix(rep(local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_x1_age[[i]]<-matrix(rep(x_age[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2_age[[i]]<-matrix(rep(x_age[[i]]*x_age[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x1t1_age[[i]]<-matrix(rep((T.y[[i]]-t)*x_age[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2t1_age[[i]]<-matrix(rep((T.y[[i]]-t)*x_age[[i]]*x_age[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2t2_age[[i]]<-matrix(rep((T.y[[i]]-t)*(T.y[[i]]-t)*x_age[[i]]*x_age[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_x1_edu[[i]]<-matrix(rep(x_edu[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2_edu[[i]]<-matrix(rep(x_edu[[i]]*x_edu[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x1t1_edu[[i]]<-matrix(rep((T.y[[i]]-t)*x_edu[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2t1_edu[[i]]<-matrix(rep((T.y[[i]]-t)*x_edu[[i]]*x_edu[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2t2_edu[[i]]<-matrix(rep((T.y[[i]]-t)*(T.y[[i]]-t)*x_edu[[i]]*x_edu[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_x1_apoe12[[i]]<-matrix(rep(x_apoe12[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2_apoe12[[i]]<-matrix(rep(x_apoe12[[i]]*x_apoe12[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x1t1_apoe12[[i]]<-matrix(rep((T.y[[i]]-t)*x_apoe12[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2t1_apoe12[[i]]<-matrix(rep((T.y[[i]]-t)*x_apoe12[[i]]*x_apoe12[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_x2t2_apoe12[[i]]<-matrix(rep((T.y[[i]]-t)*(T.y[[i]]-t)*x_apoe12[[i]]*x_apoe12[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_y[[i]]<-matrix(rep(y_mmse[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_xy_age[[i]]<-matrix(rep(y_mmse[[i]]*x_age[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_xyt_age[[i]]<-matrix(rep(y_mmse[[i]]*x_age[[i]]*(T.y[[i]]-t)*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_xy_edu[[i]]<-matrix(rep(y_mmse[[i]]*x_edu[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_xyt_edu[[i]]<-matrix(rep(y_mmse[[i]]*x_edu[[i]]*(T.y[[i]]-t)*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    
    ker_xy_apoe12[[i]]<-matrix(rep(y_mmse[[i]]*x_apoe12[[i]]*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
    ker_xyt_apoe12[[i]]<-matrix(rep(y_mmse[[i]]*x_apoe12[[i]]*(T.y[[i]]-t)*local_kernel(T.y[[i]]-t,h1_kw),count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])
  }
  
  Q_kw<-matrix(0,ncol = 10,nrow = 10)
  q_kw<-matrix(0,ncol = 1,nrow = 10)
  for (i in 1:n) {
    Q_kw[1,1]<-Q_kw[1,1]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,2]<-Q_kw[1,2]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[1,3]<-Q_kw[1,3]+colSums(ker_x1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,4]<-Q_kw[1,4]+colSums(ker_x1t1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,5]<-Q_kw[1,5]+colSums(ker_x1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,6]<-Q_kw[1,6]+colSums(ker_x1t1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,7]<-Q_kw[1,7]+colSums(ker_x1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,8]<-Q_kw[1,8]+colSums(ker_x1t1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[1,9]<-Q_kw[1,9]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[1,10]<-Q_kw[1,10]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[2,2]<-Q_kw[2,2]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t)*(T.z[[i]]-t))))/n
    Q_kw[2,3]<-Q_kw[2,3]+colSums(ker_x1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[2,4]<-Q_kw[2,4]+colSums(ker_x1t1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[2,5]<-Q_kw[2,5]+colSums(ker_x1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[2,6]<-Q_kw[2,6]+colSums(ker_x1t1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[2,7]<-Q_kw[2,7]+colSums(ker_x1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[2,8]<-Q_kw[2,8]+colSums(ker_x1t1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    Q_kw[2,9]<-Q_kw[2,9]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t)*z_fa[[i]])))/n
    Q_kw[2,10]<-Q_kw[2,10]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t)*z_fa[[i]]*(T.z[[i]]-t))))/n
    
    Q_kw[3,3]<-Q_kw[3,3]+colSums((ker_x1_age[[i]]*x_age[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[3,4]<-Q_kw[3,4]+colSums((ker_x1_age[[i]]*x_age[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[3,5]<-Q_kw[3,5]+colSums((ker_x1_age[[i]]*x_edu[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[3,6]<-Q_kw[3,6]+colSums((ker_x1_age[[i]]*x_edu[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[3,7]<-Q_kw[3,7]+colSums((ker_x1_age[[i]]*x_apoe12[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[3,8]<-Q_kw[3,8]+colSums((ker_x1_age[[i]]*x_apoe12[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[3,9]<-Q_kw[3,9]+colSums(ker_x1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[3,10]<-Q_kw[3,10]+colSums(ker_x1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[4,4]<-Q_kw[4,4]+colSums((ker_x1t1_age[[i]]*x_age[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[4,5]<-Q_kw[4,5]+colSums((ker_x1t1_age[[i]]*x_edu[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[4,6]<-Q_kw[4,6]+colSums((ker_x1t1_age[[i]]*x_edu[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[4,7]<-Q_kw[4,7]+colSums((ker_x1t1_age[[i]]*x_apoe12[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[4,8]<-Q_kw[4,8]+colSums((ker_x1t1_age[[i]]*x_apoe12[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[4,9]<-Q_kw[4,9]+colSums(ker_x1t1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[4,10]<-Q_kw[4,10]+colSums(ker_x1t1_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[5,5]<-Q_kw[5,5]+colSums((ker_x1_edu[[i]]*x_edu[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[5,6]<-Q_kw[5,6]+colSums((ker_x1_edu[[i]]*x_edu[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[5,7]<-Q_kw[5,7]+colSums((ker_x1_edu[[i]]*x_apoe12[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[5,8]<-Q_kw[5,8]+colSums((ker_x1_edu[[i]]*x_apoe12[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[5,9]<-Q_kw[5,9]+colSums(ker_x1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[5,10]<-Q_kw[7,10]+colSums(ker_x1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[6,6]<-Q_kw[6,6]+colSums((ker_x1t1_edu[[i]]*x_edu[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[6,7]<-Q_kw[6,7]+colSums((ker_x1t1_edu[[i]]*x_apoe12[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[6,8]<-Q_kw[6,8]+colSums((ker_x1t1_edu[[i]]*x_apoe12[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[6,9]<-Q_kw[6,9]+colSums(ker_x1t1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[6,10]<-Q_kw[6,10]+colSums(ker_x1t1_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[7,7]<-Q_kw[7,7]+colSums((ker_x1_apoe12[[i]]*x_apoe12[[i]])%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[7,8]<-Q_kw[7,8]+colSums((ker_x1_apoe12[[i]]*x_apoe12[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[7,9]<-Q_kw[7,9]+colSums(ker_x1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[7,10]<-Q_kw[7,10]+colSums(ker_x1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[8,8]<-Q_kw[8,8]+colSums((ker_x1t1_apoe12[[i]]*x_apoe12[[i]]*(T.y[[i]]-t))%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    Q_kw[8,9]<-Q_kw[8,9]+colSums(ker_x1t1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]]))))/n
    Q_kw[8,10]<-Q_kw[8,10]+colSums(ker_x1t1_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(z_fa[[i]])*(T.z[[i]]-t))))/n
    
    Q_kw[9,9]<-Q_kw[9,9]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]*z_fa[[i]])))/n
    Q_kw[9,10]<-Q_kw[9,10]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]*z_fa[[i]]*(T.z[[i]]-t))))/n
    
    Q_kw[10,10]<-Q_kw[10,10]+colSums(ker[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]*z_fa[[i]]*(T.z[[i]]-t)*(T.z[[i]]-t))))/n
    
    q_kw[1,]<-q_kw[1,]+colSums(ker_y[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[2,]<-q_kw[2,]+colSums(ker_y[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t))))/n
    q_kw[3,]<-q_kw[3,]+colSums(ker_xy_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[4,]<-q_kw[4,]+colSums(ker_xyt_age[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[5,]<-q_kw[5,]+colSums(ker_xy_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[6,]<-q_kw[6,]+colSums(ker_xyt_edu[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[7,]<-q_kw[7,]+colSums(ker_xy_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[8,]<-q_kw[8,]+colSums(ker_xyt_apoe12[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw))))/n
    q_kw[9,]<-q_kw[9,]+colSums(ker_y[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]])))/n
    q_kw[10,]<-q_kw[10,]+colSums(ker_y[[i]]%*%t(t(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]*(T.z[[i]]-t))))/n
    
  }
  Q_kw[2,1]<-Q_kw[1,2]
  Q_kw[3,1]<-Q_kw[1,3]
  Q_kw[4,1]<-Q_kw[1,4]
  Q_kw[5,1]<-Q_kw[1,5]
  Q_kw[6,1]<-Q_kw[1,6]
  Q_kw[7,1]<-Q_kw[1,7]
  Q_kw[8,1]<-Q_kw[1,8]
  Q_kw[9,1]<-Q_kw[1,9]
  Q_kw[10,1]<-Q_kw[1,10]
  
  Q_kw[3,2]<-Q_kw[2,3]
  Q_kw[4,2]<-Q_kw[2,4]
  Q_kw[5,2]<-Q_kw[2,5]
  Q_kw[6,2]<-Q_kw[2,6]
  Q_kw[7,2]<-Q_kw[2,7]
  Q_kw[8,2]<-Q_kw[2,8]
  Q_kw[9,2]<-Q_kw[2,9]
  Q_kw[10,2]<-Q_kw[2,10]
  
  Q_kw[4,3]<-Q_kw[3,4]
  Q_kw[5,3]<-Q_kw[3,5]
  Q_kw[6,3]<-Q_kw[3,6]
  Q_kw[7,3]<-Q_kw[3,7]
  Q_kw[8,3]<-Q_kw[3,8]
  Q_kw[9,3]<-Q_kw[3,9]
  Q_kw[10,3]<-Q_kw[3,10]
  
  Q_kw[5,4]<-Q_kw[4,5]
  Q_kw[6,4]<-Q_kw[4,6]
  Q_kw[7,4]<-Q_kw[4,7]
  Q_kw[8,4]<-Q_kw[4,8]
  Q_kw[9,4]<-Q_kw[4,9]
  Q_kw[10,4]<-Q_kw[4,10]
  
  Q_kw[6,5]<-Q_kw[5,6]
  Q_kw[7,5]<-Q_kw[5,7]
  Q_kw[8,5]<-Q_kw[5,8]
  Q_kw[9,5]<-Q_kw[5,9]
  Q_kw[10,5]<-Q_kw[5,10]
  
  Q_kw[7,6]<-Q_kw[6,7]
  Q_kw[8,6]<-Q_kw[6,8]
  Q_kw[9,6]<-Q_kw[6,9]
  Q_kw[10,6]<-Q_kw[6,10]
  
  Q_kw[8,7]<-Q_kw[7,8]
  Q_kw[9,7]<-Q_kw[7,9]
  Q_kw[10,7]<-Q_kw[7,10]
  
  Q_kw[9,8]<-Q_kw[8,9]
  Q_kw[10,8]<-Q_kw[8,10]
  
  Q_kw[10,9]<-Q_kw[9,10]
  
  theta_kw[t_ind,]<-t(ginv(Q_kw)%*%q_kw)
  
  
  #kw方法的方差估计
  p_kw<-matrix(0:0,ncol = 10,nrow = 10)
  p_kw_alpha<-rep(0,n)
  p_kw_age<-rep(0,n)
  p_kw_edu<-rep(0,n)
  p_kw_apoe12<-rep(0,n)
  p_kw_fa<-rep(0,n)
  
  for (i in 1:n) {
    ppp=y_mmse[[i]]-x_age[[i]]*(theta_kw[t_ind,3]+theta_kw[t_ind,4]*(T.y[[i]]-t))-x_edu[[i]]*(theta_kw[t_ind,5]+theta_kw[t_ind,6]*(T.y[[i]]-t))-x_apoe12[[i]]*(theta_kw[t_ind,7]+theta_kw[t_ind,8]*(T.y[[i]]-t))
    p11<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p11<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p1<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p22<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p22<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)*(T.z[[i]]-t)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p2<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p33_age<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p33_age<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw)*x_age[[i]])%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p3_age<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p33_age.t<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p33_age.t<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw)*x_age[[i]]*(T.y[[i]]-t))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p3_age.t<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p66_edu<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p66_edu<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw)*x_edu[[i]])%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p6_edu<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p66_edu.t<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p66_edu.t<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw)*x_edu[[i]]*(T.y[[i]]-t))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p6_edu.t<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p77_apoe12<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p77_apoe12<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw)*x_apoe12[[i]])%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p7_apoe12<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p77_apoe12.t<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p77_apoe12.t<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw)*x_apoe12[[i]]*(T.y[[i]]-t))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p7_apoe12.t<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    
    p99_fa<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p99_fa<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p9_fa<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    p99_fa.t<-array(0,c(count_zz[i],count_yy[i],count_yy[i],count_zz[i]))
    p99_fa.t<-t(as.vector(local_kernel(T.y[[i]]-t,h1_kw))%o%as.vector(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]*(T.z[[i]]-t)))%o%(matrix(rep(ppp,count_zz[i]),ncol = count_zz[i],nrow = count_yy[i])-matrix(rep(theta_kw[t_ind,1]+theta_kw[t_ind,2]*(T.z[[i]]-t)+z_fa[[i]]*(theta_kw[t_ind,9]+theta_kw[t_ind,10]*(T.z[[i]]-t)),count_yy[i]),ncol = count_zz[i],nrow = count_yy[i],byrow = T))
    p9_fa.t<-matrix(0, ncol=count_yy[i], nrow=count_zz[i])
    
    
    for (j in 1:count_yy[i]) {
      for (k in 1:count_zz[i]) {
        if(is.matrix(p11[,,j,k]==TRUE)) {p1[k,j]<-p11[,,j,k][k,j]}
        else {p1[k,j]<-p11[,,j,k][k]}
        
        if(is.matrix(p22[,,j,k]==TRUE)) {p2[k,j]<-p22[,,j,k][k,j]}
        else {p2[k,j]<-p22[,,j,k][k]}
        
        if(is.matrix(p33_age[,,j,k]==TRUE)) {p3_age[k,j]<-p33_age[,,j,k][k,j]}
        else {p3_age[k,j]<-p33_age[,,j,k][k]}
        if(is.matrix(p33_age.t[,,j,k]==TRUE)) {p3_age.t[k,j]<-p33_age.t[,,j,k][k,j]}
        else {p3_age.t[k,j]<-p33_age.t[,,j,k][k]}
        
        
        if(is.matrix(p66_edu[,,j,k]==TRUE)) {p6_edu[k,j]<-p66_edu[,,j,k][k,j]}
        else {p6_edu[k,j]<-p66_edu[,,j,k][k]}
        if(is.matrix(p66_edu.t[,,j,k]==TRUE)) {p6_edu.t[k,j]<-p66_edu.t[,,j,k][k,j]}
        else {p6_edu.t[k,j]<-p66_edu.t[,,j,k][k]}
        
        if(is.matrix(p77_apoe12[,,j,k]==TRUE)) {p7_apoe12[k,j]<-p77_apoe12[,,j,k][k,j]}
        else {p7_apoe12[k,j]<-p77_apoe12[,,j,k][k]}
        if(is.matrix(p77_apoe12.t[,,j,k]==TRUE)) {p7_apoe12.t[k,j]<-p77_apoe12.t[,,j,k][k,j]}
        else {p7_apoe12.t[k,j]<-p77_apoe12.t[,,j,k][k]}
        
        if(is.matrix(p99_fa[,,j,k]==TRUE)) {p9_fa[k,j]<-p99_fa[,,j,k][k,j]}
        else {p9_fa[k,j]<-p99_fa[,,j,k][k]}
        if(is.matrix(p99_fa.t[,,j,k]==TRUE)) {p9_fa.t[k,j]<-p99_fa.t[,,j,k][k,j]}
        else {p9_fa.t[k,j]<-p99_fa.t[,,j,k][k]}
        
      }
    }
    p_kw[1,1]<-p_kw[1,1]+sum(p1%o%p1)/n
    p_kw[1,2]<-p_kw[1,2]+sum(p1%o%p2)/n
    p_kw[1,3]<-p_kw[1,3]+sum(p1%o%p3_age)/n
    p_kw[1,4]<-p_kw[1,4]+sum(p1%o%p3_age.t)/n
    p_kw[1,5]<-p_kw[1,5]+sum(p1%o%p6_edu)/n
    p_kw[1,6]<-p_kw[1,6]+sum(p1%o%p6_edu.t)/n
    p_kw[1,7]<-p_kw[1,7]+sum(p1%o%p7_apoe12)/n
    p_kw[1,8]<-p_kw[1,8]+sum(p1%o%p7_apoe12.t)/n
    p_kw[1,9]<-p_kw[1,9]+sum(p1%o%p9_fa)/n
    p_kw[1,10]<-p_kw[1,10]+sum(p1%o%p9_fa.t)/n
    
    p_kw[2,2]<-p_kw[2,2]+sum(p2%o%p2)/n
    p_kw[2,3]<-p_kw[2,3]+sum(p2%o%p3_age)/n
    p_kw[2,4]<-p_kw[2,4]+sum(p2%o%p3_age.t)/n
    p_kw[2,5]<-p_kw[2,5]+sum(p2%o%p6_edu)/n
    p_kw[2,6]<-p_kw[2,6]+sum(p2%o%p6_edu.t)/n
    p_kw[2,7]<-p_kw[2,7]+sum(p2%o%p7_apoe12)/n
    p_kw[2,8]<-p_kw[2,8]+sum(p2%o%p7_apoe12.t)/n
    p_kw[2,9]<-p_kw[2,9]+sum(p2%o%p9_fa)/n
    p_kw[2,10]<-p_kw[2,10]+sum(p2%o%p9_fa.t)/n
    
    p_kw[3,3]<-p_kw[3,3]+sum(p3_age%o%p3_age)/n
    p_kw[3,4]<-p_kw[3,4]+sum(p3_age%o%p3_age.t)/n
    p_kw[3,5]<-p_kw[3,5]+sum(p3_age%o%p6_edu)/n
    p_kw[3,6]<-p_kw[3,6]+sum(p3_age%o%p6_edu.t)/n
    p_kw[3,7]<-p_kw[3,7]+sum(p3_age%o%p7_apoe12)/n
    p_kw[3,8]<-p_kw[3,8]+sum(p3_age%o%p7_apoe12.t)/n
    p_kw[3,9]<-p_kw[3,9]+sum(p3_age%o%p9_fa)/n
    p_kw[3,10]<-p_kw[3,10]+sum(p3_age%o%p9_fa.t)/n
    
    p_kw[4,4]<-p_kw[4,4]+sum(p3_age.t%o%p3_age.t)/n
    p_kw[4,5]<-p_kw[4,5]+sum(p3_age.t%o%p6_edu)/n
    p_kw[4,6]<-p_kw[4,6]+sum(p3_age.t%o%p6_edu.t)/n
    p_kw[4,7]<-p_kw[4,7]+sum(p3_age.t%o%p7_apoe12)/n
    p_kw[4,8]<-p_kw[4,8]+sum(p3_age.t%o%p7_apoe12.t)/n
    p_kw[4,9]<-p_kw[4,9]+sum(p3_age.t%o%p9_fa)/n
    p_kw[4,10]<-p_kw[4,10]+sum(p3_age.t%o%p9_fa.t)/n
    
    p_kw[5,5]<-p_kw[5,5]+sum(p6_edu%o%p6_edu)/n
    p_kw[5,6]<-p_kw[5,6]+sum(p6_edu%o%p6_edu.t)/n
    p_kw[5,7]<-p_kw[5,7]+sum(p6_edu%o%p7_apoe12)/n
    p_kw[5,8]<-p_kw[5,8]+sum(p6_edu%o%p7_apoe12.t)/n
    p_kw[5,9]<-p_kw[5,9]+sum(p6_edu%o%p9_fa)/n
    p_kw[5,10]<-p_kw[5,10]+sum(p6_edu%o%p9_fa.t)/n
    
    p_kw[6,6]<-p_kw[6,6]+sum(p6_edu.t%o%p6_edu.t)/n
    p_kw[6,7]<-p_kw[6,7]+sum(p6_edu.t%o%p7_apoe12)/n
    p_kw[6,8]<-p_kw[6,8]+sum(p6_edu.t%o%p7_apoe12.t)/n
    p_kw[6,9]<-p_kw[6,9]+sum(p6_edu.t%o%p9_fa)/n
    p_kw[6,10]<-p_kw[6,10]+sum(p6_edu.t%o%p9_fa.t)/n
    
    p_kw[7,7]<-p_kw[7,7]+sum(p7_apoe12%o%p7_apoe12)/n
    p_kw[7,8]<-p_kw[7,8]+sum(p7_apoe12%o%p7_apoe12.t)/n
    p_kw[7,9]<-p_kw[7,9]+sum(p7_apoe12%o%p9_fa)/n
    p_kw[7,10]<-p_kw[7,10]+sum(p7_apoe12%o%p9_fa.t)/n
    
    p_kw[8,8]<-p_kw[8,8]+sum(p7_apoe12.t%o%p7_apoe12.t)/n
    p_kw[8,9]<-p_kw[8,9]+sum(p7_apoe12.t%o%p9_fa)/n
    p_kw[8,10]<-p_kw[8,10]+sum(p7_apoe12.t%o%p9_fa.t)/n
    
    p_kw[9,9]<-p_kw[9,9]+sum(p9_fa%o%p9_fa)/n
    p_kw[9,10]<-p_kw[9,10]+sum(p9_fa%o%p9_fa.t)/n
    
    p_kw[10,10]<-p_kw[10,10]+sum(p9_fa.t%o%p9_fa.t)/n
    
    p_kw_alpha[i]<-sum(p1)
    p_kw_age[i]<-sum(p3_age)
    p_kw_edu[i]<-sum(p6_edu)
    p_kw_apoe12[i]<-sum(p7_apoe12)
    p_kw_fa[i]<-sum(p9_fa)
  }
  
  p_kw[2,1]<-p_kw[1,2]
  p_kw[3,1]<-p_kw[1,3]
  p_kw[4,1]<-p_kw[1,4]
  p_kw[5,1]<-p_kw[1,5]
  p_kw[6,1]<-p_kw[1,6]
  p_kw[7,1]<-p_kw[1,7]
  p_kw[8,1]<-p_kw[1,8]
  p_kw[9,1]<-p_kw[1,9]
  p_kw[10,1]<-p_kw[1,10]
  
  
  p_kw[3,2]<-p_kw[2,3]
  p_kw[4,2]<-p_kw[2,4]
  p_kw[5,2]<-p_kw[2,5]
  p_kw[6,2]<-p_kw[2,6]
  p_kw[7,2]<-p_kw[2,7]
  p_kw[8,2]<-p_kw[2,8]
  p_kw[9,2]<-p_kw[2,9]
  p_kw[10,2]<-p_kw[2,10]
  
  
  p_kw[4,3]<-p_kw[3,4]
  p_kw[5,3]<-p_kw[3,5]
  p_kw[6,3]<-p_kw[3,6]
  p_kw[7,3]<-p_kw[3,7]
  p_kw[8,3]<-p_kw[3,8]
  p_kw[9,3]<-p_kw[3,9]
  p_kw[10,3]<-p_kw[3,10]
  
  
  p_kw[5,4]<-p_kw[4,5]
  p_kw[6,4]<-p_kw[4,6]
  p_kw[7,4]<-p_kw[4,7]
  p_kw[8,4]<-p_kw[4,8]
  p_kw[9,4]<-p_kw[4,9]
  p_kw[10,4]<-p_kw[4,10]
  
  
  p_kw[6,5]<-p_kw[5,6]
  p_kw[7,5]<-p_kw[5,7]
  p_kw[8,5]<-p_kw[5,8]
  p_kw[9,5]<-p_kw[5,9]
  p_kw[10,5]<-p_kw[5,10]
  
  
  p_kw[7,6]<-p_kw[6,7]
  p_kw[8,6]<-p_kw[6,8]
  p_kw[9,6]<-p_kw[6,9]
  p_kw[10,6]<-p_kw[6,10]
  
  
  p_kw[8,7]<-p_kw[7,8]
  p_kw[9,7]<-p_kw[7,9]
  p_kw[10,7]<-p_kw[7,10]
  
  p_kw[9,8]<-p_kw[8,9]
  p_kw[10,8]<-p_kw[8,10]
  
  
  p_kw[10,9]<-p_kw[9,10]
  
  
  hatvar_kw[,,t_ind]<-ginv(Q_kw)%*%p_kw%*%ginv(Q_kw)
  diag_hatvar_theta[t_ind,]<-diag(ginv(Q_kw)%*%p_kw%*%ginv(Q_kw))
  
  
  L_alpha_kw[t_ind]<-theta_kw[t_ind,1]-1.96*sqrt(hatvar_kw[1,1,t_ind])/sqrt(n)
  R_alpha_kw[t_ind]<-theta_kw[t_ind,1]+1.96*sqrt(hatvar_kw[1,1,t_ind])/sqrt(n)
  
  L_age_kw[t_ind]<-theta_kw[t_ind,3]-1.96*sqrt(hatvar_kw[3,3,t_ind])/sqrt(n)
  R_age_kw[t_ind]<-theta_kw[t_ind,3]+1.96*sqrt(hatvar_kw[3,3,t_ind])/sqrt(n)
  
  L_edu_kw[t_ind]<-theta_kw[t_ind,5]-1.96*sqrt(hatvar_kw[5,5,t_ind])/sqrt(n)
  R_edu_kw[t_ind]<-theta_kw[t_ind,5]+1.96*sqrt(hatvar_kw[5,5,t_ind])/sqrt(n)
  
  L_apoe12_kw[t_ind]<-theta_kw[t_ind,7]-1.96*sqrt(hatvar_kw[7,7,t_ind])/sqrt(n)
  R_apoe12_kw[t_ind]<-theta_kw[t_ind,7]+1.96*sqrt(hatvar_kw[7,7,t_ind])/sqrt(n)
  
  L_fa_kw[t_ind]<-theta_kw[t_ind,9]-1.96*sqrt(hatvar_kw[9,9,t_ind])/sqrt(n)
  R_fa_kw[t_ind]<-theta_kw[t_ind,9]+1.96*sqrt(hatvar_kw[9,9,t_ind])/sqrt(n)
  
  se_kw[t_ind,]<-sqrt(diag_hatvar_theta[t_ind,])/sqrt(n)
  
  
  for (k_simu in 1:K_simu) {
    #tau<-rbinom(n,1,0.5)
    #tau[which(tau==0)]=-1
    tau<-rnorm(n,0,1)
    G_kw_alpha[k_simu,t_ind]<-ginv(Q_kw[1,1])*sum(tau*p_kw_alpha)/n
    G_kw_age[k_simu,t_ind]<-ginv(Q_kw[3,3])*sum(tau*p_kw_age)/n
    G_kw_edu[k_simu,t_ind]<-ginv(Q_kw[5,5])*sum(tau*p_kw_edu)/n
    G_kw_apoe12[k_simu,t_ind]<-ginv(Q_kw[7,7])*sum(tau*p_kw_apoe12)/n
    G_kw_fa[k_simu,t_ind]<-ginv(Q_kw[9,9])*sum(tau*p_kw_fa)/n
  }
}


abs_G_kw_alpha<-rep(0,K_simu)
abs_G_kw_age<-rep(0,K_simu)
abs_G_kw_edu<-rep(0,K_simu)
abs_G_kw_apoe12<-rep(0,K_simu)
abs_G_kw_fa<-rep(0,K_simu)

for (k_simu in 1:K_simu) {
  abs_G_kw_alpha[k_simu]<-max(abs(G_kw_alpha[k_simu,]))
  abs_G_kw_age[k_simu]<-max(abs(G_kw_age[k_simu,]))
  abs_G_kw_edu[k_simu]<-max(abs(G_kw_edu[k_simu,]))
  abs_G_kw_apoe12[k_simu]<-max(abs(G_kw_apoe12[k_simu,]))
  abs_G_kw_fa[k_simu]<-max(abs(G_kw_fa[k_simu,]))
}

quantile(abs_G_kw_alpha,0.95)
quantile(abs_G_kw_age,0.95)
quantile(abs_G_kw_edu,0.95)
quantile(abs_G_kw_apoe12,0.95)
quantile(abs_G_kw_fa,0.95)

L_scb_kw_alpha<-rep(0,length(o))
R_scb_kw_alpha<-rep(0,length(o))
L_scb_kw_age<-rep(0,length(o))
R_scb_kw_age<-rep(0,length(o))
L_scb_kw_apoe12<-rep(0,length(o))
R_scb_kw_apoe12<-rep(0,length(o))
L_scb_kw_edu<-rep(0,length(o))
R_scb_kw_edu<-rep(0,length(o))
L_scb_kw_fa<-rep(0,length(o))
R_scb_kw_fa<-rep(0,length(o))


for (t_ind in 1:length(o)) {
  L_scb_kw_alpha[t_ind]<-theta_kw[t_ind,1]-quantile(abs_G_kw_alpha,0.95)
  R_scb_kw_alpha[t_ind]<-theta_kw[t_ind,1]+quantile(abs_G_kw_alpha,0.95)
  L_scb_kw_age[t_ind]<-theta_kw[t_ind,3]-quantile(abs_G_kw_age,0.95)
  R_scb_kw_age[t_ind]<-theta_kw[t_ind,3]+quantile(abs_G_kw_age,0.95)
  L_scb_kw_edu[t_ind]<-theta_kw[t_ind,5]-quantile(abs_G_kw_edu,0.95)
  R_scb_kw_edu[t_ind]<-theta_kw[t_ind,5]+quantile(abs_G_kw_edu,0.95)
  L_scb_kw_apoe12[t_ind]<-theta_kw[t_ind,7]-quantile(abs_G_kw_apoe12,0.95)
  R_scb_kw_apoe12[t_ind]<-theta_kw[t_ind,7]+quantile(abs_G_kw_apoe12,0.95)
  L_scb_kw_fa[t_ind]<-theta_kw[t_ind,9]-quantile(abs_G_kw_fa,0.95)
  R_scb_kw_fa[t_ind]<-theta_kw[t_ind,9]+quantile(abs_G_kw_fa,0.95)
}



