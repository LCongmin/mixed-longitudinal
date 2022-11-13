#model: mmse=alpha+beta1*age+beta2*ad+beta3*edu+gamma*fa
rm(list = ls (all = TRUE))
library(MASS)

dataset2<-read.csv(choose.files())#FA_hazard_WB
fa<-dataset2[,38]#u=0.3

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


h=0.109
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




#center
#中心化过程,只针对mmse,edu,apoe
h=0.109
#x.tilde_edu<-list()
x.tilde_apoe<-list()
y.tilde<-list()

Ey<-list()
#Eedu=list()
Eapoe=list()

#中心化
for (i in 1:n) {
  Ey[[i]]<-rep(0,count_yy[i])
  #Eedu[[i]]<-rep(0,count_yy[i])
  Eapoe[[i]]<-rep(0,count_yy[i])
  
  x.tilde_apoe[[i]]<-rep(0,count_yy[i])
  #x.tilde_edu[[i]]<-rep(0,count_yy[i])
  y.tilde[[i]]<-rep(0,count_yy[i])
  for (j in 1:count_yy[i]) {
    t=T.y[[i]][j]
    #Eedu[[i]][j]=sum(unlist(x_edu)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    Eapoe[[i]][j]=sum(unlist(x_apoe12)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    Ey[[i]][j]=sum(unlist(y_mmse)*local_kernel(unlist(T.y)-t,h))/sum(local_kernel(unlist(T.y)-t,h))
    
    x.tilde_apoe[[i]][j]<-(x_apoe12[[i]][j]-Eapoe[[i]][j])
    #x.tilde_edu[[i]][j]<-(x_edu[[i]][j]-Eedu[[i]][j])
    y.tilde[[i]][j]<-(y_mmse[[i]][j]-Ey[[i]][j])
    
  }
} 

Q<-matrix(0,ncol=1,nrow=1)
q<-matrix(0, ncol=1,nrow=1)
for (i in 1:n) {
  #Q[1,1]<-Q[1,1]+sum(x.tilde_edu[[i]]*x.tilde_edu[[i]])/n
  #Q[1,2]<-Q[1,2]+sum(x.tilde_edu[[i]]*x_apoe12[[i]])/n
  
  Q[1,1]<-Q[1,1]+sum(x.tilde_apoe[[i]]*x.tilde_apoe[[i]])/n
  
  #q[1,1]<-q[1,1]+sum(y.tilde[[i]]*x.tilde_edu[[i]])/n
  q[1,1]<-q[1,1]+sum(y.tilde[[i]]*x.tilde_apoe[[i]])/n
}
#Q[2,1]<-Q[1,2]

theta_center<-ginv(Q)%*%q
as.numeric(theta_center)
#h0.109 [1] -0.4713756


sig_beta_p<-matrix(0:0,ncol = 1,nrow = 1)
for (i in 1:n) {
  pp=(y.tilde[[i]]-theta_center[1,1]*x.tilde_apoe[[i]])
  #sig_beta_p[1,1]<-sig_beta_p[1,1]+sum((x.tilde_edu[[i]]*pp)%o%(x.tilde_edu[[i]]*pp))/n
  #sig_beta_p[1,2]<-sig_beta_p[1,2]+sum((x.tilde_edu[[i]]*pp)%o%(x.tilde_apoe[[i]]*pp))/n
  
  sig_beta_p[1,1]<-sig_beta_p[1,1]+sum((x.tilde_apoe[[i]]*pp)%o%(x.tilde_apoe[[i]]*pp))/n
}
#sig_beta_p[2,1]<-sig_beta_p[1,2]

hatvar_center<-ginv(Q)%*%sig_beta_p%*%ginv(Q)
diag(hatvar_center)
se_center=sqrt(diag(hatvar_center))/sqrt(n)
se_center
#h0.109 [1] 0.07980242

z_value_center=as.numeric(theta_center)/se_center
z_value_center
2*pnorm(-abs(z_value_center))
#h0.109[1] 3.488519e-09




#center+kw(alpha+fa,alpha用profile,然后估计fa)
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


h=0.109
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



h1_kw<-0.144
h2_kw<-0.144
q1=list()
q2=list()

for (i1 in 1:n) {
  q1[[i1]]=rep(0,count_zz[i1])
  q2[[i1]]=rep(0,count_zz[i1])
  for (j1 in 1:count_zz[i1]) {
    t=T.z[[i1]][j1]
    for (i in 1:n) {
      q1[[i1]][j1]=q1[[i1]][j1]+sum((local_kernel(T.y[[i]]-t,h1_kw))%o%((T.z[[i]]-t)*local_kernel(T.z[[i]]-t,h2_kw)))
      q2[[i1]][j1]=q2[[i1]][j1]+sum((local_kernel(T.y[[i]]-t,h1_kw))%o%((T.z[[i]]-t)^2*local_kernel(T.z[[i]]-t,h2_kw)))
      
    }
  }
}

res<-list()
for (i in 1:n) {
  res[[i]]<-y_mmse[[i]]-theta_center[1,1]*x_apoe12[[i]]
}

w=list()
s11=list()
s44=list()
for (i1 in 1:n) {
  w[[i1]]=rep(0,count_zz[i1])
  s11[[i1]]=rep(0,count_zz[i1])
  s44[[i1]]=rep(0,count_zz[i1])
  for (j1 in 1:count_zz[i1]) {
    t=T.z[[i1]][j1]
    for (i in 1:n) {
      w[[i1]][j1]=w[[i1]][j1]+sum((local_kernel(T.y[[i]]-t,h1_kw))%o%(local_kernel(T.z[[i]]-t,h2_kw)*(q2[[i]]-q1[[i]]*(T.z[[i]]-t))))
      s11[[i1]][j1]=s11[[i1]][j1]+sum((res[[i]]*local_kernel(T.y[[i]]-t,h1_kw))%o%(local_kernel(T.z[[i]]-t,h2_kw)*(q2[[i]]-q1[[i]]*(T.z[[i]]-t))))
      s44[[i1]][j1]=s44[[i1]][j1]+sum((local_kernel(T.y[[i]]-t,h1_kw))%o%(local_kernel(T.z[[i]]-t,h2_kw)*z_fa[[i]]*(q2[[i]]-q1[[i]]*(T.z[[i]]-t))))
    }
  }
}

s1=list()
s4=list()
for (i in 1:n) {
  s1[[i]]=s11[[i]]/w[[i]]
  s4[[i]]=s44[[i]]/w[[i]]
}

#y=list()
x3=list()
for (i in 1:n) {
  #y[[i]]=res[[i]]-s1[[i]]
  x3[[i]]=z_fa[[i]]-s4[[i]]
}

h=0.144
ker_33<-list()
ker_y<-list()
ker_y2<-list()

for (i in 1:n) {
  ker_33[[i]]<-matrix(0, nrow = count_yy[i], ncol = count_zz[i])
  
  ker_y[[i]]<-matrix(0, nrow = count_yy[i], ncol = count_zz[i])
  ker_y2[[i]]<-matrix(0, nrow = count_yy[i], ncol = count_zz[i])
  for (j in 1:count_yy[i]) {
    for (k in 1:count_zz[i]) {
      ker_33[[i]][j,k]<-local_kernel(T.y[[i]][j]-T.z[[i]][k],h)*x3[[i]][k]*x3[[i]][k]
      
      ker_y[[i]][j,k]<-local_kernel(T.y[[i]][j]-T.z[[i]][k],h)*(res[[i]][j])*(z_fa[[i]][k]-s4[[i]][k])
      ker_y2[[i]][j,k]<-local_kernel(T.y[[i]][j]-T.z[[i]][k],h)*s1[[i]][k]*(z_fa[[i]][k]-s4[[i]][k])
    }
  }
}



Q<-0
q<-0
for (i in 1:n) {
  Q<-Q+sum(ker_33[[i]])/n
  q<-q+sum(ker_y[[i]]-ker_y2[[i]])/n
}
theta_kw<-ginv(Q)%*%q
as.numeric(theta_kw)
#h0.144[1] 0.1214683


res_yx1<-list()
for (i in 1:n) {
  res_yx1[[i]]<-matrix(0, nrow = count_yy[i], ncol = count_zz[i])
  
  for (j in 1:count_yy[i]) {
    for (k in 1:count_zz[i]) {
      res_yx1[[i]][j,k]<-local_kernel(T.y[[i]][j]-T.z[[i]][k],h)*x3[[i]][k]*(res[[i]][j]-theta_kw[1,1]*x3[[i]][k]-s1[[i]][k])
    }
  }
}

sig<-matrix(0,ncol = 1,nrow = 1)
for (i in 1:n) {
  sig[1,1]<-sig[1,1]+sum(res_yx1[[i]]%o%res_yx1[[i]])/n
}

hatvar_kw<-ginv(Q)%*%sig%*%ginv(Q)
diag(hatvar_kw)

se_kw<-sqrt(diag(hatvar_kw))/sqrt(n)
se_kw
#h0.144[1] 0.05015524

z_kw<-as.numeric(theta_kw)/se_kw
2*pnorm(-abs(z_kw))
#[1] 0.01544184

