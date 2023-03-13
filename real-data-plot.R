#model: mmse=alpha+beta1*age+beta2*mci+beta3*ad+beta4*edu+beta5*apoe1+beta6*apoe2+gamma*fa
rm(list = ls (all = TRUE))
library(MASS)

dataset2<-read.csv(choose.files())#FA_hazard_WB
fa<-dataset2[,38]#u=0.35

dataset1<-read.csv(choose.files())#z_covariate
#作为同步的X
z_id<-dataset1[,1]
z_time<-dataset1[,2]
mci<-dataset1[,4]
ad<-dataset1[,5]
edu<-dataset1[,7]
apoe1<-dataset1[,8]
apoe2<-dataset1[,9]


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

x_apoe1<-list()
x_apoe1[[1]]<-rep(apoe1[1],count_yy[1])

x_apoe2<-list()
x_apoe2[[1]]<-rep(apoe2[1],count_yy[1])


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
  x_apoe1[[i]]<-rep(apoe1[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  x_apoe2[[i]]<-rep(apoe2[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))][1],count_yy[i])
  
  y_mmse[[i]]<-y_value[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  
  T.y[[i]]<-y_time[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]
  T.z[[i]]<-z_time[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))]
}


#所有个体，横坐标个体，纵坐标观测时间
id1=rep(0,sum(count_yy))
id1[1:count_yy[1]]=1
for (i in 2:n) {
  id1[(1+sum(count_yy[1:(i-1)])):(sum(count_yy[1:i]))]=i 
}
id2=rep(0,sum(count_zz))
id2[1:count_zz[1]]=1
for (i in 2:n) {
  id2[(1+sum(count_zz[1:(i-1)])):(sum(count_zz[1:i]))]=i 
}

plot(id1,y_time,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(1,256),ylim = c(0,1),xlab = "Patient",ylab = "Visit time",main="",cex.main=0.8,cex.lab=0.75)
axis(1,seq(1,256,30), seq(1,256,30),cex.axis=0.75)
axis(2,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
lines(id2,z_time,type="p",pch=21,bg="green")
legend("topright", pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("Log hazard of FA","MMSE"),
       cex = 0.8,bty = "n",x.intersp=0.4,y.intersp=1,seg.len=0.9)
#legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
#legend=c("MMSE","Log hazard of FA"),
#cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

#挑出两个人画MMSE,FA
par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(1,2))
#par()
x = T.y[[94]]
z = T.z[[94]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(x,yy1,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Visit time",ylab = "",main="Patient 94",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(z,yy2,type="p",pch=21,bg="green")
legend("topright", pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("Log hazard of FA","MMSE"),
       cex = 0.8,bty = "n",x.intersp=0.4,y.intersp=0.9,seg.len=1)


x = T.y[[256]]
z = T.z[[256]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(x,yy1,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Visit time",ylab = "",main="Patient 183",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(z,yy2,type="p",pch=21,bg="green")
legend("topright", pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("Log hazard of FA","MMSE"),
       cex = 0.8,bty = "n",x.intersp=0.4,y.intersp=0.9,seg.len=1)


#挑出三四个个体画age,mmse,fa生长曲线
#id9,103,156,256
count_yy[9];count_zz[9]
count_yy[103];count_zz[103]
count_yy[156];count_zz[156]
count_yy[256];count_zz[256]
z_fa[[9]];y_mmse[[9]];x_age[[9]]
z_fa[[103]];y_mmse[[103]];x_age[[103]]
z_fa[[156]];y_mmse[[156]];x_age[[156]]
z_fa[[256]];y_mmse[[256]];x_age[[256]]
T.y[[9]];T.y[[103]];T.y[[156]];T.y[[256]]
T.z[[9]];T.z[[103]];T.z[[156]];T.z[[256]]

par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(3,3))
#mmse
plot(T.y[[9]],y_mmse[[9]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(0,35),xlab = "Visit time", ylab = "MMSE",main = "Patient 9")
plot(T.y[[156]],y_mmse[[156]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(0,35),xlab = "Visit time", ylab = "MMSE",main = "Patient 156")
plot(T.y[[256]],y_mmse[[256]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(0,35),xlab = "Visit time", ylab = "MMSE",main = "Patient 256")
#age
plot(T.y[[9]],x_age[[9]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(60,85),xlab = "Visit time", ylab = "Age",main = "Patient 9")
plot(T.y[[156]],x_age[[156]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(60,85),xlab = "Visit time", ylab = "Age",main = "Patient 156")
plot(T.y[[256]],x_age[[256]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(60,85),xlab = "Visit time", ylab = "Age",main = "Patient 256")
#fa
plot(T.z[[9]],z_fa[[9]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(0,2),xlab = "Visit time", ylab = "FA",main = "Patient 9")
plot(T.z[[156]],z_fa[[156]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(0,2),xlab = "Visit time", ylab = "FA",main = "Patient 156")
plot(T.z[[256]],z_fa[[256]],type='o',pch=1,lwd=1,col="red",xlim = c(0.3,0.7),ylim = c(0,2),xlab = "Visit time", ylab = "FA",main = "Patient 256")





#group.factor=T.y[[9]]
#group.value=x_age[[9]]
#dataset <- data.frame(value = group.value, group = group.factor)
#library(ggplot2)
#sp=ggplot(dataset,aes(x=group.factor,y=group.value))+geom_point(shape=1)
#sp



#install.packages("ggprism")
library(ggplot2)
library(ggprism)
library(cowplot)
df=read.csv(choose.files())#age9-156-256
p1=ggplot(df)+
  geom_line(aes(tx1,patient9),size=0.8,color="red")+
  scale_x_continuous(breaks = seq(0.3,0.7,0.1))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())+
  labs(title = "Patient 9",
       x=NULL,
       y="Age")
p2=ggplot(df)+
  geom_line(aes(tx2,patient156),size=0.8,color="green")+
  #theme_prism(palette="candy_soft",
   #           base_footaface="plain",
   #           base_size = 16,
   #           base_line_size = 0.8,
   #           axis_text_angle = 45)+
  scale_x_continuous(breaks = seq(0.3,0.7,0.1))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())+
  labs(title = "Patient 156",
       x=NULL,
       y="Age")
p3=ggplot(df)+
  geom_line(aes(tx3,patient256),size=0.8,color="blue")+
  #theme_prism(palette="candy_soft",
  #            base_footaface="plain",
  #            base_size = 16,
  #            base_line_size = 0.8,
   #           axis_text_angle = 45)+
  scale_x_continuous(breaks = seq(0.3,0.7,0.1))+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank())+
  labs(title = "Patient 256",
       x=NULL,
       y="Age")
#install.packages("cowplot")
plot_grid(p1,p2,p3,ncol=1)

plot(value ~ group, dataset,ylim = c(68,84),xlim=c(0.3,0.7),type='b', pch=1,lwd=2, ylab = "Age",xlab="Visit time")
lines(T.y[[156]],x_age[[156]],type='b', pch=2,lwd=2,col="red")
lines(T.y[[256]],x_age[[256]],type='b', pch=3,lwd=2,col="blue")




x = T.y[[24]]
z = T.z[[24]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 24",cex.main=0.8,cex.lab=0.75)
axis(1,seq(1,256,15), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend("topright", pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.4,y.intersp=0.6,seg.len=0.8)



x = T.y[[183]]
z = T.z[[183]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 183",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend("topright", pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.4,y.intersp=0.6,seg.len=0.8)


x = T.y[[94]]
z = T.z[[94]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 94",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend("topright", pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.4,y.intersp=0.6,seg.len=0.8)


x = T.y[[219]]
z = T.z[[219]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 219",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)


x = T.y[[221]]
z = T.z[[221]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 221",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[236]]
z = T.z[[236]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 236",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)












x = T.y[[9]]
z = T.z[[9]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 9")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")

x = T.y[[24]]
z = T.z[[24]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 24")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[169]]
z = T.z[[169]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 169")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[102]]
z = T.z[[102]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 102")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")



x = T.y[[225]]
z = T.z[[225]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 255")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[183]]
z = T.z[[183]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 183")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")












#挑个体，chenli
#par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(2,3))

#0
x = T.y[[169]]
z = T.z[[169]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 169",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,5.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)




#0.2
x = T.y[[173]]
z = T.z[[173]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 173",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)



#0.3
x = T.y[[147]]
z = T.z[[147]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 147",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[152]]
z = T.z[[152]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 152",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[182]]
z = T.z[[182]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 182",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[215]]
z = T.z[[215]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 215",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)




#0.4
x = T.y[[4]]
z = T.z[[4]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 4",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[186]]
z = T.z[[186]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 186",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[246]]
z = T.z[[246]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 246",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[247]]
z = T.z[[247]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 247",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[251]]
z = T.z[[251]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 251",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)



#0.5
x = T.y[[10]]
z = T.z[[10]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 10",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[170]]
z = T.z[[170]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 170",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[176]]
z = T.z[[176]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 176",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[228]]
z = T.z[[228]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 228",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[229]]
z = T.z[[229]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 229",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)




#0.6
x = T.y[[3]]
z = T.z[[3]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 3",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[5]]
z = T.z[[5]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 5",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[25]]
z = T.z[[25]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 25",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[26]]
z = T.z[[26]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 26",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[99]]
z = T.z[[99]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 99",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[187]]
z = T.z[[187]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 187",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[188]]
z = T.z[[188]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 188",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)



#0.7
x = T.y[[236]]
z = T.z[[236]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 236",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)



#0.8
x = T.y[[50]]
z = T.z[[50]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 50",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[68]]
z = T.z[[68]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 68",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)






#挑选6个
par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(2,3))
#par()
x = T.y[[169]]
z = T.z[[169]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 169",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[215]]
z = T.z[[215]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 215",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[247]]
z = T.z[[247]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 247",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[26]]
z = T.z[[26]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 26",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[236]]
z = T.z[[236]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 236",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[68]]
z = T.z[[68]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 68",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)






#挑选4个
par(mar=c(5,5,4,2),oma=c(0,0,0,0),mfrow=c(2,2))
x = T.y[[169]]
z = T.z[[169]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 169",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[215]]
z = T.z[[215]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 215",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[26]]
z = T.z[[26]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 26",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)

x = T.y[[236]]
z = T.z[[236]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 236",cex.main=1,cex.lab=1)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=1)
axis(2,seq(2,4,2), c("FA","MMSE"),cex.axis=1)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","FA"),
       cex = 1,bty = "n",x.intersp=0.8,y.intersp=1.6,seg.len=1)














x = T.y[[170]]
z = T.z[[170]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))
plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 170",cex.main=0.8,cex.lab=0.75)
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.75)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.75)
lines(x,yy1,type="p",pch=21,bg="green")
legend(0.65,6.7, pch=c(21,24),pt.bg=c("green","red"), col=c("black","black"),
       legend=c("MMSE","Log hazard of FA"),
       cex = 0.75,bty = "n",x.intersp=0.2,y.intersp=0.4,seg.len=0.6)

x = T.y[[199]]
z = T.z[[199]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 199")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[217]]
z = T.z[[217]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 217")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[144]]
z = T.z[[144]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 144")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[92]]
z = T.z[[92]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 92")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")


x = T.y[[81]]
z = T.z[[81]]
yy1 = rep(4, length(x))
yy2 = rep(2, length(z))

plot(z,yy2,xaxt = "n", yaxt = "n",type="p",pch=24,bg="red",xlim=c(0,1),ylim = c(0,6),xlab = "Observation time",ylab = "",main="Patient 81")
axis(1,seq(0,1,0.2), seq(0,1,0.2),cex.axis=0.8)
axis(2,seq(2,4,2), c("Log hazard of FA","MMSE"),cex.axis=0.8)
lines(x,yy1,type="p",pch=21,bg="green")




