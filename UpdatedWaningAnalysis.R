rm(list=ls())
setwd("~/Desktop/berkeley/Boots Lab/WaningVaccines")

library(ggplot2)

get_eq <- function(w,v,mu,n,i) {
  m = n-1
  denom = mu+v+(m*w)
  if(i==0){
    return(mu/(mu+v))
  } else if(i==1){
    return((v/(mu+v)*(((m*w)/denom)^m)))
  } else {
    return((((m*w)/denom)^(n-i))*(v/denom))
  }
}

accel_wane <- function(epsilon0,p,n,i){
  return(epsilon0*(((n-i)/(n-1))^p))
}
hill_wane <- function(epsilon0,p,k,n,i){
  return(epsilon0*(((n-i)^p)/(k^p+(n-i)^p))*((k^p+(n-1)^p)/(n-1)^p))
}

get_susceptibility_accel <- function(w,v,mu,epsilon0,p,n){
  total_s <- 0
  for(i in 0:n){
    total_s <- total_s + accel_wane(epsilon0,p,n,i)*get_eq(w,v,mu,n,i)
  }
  return(total_s)
}

get_susceptibility_hill <- function(w,v,mu,epsilon0,k,p,n){
  total_s <- 0
  for(i in 0:n){
    total_s <- total_s + hill_wane(epsilon0,p,k,n,i)*get_eq(w,v,mu,n,i)
  }
  return(total_s)
}
##############
##FIgure 1####
##############
p <- c(0.5,1,2)
shape <- c(50)
i <- seq(1,100,1)
n <- 100
pdf(file="AccelShapes.pdf",width=4,height=4)
plot(NULL, xlim=c(0,1), ylim=c(0,1), cex.lab=1.5,
     ylab="Susceptibility", xlab="Fraction of Immunity Lost",
     frame=F, main="Early/Late Waning",yaxt='n')
axis(2, at = c(0,1),
     labels = c(0,expression(epsilon[0])))
temp <- c()
type = 1
for (par in p){
  #for (k in shape){
  for (j in i) {
    temp[j] <- accel_wane(1,par,n,j)
  }
  lines(1-i/100,temp,lty=type)
  type <- type + 1
  #}
}
legend(0,1,legend=c("p=0.5", "p=1", "p=2"), lty=1:3,cex=0.75)
dev.off()

p <- c(1,3,5)
shape <- c(50)
i <- seq(1,100,1)
n <- 100
pdf(file="HillShapes1.pdf",width=4,height=4)
plot(NULL, xlim=c(0,1), ylim=c(0,1), cex.lab=1.5,
     ylab="Susceptibility", xlab="Fraction of Immunity Lost",
     frame=F, main="Hill Function Waning (fixed k)",yaxt='n')
axis(2, at = c(0,1),
     labels = c(0,expression(epsilon[0])))
temp <- c()
type <- 1
for (par in p){
  #for (k in shape){
  for (j in i) {
    temp[j] <- hill_wane(1,par,shape,n,j)
  }
  lines(1-i/100,temp,lty=type)
  type <- type + 1
  #}
}
legend(0,1,legend=c("p=1", "p=3", "p=5"), lty=1:3,cex=0.75)
dev.off()

p <- c(3)
shape <- c(25,50,75)
i <- seq(1,100,1)
n <- 100
pdf(file="HillShapes2.pdf",width=4,height=4)
plot(NULL, xlim=c(0,1), ylim=c(0,1), cex.lab=1.5,
     ylab="Susceptibility", xlab="Fraction of Immunity Lost",
     frame=F, main="Hill Function Waning (fixed p)",yaxt='n')
axis(2, at = c(0,1),
     labels = c(0,expression(epsilon[0])))
temp <- c()
type <- 1
for (par in p){
  for (k in shape){
    for (j in i) {
      temp[j] <- hill_wane(1,par,k,n,j)
    }
    lines(1-i/100,temp,lty=type)
    type <- type + 1
  }
}
legend(0,1,legend=c("k=25", "k=50", "k=75"), lty=1:3,cex=0.75)
dev.off()


################
####FIgure 2####
################


w <- 1/seq(0.01,3,0.05)
mu <- 0.02
v <- seq(0,3,0.05)
epsilon0 <- c(1)
p <- c(0.5,1,2)
epsilon0 <- c(1,0.7,0.5)
k <- c(1)
i <- seq(1,100,1)
n <- 100

grid_1 <- expand.grid(w,mu,v,epsilon0,p)
colnames(grid_1) <- c("w","mu","v","epsilon0", "p")
grid_1 <- transform(grid_1, susceptibility=get_susceptibility_accel(w,v,mu,epsilon0,p,100))
p1 <- ggplot(grid_1,aes(1/w,v))+
  geom_raster(aes(fill=1/susceptibility))+
  facet_grid(p~epsilon0, labeller = label_bquote(cols = epsilon[0] == .(epsilon0),rows=p==.(p)))+
  scale_fill_binned(breaks=c(2,4,8,16,32),low="white",high="black",name = expression(R[inv]))+
  xlab(expression(paste("Average duration of immunity (",1/omega,")")))+
  ylab(expression(paste("Vaccination rate (",nu,")")))+
  theme_classic()+
  theme(text = element_text(size=20))
ggsave(filename = "Figure2.pdf",width=7,height=6)

##############
###Figure 3###
##############

w <- 1/seq(0.01,3,0.05)
mu <- 0.02
v <- seq(0,3,0.05)
epsilon0 <- c(1,0.7,0.5)
#p <- c(0.5,1,3)
p <- c(1,3,5)
k <- c(50)
#k <- c(25,50,75)
i <- seq(1,100,1)
n <- 100

p.labs <- c("p=1", "p=3", "p=5")
names(p.labs) <- c("1", "3", "5")

grid_h <- expand.grid(w,mu,v,epsilon0,p)
colnames(grid_h) <- c("w","mu","v","epsilon0", "p")
grid_h <- transform(grid_h, susceptibility=get_susceptibility_hill(w,v,mu,epsilon0,k,p,100))
p2 <- ggplot(grid_h,aes(1/w,v))+
  geom_raster(aes(fill=1/susceptibility))+
  facet_grid(p~epsilon0, labeller = label_bquote(cols = epsilon[0] == .(epsilon0),rows=p==.(p)))+
  scale_fill_binned(breaks=c(2,4,8,16,32),low="white",high="black",name = expression(R[inv]))+
  xlab(expression(paste("Average duration of immunity (",1/omega,")")))+
  ylab(expression(paste("Vaccination rate (",nu,")")))+
  theme_classic()+
  theme(text = element_text(size=20))

w <- 1/seq(0.01,3,0.05)
mu <- 0.02
v <- seq(0,3,0.05)
epsilon0 <- c(1,0.7,0.5)
p <- c(1,3,5)
k <- c(25,50,75)
i <- seq(1,100,1)
n <- 100

p.labs <- c("k=25", "k=50", "k=75")
names(p.labs) <- c("25", "50", "75")

grid_h <- expand.grid(w,mu,v,epsilon0,k,p)
colnames(grid_h) <- c("w","mu","v","epsilon0", "k", "p")
grid_h <- transform(grid_h, susceptibility=get_susceptibility_hill(w,v,mu,epsilon0,k,p,100))
p21 <- ggplot(grid_h,aes(1/w,v))+
  geom_raster(aes(fill=1/susceptibility))+
  facet_grid(k~epsilon0, labeller = label_bquote(cols = epsilon[0] == .(epsilon0),rows=k==.(k)))+
  scale_fill_binned(breaks=c(2,4,8,16,32),low="white",high="black",name = expression(R[inv]))+
  xlab(expression(paste("Average duration of immunity (",1/omega,")")))+
  ylab(expression(paste("Vaccination rate (",nu,")")))+
  theme_classic()+
  theme(text = element_text(size=20))

library(cowplot)
legend <- get_legend(p21)

p2 <- p2+theme(legend.position = 'none')
p21 <- p21+theme(legend.position = 'none')

library(cowplot)
plot_grid(p2,p21,legend,nrow=1,rel_widths = c(1,1,0.2),labels=c("A)","B)"))

ggsave(filename = "Figure3.pdf",width=10,height=5)

##############
###FIgure 4###
##############



p <- seq(0.1,2,0.005)
epsilon0 <- seq(0.5,1,0.005)
w <- 1
mu <- 0.02
v <- 1

diff_grid <- expand.grid(w,mu,v,epsilon0,p)
colnames(diff_grid) <- c("w","mu","v","epsilon0", "p")
diff_grid <- transform(diff_grid, susceptibility=get_susceptibility_accel(w,v,mu,epsilon0,p,100))
diff_grid <- transform(diff_grid, ratio=((1/susceptibility) - (1/get_susceptibility_accel(w,v,mu,1,1,100)))/(1/get_susceptibility_accel(w,v,mu,1,1,100)))
p4 <- ggplot(diff_grid,aes(epsilon0,p))+
  geom_raster(aes(fill=ratio))+
  scale_fill_binned(breaks=c(0,.2,.4,.6,.8,1),low="white",high="black",name=expression(paste( R[inv]," fold change")))+
  xlab(expression(paste("Relative susceptibility after waning (",epsilon[0],")")))+ylab("Waning shape parameter (p)")+
  theme_classic()+
  theme(text = element_text(size=20))+
  ggtitle("Early/late function")
#print(p4)


ph <- seq(1,10,0.01)
epsilon0 <- seq(0.5,1,0.005)
w <- 1
mu <- 0.02
v <- 1

diff_grid_hill <- expand.grid(w,mu,v,epsilon0,ph)
colnames(diff_grid_hill) <- c("w","mu","v","epsilon0", "p")
diff_grid_hill <- transform(diff_grid_hill, susceptibility=get_susceptibility_hill(w,v,mu,epsilon0,50,p,100))
diff_grid_hill <- transform(diff_grid_hill, ratio=((1/susceptibility) - (1/get_susceptibility_hill(w,v,mu,1,50,3,100)))/(1/get_susceptibility_hill(w,v,mu,1,50,3,100)))
p5 <- ggplot(diff_grid_hill,aes(epsilon0,p))+
  geom_raster(aes(fill=ratio))+
  scale_fill_binned(breaks=c(0,.2,.4,.6,.8,1),low="white",high="black",name=expression(paste( R[inv]," fold change")))+
  xlab(expression(paste("Relative susceptibility after waning (",epsilon[0],")")))+ylab("Waning shape parameter (p)")+
  theme_classic()+
  theme(text = element_text(size=20))+
  ggtitle("Hill-like function (fixed k=50)")
#print(p5)

kh <- seq(1,100,1)
epsilon0 <- seq(0.5,1,0.005)
w <- 1
mu <- 0.02
v <- 1

diff_grid_hill <- expand.grid(w,mu,v,epsilon0,kh)
colnames(diff_grid_hill) <- c("w","mu","v","epsilon0", "k")
diff_grid_hill <- transform(diff_grid_hill, susceptibility=get_susceptibility_hill(w,v,mu,epsilon0,k,3,100))
diff_grid_hill <- transform(diff_grid_hill, ratio=((1/susceptibility) - (1/get_susceptibility_hill(w,v,mu,1,50,5,100)))/(1/get_susceptibility_hill(w,v,mu,1,50,5,100)))
p6 <- ggplot(diff_grid_hill,aes(epsilon0,k))+
  geom_raster(aes(fill=ratio))+
  scale_fill_binned(breaks=c(0,.2,.4,.6,.8,1),low="white",high="black",name=expression(paste( R[inv]," fold change")))+
  xlab(expression(paste("Relative susceptibility after waning (",epsilon[0],")")))+ylab("Half susceptibility constant (k)")+
  theme_classic()+
  theme(text = element_text(size=20),legend.key.size = unit(2, 'cm'),
        legend.text = element_text(size=25))+
  ggtitle("Hill-like function (fixed p=3)")
#print(p6)

library(cowplot)
legend <- get_legend(p6)

p4 <- p4+theme(legend.position = 'none')
p5 <- p5+theme(legend.position = 'none')
p6 <- p6+theme(legend.position = 'none')

library(cowplot)
plot_grid(p4,legend,p5,p6,nrow=2,rel_widths = c(.8,.8),labels=c("A)"," ", "B)","C)"))

ggsave(filename = "Figure4.pdf",width=12,height=12)


#################

