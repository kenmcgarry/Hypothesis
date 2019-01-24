# peter_hoff.R
# https://www.stat.washington.edu/people/pdhoff/book.php
# chapter 1
####
a<-2 ; b<-20
a/(a+b)
(a-1)/(a-1+b-1)
pbeta(.20,a,b) - pbeta(.05,a,b)
pbeta(.10,a,b)


pdf("fig1_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

dbinom(0,20,.05)
n<-20
x<-0:n
del<-.25
plot( range(x-del), c(0,.4),xlab="number infected in the sample",
      ylab="probability",type="n")

points( x-del,dbinom(x,n,.05),type="h",col=gray(.75),lwd=3)
points( x,dbinom(x,n,.10),type="h",col=gray(.5),lwd=3)
points( x+del,dbinom(x,n,.20),type="h",col=gray(0),lwd=3)
legend(10,.35,legend=c(
  expression(paste(theta,"=0.05",sep="")), 
  expression(paste(theta,"=0.10",sep="")),
  expression(paste(theta,"=0.20",sep="")) ),
  lwd=c(3,3,3), 
  col=gray(c(.75,.5,0)) ,bty="n") 

a<-2 ; b<-20
y<-0 ; n<-20


(a+y)/(a+b+n)
(a+y-1)/(a-1+b+n-1)
pbeta(.20,a+y,b+n-y) - pbeta(.05,a+y,b+n-y)
pbeta(.10,a+y,b+n-y)




theta<-seq(0,1,length=500)
plot(theta, dbeta(theta,a+y,b+n-y),
     type="l",
     xlab="percentage infected in the population",
     ylab="", lwd=2, ylim=c(0,16)
)
lines(theta, dbeta(theta,a,b),col="gray",lwd=2)
legend(.5,14,legend=c( expression(paste(italic("p"),"(",theta,")",sep="")), 
                       expression(paste(italic("p"),"(",theta,"|",italic("y"),")",sep=""))  ), 
       bty="n", lwd=c(2,2),col=c("gray","black"))


dev.off()

(a+y)/(b+n-y)
(a+y-1)/(a+y-1+b+n-y-1)
pbeta(.20,a+y,b+n-y) - pbeta(.05,a+y,b+n-y)
pbeta(.10,a+y,b+n-y)


#### sensitivity analysis


pdf("fig1_2.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
g<-50
th0<-seq(0.01,.5,length=g)
nu0<-seq(1,25,length=g) 

PP10<-PM<-PLQ<-PUQ<-matrix(0,g,g)
for(i in 1:g) {for(j in 1:g) {
  a<-nu0[i]*th0[j]
  b<-nu0[i]*(1-th0[j]) 
  
  PM[i,j]<- (a+y)/(a+y+b+n-y) 
  PP10[i,j]<- pbeta(.10,a+y,b+n-y)
  PLQ[i,j]<- qbeta(.05,a+y,b+n-y) 
  PUQ[i,j]<- qbeta(.95,a+y,b+n-y) 
}}

contour(nu0,th0,PM,xlab=expression(italic(w)), ylab=expression(theta[0]))
contour(nu0,th0,PP10,xlab=expression(italic(w)), 
        levels=c(0.1,0.3,.5,.70,.90,.975) )

dev.off()

a<-1 ; b<-1
(a+y)/(b+n-y)
(a+y-1)/(a+y-1+b+n-y-1)
pbeta(.10,a+y,b+n-y) 





#### adj wald interval 
a<-2 ; b<-2 
th<-  (y+a)/(n+a+b)
th+c(-1,1)*1.96*sqrt(th*(1-th)/n)

qbeta(c(.025,.975),a+y,b+n-y)
###########

########## diabetes example

library(lars) ; data(diabetes) ;source("regression_gprior.r")
yf<-diabetes$y
yf<-(yf-mean(yf))/sd(yf)

Xf<-diabetes[[3]]
Xf<-t( (t(Xf)-apply(Xf,2,mean))/apply(Xf,2,sd))

###
n<-length(yf)
set.seed(1)

i.te<-sample(1:n,100)
i.tr<-(1:n)[-i.te]

y<-yf[i.tr] ; y.te<-yf[i.te]
X<-Xf[i.tr,]; X.te<-Xf[i.te,]
#####

#####
p<-dim(X)[2]
S<-10000

if(2==3) {
  BETA<-Z<-matrix(NA,S,p)
  z<-rep(1,dim(X)[2] )
  lpy.c<-lpy.X(y,X[,z==1,drop=FALSE])
  for(s in 1:S)
  {
    for(j in sample(1:p))
    {
      zp<-z ; zp[j]<-1-zp[j]
      lpy.p<-lpy.X(y,X[,zp==1,drop=FALSE])
      r<- (lpy.p - lpy.c)*(-1)^(zp[j]==0)
      z[j]<-rbinom(1,1,1/(1+exp(-r)))
      if(z[j]==zp[j]) {lpy.c<-lpy.p}
    }
    beta<-z;if(sum(z)>0){beta[z==1]<-lm.gprior(y,X[,z==1,drop=FALSE],S=1)$beta }
    Z[s,]<-z
    BETA[s,]<-beta
    if(s%%10==0) {
      dput(list(BETA=BETA,Z=Z), file="diabetes.bma") 
    }
    
  } 
}
#####

####
bfit<-dget("diabetes.bma")
ZPM<-apply(bfit$Z,2,mean,na.rm=TRUE)
beta.bma<-apply(bfit$BETA,2,mean,na.rm=TRUE)

pdf("fig1_3.pdf",family="Times",height=3.5, width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(ZPM,xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(beta[j] != 0),"|",italic(y),",X)",sep="")),type="h",lwd=2)
dev.off()
####

####
y.te.bma<-X.te%*%beta.bma
pdf("fig1_4.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),mfrow=c(1,2))

y.te.bma<-X.te%*%beta.bma
beta.ols<-lm(y~-1+X)$coef
y.te.ols<-X.te%*%beta.ols

plot(y.te,y.te.bma,xlab=expression(italic(y)[test]),
     ylim=range(c(y.te.bma,y.te.ols,y.te)),
     xlim=range(c(y.te.bma,y.te.ols,y.te)), 
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)

plot(y.te,y.te.ols,xlab=expression(italic(y)[test]), 
     ylim=range(c(y.te.bma,y.te.ols,y.te)),
     xlim=range(c(y.te.bma,y.te.ols,y.te)), 
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)
dev.off()
mean( (y.te-y.te.ols)^2 )
mean( (y.te-y.te.bma)^2 )



#### chapter 2
###### 
pdf("fig2_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

x<-0:10
plot(x,dpois(x,2.1), type="h",lwd=1, xlab=expression(italic(y)),ylab=
       #    expression(paste(italic("p(y"),"|",theta==2.1,")")) )
       expression(paste(italic("p"),"(",italic("y"),"|",theta==2.1,")")) )

x<-0:100
plot(x,dpois(x,21), type="h",lwd=1, xlab=expression(italic(y)),ylab=
       expression(paste(italic("p"),"(",italic("y"),"|",theta==21,")")) )
dev.off()
######

######
pdf("fig2_2.pdf",height=3.5,width=7,family="Times")

mu<-10.75
sig<- .8

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

x<- seq(7.9,13.9,length=500)
plot(x,pnorm(x,mu,sig),type="l",ylab=expression(paste(italic("F"),"(",italic("y"),")")),xlab=
       expression(italic(y)),lwd=1)
abline(h=c(0,.5,1),col="gray")
plot(x,dnorm(x,mu,sig),type="l",ylab=expression(paste(italic("p"),"(",italic("y"),")")),
     xlab=
       expression(italic(y)),lwd=1)
abline(v=mu,col="gray")
dev.off()
######

######
pdf("fig2_3.pdf",height=3.5,width=7,family="Times")
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
x<-seq(7.75,13.75,length=100)
mu<-10.75 ; sig<-.8

par(mfrow=c(1,2))
plot(x, dnorm(x,mu,sig),type="l",xlab=expression(italic(y)),
     ylab= expression(paste(italic("p"),"(",italic("y"),")")) )
abline(v=mu,lty=1,col=gray(0))
abline(v=mu,lty=2,col=gray(.33))
abline(v=mu,lty=4,col=gray(.66))

#a<-6 ; b<-1.5
#plot(x, dgamma(x,6,1.5),type="l", xlab=expression(italic(y)),
#    ylab= expression(italic("p(y)")))
#abline(v=(a-1)/b,col="blue")
#abline( v=qgamma(.5,a,b),col="green",lty=2)
#abline(v=a/b,col="red",lty=4)
#legend(4.5,.25,c("mode","median","mean"),lty=c(1,2,4),col=c("blue","green","red"),
#       bty="n",cex=.85)

x<-seq( 0,300000,length=200)
mu<-10.75 ; sig<-.8
plot(x, dlnorm(x,mu,sig)*1e5,type="l", xlab=expression(italic(y)),
     ylab= expression( 10^5*paste(italic("p"),"(",italic("y"),")")) )
abline(v=24600,col=gray(0))
abline( v=qlnorm(.5,mu,sig),col=gray(.3),lty=2)
abline(v=exp(mu+.5*sig^2) , col=gray(.7),lty=4)
legend(150000,1.0,c("mode","median","mean"),
       lty=c(1,2,4),col=gray(c(0,.33,.66)),
       bty="n",cex=.85)

dev.off()

## exchange : 
x<-c(467,667,113,25)
xc<-c(x[1]+x[2],x[3]+x[4])
xc[1]/sum(xc)
## .8915

# chapter 3
pdf("fig3_1.pdf",family="Times",height=4.4,width=7)

#source("../Data/alldata.r")
# get from http://www.stat.washington.edu/~hoff/Book/Data/data/alldata
load("alldata") 

y<-Y[Y$YEAR==1998 & Y$AGE>=65 & Y$FEMALE==1, ]$HAPUNHAP
y[y>4]<-NA
y[y<=2]<-1 
y[y>2]<-0
y<-y[!is.na(y) ]
sy<-sum(y)
n<-length(y)

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(2,1))

theta<-seq(0,1,length=200)
#n<-50 ; y<-44
plot(theta,10^17*theta^sy*(1-theta)^(n-sy),type="l",ylab=
       expression(paste(10^27, paste(italic("p"),"(",italic(y[1]),",...,",italic(y[129]),sep=""), 
                        paste("|",theta,")",sep="")), sep=""),
     xlab=expression(theta))
#mtext(expression(italic( paste("n=50   ",sum(y[i]),"=118",sep=""))),side=3)

plot(theta,dbeta(theta,sy+1,n-sy+1),type="l",ylab=
       expression(paste(italic("p"),"(",theta,"|",italic(y[1]),",...,",italic(y[129]),")",sep="")),
     xlab=expression(theta))
#mtext( expression(italic( paste("n=129   ",sum(y[i]),"=118",sep=""))), side=3)
abline(h=1,col="gray")

ap<- sy+1
bp<- n-sy+1

pmd<-(ap-1)/(ap-1 + bp-1) 
pmn<-ap/(ap+bp)
pvr<- pmn*(1-pmn)/(ap+bp+1)

dev.off()
######

######
pdf("fig3_2.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

n<-10
theta<-.2
plot(0:n,dbinom(0:n,n,theta), type="h",lwd=2,xlab=expression(italic(y)),
     ylab=expression(paste("Pr(",italic("Y=y"),"|",theta==.2,italic(", n="),"10)",sep="")))
#MTEXT(EXpression(
#   italic(paste("n=",10,", ",theta==0.2))),side=3,cex=.8)
n<-10
theta<-.8
plot(0:n,dbinom(0:n,n,theta), type="h",lwd=2,xlab=expression(italic(y)),
     ylab=expression(paste("Pr(",italic("Y=y"),"|",theta==.8,italic(", n="),"10)",sep="")))
#mtext(expression(
#   italic(paste("n=",10,", ",theta==0.8))),side=3,cex=.8)
dev.off()
######

######
pdf("fig3_3.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

n<-100
theta<-.2
plot(0:n,dbinom(0:n,n,theta), type="h",lwd=2,xlab=expression(italic(y)),
     ylab=expression(paste("Pr(",italic("Y=y"),"|",theta==.2,italic(", n="),"100)",sep="")))

n<-100
theta<-.8
plot(0:n,dbinom(0:n,n,theta), type="h",lwd=2,xlab=expression(italic(y)),
     ylab=expression(paste("Pr(",italic("Y=y"),"|",theta==.8,italic(", n="),"100)",sep="")))
dev.off()
######

######
pdf("fig3_4.pdf",family="Times",height=4.4,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0),oma=c(0,0,.5,0))
par(mfrow=c(2,2))
theta<-seq(0,1,length=100)
a<-1; b<-1
n<-5 ; y<-1
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
mtext(expression(paste("beta(1,1) prior,  ", italic("n"),"=5  ",italic(sum(y[i])),"=1",sep="")), side=3,line=.1)
#abline(v=c((a+y-1)/(a+b+n-2),(a+y)/(a+b+n)),col=c("black","gray"),lty=c(2,2))
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)
legend(.45,2.4,legend=c("prior","posterior"),lwd=c(2,2),col=c("gray","black"), bty="n")

a<-3; b<-2
n<-5 ; y<-1
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
#   expression(italic(paste("p(",theta,"|y)",sep=""))), xlab=expression(theta),lwd=2)
mtext(expression(paste("beta(3,2) prior,  ", italic("n"),"=5  ",italic(sum(y[i])),"=1",sep="")), side=3,line=.1)
#abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)

a<-1 ; b<-1
n<-100; y<-20
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
#    expression(italic(paste("p(",theta,"|y)",sep=""))),    xlab=expression(theta),lwd=2)
mtext(expression(paste("beta(1,1) prior,  ", italic("n"),"=100  ",italic(sum(y[i])),"=20",sep="")), side=3,line=.1)
#abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)

a<-3 ; b<-2
n<-100; y<-20
plot(theta,dbeta(theta,a+y,b+n-y),type="l",ylab=
       expression(paste(italic("p("),theta,"|y)",sep="")), xlab=expression(theta), 
     lwd=2)
#    expression(italic(paste("p(",theta,"|y)",sep=""))),xlab=expression(theta),
#    lwd=2)
mtext(expression(paste("beta(3,2) prior,  ", italic("n"),"=100  ",italic(sum(y[i])),"=20",sep="")), side=3,line=.1)
#abline(v=c((a+y-1)/(a+b+n-2), (a+y)/(a+b+n)) , col=c("green","red") )
lines(theta,dbeta(theta,a,b),type="l",col="gray",lwd=2)

dev.off()
######



###### Alternative code, provided by Eric Vanhove

# Code for figure 3.4

a = c(1,3,1,3)
b = c(1,2,1,2)
n = c(5,5,100,100)
y = c(1,1,20,20)

# set up 2x2 panel plot
par(mfrow=c(2,2),
    mar = c(4,4,3,1) + 0.1)

aPlot <- function(a, b, n, y) {
  curve(expr = dbeta(x = x,
                     shape1 = a + y,
                     shape2 = b + n - y),
        col = "black",
        lwd = 1,
        xlab = expression(theta),
        main = bquote(~"beta("~.(a)~","~.(b)~") prior, "~ n==.(n)~"
                      "~Sigma~""~y[i]==.(y)),
        ylab = expression(p(theta ~"|"~ y))
        )
  curve(expr = dbeta(x = x, shape1 = a, shape2 = b), add = TRUE, col =
          "grey", lwd = 1)
}

for (i in 1:4){
  aPlot(a[i],b[i],n[i],y[i])
}

######








######
pdf("fig3_5.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

a<-1  ; b<-1   #prior
n<-10 ; y<-2   #data
theta.support<-seq(0,1,length=100)
plot(theta.support, dbeta(theta.support, a+y, b+n-y), type="l",
     xlab=expression(theta),ylab=expression(paste(italic("p("),theta,"|y)"))) 
abline(v=qbeta( c(.025,.975), a+y,b+n-y))
dev.off()
######

######
pdf("fig3_6.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

theta.support<-seq(0,1,length=5000)
plot(theta.support, dbeta(theta.support, a+y, b+n-y), type="l",
     xlab=expression(theta),ylab=expression(paste(italic("p("),theta,"|y)"))) 
pth<-dbeta(theta.support, a+y, b+n-y)
pth<-pth
ord<- order(-pth)
xpx<-cbind(theta.support[ord], pth[ord])
xpx<-cbind(xpx,cumsum(xpx[,2])/sum(xpx[,2]))

hpd<-function(x,dx,p){
  md<-x[dx==max(dx)]
  px<-dx/sum(dx)
  pxs<--sort(-px)
  ct<-min(pxs[cumsum(pxs)< p])
  list(hpdr=range(x[px>=ct]),mode=md) }

tmp<-hpd(xpx[,1],xpx[,2],.5)$hpdr
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(.75),lwd=2   )
tmp<-hpd(xpx[,1],xpx[,2],.75)$hpdr
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(.5),lwd=2   )
tmp<-hpd(xpx[,1],xpx[,2],.95)$hpdr
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(0),lwd=2   )

tmp<-qbeta( c(.025,.975), a+y,b+n-y)
lines( x=c(tmp[1],tmp[1],tmp[2],tmp[2]),
       y=dbeta(c(0,tmp[1],tmp[2],0),a+y,b+n-y)  ,col=gray(0),lwd=2 ,lty=2  )


legend(.5, 2.75, c("50% HPD","75% HPD","95% HPD","95% quantile-based"), 
       col=c(gray(.75),gray(.5),
             gray(0),gray(0)),lty=c(1,1,1,2),lwd=c(2,2,2,2),
       bty="n")

dev.off()
######

######
pdf("fig3_7.pdf",family="Times",height=3.5,width=7)
#source("/Users/hoff/Work/Datasets/GSS/alldata.r")

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

#CHILDS<-Y$CHILDS[Y$FEMALE==1 &  Y$YEAR==1998 & Y$AGE>=40]
#CHILDS<-Y$CHILDS[Y$FEMALE==1 &  Y$YEAR==1996 & Y$AGE==40 & !is.na(Y$DEG)]
CHILDS<-Y$CHILDS[Y$FEMALE==1&Y$YEAR==1998 & Y$AGE>=40 & Y$AGE<50 & !is.na(Y$DEG)]
CHILDS<-Y$CHILDS[Y$FEMALE==1 & Y$YEAR>=1990 & Y$AGE==40 & !is.na(Y$DEG)]
CHILDS<-CHILDS[!is.na(CHILDS)]
ecdf<-(table(c(CHILDS,0:8))-1 )/sum(table(CHILDS))
plot(0:8+.1,ecdf,type="h",lwd=5,xlab="number of children", 
     ylab=expression(paste("Pr(",italic(Y[i]==y[i]),")",sep="")),col="gray")
points(0:8-.1, dpois(0:8,mean(CHILDS,na.rm=T)),lwd=5,col="black",type="h")

legend(2.25,.31,
       legend=c("Poisson model","empirical distribution"),lwd=c(2,2),col=
         c("black","gray"),bty="n",cex=.75)

#ys<-NULL
#for(ns in 1:10000) { ys<-c(ys,sum(sample(CHILDS,10,replace=T),na.rm=T)) }
#plot(0:52-.1, (table(c(ys,0:52))-1)/10000,type="h",lwd=3 )
plot(0:52,dpois(0:52,10*mean(CHILDS)),lwd=3,col="black",type="h",
     xlab="number of children", 
     ylab=expression(paste("Pr(",italic(sum(Y[i])==y),"|",theta==1.83,")",sep="")))

dev.off()


######
pdf("fig3_8.pdf",family="Times",height=4,width=6)  
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(2,3))

a<-1 ; b<-1
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",1," b=",1,sep=""))),side=3,line=.12,cex=.8)

a<-2 ; b<-2
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",2," b=",2,sep=""))),side=3,line=.12,cex=.8)

a<-4 ; b<-4
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",4," b=",4,sep=""))),side=3,line=.12,cex=.8)

a<-2 ; b<-1
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",2," b=",1,sep=""))),side=3,line=.12,cex=.8)

a<-8 ; b<-4
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",8," b=",4,sep=""))),side=3,line=.12,cex=.8)

a<-32 ; b<-16
x<-seq(.001,10,length=100)
plot(x, dgamma(x,a,b),type="l",
     xlab=expression(theta), ylab=expression(italic(paste("p(",theta,")",sep=""))))
mtext(expression(italic(paste("a=",32," b=",16,sep=""))),side=3,line=.12,cex=.8)

dev.off()
######

###### 
y2<-Y$CHILDS[Y$FEMALE==1 &  Y$YEAR>=1990  & Y$AGE==40 & Y$DEG>=3 ]
y1<-Y$CHILDS[Y$FEMALE==1 &  Y$YEAR>=1990  & Y$AGE==40 & Y$DEG<3 ]



y2<-y2[!is.na(y2)]
y1<-y1[!is.na(y1)]

pdf("fig3_9.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

set.seed(1) 
n1<-length(y1) ; n2<-length(y2)
s1<-sum(y1)
s2<-sum(y2)

par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(table(y1), type="h",xlab=expression(italic(y)),ylab=expression(italic(n[1](y))),col=gray(.5) ,lwd=3)
mtext("Less than bachelor's",side=3)
plot(table(y2), type="h",xlab=expression(italic(y)),ylab=expression(italic(n[2](y))),col=gray(0),lwd=3)
mtext("Bachelor's or higher",side=3,lwd=3)
dev.off()
######

######
pdf("fig3_10.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
a<-2
b<-1
xtheta<-seq(0,5,length=1000)
plot(xtheta,dgamma(xtheta,a+s1,b+n1),type="l",col=gray(.5),xlab=expression(theta),
     ylab=expression(paste(italic("p("),theta,"|",y[1],"...",y[n],")",sep="")))
lines(xtheta,dgamma(xtheta,a+s2,b+n2),col=gray(0),lwd=2)
lines(xtheta,dgamma(xtheta,a,b),type="l",lty=2,lwd=2)
abline(h=0,col="black")
#dev.off()
#pdf("fig4_6.pdf",family="Times",height=3,width=6)
#par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
#par(mfrow=c(1,2))
y<-(0:12)
plot(y-.1, dnbinom(y, size=(a+s1), mu=(a+s1)/(b+n1)) , col=gray(.5) ,type="h",
     ylab=expression(paste(italic("p("),y[n+1],"|",y[1],"...",y[n],")",sep="")), 
     xlab=expression(italic(y[n+1])),ylim=c(0,.35),lwd=3)
points(y+.1, dnbinom(y, size=(a+s2), mu=(a+s2)/(b+n2)) , col=gray(0) ,type="h",lwd=3)
legend(1,.375,legend=c("Less than bachelor's","Bachelor's or higher"),bty="n",
       lwd=c(3,3),col=c(gray(.5),gray(0)))
dev.off()
######





a<-2 ; b<-1
n1<-length(y1) ; s1<-sum(y1)
n2<-length(y2) ; s2<-sum(y2)


a<-2 ; b<-1          # prior parameters
n1<-111 ; s1<-217    # data in group 1
n2<-44  ; s2<-66     # data in group 2


(a+s1)/(b+n1)        # posterior mean 
(a+s1-1)/(b+n1)      # posterior mode
qgamma( c(.025,.975),a+s1,b+n1)   # posterior 95% CI

(a+s2)/(b+n2)
(a+s2-1)/(b+n2)
qgamma( c(.025,.975),a+s2,b+n2)


th1_mc<-rgamma(100000,a+s1,b+n1)

th2_mc<-rgamma(100000,a+s2,b+n2)

mean(th1_mc>th2_mc)

y1_mc<-rpois(1000000,th1_mc)
y2_mc<-rpois(1000000,th2_mc)
mean(y1_mc>y2_mc)
mean(y1_mc>=y2_mc)
mean(y1_mc==y2_mc)


options(width=60)

y<- 0:10

dnbinom(y, size=(a+s1), mu=(a+s1)/(b+n1))

dnbinom(y, size=(a+s2), mu=(a+s2)/(b+n2))



# chapter 4
#####
pdf("fig4_1.pdf",family="Times",height=3,width=6)
par(mar=c(3,3,.25,1),mgp=c(1.75,.75,0))
par(mfrow=c(2,3))
set.seed(1)
a<-68 ; b<-45
set.seed(1)
theta.support<-seq(0,3,length=100)
theta.sim10<-rgamma(10,a,b)
theta.sim100<-rgamma(100,a,b)
theta.sim1000<-rgamma(1000,a,b)

xlim<-c(.75,2.25)
ylim=c(0,2.5)
lty=1

hist( theta.sim10, prob=T,xlim=xlim,ylim=ylim,xlab="",main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
text(2.1,2.25,expression(paste(italic(S),"=10",sep="")))

hist( theta.sim100, prob=T,xlim=xlim,ylim=ylim,xlab="",main="" ,ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
text(2.1,2.25,expression(paste(italic(S),"=100",sep="")))



hist( theta.sim1000, prob=T,xlim=xlim,ylim=ylim,xlab="",main="" ,ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
text(2.1,2.25,expression(paste(italic(S),"=1000",sep="")))


plot(density(theta.sim10),xlim=xlim,ylim=ylim,xlab=expression(theta),main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)

plot(density(theta.sim100),xlim=xlim,ylim=ylim,xlab=expression(theta),main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)

plot(density(theta.sim1000),xlim=xlim,ylim=ylim,xlab=expression(theta),main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
dev.off()
#####

set.seed(1)
a<-2  ; b<-1
sy<-66; n<-44

theta.sim10<-rgamma(10,a+sy,b+n)
theta.sim100<-rgamma(100,a+sy,b+n)
theta.sim1000<-rgamma(1000,a+sy,b+n)

(a+sy)/(b+n) 

mean(theta.sim10)
mean(theta.sim100)
mean(theta.sim1000)

pgamma(1,75,a+sy,b+n)

mean( theta.sim10<1.75)
mean( theta.sim100<1.75)
mean( theta.sim1000<1.75)

qgamma(c(.025,.975),a+sy,b+n)
quantile( theta.sim10, c(.025,.975))
quantile( theta.sim100, c(.025,.975))
quantile( theta.sim1000, c(.025,.975))



######
pdf("fig4_2.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))

set.seed(1)
a<-2   ; b<-1
sy<-66 ; n<-44

nsim<-1000
theta.sim<-rgamma(nsim,a+sy,b+n)

#cumulative mean

cmean<-cumsum(theta.sim)/(1:nsim)
cvar<- cumsum(theta.sim^2)/(1:nsim) - cmean^2
ccdf<- cumsum(theta.sim<1.75)/ (1:nsim)
cq<-NULL
for(j in 1:nsim){ cq<-c(cq,quantile(theta.sim[1:j],probs=0.975)) }

sseq<- c(1,(1:100)*(nsim/100))
cmean<-cmean[sseq] 
cq<-cq[sseq] 
ccdf<-ccdf[sseq] 

plot(sseq,cmean,type="l",xlab="# of Monte Carlo samples",ylab="cumulative mean",
     col="black")
abline(h= (a+sy)/(b+n),col="gray",lwd=2)

plot(sseq,ccdf,type="l",xlab="# of Monte Carlo samples",ylab="cumulative cdf at 1.75",col="black")
abline(h= pgamma(1.75,a+sy,b+n),col="gray",lwd=2)

plot(sseq,cq,type="l",xlab="# of Monte Carlo samples",ylab="cumulative 97.5% quantile",col="black")
abline(h= qgamma(.975,a+sy,b+n),col="gray",lwd=2)

dev.off()
######




#####
#source("../Data/alldata.r")
# get from http://www.stat.washington.edu/~hoff/Book/Data/data/alldata
load("alldata")

table(Y$DEG[Y$YEAR==1998])
y1<-Y$PRAYER[Y$YEAR==1998 & Y$RELIG==1 ]
y1<-1*(y1==1)
y1<-y1[!is.na(y1) ]
sy1<-sum(y1)
n1<-length(y1)
sy1/n1

y2<-Y$PRAYER[Y$YEAR==1998 & Y$RELIG!=1 ]
y2<-1*(y2==1)
y2<-y2[!is.na(y2) ]
sy2<-sum(y2)
n2<-length(y2)
sy2/n2


table(Y$FEMALE[Y$YEAR==1998])
y<-Y$FEMALE[Y$YEAR==1998]
y<-1*(y==1)
y<-y[!is.na(y) ]
sy<-sum(y)
n<-length(y)
sy/n


sy<-sy2
n<-n2


###



set.seed(1)

a<-1 ; b<-1 
theta.prior.sim<-rbeta(10000,a,b)
gamma.prior.sim<- log( theta.prior.sim/(1-theta.prior.sim) )

n0<-860-441 ; n1<-441
theta.post.sim<-rbeta(10000,a+n1,b+n0)
gamma.post.sim<- log( theta.post.sim/(1-theta.post.sim) )

pdf("fig4_3.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(2,3))
par(cex=.8)

par(mfrow=c(1,2),mar=c(3,3,1,1), mgp=c(1.75,.75,.0))
plot(density(gamma.prior.sim,adj=2),xlim=c(-5,5),main="", xlab=expression(gamma),
     ylab=expression(italic(p(gamma))),col="gray")
plot(density(gamma.post.sim,adj=2),xlim=c(-5,5),main="",xlab=expression(gamma),
     #  ylab=expression(italic(paste("p(",gamma,"|",y[1],"...",y[n],")")))) 
     ylab=expression(paste(italic("p("),gamma,"|",y[1],"...",y[n],")",
                           sep="")) )
lines(density(gamma.prior.sim,adj=2),col="gray")

dev.off()
#####


#####
set.seed(1)
a<-2 ; b<-1
sy1<-217 ;  n1<-111
sy2<-66  ;  n2<-44

theta1.mc<-rgamma(10000,a+sy1, b+n1)
theta2.mc<-rgamma(10000,a+sy2, b+n2)

y1.mc<-rpois(10000,theta1.mc)
y2.mc<-rpois(10000,theta2.mc)


mean(theta1.mc>theta2.mc)

mean(y1.mc>y2.mc)


pdf("fig4_4.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))
plot(density(theta1.mc/theta2.mc,adj=2),main="",xlim=c(.75,2.25),
     xlab=expression(gamma==theta[1]/theta[2]),
     #ylab=expression(italic(paste("p(",gamma,"|",bold(y[1]),",",bold(y[2]),")",
     #   sep="")) ))
     ylab=expression(paste(italic("p("),gamma,"|",bold(y[1]),",",bold(y[2]),")",
                           sep="")) )

dev.off()

pdf("fig4_5.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))

diff.mc<-y1.mc-y2.mc

ds<- -11:11
plot(ds,(table(c(diff.mc,ds))-1)/length(diff), type="h",lwd=3,
     xlab=expression(italic(D==tilde(Y)[1]-tilde(Y)[2])),
     ylab=expression(paste(italic("p(D"),"|",bold(y[1]),",",bold(y[2]),")",sep="")))
dev.off()
#####


##### 


######
y1<-Y$CHILDS[Y$FEMALE==1 &  Y$YEAR>=1990  & Y$AGE==40 & Y$DEG<3 ]
y1<-y1[!is.na(y1)]

##
set.seed(1)

a<-2 ; b<-1
t.mc<-NULL
for(s in 1:10000) {
  theta1<-rgamma(1,a+sum(y1), b+length(y1))
  y1.mc<-rpois(length(y1),theta1)
  t.mc<-c(t.mc,sum(y1.mc==2)/sum(y1.mc==1))
}

t.obs<-sum(y1==2)/sum(y1==1)
mean(t.mc>=t.obs)


#a<-2 ; b<-1
#t.mc<-NULL
#for(s in 1:10000) {
#theta1<-rgamma(1,a+sum(y1), b+length(y1))
#y1.mc<-rpois(length(y1),theta1)
#t.mc<-rbind(t.mc,c(mean(y1.mc),var(y1.mc) ))
#                    }





##
pdf("fig4_6.pdf",family="Times",height=3.5,width=7)

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

ecdf<-(table(c(y1,0:9))-1 )/sum(table(y1))
#ecdf.mc<-(table(c(y1.mc,0:9))-1 )/sum(table(y1.mc))
ecdf.mc<- dnbinom(0:9,size=a+sum(y1),mu=(a+sum(y1))/(b+length(y1)))
plot(0:9+.1,ecdf.mc,type="h",lwd=5,xlab="number of children",
     ylab=expression(paste("Pr(",italic(Y[i]==y[i]),")",sep="")),col="gray",
     ylim=c(0,.35))
points(0:9-.1, ecdf,lwd=5,col="black",type="h")

legend(1.8,.35,
       legend=c("empirical distribution","predictive distribution"),
       lwd=c(2,2),col=
         c("black","gray"),bty="n",cex=.8)


hist(t.mc,prob=T,main="",ylab="",xlab=expression(t(tilde(Y))) )

segments(t.obs,0,t.obs,.25,col="black",lwd=3)
dev.off()


# chapter 5
#####
pdf("fig5_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

par(mfrow=c(1,1))
theta<-2 ; sigma<- .5
x<-seq(0,10,length=100)
plot(x,dnorm(x,theta,sigma),type="l",
     xlab=expression(italic(y)),
     ylab=expression(paste(italic("p(y"),"|",theta,",",sigma^2,")",sep="")),col=gray(.15),lwd=2)

theta<-5 ; sigma<-2
lines(x, dnorm(x,theta,sigma),col=gray(.5),lwd=2)

theta<-7 ; sigma<-1
lines(x, dnorm(x,theta,sigma),col=gray(.85),lwd=2)
legend( 5,.7, 
        c(expression(paste(theta==2,",",sigma^2==.25,sep="")),
          expression(paste(theta==5,",",sigma^2==4,sep="")),
          expression(paste(theta==7,",",sigma^2==1,sep="")) 
        ), 
        col=gray(c(.15,.5,.85)), lwd=c(2,2,2),lty=c(1,1,1),bty="n")

dev.off()
#####

#####
pdf("fig5_2.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,.25,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))

library(alr3)
data(heights)
y<-heights[,2]
#y<-round(heights[,2])
hist(y,main="",ylab="",
     #     xlab=expression(italic(y)),
     xlab="height in inches",
     prob=TRUE,nclass=15,
     col="gray")

y.av<-mean(y) ; y.sd<-sd(y)
ys<-seq(min(y)*.9,max(y)*1.1,length=100)
lines(ys,dnorm(ys,y.av,y.sd),lwd=2)
dev.off()
#####


#####
pdf("fig5_3.pdf",family="Times",height=3.5,width=7)

library(Flury)
data(midge)  
y<-midge[midge[,1]=="Af",3]

mu0<-1.9
tau0<-(.5*1.9)
t02<-tau0^2

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

ybar<- mean(y)
s2<- var(y)
n<- length(y)
sigma<-sqrt(s2)

mun<-( mu0/(t02) + n*ybar/s2)/( 1/t02+n/s2)
t2n<-1/(1/t02 +n/s2)

qnorm(c(.025,.975),mun,sqrt(t2n))

ys<-seq(0,mu0*2,length=500)
plot(ys,dnorm(ys,mun,sqrt(t2n)),type="l",col="black",xlab=expression(theta),
     ylab=expression(paste(italic("p("),theta,"|",italic(y[1]),"...",italic(y[n]),",",
                           sigma^2==0.017,")",sep="")),lwd=2)  
lines(ys,dnorm(ys,mu0,sqrt(t02)),type="l",col="gray",lwd=2)
dev.off()
#####


#####
# prior
mu0<-1.9  ; k0<-1
s20<-.01 ; nu0<-1

# data
y<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n<-length(y) ; ybar<-mean(y) ; s2<-var(y)

# posterior inference
kn<-k0+n ; nun<-nu0+n
mun<- (k0*mu0 + n*ybar)/kn  
s2n<- (nu0*s20 +(n-1)*s2 +k0*n*(ybar-mu0)^2/(kn))/(nun)
mun
s2n

dinvgamma<-function(x,a,b) {
  ld<- a*log(b) -lgamma(a) -(a+1)*log(x)  -b/x 
  exp(ld)
}

gs<-100
theta<-seq(1.6,2.0,length=gs)
is2<-seq(15,160 ,length=gs)
s2g<-seq(.001,.045,length=gs)

ld.th.is2<-ld.th.s2<-matrix(0,gs,gs)
for(i in 1:gs) { for(j in 1:gs) {
  ld.th.is2[i,j]<- dnorm(theta[i],mun,1/sqrt( is2[j] *10) ,log=TRUE) +
    dgamma(is2[j],10/2,10*s2n/2, log=TRUE )
  ld.th.s2[i,j]<- dnorm(theta[i],mun,sqrt(s2g[j]/10) ,log=TRUE) +
    log( dinvgamma(s2g[j],10/2,10*s2n/2 ))
}} 
pdf("fig5_4.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
grays<- gray( (10:0)/10)
image(theta,is2,exp(ld.th.is2),col=grays,xlab=expression(theta),
      ylab=expression(tilde(sigma)^2) ) 
image(theta,s2g,exp(ld.th.s2), col=grays,xlab=expression(theta),
      ylab=expression(sigma^2) )
dev.off()


# monte carlo sampling
pdf("fig5_5.pdf",family="Times",height=7,width=7)
set.seed(1)
S<-10000
s2.postsample<-1/rgamma(S,  (nu0+n)/2, s2n*(nu0+n)/2 )
theta.postsample<-rnorm(S, mun, sqrt(s2.postsample/(k0+n)))
quantile(theta.postsample, c(.025,.975))
#plot(density(theta.postsample,adjust=3) )
#x<-seq(1.6,2.0,length=200) ; lines(x,dnorm(x,mun,sqrt(t2n)),col="red")
###
layout(matrix(c(1,1,2,3),2,2,byrow=T))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
image(theta,s2g,exp(ld.th.s2), col=grays,xlab=expression(theta),
      ylab=expression(sigma^2),xlim=c(1.60,2.0),ylim=c(.001,.07) )
points(theta.postsample[1:5000], s2.postsample[1:5000],pch=".",
       xlab=expression(theta),ylab=expression(sigma^2),xlim=c(1.65,1.95),ylim=c(.005,
                                                                                .07) )
plot(density(s2.postsample,adjust=3),main="",xlab=expression(sigma^2), 
     xlim=c(0,.075),ylab=expression( paste(italic("p("),
                                           sigma^2,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 
# abline(v=s2n)
plot(density(theta.postsample,adjust=3),main="",xlab=expression(theta),xlim=c(1.60,2.0),ylab=expression( paste(italic("p("),
                                                                                                               theta,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 
#abline(v=mun)
quantile( theta.postsample,c(.025,.975))
abline(v=quantile( theta.postsample,c(.025,.975)),col="gray",lwd=2)
### t-test based confidence interval
n<-length(y) ; ybar<-mean(y) ; s2<-var(y)

ybar+qt( c(.025,.975), n-1) *sqrt(s2/n)
abline( v= ybar+qt( c(.025,.975), n-1) *sqrt(s2/n), col="black",lwd=2)
#legend(1.9,6,legend=c("t interval", "posterior interval"),lwd=c(2,2),
#   col=c("black","gray"),bty="y")

dev.off()
#####

#####
pdf("fig5_6.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

b<-(100-112)^2
s2<-13^2
n<-1:50 

k<-1 ; brk1<- (n/(k+n))^2 + n*(k/(k+n))^2*b/s2
k<-2 ; brk2<- (n/(k+n))^2 + n*(k/(k+n))^2*b/s2
k<-3 ; brk3<- (n/(k+n))^2 + n*(k/(k+n))^2*b/s2

plot(range(n),c(0.4,1.1),type="n",xlab="sample size", ylab="relative MSE")
abline(h=1,lty=2,lwd=2)
lines(n, brk1,col=gray(.25),lwd=2)
lines(n, brk2,col=gray(.5),lwd=2)
lines(n, brk3,col=gray(.75),lwd=2)
legend(20,.8,legend=c(expression(kappa[0]==0),expression(kappa[0]==1), expression(kappa[0]==2), 
                      expression(kappa[0]==3) ),lwd=c(2,2,2),lty=c(2,1,1,1),col=c(gray(c(0,.25,.5,.75))),bty="n")

####
theta0<-112
mu0<-100
n<-10 
s2m<-s2/n
x<-seq(theta0-4*sqrt(s2m),theta0+4*sqrt(s2m), length=100)
plot(x,dnorm(x,theta0,sqrt(s2m)),type="l",lwd=2,ylim=c(0,.13),lty=2, xlab="IQ",
     ylab="")
abline(v=theta0)
for(k in 1:3) {
  w<- n/(n+k) 
  lines(x,dnorm(x,w*theta0+(1-w)*mu0,sqrt(w^2*s2m)),type="l",col=gray(k/4),lwd=2) 
} 
dev.off()
####################


#source("../Data/alldata.r")
# get from http://www.stat.washington.edu/~hoff/Book/Data/data/alldata
load("alldata")

pdf("fig5_7.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))

###
#hist(Y$INCOME98[Y$YEAR==1998] )
#plot(table( Y$INCOME98[Y$YEAR==1998]  ))
#hist( Y$AGE[Y$YEAR==1998]  )
CHILDS<-Y$CHILDS[Y$FEMALE==1&Y$YEAR==1998 & Y$AGE>=40  ]
CHILDS<-CHILDS[!is.na(CHILDS) ]
###

set.seed(1)
NY<-NULL
N<-c(5,15,45)
for(n in N) {
  for(sim in 1:5000) {
    y<-sample(CHILDS,n)
    NY<-rbind(NY, c(n,mean(y),var(y)) )
  } }

plot(table(CHILDS)/sum(table(CHILDS)),type="h",
     xlab=expression(paste(italic(y),"=number of children",sep="" )), 
     ylab=expression(italic(p(y)) ))

x<-seq(0,6,length=200)
plot( range(NY[,2]),c(0,1.7),type="n",xlab="number of children", 
      ylab=expression( italic( p(bar(y))) ) )
for( n in N) {
  yb<- NY[NY[,1]==n,2]
  lines(density( yb,adj=2) ,col=gray(1-( n/5)/9 ), lwd=2)
  cat(n,mean(yb),sd(yb),"\n")
}
abline(v=mean(CHILDS))
legend( 2.35,1.8,legend=c("n=5","n=15","n=45"),lwd=c(2,2,2),col=
          gray(1-( N/5)/9 ), bty="n")


# plot( NY[NY[,1]==45,2:3],xlab=expression(italic(bar(y))),
#    ylab=expression(italic(s^2)) )

source("hdr_2d.r")
plot.hdr2d(NY[NY[,1]==45,2:3],xlab=expression(italic(bar(y))), 
           ylab=expression(italic(s^2)) ,ylim=c(1,6) )


dev.off()



####
mean(CHILDS)
pq<-.5
quantile(CHILDS,pq)
mean( ( CHILDS-mean(CHILDS))^3) /( mean( ( CHILDS-mean(CHILDS))^2)^(3/2)) 

###
k0<-1; mu0<-2 ; nu0<-1 ; s20<-1

n<-45
set.seed(1)
y<-sample(CHILDS,n)
ybar<-mean(y) ; s2<-var(y)

kn<-k0+n ; nun<-nu0+n
mun<- (k0*mu0 + n*ybar)/kn    
s2n<- (nu0*s20 +(n-1)*s2 +k0*n*(ybar-mu0)^2/(kn))/(nun)
###

S<-10000
s2.postsample <- 1/rgamma(S, (nu0+n)/2, s2n*(nu0+n)/2 )
theta.postsample <- rnorm(S, mun, sqrt(s2.postsample/(k0+n)))


#####



hist(theta.postsample) ; abline(v=mean(CHILDS))

hist( qnorm( pq, theta.postsample, sqrt(s2.postsample))) 
abline(v=quantile(CHILDS,pq))

hist( pnorm( 2, theta.postsample, sqrt(s2.postsample)))   
abline(v=quantile(CHILDS,pq))

1-pnorm(7, mean(CHILDS), sd(CHILDS))
1-pnorm(6.5, mean(CHILDS), sd(CHILDS))
mean( CHILDS>= 7 )



# chapter 6
# priors
mu0<-1.9  ; t20<-0.95^2
s20<-.01 ; nu0<-1

#data
y<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n<-length(y) ; mean.y<-mean(y) ; var.y<-var(y)


####
G<-100 ; H<-100

mean.grid<-seq(1.505,2.00,length=G) 
prec.grid<-seq(1.75,175,length=H) 

post.grid<-matrix(nrow=G,ncol=H)

for(g in 1:G) {
  for(h in 1:H) { 
    
    post.grid[g,h]<- dnorm(mean.grid[g], mu0, sqrt(t20)) *
      dgamma(prec.grid[h], nu0/2, s20*nu0/2 ) *
      prod( dnorm(y,mean.grid[g],1/sqrt(prec.grid[h])) )
  }
}

post.grid<-post.grid/sum(post.grid)

pdf("fig6_1.pdf",height=1.75,width=5,family="Times")
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
       xlab=expression(theta), ylab=expression(tilde(sigma)^2) )

mean.post<- apply(post.grid,1,sum)
plot(mean.grid,mean.post,type="l",xlab=expression(theta),
     ylab=expression( paste(italic("p("),
                            theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))

prec.post<-apply(post.grid,2,sum)
plot(prec.grid,prec.post,type="l",xlab=expression(tilde(sigma)^2),
     ylab=expression( paste(italic("p("),
                            tilde(sigma)^2,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 

dev.off()
#####







###
set.seed(1)
S<-1000
PHI<-matrix(nrow=S,ncol=2)
PHI[1,]<-phi<-c( mean.y, 1/var.y)

### Gibbs sampling
for(s in 2:S) {
  
  # generate a new theta value from its full conditional
  mun<-  ( mu0/t20 + n*mean.y*phi[2] ) / ( 1/t20 + n*phi[2] )
  t2n<- 1/( 1/t20 + n*phi[2] )
  phi[1]<-rnorm(1, mun, sqrt(t2n) )
  
  # generate a new sigma^2 value from its full conditional
  nun<- nu0+n
  s2n<- (nu0*s20 + (n-1)*var.y + n*(mean.y-phi[1])^2 ) /nun
  phi[2]<- rgamma(1, nun/2, nun*s2n/2)
  
  PHI[s,]<-phi         }
###



pdf("fig6_2.pdf",height=1.75,width=5,family="Times")
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
m1<-5
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-15
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-100
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )
dev.off()
#####


#####
pdf("fig6_3.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
sseq<-1:1000


image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
       xlab=expression(theta), ylab=expression(tilde(sigma)^2) ,
       xlim=range(PHI[,1]),ylim=range(PHI[,2]) )
points(PHI[sseq,1],PHI[sseq,2],pch=".",cex=1.25 )

plot(density(PHI[,1],adj=2),  xlab=expression(theta),main="",
     xlim=c(1.55,2.05),
     ylab=expression( paste(italic("p("),
                            theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))
abline(v=quantile(PHI[,1],prob=c(.025,.975)),lwd=2,col="gray")
### t-test based confidence interval
n<-length(y) ; ybar<-mean(y) ; s2<-var(y)
ybar+qt( c(.025,.975), n-1) *sqrt(s2/n)
abline( v= ybar+qt( c(.025,.975), n-1) *sqrt(s2/n), col="black",lwd=1)


plot(density(PHI[,2],adj=2), xlab=expression(tilde(sigma)^2),main="",
     ylab=expression( paste(italic("p("),
                            tilde(sigma)^2,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 



dev.off()
#####

quantile(PHI[,1],c(.025,.5,.975))
quantile(PHI[,2],c(.025,.5, .975))
quantile(1/sqrt(PHI[,2]),c(.025,.5, .975))



##### 
pdf("fig6_5.pdf",family="Times",height=4,width=4)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
mu0<-10 ; tau0<-5
s20<-2  ; nu0<-1  

K<-50
L<-50
theta.grid<-seq( 12,18, length=K)
prec.grid<-seq(.05,.25, length=L)
post.grid<-matrix(nrow=K,ncol=L)

plot( 0,0,xlim=range(theta.grid), ylim=range(prec.grid),
      xlab=expression(theta), ylab=expression(1/sigma^2) )
abline( v=theta.grid, col="gray")
abline( h=prec.grid, col="gray")
dev.off()
#####


#####
pdf("fig6_6.pdf",family="Times",height=4,width=8)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

set.seed(1)
y<-rnorm(25,15,3)
for(k in 1:K) {
  for(l in 1:L) {
    post.grid[k,l]<- dnorm(theta.grid[k],mu0,tau0) *
      dgamma(prec.grid[l], nu0/2, s20*nu0/2 ) *
      prod( dnorm(y,theta.grid[k],1/sqrt(prec.grid[l])) )
  }}
post.grid<-post.grid/sum(post.grid)

par(mfrow=c(1,2))
image( theta.grid,prec.grid,post.grid,col=terrain.colors(15),
       xlab=expression(theta), ylab=expression(1/sigma^2))
abline(v=mean(y),lty=2)
abline(h=1/var(y), lty=2)

#abline(v=mu0,lty=2,col="blue")
#abline(h=1/tau0, lty=2,col="blue")

contour( theta.grid,prec.grid,post.grid,
         xlab=expression(theta), ylab=expression(1/sigma^2))
abline(v=mean(y),lty=2)
abline(h=1/var(y), lty=2)
dev.off()
#####




###### Intro to MCMC diagnostics
mu<-c(-3,0,3)
s2<-c(.33,.33,.33)
w<-c(.45,.1,.45)

ths<-seq(-5,5,length=100)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" )



#### MC Sampling
set.seed(1)
S<-2000
d<-sample(1:3,S, prob=w,replace=TRUE)
th<-rnorm(S,mu[d],sqrt(s2[d]))
THD.MC<-cbind(th,d)
####

#### MCMC sampling
th<-0
THD.MCMC<-NULL
S<-10000
set.seed(1)
for(s in 1:S) {
  d<-sample(1:3 ,1,prob= w*dnorm(th,mu,sqrt(s2))   )
  th<-rnorm(1,mu[d],sqrt(s2[d]) )
  THD.MCMC<-rbind(THD.MCMC,c(th,d) )
}

plot(THD.MCMC[,1])
lines( mu[THD.MCMC[,2]])

###
pdf("fig6_4.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
ths<-seq(-6,6,length=1000)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),ylab=
       expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MC[,1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
dev.off()
###
###
pdf("fig6_5.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
Smax<-1000
ths<-seq(-6,6,length=1000)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),
     ylab=expression( paste( italic("p("),theta,")",sep="") ),lwd=2 ,ylim=c(0,.40))
hist(THD.MCMC[1:Smax,1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
plot(THD.MCMC[1:Smax,1],xlab="iteration",ylab=expression(theta))
dev.off()
####
####
pdf("fig6_6.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
Smax<-10000
ths<-seq(-6,6,length=1000)
plot(ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
       w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
       w[3]*dnorm(ths,mu[3],sqrt(s2[3])) ,type="l" , xlab=expression(theta),
     ylab=expression( paste( italic("p("),theta,")",sep="") ),lwd=2,ylim=c(0,.40) )
hist(THD.MCMC[1:Smax,1],add=TRUE,prob=TRUE,nclass=20,col="gray")
lines( ths, w[1]*dnorm(ths,mu[1],sqrt(s2[1])) +
         w[2]*dnorm(ths,mu[2],sqrt(s2[2])) +
         w[3]*dnorm(ths,mu[3],sqrt(s2[3])),lwd=2 )
plot(THD.MCMC[1:Smax,1],xlab="iteration",ylab=expression(theta))
dev.off()
#####


acf(THD.MCMC[,1],lag.max=50)

library(coda)
effectiveSize(THD.MCMC[,1])

##### MCMC diagnostics for semiconjugate normal analysis

stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }



#####
pdf("fig6_7.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(PHI[,1],xlab="iteration",ylab=expression(theta))  
plot(1/PHI[,2],xlab="iteration",ylab=expression(sigma^2))  

dev.off()





acf(PHI[,1]) -> tmp1
acf(1/PHI[,2]) -> tmp2

effectiveSize( PHI[,1] )
effectiveSize(1/PHI[,2] )


# chapter 7
library(MASS) ; source("hdr_2d.r")
###
rmvnorm<-
  function(n,mu,Sigma) {
    p<-length(mu)
    res<-matrix(0,nrow=n,ncol=p)
    if( n>0 & p>0 ) {
      E<-matrix(rnorm(n*p),n,p)
      res<-t(  t(E%*%chol(Sigma)) +c(mu))
    }
    res
  }
###


###
rinvwish<-function(n,nu0,iS0) 
{
  sL0 <- chol(iS0) 
  S<-array( dim=c( dim(L0),n ) )
  for(i in 1:n) 
  {
    Z <- matrix(rnorm(nu0 * dim(L0)[1]), nu0, dim(iS0)[1]) %*% sL0  
    S[,,i]<- solve(t(Z)%*%Z)
  }     
  S[,,1:n]
}
###

###
ldmvnorm<-function(y,mu,Sig){  # log mvn density
  c(  -(length(mu)/2)*log(2*pi) -.5*log(det(Sig)) -.5*
        t(y-mu)%*%solve(Sig)%*%(y-mu)   )  
}
####

### sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}
###



####
pdf("fig7_1.pdf",family="Times",height=1.75,width=5)

mu<-c(50,50)
sig<-c(8,12)
ng<-50
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
y<-seq(20,80,length=ng)
set.seed(1)

for( rho in c(-.5,0,.5)) {
  
  Sigma<- matrix( c(1,rho,rho,1),2,2) * outer(sig,sig)
  
  Y<-rmvnorm(30,mu,Sigma)
  
  LDY<-matrix(nrow=ng,ncol=ng)
  for(i in 1:ng){
    for(j in 1:ng){ LDY[i,j]<-ldmvnorm( c(y[i],y[j]),mu,Sigma) }}
  
  plot(range(y),range(y),type="n", 
       xlab=expression(italic(y[1])),ylab=expression(italic(y[2])) )
  filledcontour(y,y,exp(LDY))
  points(Y,pch="x")
}
dev.off()
#####

#####
Y<-dget("readpp.dat")

mu0<-c(50,50)
L0<-matrix( c(625,312.5,312.5,625),nrow=2,ncol=2)

nu0<-4
S0<-matrix( c(625,312.5,312.5,625),nrow=2,ncol=2)

n<-dim(Y)[1] ; ybar<-apply(Y,2,mean)
Sigma<-cov(Y) ; THETA<-SIGMA<-NULL
YS<-NULL
set.seed(1)

for(s in 1:5000) 
{
  
  ###update theta
  Ln<-solve( solve(L0) + n*solve(Sigma) )
  mun<-Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
  theta<-rmvnorm(1,mun,Ln)  
  ### 
  
  ###update Sigma
  Sn<- S0 + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta) ) 
  #  Sigma<-rinvwish(1,nu0+n,solve(Sn))
  Sigma<-solve( rwish(1, nu0+n, solve(Sn)) )
  ###
  
  ###
  YS<-rbind(YS,rmvnorm(1,theta,Sigma)) 
  ###
  
  ### save results 
  THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
  ###
  cat(s,round(theta,2),round(c(Sigma),2),"\n")
}

quantile(  SIGMA[,2]/sqrt(SIGMA[,1]*SIGMA[,4]), prob=c(.025,.5,.975) )
quantile(   THETA[,2]-THETA[,1], prob=c(.025,.5,.975) )
mean( THETA[,2]-THETA[,1])
mean( THETA[,2]>THETA[,1]) 
mean(YS[,2]>YS[,1])


pdf("fig7_2.pdf",family="Times",height=3.5,width=7)

par(mfrow=c(1,2),mgp=c(1.75,.75,0),mar=c(3,3,1,1))

plot.hdr2d(THETA,xlab=expression(theta[1]),ylab=expression(theta[2]) )
abline(0,1)

plot.hdr2d(YS,xlab=expression(italic(y[1])),ylab=expression(italic(y[2])), 
           xlim=c(0,100),ylim=c(0,100) )
points(Y[,1],Y[,2],pch=16,cex=.7)
abline(0,1)

dev.off()


##### get data
library(MASS) ; data(Pima.tr)
Y0<-Pima.tr[,2:5]
Y<-Y0
n<-dim(Y)[1]
p<-dim(Y)[2]

set.seed(1)
O<-matrix(rbinom(n*p,1,.9),n,p)
Y[O==0]<-NA
#####

#####
pdf("fig7_3.pdf",family="Times", height=6,width=6)
par(mar=c(1,1,.5,.5)*1.75,mfrow=c(p,p),mgp=c(1.75,.75,0))
for(j1 in 1:p) {
  for(j2 in 1:p) {
    if(j1==j2){hist(Y[,j1],main="");mtext(colnames(Y)[j1],side=3,line=-.1,cex=.7)}
    if(j1!=j2) { plot(Y[,j1],Y[,j2],xlab="",ylab="",pch=16,cex=.7)} 
  }}
dev.off()
#####

############################


### prior parameters
p<-dim(Y)[2]
mu0<-c(120,64,26,26)
sd0<-(mu0/2)
L0<-matrix(.1,p,p) ; diag(L0)<-1 ; L0<-L0*outer(sd0,sd0)
nu0<-p+2 ; S0<-L0
###

### starting values
Sigma<-S0
Y.full<-Y
for(j in 1:p)
{
  Y.full[is.na(Y.full[,j]),j]<-mean(Y.full[,j],na.rm=TRUE)
}
###


### Gibbs sampler
THETA<-SIGMA<-Y.MISS<-NULL
set.seed(1)

d45<-dget("data.f7_4.f7_5")
SIGMA<-d45$SIGMA ; THETA<-d45$THETA; Y.MISS<-d45$Y.MISS

if(!exists("d45")) {
  for(s in 1:1000)
  {
    
    ###update theta
    ybar<-apply(Y.full,2,mean)
    Ln<-solve( solve(L0) + n*solve(Sigma) )
    mun<-Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
    theta<-rmvnorm(1,mun,Ln)
    ###
    
    ###update Sigma
    Sn<- S0 + ( t(Y.full)-c(theta) )%*%t( t(Y.full)-c(theta) )
    #  Sigma<-rinvwish(1,nu0+n,solve(Sn))
    Sigma<-solve( rwish(1, nu0+n, solve(Sn)) )
    ###
    
    ###update missing data
    for(i in 1:n)
    { 
      b <- ( O[i,]==0 )
      a <- ( O[i,]==1 )
      iSa<- solve(Sigma[a,a])
      beta.j <- Sigma[b,a]%*%iSa
      s2.j   <- Sigma[b,b] - Sigma[b,a]%*%iSa%*%Sigma[a,b]
      theta.j<- theta[b] + beta.j%*%(t(Y.full[i,a])-theta[a])
      Y.full[i,b] <- rmvnorm(1,theta.j,s2.j )
    }
    
    ### save results
    THETA<-rbind(THETA,theta) ; SIGMA<-rbind(SIGMA,c(Sigma))
    Y.MISS<-rbind(Y.MISS, Y.full[O==0] )
    ###
    
    cat(s,theta,"\n")
    
  }
}
#############

apply(THETA,2,mean)

COR <- array( dim=c(p,p,1000) )
for(s in 1:1000)
{
  Sig<-matrix( SIGMA[s,] ,nrow=p,ncol=p)
  COR[,,s] <- Sig/sqrt( outer( diag(Sig),diag(Sig) ) )
}


apply(COR,c(1,2),mean)


##
##

#####
pdf("fig7_5.pdf",family="Times", height=7,width=7)
Y.true<-Y0
V<-matrix(1:p,nrow=n,ncol=p,byrow=TRUE)

v.miss<-V[O==0]
y.pred<-apply(Y.MISS,2,mean)
y.true<-Y.true[O==0]
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
for(j in 1:p){ plot(y.true[v.miss==j], y.pred[v.miss==j], 
                    xlab=paste("true", colnames(Y.true)[j]), 
                    ylab=paste("predictied", colnames(Y.true)[j]),pch=16 )
  abline(0,1)
  cat(j, mean( (y.true[v.miss==j]- y.pred[v.miss==j])^2),"\n") }
dev.off()
#####

##### convert SIGMA to an array of correlation parameters
library(sbgcop)
COR<-array(dim=c(p,p,dim(SIGMA)[1]) )
for(s in 1:dim(SIGMA)[1]) {
  Sig<-matrix( SIGMA[s,] ,nrow=p,ncol=p)
  COR[,,s] <- Sig/sqrt(outer( diag(Sig),diag(Sig)))
}
####
colnames(COR)<-rownames(COR)<-colnames(Y)

pdf("fig7_4.pdf",height=6,width=6,family="Times")

par(mfcol=c(4,2),mar=c(1,2.75,1,1),mgp=c(1.75,.75,0),oma=c(1.5,0,0,0))
plotci.sA(COR)

REG<-sR.sC(COR)
plotci.sA(REG)
dev.off()

CQ<-apply(COR, c(1,2), quantile,prob=c(.025,.5,.975) )

round(CQ[1,,],2)
round(CQ[2,,],2)
round(CQ[3,,],2)

round(apply(COR,c(1,2),mean),2)



# chpater 8
#####
dat<-dget("../Data/nels_2002_data")

XSCH<- dat[,2:13]
XSCH<- dat[,2:5]


MSCH<-NULL
for(j in 1:dim(XSCH)[2]) {
  MSCH<-cbind(MSCH, tapply(XSCH[,j],dat[,1],mean,na.rm=TRUE) )
}

SSCH<- t(  (t(MSCH)-apply(MSCH,2,mean))/apply(MSCH,2,sd) )

dsch<-as.matrix( dist( SSCH) )

nsch<-table(dat[,1])

MSCH<-cbind(sort(unique(dat[,1])) , MSCH)

id<-MSCH[,1]

i1<-1
i2<- which.min( dsch[1,-1])+1

y1<- dat[ dat[,1]==id[i1],14]
y2<- dat[ dat[,1]==id[i2],14]
###
if(22==3){
  idx<-id[  nsch < 15  & MSCH[,2]==3 & MSCH[,4]==1 & MSCH[,5]==2 ]
  PV<-matrix(0,length(idx),length(idx))
  for(i in 1:length(idx)) { for(j in 1:length(idx)) {
    #y1<- dat[ dat[,1]==idx[i],14]
    #y2<- dat[ dat[,1]==idx[j],14]
    PV[i,j]<- t.test(y1,y2,var.equal=T)$p.val
  }}
  round(PV,3)
  idx1<-idx[2]
  idx2<-idx[3]
  idx1<-idx[3]
  idx2<-idx[6]
  y1<- dat[ dat[,1]==idx1,14]
  y2<- dat[ dat[,1]==idx2,14]
  
}

###

pdf("fig8_1.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
boxplot(list(y1,y2),range=0,ylab="score",names=c("school 1","school 2"))

n1<-length(y1)
n2<-length(y2)
mean(y1)
mean(y2)
sd(c(y1,y2))
s2p<- ( var(y1)*(n1-1) + var(y2)*(n2-1) )/(n1+n2-2 )
tstat<- ( mean(y1)-mean(y2) ) /
  sqrt( s2p*(1/length(y1)+1/length(y2)))
t.test(y1,y2)


ts<-seq(-4,4,length=100)
plot(ts,dt(ts,n1+n2-1),type="l",xlab=expression(italic(t)),ylab="density")
abline(v=tstat,lwd=2,col="gray")
dev.off()

####





#####

##### data 
n1<-length(y1) ; n2<-length(y2)

##### prior parameters
mu0<-50 ; g02<-625
del0<-0 ; t02<-625
s20<-100; nu0<-1
#####

##### starting values
mu<- ( mean(y1) + mean(y2) )/2
del<- ( mean(y1) - mean(y2) )/2
#####

##### Gibbs sampler
MU<-DEL<-S2<-NULL
Y12<-NULL
set.seed(1)
for(s in 1:5000) 
{
  
  ##update s2
  s2<-1/rgamma(1,(nu0+n1+n2)/2, 
               (nu0*s20+sum((y1-mu-del)^2)+sum((y2-mu+del)^2) )/2)
  ##
  
  ##update mu
  var.mu<-  1/(1/g02+ (n1+n2)/s2 )
  mean.mu<- var.mu*( mu0/g02 + sum(y1-del)/s2 + sum(y2+del)/s2 )
  mu<-rnorm(1,mean.mu,sqrt(var.mu))
  ##
  
  ##update del
  var.del<-  1/(1/t02+ (n1+n2)/s2 )
  mean.del<- var.del*( del0/t02 + sum(y1-mu)/s2 - sum(y2-mu)/s2 )
  del<-rnorm(1,mean.del,sqrt(var.del))
  ##
  
  ##save parameter values
  MU<-c(MU,mu) ; DEL<-c(DEL,del) ; S2<-c(S2,s2) 
  Y12<-rbind(Y12,c(rnorm(2,mu+c(1,-1)*del,sqrt(s2) ) ) )
}                 
#####




plot(MU+DEL,MU-DEL)
cor(MU+DEL,MU-DEL)

pdf("fig8_2.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

plot( density(MU,adj=2),xlim=c(mu0-sqrt(g02),mu0+sqrt(g02)), 
      main="",xlab=expression(mu),ylab="density",lwd=2 )
ds<-seq(mu0-sqrt(g02),mu0+sqrt(g02),length=100)
lines(ds,dnorm(ds,mu0,sqrt(g02)),lwd=2,col="gray" )
legend(22,.27,legend=c("prior","posterior"),lwd=c(2,2),col=c("black","gray"),
       bty="n")


plot( density(DEL,adj=2),xlim=c(-sqrt(t02),sqrt(t02)),
      main="",xlab=expression(delta),ylab="density",lwd=2 )
ds<-seq(-sqrt(t02),sqrt(t02),length=100)
lines(ds,dnorm(ds,0,sqrt(t02)),lwd=2,col="gray" )
legend(-28,.27,legend=c("prior","posterior"),lwd=c(2,2),col=c("black","gray"),
       bty="n")


dev.off()

## need funciton for dinvgamma
quantile(DEL,c(.025,.5,.975))
quantile(DEL*2,c(.025,.5,.975))
mean(DEL>0)
mean(Y12[,1]>Y12[,2])


#####




stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }



#####
#hierarchical data

pdf("fig8_4.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

ids<-MSCH[ MSCH[,2]>=5  &   MSCH[,4]==1 & MSCH[,5]==1,1 ]

Y<-list()
YM<-NULL
J<-length(ids)
n<-ybar<-ymed<-s2<-rep(0,J)
for(j in 1:J) {
  Y[[j]]<- dat[ dat[,1]==ids[j],14]
  ybar[j]<-mean(Y[[j]])
  ymed[j]<-median(Y[[j]])
  n[j]<-length(Y[[j]])
  s2[j]<-var(Y[[j]])
  YM<-rbind( YM, cbind( rep(j,n[j]), Y[[j]] ))
}

colnames(YM)<-c("school","mathscore")



par(mfrow=c(1,1))
plot(c(1,J),range(Y) ,type="n",ylab="math score",xlab="rank of  school-specific math score  average")

for(l in 1:J)  {
  j<-order(ybar)[l]
  points( rep(l,n[j]), Y[[j]],pch=16,cex=.6 )
  segments( l,min(Y[[j]]),l,max(Y[[j]]))
}

abline(h=mean(ybar))

dev.off()
#####


#####
pdf("fig8_5.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
hist(ybar,main="",xlab="sample mean")
plot(n,ybar,xlab="sample size",ylab="sample mean")
dev.off()

#####
pdf("fig8_3.pdf",family="Times",height=4,width=8)
par(mar=c(0,0,0,0),mgp=c(0,0,0))
plot(c(.15,.85),c(.15,.85),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")

text( .5,.8 , expression(paste( mu ,",", tau^2 )) )

text(  c(.2,.3,.5,.7,.8),  rep(.6,5),
       expression( theta[1],theta[2],".....",theta[m-1],theta[m]) )

segments(  rep(.5,4),rep(.75,4), c(.2,.3,.7,.8),rep(.65,4) )

text( c(.2,.3,.5,.7,.8), rep(.4,5),
      expression( bar(y)[1], bar(y)[2],".....",bar(y)[m-1],bar(y)[m]) )

segments(  c(.2,.3,.7,.8),rep(.55,4), c(.2,.3,.7,.8),rep(.45,4) )


segments(  rep(.5,4),rep(.25,4), c(.2,.3,.7,.8),rep(.35,4) )

text( .5,.2,expression(sigma^2))
dev.off()
#####


##### MCMC analysis for school data

### weakly informative priors
nu0<-1  ; s20<-100
eta0<-1 ; t20<-100
mu0<-50 ; g20<-25
###

### starting values
m<-length(Y) 
n<-sv<-ybar<-rep(NA,m) 
for(j in 1:m) 
{ 
  ybar[j]<-mean(Y[[j]])
  sv[j]<-var(Y[[j]])
  n[j]<-length(Y[[j]]) 
}
theta<-ybar
sigma2<-mean(sv)
mu<-mean(theta)
tau2<-var(theta)
###

### setup MCMC
set.seed(1)
S<-5000
THETA<-matrix( nrow=S,ncol=m)
MST<-matrix( nrow=S,ncol=3)
###

### MCMC algorithm
for(s in 1:S) 
{
  
  # sample new values of the thetas
  for(j in 1:m) 
  {
    vtheta<-1/(n[j]/sigma2+1/tau2)
    etheta<-vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j]<-rnorm(1,etheta,sqrt(vtheta))
  }
  
  #sample new value of sigma2
  nun<-nu0+sum(n)
  ss<-nu0*s20;for(j in 1:m){ss<-ss+sum((Y[[j]]-theta[j])^2)}
  sigma2<-1/rgamma(1,nun/2,ss/2)
  
  #sample a new value of mu
  vmu<- 1/(m/tau2+1/g20)
  emu<- vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu<-rnorm(1,emu,sqrt(vmu)) 
  
  # sample a new value of tau2
  etam<-eta0+m
  ss<- eta0*t20 + sum( (theta-mu)^2 )
  tau2<-1/rgamma(1,etam/2,ss/2)
  
  #store results
  THETA[s,]<-theta
  MST[s,]<-c(mu,sigma2,tau2)
  
} 
###

mcmc1<-list(THETA=THETA,MST=MST)

pdf("fig8_6.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))


stationarity.plot(MST[,1],xlab="iteration",ylab=expression(mu))
stationarity.plot(MST[,2],xlab="iteration",ylab=expression(sigma^2))
stationarity.plot(MST[,3],xlab="iteration",ylab=expression(tau^2))

dev.off()

library(coda)
effectiveSize(MST)
par(mfrow=c(1,3))
acf(MST[,1]) -> a1
acf(MST[,2]) -> a2
acf(MST[,3]) -> a3

effectiveSize(THETA) -> esTHETA


MCERR<-  apply(MST,2,sd)/sqrt( effectiveSize(MST) )
apply(MST,2,mean)

TMCERR<-  apply(THETA,2,sd)/sqrt( effectiveSize(THETA) )

#####

#####
pdf("fig8_7.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(density(MST[,1],adj=2),xlab=expression(mu),main="",lwd=2,
     ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,1],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
plot(density(MST[,2],adj=2),xlab=expression(sigma^2),main="", lwd=2,
     ylab=expression(paste(italic("p("),sigma^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,2],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
plot(density(MST[,3],adj=2),xlab=expression(tau^2),main="",lwd=2,
     ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,3],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
dev.off()
#####

mean((MST[,1]))
mean(sqrt(MST[,2]))
mean(sqrt(MST[,3]))



#####
pdf("fig8_8.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

theta.hat<-apply(THETA,2,mean)
plot(ybar,theta.hat,xlab=expression(bar(italic(y))),ylab=expression(hat(theta)))
abline(0,1)
plot(n,ybar-theta.hat,ylab=expression( bar(italic(y))-hat(theta) ),xlab="sample size")
abline(h=0)
dev.off()
#####


#####
pdf("fig8_9.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

theta.order<-order(theta.hat)
theta.order[1:20]

idx<-c(46,82)


ybar.order<-order(ybar)
ybar.order[1:20]

ybar[c(46,82)]
n[c(46,82)]
theta.hat[c(46,82)]


mean(THETA[,46]<THETA[,82])

par(mfrow=c(1,1))
plot(density(THETA[,46],adj=2),col="black",xlim=
       range(c(Y[[46]],Y[[82]],THETA[,c(46,82)])),lwd=2,
     main="",xlab="math score",ylim=c(-.05,.2),ylab="",yaxt="n")
axis(side=2,at=c(0,0.10,0.20) )
lines(density(THETA[,82],adj=2),col="gray",lwd=2)
abline(h=0)

points( Y[[46]],rep(-0.01666,n[46]), col="black",pch=16)
points( ybar[46],-.01666,col="black",pch=16 ,cex=1.5)
abline( h=-.01666,col="black")

points( Y[[82]],rep(-0.0333,n[82]), col="gray",pch=16)
points( ybar[82],-.0333,col="gray",pch=16 ,cex=1.5)
abline( h=-.0333,col="gray")

segments(mean(MST[,1]), 0,mean(MST[,1]),1,lwd=2,lty=2 )

legend(52.5,.15,legend=c("school 46","school 82",
                         expression(paste("E[", mu,"|",italic(y[1]),"...",italic(y[m]),"]"))),
       lwd=c(2,2),lty=c(1,1,2),col=c("black","gray"),bty="n")

dev.off()




##### hm for mean and variance
### weakly informative priors
nu0<-1  ; s20<-100
eta0<-1 ; t20<-100
mu0<-50 ; g20<-25
a0<-1 ; b0<-1/100 ; wnu0<-1
###

### starting values
m<-length(Y)
n<-sv<-ybar<-rep(NA,m)
for(j in 1:m)
{
  ybar[j]<-mean(Y[[j]])
  sv[j]<-var(Y[[j]])
  n[j]<-length(Y[[j]])
}
theta<-ybar
sigma2<-sv 
mu<-mean(theta)
tau2<-var(theta)
s20<-1/mean(1/sv)
nu0<-10
###

### setup MCMC
set.seed(1)
S<-5000
SIGMA2<-THETA<-matrix( nrow=S,ncol=m)
MTSN<-matrix( nrow=S,ncol=4)
s2.pp<-NULL
###
nu0s<-1:5000
### MCMC algorithm
for(s in 1:S)
{
  # sample new values of the thetas
  for(j in 1:m)
  {
    vtheta<-1/(n[j]/sigma2[j]+1/tau2)
    etheta<-vtheta*(ybar[j]*n[j]/sigma2[j]+mu/tau2)
    theta[j]<-rnorm(1,etheta,sqrt(vtheta))
  }
  
  #sample new value the sigma2s
  for(j in 1:m) 
  { 
    nun<-nu0+n[j]
    ss<-nu0*s20+ sum((Y[[j]]-theta[j])^2)
    sigma2[j]<-1/rgamma(1,nun/2,ss/2)
  }
  
  #sample new s20
  s20<-rgamma(1,a0+m*nu0/2,b0+nu0*sum(1/sigma2)/2)
  
  lpnu0<- .5*nu0s*m*log(s20*nu0s/2)-m*lgamma(nu0s/2)+(nu0s/2-1)*sum(log(1/sigma2)) -
    nu0s*s20*sum(1/sigma2)/2   - wnu0*nu0s
  nu0<-sample(nu0s,1,prob=exp( lpnu0-max(lpnu0)) )
  
  #sample a new value of mu
  vmu<- 1/(m/tau2+1/g20)
  emu<- vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu<-rnorm(1,emu,sqrt(vmu))
  
  # sample a new value of tau2
  etam<-eta0+m
  ss<- eta0*t20 + sum( (theta-mu)^2 )
  tau2<-1/rgamma(1,etam/2,ss/2)
  
  #store results
  THETA[s,]<-theta
  SIGMA2[s,]<-sigma2
  MTSN[s,]<-c(mu,tau2,s20,nu0)
  s2.pp<-c(s2.pp,1/rgamma(1,nu0/2,nu0*s20/2))
  if(s %%25 ==0) { hist(1/sigma2,prob=T) ;
    x<-seq(0,20,100); lines(x,dgamma(x,nu0/2,nu0*s20/2)) }
  
}
###

mcmc2<-list(THETA=THETA,SIGMA2=SIGMA2,MTSN=MTSN,s2.pp=s2.pp)

#compare post pred of sigma to inv gamma?

par(mfrow=c(1,2))
apply(SIGMA2,2,mean) -> sigma2.hat
plot(sv,sigma2.hat)
abline(0,1)
plot(n, sv-sigma2.hat)  ; abline(h=0)

dinvgamma<-function(x,a,b) {
  ld<- a*log(b) -lgamma(a) -(a+1)*log(x)  -b/x
  exp(ld)
}

pdf("fig8_11.pdf",family="Times",height=3.5,width=6)
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(density(MST[,1],adj=2),lwd=2,main="",col="gray",xlab=expression(mu),
     ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")))
lines(density(MTSN[,1],adj=2),lwd=2)

plot(density(MST[,3],adj=2),lwd=2,main="",col="gray",xlab=expression(tau^2), 
     ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
lines(density(MTSN[,2],adj=2),lwd=2)

plot(table(MTSN[,4]),xlab=expression(nu[0]),
     ylab=expression(paste(italic("p("),nu[0],"|",italic(y[1]),"...",italic(y[m]),")")))


plot(density(MTSN[,3],adj=2),lwd=2,main="",xlab=expression(sigma[0]^2),
     ylab=expression(paste(italic("p("),sigma[0]^2,"|",italic(y[1]),"...",italic(y[m]),")")))


dev.off()

pdf("fig8_12.pdf",family="Times",height=3.5,width=7)
apply(SIGMA2,2,mean) -> sigma2.hat
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(sv,sigma2.hat,xlab=expression(s^2),ylab=expression(hat( sigma^2)) )
abline(0,1)
plot(n, sv-sigma2.hat,xlab="sample size",ylab=expression(s^2-hat(sigma^2)))  
abline(h=0)

dev.off()

bsh<-(1:m)[ abs(sv-sigma2.hat)>30 ]



hist(sv,prob=TRUE,nclass=15) ; lines( density(s2.pp,adj=2))
hist(sigma2,prob=TRUE,nclass=20) ; lines( density(s2.pp,adj=2))











##### 



stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }

# chapter 9
## simple mixed effects model data -library(nlme) ; data(PBG)
source("regression_gprior.r")
source("backselect.r")
#####
pdf("fig9_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))

x1<-c(0,0,0,0,0,0,1,1,1,1,1,1)
x2<-c(23,22,22,25,27,20,31,23,27,28,22,24)
y<-c(-0.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,0.26,2.51)

par(mfrow=c(1,1))
plot(y~x2,pch=16,xlab="age",ylab="change in maximal oxygen uptake", 
     col=c("black","gray")[x1+1])
legend(27,0,legend=c("aerobic","running"),pch=c(16,16),col=c("gray","black"))

dev.off()
#####


#####
pdf("fig9_2.pdf",family="Times",height=5.5,width=6)
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))


plot(y~x2,pch=16,col=c("black","gray")[x1+1],ylab="change in maximal oxygen uptake",xlab="",xaxt="n")
abline(h=mean(y[x1==0]),col="black") 
abline(h=mean(y[x1==1]),col="gray")
mtext(side=3,expression(paste(beta[3]==0,"  ",beta[4]==0)) )

plot(y~x2,pch=16,col=c("black","gray")[x1+1],xlab="",ylab="",xaxt="n",yaxt="n")
abline(lm(y~x2),col="black")
abline(lm((y+.5)~x2),col="gray")
mtext(side=3,expression(paste(beta[2]==0,"  ",beta[4]==0)) )

plot(y~x2,pch=16,col=c("black","gray")[x1+1],
     xlab="age",ylab="change in maximal oxygen uptake" )
fit<-lm( y~x1+x2)
abline(a=fit$coef[1],b=fit$coef[3],col="black")
abline(a=fit$coef[1]+fit$coef[2],b=fit$coef[3],col="gray")
mtext(side=3,expression(beta[4]==0)) 

plot(y~x2,pch=16,col=c("black","gray")[x1+1],
     xlab="age",ylab="",yaxt="n")
abline(lm(y[x1==0]~x2[x1==0]),col="black")
abline(lm(y[x1==1]~x2[x1==1]),col="gray")

dev.off()
#####


##### LS estimation 
n<-length(y)
X<-cbind(rep(1,n),x1,x2,x1*x2)
p<-dim(X)[2]

beta.ols<- solve(t(X)%*%X)%*%t(X)%*%y


## probability of intersection?


####
n<-length(y)
X<-cbind(rep(1,n),x1,x2,x1*x2)
p<-dim(X)[2]

fit.ls<-lm(y~-1+ X)
beta.0<-rep(0,p) ; Sigma.0<-diag(c(150,30,6,5)^2,p)
nu.0<-1 ; sigma2.0<- 15^2

beta.0<-fit.ls$coef
nu.0<-1  ; sigma2.0<-sum(fit.ls$res^2)/(n-p)
Sigma.0<- solve(t(X)%*%X)*sigma2.0*n


S<-5000
####

rmvnorm<-function(n,mu,Sigma) 
{ # samples from the multivariate normal distribution
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}
###

### some convenient quantites
n<-length(y)
p<-length(beta.0)
iSigma.0<-solve(Sigma.0)
XtX<-t(X)%*%X

### store mcmc samples in these objects
beta.post<-matrix(nrow=S,ncol=p)
sigma2.post<-rep(NA,S)

### starting value
set.seed(1)
sigma2<- var( residuals(lm(y~0+X)) )

### MCMC algorithm
for( scan in 1:S) {
  
  #update beta
  V.beta<- solve(  iSigma.0 + XtX/sigma2 )
  E.beta<- V.beta%*%( iSigma.0%*%beta.0 + t(X)%*%y/sigma2 )
  beta<-t(rmvnorm(1, E.beta,V.beta) )
  
  #update sigma2
  nu.n<- nu.0+n
  ss.n<-nu.0*sigma2.0 + sum(  (y-X%*%beta)^2 )
  sigma2<-1/rgamma(1,nu.n/2, ss.n/2)
  
  #save results of this scan
  beta.post[scan,]<-beta
  sigma2.post[scan]<-sigma2
}

#####
round( apply(beta.post,2,mean), 3)

#####
tmp<-lm.gprior(y,X )
beta.post<-tmp$beta
beta.ols<-lm(y~-1+X)$coef
g<-n ; nu0=1 ; s20<-summary( lm(y~ -1+X))$sigma^2 
beta.ols*g/(g+1)
iXX<-solve(t(X)%*%X)

mdt<-function(t,mu,sig,nu){ 
  
  gamma(.5*(nu+1))*(1+ ( (t-mu)/sig )^2/nu )^(-.5*(nu+1))/ 
    ( sqrt(nu*pi)*sig* gamma(nu/2)  )
}

#####
pdf("fig9_3.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

x<-seq(-85,130,length=200)
plot(density(beta.post[,2],adj=2),xlab=expression(beta[2]),main="",ylab="",lwd=2)
abline(v=0,col="gray")
lines(x,mdt(x,0,sqrt(n*s20*iXX[2,2]),nu0 ),col="gray")

x<-seq(-5,5,length=100)
plot(density(beta.post[,4],adj=2),xlab=expression(beta[4]),main="",ylab="",lwd=2)
abline(v=0,col="gray")
lines(x,mdt(x,0,sqrt(n*s20*iXX[4,4]),nu0 ),col="gray")


source("hdr_2d.r")
plot.hdr2d( beta.post[,c(2,4)],xlab=expression(beta[2]),
            ylab=expression(beta[4]))
abline(h=0,col="gray") ; abline(v=0,col="gray")
dev.off()
#####

BX<-NULL
for(s in 1:dim(beta.post)[1]) { 
  BX<-rbind(BX, beta.post[s,2] + (min(X[,3]):max(X[,3]))*beta.post[s,4] )
}

###
########

###
qboxplot<-function(x,at=0,width=.5,probs=c(.025,.25,.5,.75,.975))
{
  qx<-quantile(x,probs=probs)
  segments(at,qx[1],at,qx[5])
  polygon(x=c(at-width,at+width,at+width,at-width),
          y=c(qx[2],qx[2],qx[4],qx[4]) ,col="gray")
  segments(at-width,qx[3],at+width,qx[3],lwd=3)
  segments(at-width/2,qx[1],at+width/2,qx[1],lwd=1)
  segments(at-width/2,qx[5],at+width/2,qx[5],lwd=1)
} 
###



######
pdf("fig9_4.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(range(X[,3]),range(y),type="n",xlab="age",
     #   ylab="expected difference in change score")
     ylab=expression(paste( beta[2] + beta[4],"age",sep="") ) )
for( age  in  1:dim(BX)[2]  ) {
  qboxplot( BX[,age] ,at=age+19 , width=.25) }  

abline(h=0,col="gray")
dev.off() 
########



####################
library(lars) ; data(diabetes)
yf<-diabetes$y
yf<-(yf-mean(yf))/sd(yf)

Xf<-diabetes[[3]]
Xf<-t( (t(Xf)-apply(Xf,2,mean))/apply(Xf,2,sd))

###
n<-length(yf)
set.seed(1)

i.te<-sample(1:n,100)
i.tr<-(1:n)[-i.te]

y<-yf[i.tr] ; y.te<-yf[i.te]
X<-Xf[i.tr,]; X.te<-Xf[i.te,]

yperm<-sample(y)
######################


pdf("fig9_5.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.5,.5,0))
olsfit<-lm(y~-1+X)
y.te.ols<-X.te%*%olsfit$coef
plot(y.te,y.te.ols,xlab=expression(italic(y)[test]),
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)
mean( (y.te-y.te.ols )^2 )
plot(olsfit$coef,type="h",lwd=2,xlab="regressor index",ylab=expression(hat(beta)[ols]))
###

### back elim
vars<-bselect.tcrit(y,X,tcrit=1.65)
bslfit<-lm(y~-1+X[,vars$remain])
y.te.bsl<-X.te[,vars$remain]%*%bslfit$coef
mean( (y.te-y.te.bsl)^2)
plot(y.te,y.te.bsl,ylim=range( c(y.te.bsl,y.te.ols)),
     xlab=expression(italic(y)[test]),ylab=expression(hat(italic(y))[test]))
abline(0,1)
###
dev.off()




#####
### back elim with permuted data
pdf("fig9_6.pdf",family="Times",height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
fit.perm<-lm(yperm~-1+X)
t.perm<-summary(fit.perm)$coef[,3]
b.perm<-summary(fit.perm)$coef[,1]
plot(t.perm,type="h",lwd=2,xlab="regressor index",ylab="t-statistic",ylim=c(-4.8,4.8))

vars.perm<-bselect.tcrit(yperm,X,tcrit=1.65)
bslfit.perm<-lm(yperm~-1+X[,vars.perm$remain])
t.bslperm<-t.perm*0
b.bslperm<-b.perm*0
t.bslperm[vars.perm$remain]<-summary(bslfit.perm)$coef[,3]
b.bslperm[vars.perm$remain]<-summary(bslfit.perm)$coef[,1]
plot(t.bslperm,type="h",lwd=2,xlab="regressor index",ylab="t-statistic",
     ylim=c(-4.8,4.8) )
dev.off()
#####



#####
source("regression_gprior.r")
#load("diabetes.bma")
tmp<-dget("diabetes.bma")
if(2==3) {
  par(mfrow=c(1,2))
  p<-dim(X)[2]
  S<-10000
  BETA<-Z<-matrix(NA,S,p)
  z<-rep(1,dim(X)[2] )
  lpy.c<-lpy.X(y,X[,z==1,drop=FALSE])
  
  for(s in 1:S)
  {
    for(j in sample(1:p))
    {
      zp<-z ; zp[j]<-1-zp[j]
      lpy.p<-lpy.X(y,X[,zp==1,drop=FALSE])
      r<- (lpy.p - lpy.c)*(-1)^(zp[j]==0)
      z[j]<-rbinom(1,1,1/(1+exp(-r)))
      if(z[j]==zp[j]) {lpy.c<-lpy.p}
    }
    beta<-z;if(sum(z)>0){beta[z==1]<-lm.gprior(y,X[,z==1,drop=FALSE],S=1)$beta }
    Z[s,]<-z
    BETA[s,]<-beta
    if(s>1 ) {
      bpm<-apply(BETA[1:s,],2,mean) ; plot(bpm)
      cat(s,mean(z), mean( (y.te-X.te%*%bpm)^2),"\n")
      Zcp<- apply(Z[1:s,,drop=FALSE],2,cumsum)/(1:s)
      plot(c(1,s),range(Zcp),type="n") ; apply(Zcp,2,lines)
    }
    
  }
  save.image("diabetes.bma")
}


BETA<-tmp$BETA ; Z<-tmp$Z

#####

pdf("fig9_7.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

beta.bma<-apply(BETA,2,mean,na.rm=TRUE)
y.te.bma<-X.te%*%beta.bma
mean( (y.te-y.te.bma)^2)

layout( matrix(c(1,1,2),nrow=1,ncol=3) )

plot(apply(Z,2,mean,na.rm=TRUE),xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2)

plot(y.te,y.te.bma,xlab=expression(italic(y)[test]),
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)

dev.off()

# chapter 10
fdat<-read.table("female.dat",header=TRUE,na.string=".")
fss<- tapply(fdat$fpop,fdat$year,length)
fk<-fdat[fdat$year==14,c(4,5)  ]
#fk<-fdat[,c(4,5)  ]
fk<-fk[ fk[,1]<7 ,]
age<-fk[,1]
spf<-fk[,2]
age2<-age^2
plot(spf~as.factor(age),range=0)
summary(glm(spf~age,family="poisson"))
summary(glm(spf~age+age2,family="poisson"))
##########

########## the data above can be obtained via 
########## yX.sparrow<-dget("http://www.stat.washington.edu/~hoff/Book/Data/data/yX.sparrow")



pdf("fig10_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(spf~as.factor(age),range=0,xlab="age",ylab="offspring",
     col="gray")

#spf<-spfage$spf; age<-spfage$age ; age2<-age^2
summary(glm(spf~age+age2,family="poisson"))

dev.off()

p<-3
beta0<-rep(0,p)
S0<-diag( rep(100,3))
gs<-100
LPB<-array(0,dim=rep(gs,p))

beta1<-seq(.27-1.75,.27+1.75,length=gs)
beta2<-seq(.68-1.5,.68+1.5,length=gs)
beta3<-seq(-.13-.25,-.13+.25,length=gs)

beta1<-seq(.27-2.5,.27+2.5,length=gs)
beta2<-seq(.68-2,.68+2,length=gs)
beta3<-seq(-.13-.5,-.13+.5,length=gs)



for(i in 1:gs) { for(j in 1:gs) { for(k in 1:gs) {
  theta<-beta1[i]+beta2[j]*age+beta3[k]*age^2
  LPB[i,j,k]<-dnorm(beta1[i],beta0[1],sqrt(S0[1,1]),log=TRUE)  +
    dnorm(beta2[j],beta0[2],sqrt(S0[2,2]),log=TRUE)  +
    dnorm(beta3[k],beta0[3],sqrt(S0[3,3]),log=TRUE)  +
    sum( dpois(spf,exp(theta),log=TRUE )  )
}} 
  cat(i,"\n") }

PB<-exp( LPB - max(LPB) )
PB<-PB/sum(PB)

PB1<-apply(PB,1,sum)
PB2<-apply(PB,2,sum)
PB3<-apply(PB,3,sum)
PB23<-apply(PB,c(2,3),sum)

S<-50000
BETAg<-matrix(nrow=S,ncol=3)
for(s in 1:S) {
  i<-sample(1:gs,1,prob=PB2)
  j<-sample(1:gs,1,prob=PB23[i,] )
  k<-sample(1:gs,1,prob=PB[,i,j] )
  BETAg[s,]<-c(beta1[k],beta2[i],beta3[j])  }


pdf("fig10_2.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(beta2,PB2*length(beta2)/(max(beta2)-min(beta2)) ,type="l",xlab=expression(beta[2]),ylab=expression(paste(italic("p("),beta[2],"|",italic("y)"),sep="") ) )
plot(beta3,PB3*length(beta3)/(max(beta3)-min(beta3)),type="l",xlab=expression(beta[3]),ylab=expression(paste(italic("p("),beta[3],"|",italic("y)"),sep="") ))

Xs<-cbind(rep(1,6),1:6,(1:6)^2)
eXB.post<- exp(t(Xs%*%t(BETAg )) )
qE<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

source("hdr_2d.r")
plot.hdr2d(BETAg[,2:3],bw=c(15,15),xlab=expression(beta[2]),
           ylab=expression(beta[3]))

#plot( c(1,6),range(c(0,qE)),type="n",xlab="age",
#   ylab="expected number of offspring")
#lines( qE[1,],col="black",lwd=1)
#lines( qE[2,],col="black",lwd=2)
#lines( qE[3,],col="black",lwd=1)
#legend(1.25,1.15,legend=c("97.5%","50%","2.5%"),lwd=c(1,2,1),bty="n")
dev.off()


Xs<-cbind(rep(1,6),1:6,(1:6)^2)
eXB.post<- exp(t(Xs%*%t(BETAg )) )
qEg<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

plot( c(1,6),range(c(0,qE)),type="n",xlab="age",
      ylab="expected number of offspring")
lines( qE[1,],col="black",lwd=1)
lines( qE[2,],col="black",lwd=2)
lines( qE[3,],col="black",lwd=1)



# MH algorithm for one-sample normal problem with 
# known variance

s2<-1 
t2<-10 ; mu<-5

set.seed(1)
n<-5
y<-round(rnorm(n,10,1),2)

mu.n<-( mean(y)*n/s2 + mu/t2 )/( n/s2+1/t2) 
t2.n<-1/(n/s2+1/t2)


#####
s2<-1 ; t2<-10 ; mu<-5 
y<-c(9.37, 10.18, 9.16, 11.60, 10.33)
theta<-0 ; delta<-2 ; S<-10000 ; THETA<-NULL ; set.seed(1)

for(s in 1:S)
{
  
  theta.star<-rnorm(1,theta,sqrt(delta))
  
  log.r<-( sum(dnorm(y,theta.star,sqrt(s2),log=TRUE)) +
             dnorm(theta.star,mu,sqrt(t2),log=TRUE) )  -
    ( sum(dnorm(y,theta,sqrt(s2),log=TRUE)) +
        dnorm(theta,mu,sqrt(t2),log=TRUE) ) 
  
  if(log(runif(1))<log.r) { theta<-theta.star }
  
  THETA<-c(THETA,theta)
  
}
#####

pdf("fig10_3.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

skeep<-seq(10,S,by=10)
plot(skeep,THETA[skeep],type="l",xlab="iteration",ylab=expression(theta))

hist(THETA[-(1:50)],prob=TRUE,main="",xlab=expression(theta),ylab="density")
th<-seq(min(THETA),max(THETA),length=100)
lines(th,dnorm(th,mu.n,sqrt(t2.n)) )
dev.off()


#### MH
par(mfrow=c(4,3))
ACR<-ACF<-NULL
THETAA<-NULL
for(delta2 in 2^c(-5,-1,1,5,7) ) {
  set.seed(1)
  THETA<-NULL
  S<-10000
  theta<-0
  acs<-0
  delta<-2
  
  for(s in 1:S) 
  {
    
    theta.star<-rnorm(1,theta,sqrt(delta2))
    log.r<-sum( dnorm(y,theta.star,sqrt(s2),log=TRUE)-
                  dnorm(y,theta,sqrt(s2),log=TRUE)  )  +
      dnorm(theta.star,mu,sqrt(t2),log=TRUE)-dnorm(theta,mu,sqrt(t2),log=TRUE) 
    
    if(log(runif(1))<log.r)  { theta<-theta.star ; acs<-acs+1 }
    THETA<-c(THETA,theta) 
    
  }
  plot(THETA[1:1000])
  
  ACR<-c(ACR,acs/s) 
  ACF<-c(ACF,acf(THETA,plot=FALSE)$acf[2]  )
  THETAA<-cbind(THETAA,THETA)
}
plot(ACR,ACF) ; lines(ACR,ACF)
#####

pdf("fig10_4.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
laby<-c(expression(theta),"","","","")

for(k in c(1,3,5)) {
  plot(THETAA[1:500,k],type="l",xlab="iteration",ylab=laby[k], 
       ylim=range(THETAA) )
  abline(h=mu.n,lty=2)
}
dev.off()

THCM<-apply(THETAA,2,cumsum)
THCM<- THCM/(1:dim(THCM)[1]) 


###################


#####
fit.mle<-glm(spf~age+age2,family="poisson")
summary(fit.mle)

y<-spf ; X<-cbind(rep(1,length(y)),age,age^2)
yX<-cbind(y,X)
colnames(yX)<-c("fledged","intercept","age","age2") 
dput(yX,"yX.sparrow")

source("~hoff/USBWork/rfunctions.r")
n<-length(y) ; p<-dim(X)[2]

pmn.beta<-rep(0,p)
psd.beta<-rep(10,p)

var.prop<- var(log(y+1/2))*solve( t(X)%*%X )
beta<-rep(0,p)
S<-10000
BETA<-matrix(0,nrow=S,ncol=p)
ac<-0
set.seed(1)

for(s in 1:S) {
  
  #propose a new beta
  
  beta.p<- t(rmvnorm(1, beta, var.prop ))
  
  lhr<- sum(dpois(y,exp(X%*%beta.p),log=T)) -
    sum(dpois(y,exp(X%*%beta),log=T)) +
    sum(dnorm(beta.p,pmn.beta,psd.beta,log=T)) -
    sum(dnorm(beta,pmn.beta,psd.beta,log=T))
  
  if( log(runif(1))< lhr ) { beta<-beta.p ; ac<-ac+1 }
  
  BETA[s,]<-beta
}
cat(ac/S,"\n")

#######

library(coda)
apply(BETA,2,effectiveSize)



####
pdf("fig10_5.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
blabs<-c(expression(beta[1]),expression(beta[2]),expression(beta[3]))
thin<-c(1,(1:1000)*(S/1000))
j<-3
plot(thin,BETA[thin,j],type="l",xlab="iteration",ylab=blabs[j])
abline(h=mean(BETA[,j]) )

acf(BETA[,j],ci.col="gray",xlab="lag")
acf(BETA[thin,j],xlab="lag/10",ci.col="gray")
dev.off()
####


####
pdf("fig10_6.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot(beta2,PB2*length(beta2)/(max(beta2)-min(beta2)) ,type="l",xlab=expression(beta[2]),ylab=expression(paste(italic("p("),beta[2],"|",italic("y)"),sep="") ) ,lwd=2,lty=2,col="gray")
lines(density(BETA[,2],adj=2),lwd=2)

plot(beta3,PB3*length(beta3)/(max(beta3)-min(beta3)),type="l",xlab=expression(beta[3]),ylab=expression(paste(italic("p("),beta[3],"|",italic("y)"),sep="") ),lwd=2,col="gray",lty=2)
lines(density(BETA[,3],adj=2),lwd=2)

Xs<-cbind(rep(1,6),1:6,(1:6)^2) 
eXB.post<- exp(t(Xs%*%t(BETA )) )
qE<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

plot( c(1,6),range(c(0,qE)),type="n",xlab="age",
      ylab="number of offspring")
lines( qE[1,],col="black",lwd=1)
lines( qE[2,],col="black",lwd=2)
lines( qE[3,],col="black",lwd=1)


dev.off()
####

#############

dtmp<-as.matrix(read.table("vostok.1999.temp.dat",header=TRUE))
dco2<-as.matrix(read.table("vostok.icecore.co2.dat",header=TRUE))
dtmp[,2]<- -dtmp[,2]
dco2[,2]<- -dco2[,2]
library(nlme)

#### get evenly spaced temperature points
ymin<-max( c(min(dtmp[,2]),min(dco2[,2])))
ymax<-min( c(max(dtmp[,2]),max(dco2[,2])))
n<-200
syear<-seq(ymin,ymax,length=n)
dat<-NULL
for(i in 1:n) {
  tmp<-dtmp[ dtmp[,2]>=syear[i] ,]
  dat<-rbind(dat,  tmp[dim(tmp)[1],c(2,4)] )
}
dat<-as.matrix(dat)
####

####
dct<-NULL
for(i in 1:n) {
  xc<-dco2[ dco2[,2] < dat[i,1] ,,drop=FALSE]
  xc<-xc[ 1, ]
  dct<-rbind(dct, c( xc[c(2,4)], dat[i,] ) )
}

mean( dct[,3]-dct[,1])


dct<-dct[,c(3,2,4)]
colnames(dct)<-c("year","co2","tmp")
rownames(dct)<-NULL
dct<-as.data.frame(dct)


########
pdf("fig10_7.pdf",family="Times",height=1.75,width=5)

par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
layout(matrix( c(1,1,2),nrow=1,ncol=3) )

#plot(dct[,1],qnorm( rank(dct[,3])/(length(dct[,3])+1 )) ,
plot(dct[,1],  (dct[,3]-mean(dct[,3]))/sd(dct[,3]) ,
     type="l",col="black",
     xlab="year",ylab="standardized measurement",ylim=c(-2.5,3))
legend(-115000,3.2,legend=c("temp",expression(CO[2])),bty="n",
       lwd=c(2,2),col=c("black","gray"))
lines(dct[,1],  (dct[,2]-mean(dct[,2]))/sd(dct[,2]),
      #lines(dct[,1],qnorm( rank(dct[,2])/(length(dct[,2])+1 )),
      type="l",col="gray")

plot(dct[,2], dct[,3],xlab=expression(paste(CO[2],"(ppmv)")),ylab="temperature difference (deg C)")

dev.off()
########



########
pdf("fig10_8.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

lmfit<-lm(dct$tmp~dct$co2)
hist(lmfit$res,main="",xlab="residual",ylab="frequency")
#plot(dct$year, lmfit$res,xlab="year",ylab="residual",type="l" ); abline(h=0)
acf(lmfit$res,ci.col="gray",xlab="lag")
dev.off()
########



######## starting values
n<-dim(dct)[1]
y<-dct[,3]
X<-cbind(rep(1,n),dct[,2])
DY<-abs(outer( (1:n),(1:n) ,"-"))

lmfit<-lm(y~-1+X)
fit.gls <- gls(y~X[,2], correlation=corARMA(p=1), method="ML")
beta<-lmfit$coef
s2<-summary(lmfit)$sigma^2
phi<-acf(lmfit$res,plot=FALSE)$acf[2]
nu0<-1 ; s20<-1 ; T0<-diag(1/1000,nrow=2)
###
set.seed(1)
S<-1000 ; odens<-S/1000
OUT<-NULL ; ac<-0 ; par(mfrow=c(1,2))
for(s in 1:S)
{
  
  Cor<-phi^DY  ; iCor<-solve(Cor)
  V.beta<- solve( t(X)%*%iCor%*%X/s2 + T0)
  E.beta<- V.beta%*%( t(X)%*%iCor%*%y/s2  )
  beta<-t(rmvnorm(1,E.beta,V.beta)  )
  
  s2<-1/rgamma(1,(nu0+n)/2,(nu0*s20+t(y-X%*%beta)%*%iCor%*%(y-X%*%beta)) /2 )
  
  phi.p<-abs(runif(1,phi-.1,phi+.1))
  phi.p<- min( phi.p, 2-phi.p)
  lr<- -.5*( determinant(phi.p^DY,log=TRUE)$mod -
               determinant(phi^DY,log=TRUE)$mod  +
               tr( (y-X%*%beta)%*%t(y-X%*%beta)%*%(solve(phi.p^DY) -solve(phi^DY)) )/s2 )
  
  if( log(runif(1)) < lr ) { phi<-phi.p ; ac<-ac+1 }
  
  if(s%%odens==0)
  {
    cat(s,ac/s,beta,s2,phi,"\n") ; OUT<-rbind(OUT,c(beta,s2,phi))
    #      par(mfrow=c(2,2))
    #      plot(OUT[,1]) ; abline(h=fit.gls$coef[1])
    #      plot(OUT[,2]) ; abline(h=fit.gls$coef[2])
    #      plot(OUT[,3]) ; abline(h=fit.gls$sigma^2)
    #      plot(OUT[,4]) ; abline(h=.8284)
    
  }
}
#####

OUT.1000<-OUT
library(coda)
apply(OUT,2,effectiveSize )


OUT.25000<-dget("data.f10_10.f10_11")
apply(OUT.25000,2,effectiveSize )


pdf("fig10_9.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(OUT.1000[,4],xlab="scan",ylab=expression(rho),type="l")
acf(OUT.1000[,4],ci.col="gray",xlab="lag")
dev.off()


pdf("fig10_10.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(OUT.25000[,4],xlab="scan/25",ylab=expression(rho),type="l")
acf(OUT.25000[,4],ci.col="gray",xlab="lag/25")
dev.off()

pdf("fig10_11.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

plot(density(OUT.25000[,2],adj=2),xlab=expression(beta[2]),
     ylab="posterior marginal density",main="")

plot(y~X[,2],xlab=expression(CO[2]),ylab="temperature")
abline(mean(OUT.25000[,1]),mean(OUT.25000[,2]),lwd=2)
abline(lmfit$coef,col="gray",lwd=2)
legend(180,2.5,legend=c("GLS estimate","OLS estimate"),bty="n",
       lwd=c(2,2),col=c("black","gray"))
dev.off()




quantile(OUT.25000[,2],probs=c(.025,.975) )



plot(X[,2],y,type="l")
points(X[,2],y,cex=2,pch=19)
points(X[,2],y,cex=1.9,pch=19,col="white")
text(X[,2],y,1:n)

iC<-solve( mean(OUT[,4])^DY )
Lev.gls<-solve(t(X)%*%iC%*%X)%*%t(X)%*%iC
Lev.ols<-solve(t(X)%*%X)%*%t(X)

plot(y,Lev.ols[2,] )
plot(y,Lev.gls[2,] )

# chapter 11
########
odat<-dget("../Data/nels_2002_data")
colnames(odat)<-c("sch_id","sch_enroll","sch_freelunch","sch_cnrtl",
                  "sch_urban","mteach_deg","eteach_deg","mteach_years","eteach_years" , 
                  "stu_sex","stu_lang","stu_pared","stu_income","stu_mathscore",
                  "stu_readscore","stu_mhw","stu_ehw","stu_readhours","stu_ses")
odat<-as.data.frame(odat)

ids<-dget("ids_selectschools")
group<-odat$sch_id
indset<-apply( odat[,1,drop=F],1,is.element,ids)
dat<-odat[indset,]

mathdat<-dat[,c(1,3,19,14)]
mathdat[,3]<-(mathdat[,3]-mean(mathdat[,3]))/sd(mathdat[,3]) 
########


########

groups<-ids
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-mathdat[mathdat[,1]==ids[j], 4] 
  N[j]<- sum(dat$sch_id==ids[j])
  xj<-mathdat[mathdat[,1]==ids[j], 3] 
  xj<-(xj-mean(xj))
  X[[j]]<-cbind( rep(1,N[j]), xj  )
}
#######

S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[[j]]~-1+X[[j]] )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 
####

#####
pdf("fig11_1.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot( range(mathdat[,3]),range(mathdat[,4]),type="n",xlab="SES", 
      ylab="math score")
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }

BETA.MLS<-apply(BETA.LS,2,mean)
abline(BETA.MLS[1],BETA.MLS[2],lwd=2)

plot(N,BETA.LS[,1],xlab="sample size",ylab="intercept")
abline(h= BETA.MLS[1],col="black",lwd=2)
plot(N,BETA.LS[,2],xlab="sample size",ylab="slope")
abline(h= BETA.MLS[2],col="black",lwd=2)

dev.off()
#####

if(2==3) {
  ##### hierarchical regression model
  p<-dim(X[[1]])[2]
  theta<-mu0<-apply(BETA.LS,2,mean)
  nu0<-1 ; s2<-s20<-mean(S2.LS)
  eta0<-p+2 ; Sigma<-S0<-L0<-cov(BETA.LS) ; BETA<-BETA.LS
  THETA.b<-S2.b<-NULL
  iL0<-solve(L0) ; iSigma<-solve(Sigma)
  source("~hoff/USBWork/rfunctions.r")
  Sigma.ps<-matrix(0,p,p)
  SIGMA.PS<-NULL
  BETA.ps<-BETA*0
  BETA.pp<-NULL
  set.seed(1)
  mu0[2]+c(-1.96,1.96)*sqrt(L0[2,2])
  for(s in 1:10000) {
    ##update beta_j 
    for(j in 1:m) 
    {  
      Vj<-solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
      Ej<-Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[[j]]/s2 )
      BETA[j,]<-rmvnorm(1,Ej,Vj) 
    } 
    ##
    
    ##update theta
    Lm<-  solve( iL0 +  m*iSigma )
    mum<- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum))
    theta<-t(rmvnorm(1,mum,Lm))
    ##
    
    ##update Sigma
    mtheta<-matrix(theta,m,p,byrow=TRUE)
    iSigma<-rwish( solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) )  ,  eta0+m) 
    ##
    
    ##update s2
    RSS<-0
    for(j in 1:m) { RSS<-RSS+sum( (Y[[j]]-X[[j]]%*%BETA[j,] )^2 ) }
    s2<-1/rgamma(1,(nu0+sum(N))/2, (nu0*s20+RSS)/2 )
    ##
    ##store results
    if(s%%10==0) 
    { 
      cat(s,s2,"\n")
      S2.b<-c(S2.b,s2);THETA.b<-rbind(THETA.b,t(theta))
      Sigma.ps<-Sigma.ps+solve(iSigma) ; BETA.ps<-BETA.ps+BETA
      SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
      BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
    }
    ##
  }
  
  save.image("data.f11_3")
}

load("data.f11_3")

#####
library(coda)
effectiveSize(S2.b)
effectiveSize(THETA.b[,1])
effectiveSize(THETA.b[,2])

apply(SIGMA.PS,2,effectiveSize)

tmp<-NULL;for(j in 1:dim(SIGMA.PS)[2]) { tmp<-c(tmp,acf(SIGMA.PS[,j])$acf[2]) }


acf(S2.b)
acf(THETA.b[,1])
acf(THETA.b[,2])
#####

pdf("fig11_3.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

plot(density(THETA.b[,2],adj=2),xlim=range(BETA.pp[,2]), 
     main="",xlab="slope parameter",ylab="posterior density",lwd=2)
lines(density(BETA.pp[,2],adj=2),col="gray",lwd=2)
legend( -3 ,1.0 ,legend=c( expression(theta[2]),expression(tilde(beta)[2])), 
        lwd=c(2,2),col=c("black","gray"),bty="n") 

quantile(THETA.b[,2],prob=c(.025,.5,.975))
mean(BETA.pp[,2]<0) 

BETA.PM<-BETA.ps/1000
plot( range(mathdat[,3]),range(mathdat[,4]),type="n",xlab="SES",
      ylab="math score")
for(j in 1:m) {    abline(BETA.PM[j,1],BETA.PM[j,2],col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,2]),lwd=2 )
dev.off()

########################


xs<-seq(5,100,5)/100
pops<-c("akr","b6","f1","mlh1","mom1","msh2","rb1","rb9")
DAT<-list()

for(j in 1:length(pops))
{
  err<-i<-0
  N<-NULL
  while(err==0)
  {
    x<-scan(paste("../Data/TumorData/",pops[j],".txt",sep=""),skip=i,nlines=1,quiet=T)
    if(length(x)>0)
    {
      pos<-c(1:length(x))/length(x)
      tm<-rep(pos,x)
      cnts<-xs*NA
      for(k in 1:length(cnts)){ cnts[k]<-sum( tm<=xs[k]) }
      cnts[2:length(cnts)]<-cnts[2:length(cnts)]-cnts[1:(length(cnts)-1)]
      N<-rbind(N,cnts)
      i<-i+2
    }
    if(length(x)==0){err<-1}
  }
  DAT[[j]]<-N
}


pdeg<-3
X<-cbind(rep(1,20),poly(xs,degree=pdeg) )

par(mfrow=c(4,2))
for(j in 1:8)
{
  plot(apply(DAT[[j]],2,mean),type="h",lwd=4)
  cat(j,dim(DAT[[j]])[1],mean(DAT[[j]]),"\n")
}



k<-3
m<-dim(DAT[[k]])[1] ; n<-dim(DAT[[k]])[2] ; p<-dim(X)[2]
Y<-DAT[[k]]

pdf("fig11_4.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(c(0,1),range(Y),type="n",xlab="location",ylab="number of tumors")
for(j in 1:m) { lines(xs,Y[j,],col="gray") }
lines( xs,apply(Y,2,mean),lwd=3)

lya<-log(apply(Y,2,mean))
Xs<-cbind( rep(1,n),poly(xs,deg=4,raw=TRUE))
fit2<- lm(lya~-1+Xs[,1:3] )
fit3<- lm(lya~-1+Xs[,1:4] )
fit4<- lm(lya~-1+Xs[,1:5] )

yh2<-Xs[,1:3]%*%fit2$coef
yh3<-Xs[,1:4]%*%fit3$coef
yh4<-Xs[,1:5]%*%fit4$coef

plot(xs,lya,type="l",lwd=3,xlab="location",ylab="log average number of tumors",
     ylim=range(c(lya,yh2,yh3,yh4)) )

points(xs,yh2,pch="2",col="black")
lines(xs,yh2,col="gray")
points(xs,yh3,pch="3",col="black")
lines(xs,yh3,col="gray")
points(xs,yh4,pch="4",col="black")
lines(xs,yh4,col="gray")
dev.off()




load("tumor.k4.S50000.rawTRUE")

library(coda)
round(apply(THETA.PS,2,effectiveSize),2)
round(apply(SIGMA.PS,2,effectiveSize),2)

### compare prior and posterior variance
SIGMA.PPS<-NULL
iS0<-solve(S0)
for(s in 1:5000) {
  tmp<-solve(rwish(iS0,eta0))
  SIGMA.PPS<-rbind(SIGMA.PPS,c(tmp))
}

par(mfrow=c(2,3))
for(j in c(1,7,13,19,25)) {
  plot(density(log(SIGMA.PS[,j]),adj=2) ,type="l",xlab=expression(paste("log ",Sigma[11])),ylab="density",main="",xlim=range(c(-7,log(SIGMA.PPS[,j]))),lwd=2 )
  lines( density(log(SIGMA.PPS[,j]),adj=2),col="gray",lwd=2)
  legend(-7.6,.72,legend=c("prior","posterior"),
         lwd=c(2,2),col=c("gray","black"),bty="n")
}
###


###
eXB.post<-NULL
for(s in 1:dim(THETA.PS)[1])
{
  beta<-rmvnorm(1,THETA.PS[s,],matrix(SIGMA.PS[s,],p,p))
  eXB.post<-rbind(eXB.post,t(exp(X%*%t(beta) )) )
}

qEB<-apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

eXT.post<- exp(t(X%*%t(THETA.PS )) )
qET<-apply( eXT.post,2,quantile,probs=c(.025,.5,.975))
yXT.pp<-matrix( rpois(prod(dim(eXB.post)),eXB.post),
                dim(eXB.post)[1],dim(eXB.post)[2] )

qYP<-apply( yXT.pp,2,quantile,probs=c(.025,.5,.975))


pdf("fig11_5.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot( c(0,1),range(c(0,qET,qEB,qYP)),type="n",xlab="location",
      ylab="number of tumors")
lines(xs, qET[1,],col="black",lwd=1)
lines(xs, qET[2,],col="black",lwd=2)
lines(xs, qET[3,],col="black",lwd=1)

plot( c(0,1),range(c(0,qET,qEB,qYP)),type="n",xlab="location",
      ylab="")
lines(xs, qEB[1,],col="black",lwd=1)
lines(xs, qEB[2,],col="black",lwd=2)
lines(xs, qEB[3,],col="black",lwd=1)

plot( c(0,1),range(c(0,qET,qEB,qYP)),type="n",xlab="location",
      ylab="")
lines(xs, qYP[1,],col="black",lwd=1)
lines(xs, qYP[2,],col="black",lwd=2)
lines(xs, qYP[3,],col="black",lwd=1)

dev.off()


# chapter 12
dat<-read.table("http://lib.stat.cmu.edu/aoas/107/data.txt",header=TRUE)
####


yincc<-match(dat$INC,sort(unique(dat$INC)))
ydegr<-dat$DEGREE+1
yage<-dat$AGE
ychild<-dat$CHILD
ypdeg<-1*(dat$PDEG>2)
tmp<-lm(ydegr~ychild+ypdeg+ychild:ypdeg)


#####
pdf("fig12_1.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(table(dat$DEG+1)/sum(table(dat$DEG+1)),
     lwd=2,type="h",xlab="DEG",ylab="probability")
plot(table(dat$CHILD)/sum(table(dat$CHILD)),lwd=2,type="h",xlab="CHILD",ylab="probability" )
dev.off()
#####

#####

X<-cbind(ychild,ypdeg,ychild*ypdeg)
y<-ydegr
keep<- (1:length(y))[ !is.na( apply( cbind(X,y),1,mean) ) ]
X<-X[keep,] ; y<-y[keep]
ranks<-match(y,sort(unique(y))) ; uranks<-sort(unique(ranks))
n<-dim(X)[1] ; p<-dim(X)[2]
iXX<-solve(t(X)%*%X)  ; V<-iXX*(n/(n+1)) ; cholV<-chol(V)

###starting values
set.seed(1)
beta<-rep(0,p) 
z<-qnorm(rank(y,ties.method="random")/(n+1))
g<-rep(NA,length(uranks)-1)
K<-length(uranks)
BETA<-matrix(NA,1000,p) ; Z<-matrix(NA,1000,n) ; ac<-0
mu<-rep(0,K-1) ; sigma<-rep(100,K-1)
S<-25000
for(s in 1:S) 
{
  
  #update g 
  for(k in 1:(K-1)) 
  {
    a<-max(z[y==k])
    b<-min(z[y==k+1])
    u<-runif(1, pnorm( (a-mu[k])/sigma[k] ),
             pnorm( (b-mu[k])/sigma[k] ) )
    g[k]<- mu[k] + sigma[k]*qnorm(u)
  }
  
  #update beta
  E<- V%*%( t(X)%*%z )
  beta<- cholV%*%rnorm(p) + E
  
  #update z
  ez<-X%*%beta
  a<-c(-Inf,g)[ match( y-1, 0:K) ]
  b<-c(g,Inf)[y]  
  u<-runif(n, pnorm(a-ez),pnorm(b-ez) )
  z<- ez + qnorm(u)
  
  
  #help mixing
  c<-rnorm(1,0,n^(-1/3))  
  zp<-z+c ; gp<-g+c
  lhr<-  sum(dnorm(zp,ez,1,log=T) - dnorm(z,ez,1,log=T) ) + 
    sum(dnorm(gp,mu,sigma,log=T) - dnorm(g,mu,sigma,log=T) )
  if(log(runif(1))<lhr) { z<-zp ; g<-gp ; ac<-ac+1 }
  
  if(s%%(S/1000)==0) 
  { 
    cat(s/S,ac/s,"\n")
    BETA[s/(S/1000),]<-  beta
    Z[s/(S/1000),]<- z
  }
  
} 
#####
pdf("fig12_2.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(X[,1]+.25*(X[,2]),Z[1000,],
     pch=15+X[,2],col=c("gray","black")[X[,2]+1],
     xlab="number of children",ylab="z", ylim=range(c(-2.5,4,Z[1000,])),
     xlim=c(0,9))

beta.pm<-apply(BETA,2,mean)
ZPM<-apply(Z,2,mean)
abline(0,beta.pm[1],lwd=2 ,col="gray")
abline(beta.pm[2],beta.pm[1]+beta.pm[3],col="black",lwd=2 )
#legend(3.75,4.25,legend=c("parents without college","parents with college"),pch=c(15,16),col=c("gray","black"))
legend(5,4,legend=c("PDEG=0","PDEG=1"),pch=c(15,16),col=c("gray","black"))


plot(density(BETA[,3],adj=2),lwd=2,xlim=c(-.5,.5),main="",
     xlab=expression(beta[3]),ylab="density")
sd<-sqrt(  solve(t(X)%*%X/n)[3,3] )
x<-seq(-.7,.7,length=100)
lines(x,dnorm(x,0,sd),lwd=2,col="gray")
legend(-.5,6.5,legend=c("prior","posterior"),lwd=c(2,2),col=c("gray","black"),bty="n")
dev.off()
#####

beta.pm<-apply(BETA,2,mean)
beta.pm[1]+beta.pm[3]
quantile(BETA[,3],prob=c(.025,.0975))
quantile(BETA[,3],prob=c(0.025,0.975))

source("treg.r")
rfit<-treg(y,X)

pdf("fig12_3.pdf",family="Times",height=1.75,width=5) 
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
lab<-c(expression(beta[1]),expression(beta[2]),expression(beta[3]))
ymx<-c(12,3.1,7.25)
laby<-c("density","","")
for(j in 1:3) {
  plot(density(rfit$BETA[,j],adj=2),lwd=2,main="",
       xlab=lab[j],col="black",ylim=c(0,ymx[j]),ylab=laby[j])
  lines(density(BETA[,j],adj=2),col="gray",lwd=2)
  if(j==4) {
    legend(-.2,12,legend=c("ordered probit","rank likelihood"), 
           lwd=c(2,2),col=c("gray","black"))
  } 
}
dev.off()
#################





