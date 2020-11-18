u =  1-exp(-b/a(e^(ax)-1))
exp(-b/a(e^(ax)-1))=1-u
-b/a(e^(ax)-1)= log(1-u)
e^(ax)-1 = -a/b*log(1-u)
e^(ax)= 1 -a/b*log(1-u)
x= log(1 -a/b*log(1-u))/a


#############################
#  Histograma dos valores gerados #
#############################

# Códigos - Alessandro

library(flexsurv)
set.seed(123) #vai fornecer os mesmos números usados
rfun1 <- function(n,a,b){
  x <- runif(n)
  Z = log((a*((b/a)-log(1-x)))/b)/a  #Inversa da acumulada
  return(Z)
}
pfun1 <- function(x,a,b){
  f <- (b*(exp(1)^(a*x)))*exp((-b/a)*((exp(1)^(a*x))-1)) #Densidade
  return(f)
}


n=200
a=2
b=3
Z <- rfun1(n,a,b)
ZR = rgompertz(200, shape = 2, rate = 3)
x <- seq(0, max(Z), 0.01)


par(mfrow=c(1,2))
hist(ZR,freq = FALSE,main="",xlab="x",ylab="Density-fun",lwd=2)
lines(x,pfun1(x,a,b),lty=1,lwd=2)
box()
hist(Z,freq = FALSE,main="",xlab="x",ylab="Density-uni",lwd=2)
lines(x,pfun1(x,a,b),lty=1,lwd=2)
box()



#######################################################################
#  Código da estimação dos parâmetros a e b pelo método de máxima verossimilhança #
#######################################################################

# Códigos - Ivo
library(flexsurv)

rgomp=function(n,a,b){
  u=runif(n)
  log(1 -a/b*log(1-u))/a
}
n=100
a=1
b=1
dados=rgomp(n,a,b)
#dados2=rgompertz(100,3,5)
#par(mfrow=c(1,2))
#hist(dados)
#hist(dados2)

logdgomp=function(ab,x){
  a=ab[1]
  b=ab[2]
  -sum(log(b)+a*x-b/a*(exp(a*x)-1))
}

az=abs(log(mean(dados))/(mean(dados)*log(sd(dados))))
#achute=sum(1-exp(-az*dados))*(1-az*mean(dados))/sum(dados*exp(-az*dados))
achute=az
bchute=n*achute/(sum(1-exp(-achute*dados)))
chute=c(achute,bchute)

chute
estimativas=optim(par=chute,fn=logdgomp,x=dados)$par
ahat=estimativas[1]
bhat=estimativas[2]
A=n/bhat^2
B=mean(dados)*sum(dados*exp(-ahat*dados))/(1-ahat*mean(dados))
C=2*bhat*mean(dados)*sum(dados*exp(-ahat*dados))/(ahat*(1-ahat*mean(dados)))-bhat/ahat*sum(dados^2*exp(-ahat*dados))
varahat=A/(A*C-B^2)
varbhat=C/(A*C-B^2)
covabhat=-B/(A*C-B^2)
ICa=ahat+c(-1,1)*1.96*varahat
ICb=bhat+c(-1,1)*1.96*varbhat

matriza=matrix(NA,nrow=7,ncol=7)
matrizb=matrix(NA,nrow=7,ncol=7)
aa=c(0.01,0.1,0.5,1,3,5,10)
bb=aa
rownames(matriza)=aa
colnames(matriza)=bb
dimnames(matrizb)=dimnames(matriza)
for(i in 1:7){
  for(j in 1:7){
    a=aa[i]
    b=bb[j]
    n=1000
    dados=rgomp(n,a,b)
    az=abs(log(mean(dados))/(mean(dados)*log(sd(dados))))
    achute=az
    matriza[i,j]=achute
    matrizb[i,j]=n*achute/(sum(1-exp(-achute*dados)))
  }
}

round(matriza,2)
round(matrizb,2)

