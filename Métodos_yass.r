rgomp=function(n,a,b){
  u=runif(n)
  log(1 -a/b*log(1-u))/a #Inversa da função de distribuição acumulada
}

pfun1 <- function(x,a,b){
  f <- (b*(exp(1)^(a*x)))*exp((-b/a)*((exp(1)^(a*x))-1)) #Densidade
  return(f)
}


n=5000
a=1
b=1
dados=rgomp(n,a,b)
dados2=rgompertz(n,a,b)
x <- seq(0, max(dados), 0.01)

par(mfrow=c(1,2))

#Histograma dos dados gerados.
hist(dados2,freq = FALSE,main="",xlab="x",ylab="Density-fun",lwd=2)
lines(x,pfun1(x,a,b),lty=1,lwd=2)
box()

#Histograma dos dados obtidos pela função densidade.
hist(dados,freq = FALSE,main="",xlab="x",ylab="Density-uni",lwd=2)
lines(x,pfun1(x,a,b),lty=1,lwd=2)
box()

estimativas=optim(par=chute,fn=logdgomp,dados=dados)$par
(ahat=estimativas[1]) #Estimador pontual para a
(bhat=estimativas[2]) #Estimador pontual para b

C=n/bhat^2
B=-sum(exp(ahat*dados)*(1-ahat*dados)-1)/ahat^2
A=bhat/ahat^3*(2*(sum(exp(ahat*dados)*(1-ahat*dados))-1) -
                 ahat*sum(exp(ahat*dados)*dados*((1-ahat*dados)^2-1))
)
varahat=C/(A*C-B^2)
varbhat=A/(A*C-B^2)
covabhat=-B/(A*C-B^2)
(ICa=ahat+c(-1,1)*1.96*sqrt(varahat))#Intervalo de confiança para a
(ICb=bhat+c(-1,1)*1.96*sqrt(varbhat))



logdgomp=function(x,dados){
  a=x[1]
  b=x[2]
  -sum(log(b)+a*dados-b/a*(exp(a*dados)-1))
}


achute=1 #Chute inicial para o parâmetro a
bchute=1 #Chute inicial para o parâmetro b
(chute=c(achute,bchute))

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

n = 25
R = 5000
B = 500
a = 2
b = 3
dados=rgomp(n,a,b)
set.seed(10241226)

# Monte Carlo
foreach(i=1:R) %do%{
 Z = rgomp(n,a,b)
}

stopCluster(cl)

