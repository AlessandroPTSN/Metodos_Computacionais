
#Pacotes utilizados:
require(flexsurv)
library(doParallel)

#Fun??o para gera??o de n?meros aleat?rios da densidade da distribui??o gompertz(a,b)
rgomp=function(n,a,b){
  u=runif(n)
  log(1 -a/b*log(1-u))/a #Inversa da função de distribuição acumulada
}


# Fun??o de densidade de probabilidade. (Essa a gente usa no histograma)
pfun1 <- function(x,a,b){
  f <- (b*(exp(1)^(a*x)))*exp((-b/a)*((exp(1)^(a*x))-1)) #Densidade
  return(f)
}


# Gerando valores aleat?rios da fun??o rgomp que criamos e da fun??o rgompertz do pacote flexsurv para comparar a aproxima??o via histograma.
set.seed(1)
n=5000
a=1
b=1
dados=rgomp(n,a,b)
dados2=rgompertz(n,a,b)
x <- seq(0, max(dados), 0.01)

# Comparando as aproxima??es via histograma

par(mfrow=c(1,2))

hist(dados2,freq = FALSE,main="",xlab="x",ylab="Density-fun",lwd=2)
lines(x,pfun1(x,a,b),lty=1,lwd=2)
box()

#Histograma dos dados obtidos pela fun??o densidade.
hist(dados,freq = FALSE,main="",xlab="x",ylab="Density-uni",lwd=2)
lines(x,pfun1(x,a,b),lty=1,lwd=2)
box()

estimativas=optim(par=chute,fn=logdgomp,dados=dados)$par
(ahat=estimativas[1]) #Estimador pontual para a
(bhat=estimativas[2]) #Estimador pontual para b

# Calculando termos que s?o utilizados no c?lculo da vari?ncia da distribui??o gompertz.
C=n/bhat^2
B=-sum(exp(ahat*dados)*(1-ahat*dados)-1)/ahat^2
A=bhat/ahat^3*(2*(sum(exp(ahat*dados)*(1-ahat*dados))-1) -
                 ahat*sum(exp(ahat*dados)*dados*((1-ahat*dados)^2-1))
)
# C?lculo da vari?ncia e do intervalo de confian?a para os estimadores da distribui??o gompertz a partir dos dados que foram gerados.
varahat=C/(A*C-B^2)
varbhat=A/(A*C-B^2)
covabhat=-B/(A*C-B^2)
(ICa=ahat+c(-1,1)*1.96*sqrt(varahat))#Intervalo de confiança para a
(ICb=bhat+c(-1,1)*1.96*sqrt(varbhat))


# Fun??o log-verossimilhan?a da Gompertz:

logdgomp=function(x,dados){
  a=x[1]
  b=x[2]
  -sum(log(b)+a*dados-b/a*(exp(a*dados)-1))
}


achute=1 #Chute inicial para o par?metro a
bchute=1 #Chute inicial para o par?metro b
(chute=c(achute,bchute))


#Aloca??o de mem?ria para otimizar os la?os do monte carlo e bootstrap
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

#Par?metros 

n = 25
R = 5000
B = 500
a = 2
b = 3
dados=rgomp(n,a,b)
set.seed(10241226)

# Monte Carlo

ahat<-rep(NA,R)
bhat<-rep(NA,R)

foreach(i=1:R) %do%{
 Z = rgomp(n,a,b)
 
 #Estima??o pelo m?todo de m?xima verossimilhan?a, sei que precisamos utilizar um m?todo num?rico para encontrar esse estimador. Encontrei algo interessante no portal action: http://www.portalaction.com.br/inferencia/341-metodos-numericos-dos-estimadores-de-maxima-verossimilhanca
 
 estimados_MV #<- optim(par=c(a, b), fn = logdgomp, x=Z, method="BFGS")
 
 
 ahat_MV[k] <- estimados_MV$par[1]
 bhat_MV[k]      <- estimados_MV$par[2]
}

#Fechando a aloca??o de mem?ria
stopCluster(cl)

#Estima??o pontual pelo m?todo de monte carlo
estimate_a <- mean(ahat_MV)
estimate_b <- mean(bhat_MV)

estimate_a
estimate_b

#C?lculo do Vi?s e do Erro Quadr?tico M?dio (EQM)
vies_ahat <- mean(ahat_MV) - a
vies_bhat <- mean(bhat_MV) - b

vies_ahat
vies_bhat 


EQM_a <- vies_ahat^2 + var(ahat_MV)
EQM_b <- vies_bhat^2 + var(bhat_MV)

EQM_a
EQM_b 

#________________//________________//________________//________________//_________

