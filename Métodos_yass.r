
#Pacotes utilizados:
require(flexsurv)
library(doParallel)
#rm(list=ls()) #AVISO, VAI LIMPAR TODAS AS VARIAVEIS SALVAS
##################################################################################

#                           Funções usadas no script

##################################################################################

#Função para geração de números aleatórios da densidade da distribuição gompertz(a,b)
rgomp=function(n,a,b){
  u=runif(n)
  log(1 -a/b*log(1-u))/a #Inversa da função de distribuição acumulada
}


# Função de densidade de probabilidade. (Essa a gente usa no histograma)
pfun1 <- function(x,a,b){
  f <- (b*(exp(1)^(a*x)))*exp((-b/a)*((exp(1)^(a*x))-1)) #Densidade
  return(f)
}

# Função log-verossimilhança da Gompertz:

logdgomp=function(x,dados){
  a=x[1]
  b=x[2]
  -sum(log(b)+a*dados-b/a*(exp(a*dados)-1))
}

##################################################################################


# Gerando valores aleatórios da função rgomp que criamos e da função rgompertz do pacote flexsurv para comparar a aproximação via histograma.
set.seed(1)
n=5000
a=1
b=4
dados=rgomp(n,a,b)
dados2=rgompertz(n,a,b)
achute=1 #Chute inicial para o par?metro a
bchute=1 #Chute inicial para o par?metro b
(chute=c(achute,bchute))
x <- seq(0, max(dados), 0.01)

# Comparando as aproximaçães via histograma

par(mfrow=c(1,2))

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

# Calculando termos que s?o utilizados no c?lculo da vari?ncia da distribuição gompertz.
C=n/bhat^2
B=-sum(exp(ahat*dados)*(1-ahat*dados)-1)/ahat^2
A=bhat/ahat^3*(2*(sum(exp(ahat*dados)*(1-ahat*dados))-1) -
                 ahat*sum(exp(ahat*dados)*dados*((1-ahat*dados)^2-1))
)
# Cálculo da variância e do intervalo de confiança para os estimadores da distribuição gompertz a partir dos dados que foram gerados.
varahat=C/(A*C-B^2)
varbhat=A/(A*C-B^2)
covabhat=-B/(A*C-B^2)
(ICa=ahat+c(-1,1)*1.96*sqrt(varahat))#Intervalo de confiança para a
(ICb=bhat+c(-1,1)*1.96*sqrt(varbhat))


#_________________________________________________________________________________
#
# Estimação de a e b pelo método da máxima verossimilhança utilizando Monte Carlo
#_________________________________________________________________________________


#Parâmetros 
n = 25
R = 5000
B = 500
a = 2
b = 3
achute=2 #Chute inicial para o parâmetro a
bchute=3 #Chute inicial para o parâmetro b
chute=c(achute,bchute)
set.seed(1024122)


#Alocação de memória
ahat_MV<-rep(NA,R)
bhat_MV<-rep(NA,R)


#Alocação de clusteres para otimizar os laços do monte carlo e bootstrap
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  


foreach(k=1:R) %do%{
  Z = rgomp(n,a,b)
  
  # Utilizando o método numérico "BFGS" para estimar os parâmtros a e b da distribuição gompertz pelo método da máxima verossimilhança
  estimados_MV <- optim(par=chute, method="BFGS", fn = logdgomp, dados=Z)
  
  
  ahat_MV[k] <- estimados_MV$par[1]
  bhat_MV[k]      <- estimados_MV$par[2]
  
}

#Fechando a alocação de clusteres
stopCluster(cl)

#Estimação pontual pelo m?todo de monte carlo
estimate_a <- mean(ahat_MV)
estimate_b <- mean(bhat_MV)

estimate_a
estimate_b

#Cálculo do Viés e do Erro Quadrático Médio (EQM)
vies_ahat <- mean(ahat_MV) - a
vies_bhat <- mean(bhat_MV) - b

vies_ahat
vies_bhat 


EQM_a <- vies_ahat^2 + var(ahat_MV)
EQM_b <- vies_bhat^2 + var(bhat_MV)

EQM_a
EQM_b 

# IC 95%

LimInf_a_MV = quantile(ahat_MV,0.025)
LimSup_a_MV = quantile(ahat_MV,0.975)


LimInf_b_MV = quantile(bhat_MV,0.025)
LimSup_b_MV = quantile(bhat_MV,0.975)

#_________________________________________________________________________________
#
#Estimação de a e b pelo método da máxima verossimilhança utilizando Monte Carlo e Bootstrap
#_________________________________________________________________________________



#Alocação de memória

ahat_MV_boot<-NA
bhat_MV_boot<-NA
ahat_boot<-rep(NA,B)
bhat_boot<-rep(NA,B)
Z = rgomp(n,a,b)
estimados_MV <- optim(par=chute, method="BFGS", fn = logdgomp, dados=Z)
ahat_MV_boot <- estimados_MV$par[1]
bhat_MV_boot <- estimados_MV$par[2]


#Alocação de clusteres para otimizar os laços do monte carlo e bootstrap
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  


  
  # Utilizando o método numérico "BFGS" para estimar os parâmtros a e b da distribuição gompertz pelo método da máxima verossimilhança

  # Criando o algoritmo do bootstrap paramétrico
  foreach(i=1:B) %do% {
    Z_boot = rgomp(n,ahat_MV_boot,bhat_MV_boot)
    estimados_boot <- optim(par=c(ahat_MV_boot,bhat_MV_boot), method="BFGS", fn = logdgomp, dados=Z_boot)
    
    ahat_boot[i] <- estimados_boot$par[1]
    bhat_boot[i] <- estimados_boot$par[2]
  } 

#Fechando a alocação de clusteres
stopCluster(cl)

#Estimação pontual pelo método de monte carlo e bootstrap
estimate_a_boot <- mean(ahat_boot)
estimate_b_boot <- mean(bhat_boot)

estimate_a_boot
estimate_b_boot

#Cálculo do Viés e do Erro Quadrático Médio (EQM)
vies_ahat_boot <- mean(ahat_boot) - a
vies_bhat_boot <- mean(bhat_boot) - b

vies_ahat_boot
vies_bhat_boot 

a_corrigido <- 2*ahat_MV_boot - estimate_a_boot
b_corrigido <- 2*bhat_MV_boot - estimate_b_boot

a_corrigido
b_corrigido

EQM_a_boot <- vies_ahat^2 + var(ahat_boot)
EQM_b_boot <- vies_bhat^2 + var(bhat_boot)

EQM_a_boot
EQM_b_boot 

# IC 95%

LimInf_a_boot = quantile(ahat_boot,0.025)
LimSup_a_boot = quantile(ahat_boot,0.975) 

LimInf_b_boot = quantile(bhat_boot,0.025)
LimSup_b_boot = quantile(bhat_boot,0.975) 


#_________//_________//_________//_________//_________//_________//_________//____


# Criando tabela que apresenta os resultados:

Método<-rep(c("Máxima Verossimilhança","Bootstrap","Corrigido"),each=2)

# Monte carlo
Estimativa_mv<-as.matrix(c(estimate_a,estimate_b))

Vies_MV<-as.matrix(c(vies_ahat,vies_bhat ))

EQM_MV<-as.matrix(c(EQM_a,EQM_b))


#Bootstrap
Estimativa_boot<-as.matrix(c(estimate_a_boot,estimate_b_boot))

Vies_boot<-as.matrix(c(vies_ahat_boot,vies_bhat_boot))

EQM_boot<-as.matrix(c(EQM_a_boot,EQM_b_boot))

# Corrigido

Estimativa_cor <- as.matrix(c(a_corrigido,b_corrigido))

Vies_cor <- as.matrix(c(a_corrigido-a,b_corrigido-b))

EQM_Corrigido <- as.matrix(c(NA,NA))

# Criando a tabela
Parâmetro<-rep(c("a","b"),times=3)

Estimativa=rbind(Estimativa_mv,Estimativa_boot,Estimativa_cor)
Vies<-rbind(Vies_MV,Vies_boot,Vies_cor)
EQM<-rbind(EQM_MV,EQM_boot,EQM_Corrigido)
ICinf=rbind(LimInf_a_MV,LimInf_b_MV,LimInf_a_boot,LimInf_b_boot,NA,NA)
ICsup=rbind(LimSup_a_MV,LimSup_b_MV,LimSup_a_boot,LimSup_b_boot,NA,NA)
nrep=rep(n,6)
tabela<-data.frame(nrep,Método,Parâmetro,Estimativa,Vies,EQM,ICinf,ICsup)
row.names(tabela)=1:6
tabela




# Save/Load
save.image(file = "Ivo.RData")#Salvando a img do R
#rm(list=ls()) #AVISO, VAI LIMPAR TODAS AS VARIAVEIS SALVAS

#Carregando todas as tabelas já criadas
load("Ivo.RData")#Carregando a img do R
tabela25<-tabela
load("leo.RData")#Carregando a img do R
tabela50<-tabela
load("Fabricio.RData")#Carregando a img do R
tabela100<-tabela
load("Alessandro.RData")#Carregando a img do R
tabela200<-tabela
