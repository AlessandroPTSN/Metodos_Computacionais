#Pacotes utilizados:
require(flexsurv)
library(doParallel)
rm(list=ls()) #AVISO, VAI LIMPAR TODAS AS VARIAVEIS SALVAS
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
n = 100
R = 5000
B = 500
a = 2
b = 3
chute=c(a,b)
set.seed(1024122)
#Como estamos utilizando achute = a = 2 e bchute = b = 3 então eu coloquei que chute = c(a,b) só para facilitar nas alterações
#achute=2 #Chute inicial para o parâmetro a
#bchute=3 #Chute inicial para o parâmetro b

#Alocação de memória
ahat_MV<-rep(NA,R)
bhat_MV<-rep(NA,R)
ahat_boot<-matrix(NA,nrow=R,ncol=B)
bhat_boot<-matrix(NA,nrow=R,ncol=B)
ahat_MV_boot<-rep(NA,R)
bhat_MV_boot<-rep(NA,R)

#Alocação de clusteres para otimizar os laços do monte carlo e bootstrap
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  


foreach(k=1:R)%do%{
  Z = rgomp(n,a,b)
  
  # Utilizando o método numérico "BFGS" para estimar os parâmtros a e b da distribuição gompertz pelo método da máxima verossimilhança
  estimados_MV <- optim(par=chute, method="BFGS", fn = logdgomp, dados=Z)
  ahat_MV[k] <- estimados_MV$par[1]
  bhat_MV[k]      <- estimados_MV$par[2]
  
  
  foreach(i=1:B)%do%{
    Z_boot = rgomp(n,ahat_MV[k],bhat_MV[k])
    estimados_boot <- optim(par=c(ahat_MV[k],bhat_MV[k]), method="BFGS", fn = logdgomp, dados=Z_boot)
    
    ahat_boot[k,i] <- estimados_boot$par[1]
    bhat_boot[k,i] <- estimados_boot$par[2]
  }
  
  ahat_MV_boot[k] <- mean(ahat_boot[k,],na.rm=T)
  bhat_MV_boot[k] <- mean(bhat_boot[k,],na.rm=T)
  
  print(paste0(round(((R)/5000)*100,1),"% concluído"))
}

#Fechando a alocação de clusteres
stopCluster(cl)

#Estimação pontual pelo método de monte carlo
estimate_a <- mean(ahat_MV)
estimate_b <- mean(bhat_MV)

estimate_a
estimate_b

#Cálculo do Viés e do Erro Quadrático Médio (EQM)
vies_ahat <- mean(ahat_MV) - a
vies_bhat <- mean(bhat_MV) - b

vies_ahat
vies_bhat 

var_ahat<-var(ahat_MV)
var_bhat<-var(bhat_MV)

EQM_a <- vies_ahat^2 + var_ahat
EQM_b <- vies_bhat^2 + var_bhat

EQM_a
EQM_b 

# IC 95%

LimInf_a_MV = quantile(ahat_MV,0.025)
LimSup_a_MV = quantile(ahat_MV,0.975)


LimInf_b_MV = quantile(bhat_MV,0.025)
LimSup_b_MV = quantile(bhat_MV,0.975)



#Cálculo do Viés e do Erro Quadrático Médio (EQM)


a_corrigido <- 2*ahat_MV - ahat_MV_boot
b_corrigido <- 2*bhat_MV - bhat_MV_boot

a_corrigido
b_corrigido

estimate_a_cor = mean(a_corrigido,na.rm=T)
estimate_b_cor = mean(b_corrigido,na.rm=T)

vies_a_cor <- estimate_a_cor - a
vies_b_cor <- estimate_b_cor - b

var_a_cor<-var(a_corrigido,na.rm=T)
var_b_cor<-var(b_corrigido,na.rm=T)

EQM_a_cor <- vies_a_cor^2 + var_a_cor
EQM_b_cor <- vies_b_cor^2 + var_b_cor


# IC 95%

LimInf_a_cor = quantile(a_corrigido,0.025,na.rm=T)
LimSup_a_cor = quantile(a_corrigido,0.975,na.rm=T) 

LimInf_b_cor = quantile(b_corrigido,0.025,na.rm=T)
LimSup_b_cor = quantile(b_corrigido,0.975,na.rm=T) 

#_________//_________//_________//_________//_________//_________//_________//____


# Criando a tabela que apresenta os resultados:

Método<-rep(c("Máxima Verossimilhança","Corrigido"),each=2)

# Monte carlo
Estimativa_mv<-as.matrix(c(estimate_a,estimate_b))

Vies_MV<-as.matrix(c(vies_ahat,vies_bhat))

Var_MV<-as.matrix(c(var_ahat,var_bhat))

EQM_MV<-as.matrix(c(EQM_a,EQM_b))

# Corrigido

Estimativa_cor <- as.matrix(c(estimate_a_cor,estimate_b_cor))

Vies_cor <- as.matrix(c(vies_a_cor,vies_b_cor))

Var_cor<-as.matrix(c(var_a_cor,var_b_cor))

EQM_Corrigido <- as.matrix(c(EQM_a_cor,EQM_b_cor))

# Criando a tabela
Parâmetro<-rep(c("a","b"),times=2)
Valor_do_parametro=rep(c(a,b),times=2)
Estimativa=rbind(Estimativa_mv,Estimativa_cor)
Vies<-rbind(Vies_MV,Vies_cor)
Var<-rbind(Var_MV,Var_cor)
EQM<-rbind(EQM_MV,EQM_Corrigido)
ICinf=rbind(LimInf_a_MV,LimInf_b_MV,LimInf_a_cor,LimInf_b_cor)
ICsup=rbind(LimSup_a_MV,LimSup_b_MV,LimSup_a_cor,LimSup_b_cor)
nrep=rep(n,4)
tabela<-data.frame(nrep,Método,Parâmetro,Valor_do_parametro,Estimativa,Vies,Var,EQM,ICinf,ICsup)


row.names(tabela)=1:4
names(tabela)[c(4,7,9,10)]<-c("Valor do parâmetro","Variância","Lim. Inf. 2.5%","Lim. Sup. 97.5%")


tabela



# SALVANDO O TRABALHO REALIZADO. É IMPORTANTE SEGUIR O PASSO A PASSO ABAIXO PARA SALVAR ADEQUADAMENTE SEU TRABALHO.

#     PASSO 1- Ao rodar o algoritmo com a=2 e b=3, rode o primeiro comando write.csv2 que se encontra abaixo para salvar a tabela no seu diretório.
#     2- Retorne lá nos parâmetros e coloque a=0.5 e b=0.5 e repita o processo de criação da tabela. Feito isso, rode o segundo comando write.csv2 para salvar essa segunda tabela no seu diretório
#     3- após os passos 1 e 2 vá para o passo 3 onde você unirá as duas tabelas em uma única e será a que utilizaremos no trabalho.



# PASSO 1
#write.csv2(tabela,"Tabela100_1.csv",row.names = F)

# PASSO 2
#write.csv2(tabela,"Tabela100_2.csv",row.names = F)


# PASSO 3
#carregando as tabelas e juntando-as
tab1<-read.csv2("Tabela100_1.csv")
tab2<-read.csv2("Tabela100_2.csv")
tabela.f<-rbind(tab1,tab2)
colnames(tabela.f)=c("nrep","Método","Parâmetro","Valor do parâmetro","Estimativa","Vies","Variância","EQM","Lim. Inf. 2.5%","Lim. Sup. 97.5%")

# Escrevendo a tabela final:
#write.csv2(tabela.f,"Tabela100.csv",row.names = F)

# Lendo a tabela final: (Opcional só para ver se deu td certo)
#tabela.f<-read.csv2("Tabela100.csv")


#______________\\______________\\______________\\______________\\___

