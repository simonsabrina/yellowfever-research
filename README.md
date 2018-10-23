# yellowfever-research
Exploratory and ongoing analysis on Yellow Fever pattern in Brazil - PhD Thesis/Dissertation

---
title: "Modelagem da Febre Amarela"
author: "S. Simon"
date: "October 19, 2018"
output:
  pdf_document: default
  html_document:
    df_print: paged
theme: cerulean
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## INTRODUÇÃO

<font size="4">O documento a seguir apresenta uma análise dos dois eventos de epizootia de Febre Amarela ocorridos no Brasil entre 1999 e 2018. Após o tratamento dos dados, o R0 dos dois eventos é estimado a partir do ajuste de uma curva exponencial aos dados. Em seguida, é executada uma simulação do modelo SIR-SEI para o cenário ocorrido em 2016-2018 no Espírito Santo, utilizando os parâmetros para Febre Amarela em humanos propostos por Massad et al, 2018.

### Estimando o R0 para as epizootias de Febre Amarela de 2008 e 2017

Para estimar o R0, serão utilizados todos os registros de óbitos de pnh por febre amarela no Brasil, de 1999 a 2018, sem discriminar espécies ou estados onde ocorreram as epizootias.

```{r}
# lendo os dados:
dados <- read.csv("pnh_1999-2018.csv")
names(dados)
fa <- data.frame(dados$data, dados$municipio, dados$uf, dados$n.animais)
names(fa) <- c("data","municipio","uf","casos")
```

Como os óbitos foram registrados por data, vamos converter as datas em número de dias.

```{r, message=FALSE}
library(lubridate)
library(dplyr)
```

```{r}
# Convertendo data em vetor de valores inteiros (= contagem de dias):
fa$data <- ymd(fa$data)
fa <- arrange(fa,data) 

inicio <- fa[1,1]
nextday <- as.integer(fa[,1] - inicio)
dia <- nextday+1
fa0 <- data.frame(fa,dia)

# Nesse novo df, agrupar o número de registros por dia:
fa0 <- group_by(fa0, dia)
curve <- summarise(fa0,sum(casos))
names(curve) <- c("dia", "casos")

plot(y = curve$casos, x = curve$dia, type="l", main = "obitos de PNH entre 1999 e 2018", xlab = "Dia", ylab = "Numero de casos")
```


### Epidemia 2008 - 2009

Para conhecer o R0 é necessário selecionar um periodo na fase de crescimento, separando dados da subida da curva para o ajuste.

```{r}
# Selecionando dados
summary(curve) # pico = 153 casos no dia 3518.
curve3 <- curve[curve$dia>=3460 & curve$dia<=3515,]
names(curve3) <- c("dia","casos")
t <- c(curve3$dia-3460)
curve3 <- data.frame(curve3, t)
head(curve3,2)

plot(curve3$t, curve3$casos, type= "l", ylab = "Numero de casos por dia", xlab = "Dia", 
     main = "Numero de obitos de PNH por FA entre 2008 e 2009 no Brasil",
     sub = "Fase de ascencao da epidemia")
```

O registro de óbitos apresenta um padrão variável que não é característico de uma curva exponencial. Vamos agrupar os registros por intervalos de tempo semanal, sendo a semana 1 aquela em que começou a epizootia (não corresponde necessariamente à semanas epidemiológicas).

```{r}
sem <- c((curve3$t %/% 7)+1) # x %/% y = Parte inteira da divisão de x por y
curve4 <- data.frame(curve3$casos, sem)
curve4 <- group_by(curve4, sem)
curve4 <- summarise(curve4, sum(curve3.casos))
names(curve4) <- c("sem", "casos")
```

O ajuste da subida da curva é feito utilizando uma função exponencial Ch x exp(k x tempo) onde Ch e k sao parametros de ajuste. O valor inicial de Ch é um chute para o numero de casos inicial.

```{r}
ajuFA1 <- nls(casos ~ Ch*exp(k*sem), data = curve4, 
             start = list(Ch=1, k=1), # chute inicial
             trace = TRUE)
plot(curve4$sem, curve4$casos, ylab = "Numero de casos por semana", xlab = "Semana", 
     main = "Numero de obitos de PNH por FA entre 2008 e 2009 no Brasil",
     sub = "Fase de ascencao da epidemia", lines(predict(ajuFA1, list(x = curve4$sem))))
coef(ajuFA1)
k1 <- as.numeric(coef(ajuFA1)['k'])
k1
```

### R0 para a Epidemia 2016 - 2017

```{r}
# Selecionando dados da ultima epidemia para estimar os parametros:
curve1 <- curve[curve$dia>=6400 & curve$dia<=6475,] #dia do pico de 110 casos: 6473
names(curve1) <- c("dia","casos")
t <- c(curve1$dia-6400)
curve1 <- data.frame(curve1, t)
head(curve1,2)

plot(curve1$t, curve1$casos, type= "l", ylab = "Numero de casos por dia", xlab = "Dia", 
     main = "Numero de obitos de PNH por FA entre 2016 e 2017 no Brasil",
     sub = "Fase de ascencao da epidemia")
```

```{r}
# Intervalo de tempo semanal
# Que não corresponde necessariamente à semana epidemiológica
sem <- c((curve1$t %/% 7)+1) # x %/% y = Parte inteira da divisão de x por y
curve2 <- data.frame(curve1$casos, sem)
curve2 <- group_by(curve2, sem)
curve2 <- summarise(curve2, sum(curve1.casos))
names(curve2) <- c("sem", "casos")
```

Ajuste da subida da curva utilizando uma funcao exponencial Ih = Ch x exp(k x tempo) onde Ch e k sao parametros de ajuste.

```{r}
ajuFA2 <- nls(casos ~ Ch*exp(k*sem), data = curve2, 
              start = list(Ch=1, k=1),
              trace = TRUE)

plot(curve2$sem, curve2$casos, ylab = "Numero de casos por semana", xlab = "Semana", 
     main = "Numero de obitos de PNH por FA entre dez/2016 e dez/2017 no Brasil",
     sub = "Fase de ascencao da epidemia", lines(predict(ajuFA2, list(x = curve2$sem))))
coef(ajuFA2)
k2 <- as.numeric(coef(ajuFA2)['k'])
k2
```


## MODELAGEM

### Estimativas da população de primatas não humanos no Brasil

1. Os passos consistem em estimar o tamanho da área de floresta remanescente nos biomas onde houve registro de óbitos (Amazônia, Pantanal, Cerrado e Mata Atlântica).
2. Em seguida, multiplicar a área disponível pela densidade média da comunidade de primatas nesses biomas.

- Fonte de informação de área nativa remanescente para os biomas: <http://www.fao.org/3/a-az172e.pdf>

- Fonte de informação sobre a densidade de primatas no fragmento:
Mata Atlântica: <https://www.researchgate.net/publication/250038409_Densidade_e_tamanho_populacional_de_primatas_em_um_fragmento_florestal_no_sudeste_do_Brasil>
<https://www.researchgate.net/publication/232240217_Densidade_tamanho_populacional_e_conservacao_de_primatas_em_fragmento_de_Mata_Atlantica_no_sul_do_Estado_de_Minas_Gerais_Brasil>
Amazônia: <https://tede.ufam.edu.br/bitstream/tede/3613/1/Disserta%C3%A7%C3%A3o%20-%20L%C3%ADvia%20Rodrigues%20da%20Silva.pdf>
Pantanal e Cerrado: valores médios

Construindo um dataframe com as informações:

```{r}
biom <- c("Pantanal", "Amazon", "Cerrado", "Atlantic Forest")
area <- c(9376913, 354221815, 77929220, 22134124)
prim <- c(0.3, 0.9, 0.3, 0.3)
pnh.dens <- c(0.2246666, 0.215, 0.2246666, 0.2343333)
remain <- data.frame(biom, area, prim, pnh.dens)
remain
```

Assim, a estimativa do tamanho da população de primatas nos biomas brasileiros:

```{r}
habitat <- remain$area*remain$prim
pop <- habitat*remain$pnh.dens
sum(pop)
pop.es <-  483158*0.3*0.2343333
pop.es
```

Nem toda a área vegetada disponível é habitada pelos primatas, dado que, em diferentes intensidades, esses animais são fortemente influenciados pela borda de mata. Por isso, apenas fragmentos com uma área *core* satisfatória podem ser usados para estimar a população real. É preciso refinar os dados de área remanescente nos biomas e da instensidade da fragmentação para obter a proporção habitável, e a partir daí, estimar a população total desses animais. 
Considerando uma letalidade média entre 80% e 90%, e assumindo que não há subnotificação, o total de casos de FA em PNH no ES foi:

```{r}
# selecionando os dados:
fa.es <- fa[fa$uf=="ES",]
sum(fa.es$casos)
# considerando que os óbitos representam 80% dos casos (alpha = 0.8):
casos.totais <- (sum(fa.es$casos)*100)/80
casos.totais # assumindo que não há subnotificação
```


```{r}
# Convertendo data em vetor de valores inteiros (= contagem de dias):
fa.es$data <- ymd(fa.es$data)
fa.es <- arrange(fa.es,data) 

inicio <- fa.es[1,1]
nextday <- as.integer(fa.es[,1] - inicio)
dia <- nextday+1
fa.es0 <- data.frame(fa.es,dia)

# Nesse novo df, agrupar o número de registros por dia:
fa.es0 <- group_by(fa.es0, dia)
curve.es <- summarise(fa.es0,sum(casos))
names(curve.es) <- c("dia", "casos")

plot(y = curve.es$casos, x = curve.es$dia, type="l", main = "obitos de PNH no ES entre 2016 e 2018", xlab = "Dia", ylab = "Numero de casos")
```


### Modelo determinístico para doença transmitida por vetores (SIR/SI)

A formulação apresentada a seguir presupõe que as populações sejam constantes.

```{r, message=FALSE}
library(deSolve)
```

A proposta a seguir é reproduzir o padrão da epizootia de febre amarela ocorrida no estado do Espírito Santo entre dezembro de 2016 e março de 2018, utilizando os parâmetros estimados por Massad et al (2018) no artigo intitulado *The risk of urban yellow fever resurgence in Aedes-infested American cities*. No artigo as taxas são apresentadas por mês, e aqui foram convertidos para taxas diárias, conforme abaixo:

```{r}
a <- 0.137 # taxa de picadas
b <- 0.6062 # fracao de picadas infectantes
c <- 0.5265 # suscetibilidade do Ae. aegypti a YF
gamam <- 0.01072 # taxa de latencia em Ae. aegypti a YF
mum <- 0.0777 # taxa de mortalidade natural de mosquitos
```

Mesmo que o número de casos nos dados utilizados abranja todas as espécies de primatas não humanos acometidas, a taxa de mortalidade natural corresponde à média para os três gêneros mais afetados pela febre amarela: *Alouatta*, *Callithrix* e *Cebus* (~1/25 anos). A taxa de recuperação utilizada foi obtida do mesmo artigo citado anteriormente, referente aos hospedeiros humanos, dada a escassez de dados específicos para os PNH.

```{r}
muh <- 1.096e-4 # taxa de mortalidade natural de *PNH* (1/25 anos)
gamah <- 0.01072 # taxa de recuperacao YF em *humanos*
alfah <-  0.198 # taxa de letalidade por YF em *humanos* ~ 0.198

par.vetor <- c(a=a, b=b, c=c, muh=muh, gamah=gamah, alfah=alfah, gamam=gamam, mum=mum)
```

Estabelecendo o vetor com as variáveis e as condições iniciais:

```{r}
Sh0 <- 34000
Sm0 <- 1020000
m <- Sm0/Sh0 # Ross sugere limiar = 1.5 

R0 <- a*a*b*c*m*gamam/(mum*(mum+gamam)*(muh+alfah+gamah))
R0

state.vetor <- c(Sh=Sh0,Ih=5,Rh=0,Sm=Sm0,Lm=0,Im=0)
```
Este modelo corresponde à condição de entrada de um hospedeiro infectado no sistema, composto por uma população de suscetíveis em sua totalidade.

```{r}
# Tempo de simulacao
tsim <- 5*365
Dt <- 1
```

A função abaixo corresponde ao modelo de Ross-Macdonald, cujo termo de transmissão é frequência dependente, para mosquitos adultos. Neste modelo, sendo considerado o pressuposto de população constante, a taxa de natalidade equivale à soma das taxas de mortalidade natural incidentes sobre cada compartimento (muh x Sh + muh x Ih + muh x Rh = muh = natalidade). 

```{r}
mod.vetor <- function(t,state,parameters){
        with(as.list(c(state,parameters)),{
                
                # Equacoes para hospedeiros
                dSh <- - a*b*Im*(Sh/(Sh+Ih+Rh)) + muh*Ih + muh*Rh + alfah*Ih
                dIh <- a*b*Im*(Sh/(Sh+Ih+Rh)) - gamah*Ih - muh*Ih - alfah*Ih
                dRh <- gamah*Ih - muh*Rh
                
                # Equacoes para vetores
                dSm <- -a*c*Sm*(Ih/(Sh+Ih+Rh)) + mum*Lm + mum*Im
                dLm <- a*c*Sm*(Ih/(Sh+Ih+Rh)) - (gamam + mum)*Lm
                dIm <- gamam*Lm - mum*Im
                
                # return the output of the model
                return(list(c(dSh, dIh, dRh, dSm, dLm, dIm)))
                
        })
}

tempos <- seq(from=0,to=tsim,by=Dt)

mod <- ode(y = state.vetor, times = tempos, func = mod.vetor, parms = par.vetor, method = "ode45")

# Criando um dataframe com os resultados:
mod <- as.data.frame(mod)
names(mod) <- c("t","Sh", "Ih", "Rh", "Sm", "Lm", "Im")

# Grafico de infectados
par(mfrow=c(1,2))
plot(mod$t, mod$Ih, col="red", type="l",
     xlab="Tempo (dias)", ylab="Número de primatas infectados")
plot(mod$t, mod$Im, col="purple", type="l",
     xlab="Tempo (dias)", ylab="Número de mosquitos infectados")
```

## CONSIDERAÇÕES FINAIS

Os motivos pelos quais a simulação não reproduz exatamente o gráfico de registros de óbitos de pnh para o ES podem ser muitos:

- Em primeiro lugar, as taxas utilizadas são estimadas a partir do surto de Dengue no Rio de Janeiro em 2012. Ou seja, as espécies de vetor, hospedeiro e patógeno, assim como o ambiente, são diferentes.

- Em segundo lugar, a imprecisão dos registros pode ocultar a realidade do evento. Subnotificações e outras dificuldades mascaram aspectos importantes.

- Entretanto, a modelagem pode sugerir o que os dados não apresentam: a quantidade de indivíduos infectados durante o surto.

De forma geral, a simulação subestima a velocidade real do evento. A epizootia de 2017-2018 no ES ocorreu de forma mais súbita e com um número máximo de óbitos registrados por dia igual a 110.
