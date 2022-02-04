
## ----------------------------------------------------------------
## Topicos de Macroeconometria
## Trabajo final:  Replicacion del articulo titulado
##                 Fiscal policy in the US: Sustainable after all?
##                 Autores: Pierre Aldama y Jerome Creel
## Presentado por: Carlos Andres Zapata Q.
## Actualizacion de los datos hasta 2020
## ----------------------------------------------------------------

# Librerias usadas
library(TSA)
library(mFilter)
library(seasonal)
library(car)
library(zoo)
library(tseries)
library(readr)
library(orcutt)
library(ggplot2)
library(MSwM)
library(sandwich)
library(strucchange)

rm(list=ls())
# Carga de datos (actualizar directorio)
setwd('C:/Users/DELL/Desktop/data/')
datos <- read_delim("data_all.csv", ",", escape_double = FALSE, trim_ws = TRUE)

## ----------------------------------------------------------------
## ----------------------------------------------------------------

# Series de tiempo de datos anuales
bt <- ts(datos$b/100, start = c(1940), freq = 1)    # Deuda publica
st <- ts(datos$s/100, start = c(1940), freq = 1)        # Balance primario
gt <- ts(datos$g, start = c(1940), freq = 1)          # Gasto total
tt <- ts(datos$t, start = c(1940), freq = 1)     # Ingresos
it <-  ts(datos$i, start = c(1940), freq = 1)      # Interes neto nominal
rt <-  ts(datos$r, start = c(1940), freq = 1)         # Interes real
gdpn <- ts(log(datos$gdpn), start = c(1940), freq = 1) # Log GDP
gdp <- ts(log(datos$gdp), start = c(1940), freq = 1) # Log GDP
expend <- ts(log(datos$expend), start = c(1940), freq = 1) # Log expenditures
yt <-  ts(datos$y, start = c(1940), freq = 1)         # Tasa de crecimiento
rtadj <-  ts(datos$radj, start = c(1940), freq = 1)   # tasa interes ajustada

## ----------------------------------------------------------------
## ----------------------------------------------------------------

# Fig. 1. Surplus-debt correlations (United States, 1940-2016).
# Panel (a)

windows()
plot(bt[2:80],st[3:81],ylab ="Primary Balance/GDP", xlab="Lagged public debt/GDP",pch=20,ylim = c(-0.3,0.1),col="darkblue")
abline(lm(st[3:81]~bt[2:80]),col="darkblue")

# Figure 2. Federal debt and primary federal surplus in the US (1940-2016).

years <- 1940:2020
df1 <- data.frame(bt,years)
df2 <- data.frame(st,years)
ylim.prim <- c(0, 1.2)   
ylim.sec <- c(-0.3,0.08)  
b1 <- diff(ylim.prim)/diff(ylim.sec)
a1 <- ylim.prim[1] - b1*ylim.sec[1] 

windows()
ggplot(df1) + 
    geom_area(aes(x=years, y=bt),fill = "grey",limits=c(0,1)) + 
    geom_line(aes(x=years, y=a1+st*b1)) +
    scale_y_continuous(name = "Public Debt ", 
                       sec.axis = sec_axis(~(.-a1)/b1, name="Primary Balance")) +
    ggtitle("Fig 2. Federal debt and primary federal surplus in the US (1940-2020)")

## ----------------------------------------------------------------
## ----------------------------------------------------------------

# Filtro HP --- Brecha del producto y gasto ciclico

gdp.hp <- hpfilter(gdp, freq = 100)
gdp.gap <- gdp.hp$cycle

windows()
par(mfrow = c(1, 2), cex = 0.8)
plot.ts(gdp, ylab = "Log GDP")  # plot time series
lines(gdp.hp$trend, col = "red")  # include HP trend
legend("topleft", legend = c("Log GDP", "Trend"), lty = 1, 
       col = c("black", "red"), bty = "n")
plot.ts(gdp.hp$cycle, ylab = "", main="Brecha del producto")  # plot cycle
legend("topright", legend = c("Cycle"), lty = 1, col = c("black"), 
       bty = "n")

g.hp <- hpfilter(expend, freq = 100)
g.cycle <- g.hp$cycle

windows()
par(mfrow = c(1, 2), cex = 0.8)
plot.ts(expend, ylab = "Log Gasto")  # plot time series
lines(g.hp$trend, col = "red")  # include HP trend
legend("topleft", legend = c("Log Gasto", "Trend"), lty = 1, 
       col = c("black", "red"), bty = "n")
plot.ts(g.hp$cycle, ylab = "",main="Comp. ciclico del gasto")  # plot cycle
legend("topright", legend = c("Cycle"), lty = 1, col = c("black"), 
       bty = "n")


years <- 1940:2020
df3 <- data.frame(gdp.gap,years)
df4 <- data.frame(g.cycle,years)
ylim.prim2 <- c(-0.12, 0.12)   
ylim.sec2 <- c(-0.6,0.7)  
d <- diff(ylim.prim2)/diff(ylim.sec2)
c <- ylim.prim2[1] - d*ylim.sec2[1] 

windows()
ggplot(df3) + 
    geom_line(aes(x=years, y=gdp.gap)) + 
    geom_line(aes(x=years, y=c+g.cycle*d),linetype=2) +
    scale_y_continuous(name = "Output gap", 
                       sec.axis = sec_axis(~(.-c)/d, name="Cyclical public spending")) +
    ggtitle("Fig 3. HP filtered output gap and cyclical public spending (1940-2016).") +
    theme(legend.position="bottom")

## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Estimaciones

# Estimación FRF de Bohn (Variables de control: gdp.gap, g.cycle)
# Periodo 1942-2016

st <- ts(st[3:81])
lagbt <- ts(bt[2:80])
gdp.gapt <- ts(gdp.gap[3:81])
g.cyclet <- ts(g.cycle[3:81])
# Componentes cuadratico y cubico del modelo
lagbt2 <- lagbt^2
lagbt3 <- lagbt^3

# 1. Modelo FRF lineal
model1 <- lm(st ~ lagbt + gdp.gapt + g.cyclet) # Ec. (3)
#summary(model1)
NLLSco1 = cochrane.orcutt(model1)
summary(NLLSco1)

# 2. Modelo FRF Cuadratica
model2 <- lm(st ~ lagbt + lagbt2 + gdp.gapt + g.cyclet)
#summary(model2)
NLLSco2 = cochrane.orcutt(model2)
summary(NLLSco2)

# 3. Modelo FRF Cubica
model3 <- lm(st ~ lagbt + lagbt2 + lagbt3 + gdp.gapt + g.cyclet)
#summary(model3)
NLLSco3 = cochrane.orcutt(model3)
summary(NLLSco3)


## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Estimación del modleo de Markov-Switching 

model1.MS <- msmFit(model1, k=2, sw=c(T,T,T,T,F,T), p=1,control=list(parallel=T))
summary(model1.MS)

windows()
plotProb(model1.MS,which=1)

#slotNames(mod@Fit)
filtprob=(model1.MS@Fit@filtProb)  # Filtered probability
filtprob.res = ts(filtprob[,1],start=c(1942),end=c(2020),freq=1)
smoprob=(model1.MS@Fit@smoProb) # Smoothed probability
smoprob.res = ts(smoprob[,1],start=c(1942),end=c(2020),freq=1)


windows()
plot(filtprob.res,type="h",lwd =10,col="grey", ylab= "Probability",
     main="Fig. 4. Baseline model, probabilities of sustainable regime",cex=0.7)
lines(smoprob.res,col="black")
legend("top",legend = c("Smoothed probability","Filtered probability"),
       lty=c(1,1),cex=0.8,col=c("gray","black"))

## ------------------------------------------------------------
## ------------------------------------------------------------

# Debt-Stabilizing condition and Structural breaks

bp.ri <- breakpoints(rtadj ~ 1, h = 20)
bp.ri2 <- lm(rtadj ~ breakfactor(bp.ri, breaks = 2))
summary(bp.ri2)

years <- 1941:2020
bp.rits <- ts(fitted(bp.ri2, start = 1940))
df5 <- data.frame(rtadj[2:81],years)
df6 <- data.frame(bp.rits[2:81],years)
ylim.prim3 <- c(-0.30, 0.1)  

windows()
ggplot(df1) + 
    geom_area(aes(x=years, y=rtadj),fill = "grey",limits=c(-0.3,0.2)) + 
    geom_line(aes(x=years, y=bp.rits)) +
    ggtitle("Fig. 5. Growth-adjusted real interest rate (1941-2016)")


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

