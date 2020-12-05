


#-------------------------------------
# Base de datos
datos <- read.csv2("Insectos.csv",row.names=1)
datos = read.csv2("Insectos.csv",row.names=1)
datos <- read.csv2(file.choose(),row.names=1)

#-------------------------------------
# LIBRERÍAS REQUERIDAS
library(lattice)
library(ellipse)
require(SciViews)
require(stats)

# Estructura de la base de datos
str(datos)
datos$cuenca=as.factor(datos$cuenca)
str(datos)

summary(datos[,2:9])

# 1. Gráfica por pares
pairs(datos[,2:8])
pairs(log10(datos[,2:8]))


# 2.Figura elipses 
plotcorr(cor(datos[,2:9]))
help(plotcorr)

corr.mtcars <- cor(datos[,2:9])
ord <- order(corr.mtcars[1,])
xc <- corr.mtcars[ord, ord]
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")   
plotcorr(xc, col=colors[5*xc + 6])


# 3. Figuras de pares
pairs ((datos[,c(2:9)]),panel=function(x,y)
  {abline(lsfit(x,y)$coef,lwd=2,col=3)
  lines(lowess(x,y),lty=2,lwd=2,col=2)
  points(x,y,cex=1)})

pairs ((datos[,c(2,6,7,9)]),panel=function(x,y)
  {abline(lsfit(x,y)$coef,lwd=2,col=3)
  lines(lowess(x,y),lty=2,lwd=2,col=2)
  points(x,y,col=datos$cuenca, cex=1.4)})

# 4. 
x11()
pairs(datos[, 2:9], diag.panel = panel.hist, 
  upper.panel = panel.smooth, lower.panel = panel.cor)


# 5. Figura con tres variables (Función: coplot)
with(datos,coplot(Efem~pH|temp))

with(datos, {
  coplot(Efem~pH|temp, number = 3,
 panel = function(x, y, ...) panel.smooth(x, y, span = .8, ...))
  coplot(Efem~pH|temp,
 panel = panel.smooth)
})


# 6. Coplot 
summary(datos[,2:8])

clasetemp<-cut(datos$temp,seq(15,20,1.2),include.lowest=T)
clasetemp
clasepH<-cut(datos$pH,seq(5,8,1,include.lowest=T))
clasepH

coplot(Efem~pH | clasetemp, pch=19,
   panel = panel.lm, data=datos)

# 7. 
panel.lm = function(x, y, ...) {
  tmp<-lm(y~x,na.action=na.omit)
  abline(tmp, lwd = 1.5, col= 2)
  points(x,y, ...)}

coplot(Efem ~ pH | clasetemp, pch=19,
   panel = panel.lm, data=datos)


# Categorizando variables continuas
splom(~datos[,4:8]|clasepH,pscales=0)   
splom(~datos[,4:8]|clasepH+clasetemp,pscales=0)

# Reación por niveles del factor Cuenca
# 
datos$cuenca<-factor(datos$cuenca, 
 levels=c("cuen3","cuen4","cuen1","cuen2"))

# 7. Figuras xyplot 
xyplot(Efem~pH|cuenca,data=datos,
   panel = panel.lm, data=datos)


# 8. Histogramas  
histogram	(~Ab,data=datos, ylab="Porcentaje del Total",
   xlab="Abundancia de insectos")
x11()
histogram	(~Ab|cuenca,data=datos, ylab="Porcentaje del Total",
   xlab="Abundancia de insectos")


# 8. Figuras de densidad
densityplot(~Ab,data=datos, ylab="Porcentaje del Total",
xlab="Abundancia de insectos")  
densityplot(~Ab|cuenca,data=datos, ylab="Porcentaje del Total",
xlab="Abundancia de insectos")


# 9. qqplot
panel<-par(mfrow=c(1,2), mar=c(4,3,3,2))

# figura con datos crudos
qqnorm (datos$Ab, main="Abundancia de Insectos",
ylab="Cuantiles de la muestra",
xlab="Cuantiles teóricos") 
qqline(datos$Ab)

# 10. figura con raíz log de abndancias
Ab.log <- log10(datos$Ab+1)
qqnorm (Ab.log, main="Log de Abundancia de Insectos",
ylab="Cuantiles de la muestra",
xlab="Cuantiles teóricos") 
qqline(Ab.log)
par(panel)
panel<-par(mfrow=c(1,1))



# 11. 
plot(Efem~Plec,col=as.integer(cuenca),data=datos,ylab="",
 xlab="Plecópteros") 
legend(0,27,legend=levels(datos$cuenca),pch=19,col=1:4,cex=0.8)
lines(abline(lm(datos$Efem~datos$Plec),lwd=2,col=2, lty=2))
par(panel)



# 14. Colocar colores al grafico por cuenca
xyplot(Efem~Plec,group=cuenca,auto.key=T,data=datos)


# 15. 
plot(Efem~Plec,col=as.integer(cuenca),data=datos) 
legend(0,25,legend=levels(datos$cuenca),pch=19,col=1:4,cex=0.8)
lines(abline(lm(datos$Efem~datos$Plec),lwd=2,col=2, lty=2))



# 16. Figuras de Cajas y cinturas
datos$cuenca<-factor(datos$cuenca, 
 levels=c("cuen1","cuen2","cuen3","cuen4"))

boxplot(Ab~cuenca,data=datos,
        xlab="Cuencas",ylab="Abundancia",
        col="lightgray", cex.lab=1.3)

boxplot(Ab~cuenca,data=datos,notch=TRUE,
xlab="Cuencas",ylab="Abundancia",
col="lightgray", cex.lab=1.3)






#-------------
# PRUEBA DE HIPÓTESIS
# 1) Anova 1 vía (Abundancia vs. cuenca)

# Qué es el Anova a una vía?
# R./
# Cuál será la hipótesis nula (Ho)?
# R./

Ab.anov <- aov(Ab ~ cuenca , data=datos)
summary(Ab.anov) 
# Se acepta Ho? Por qué?
# R./


# 1.1) Supuesto de Normalidad
shapiro.test(Ab.anov$residuals)
# Se acepta el supuesto? Por qué?
# Ho Hay normalidad en los residuales del Anova
# R./  


# 1.2) Homogeneidad de varianzas
bartlett.test (Ab ~ cuenca , data=datos)
# Se acepta el supuesto? Por qué?
# R./


# 1.3) Independencia de errores - Durbin Watson 
library(car)
# Se debe generar un modelo lineal relacionado con el Anova
modelo<-lm(Ab ~ cuenca, data=datos)
names(modelo)
durbinWatsonTest(modelo)
# Se acepta el supuesto? Por qué?
# R./

# Grafico de incremento de los residuales de Durbin Watson
plot(residuals(modelo),type="l",ylab="Residuales", xlab="Indice")


# Otros tipos de Anovas

# 1) Prueba Kruskal-Wallis
# Qué es estadístico de Kruskal Wallis?
# R./
Ab.kw <- kruskal.test (Ab ~ cuenca , data=datos)
Ab.kw
# Se acepta Ho? Por qué?
# R./


# 2) Prueba de Welch (p.333 Murray)
# Qué es estadístico de Welch?
# R./
Ab.welch <- oneway.test(Ab ~ cuenca , data=datos, var.equal = F)
Ab.welch



#--------
# Resumen de la prueba de Hipótesis
# Anova, K-W y Welch diferencias entre cuencas=
# Normalidad = 
# Homogeneidad = 


# Comparación múltiple de medianas
boxplot(Ab ~ cuenca , data=datos,
col = "lightgray", notch=TRUE,
xlab= "Cuenca", ylab="Abundancia", cex.lab=1.5)


#------------
# INVESTIGAR: Comparación multiple o prueba 
# "Contraste no paramétrico de efectos relativos" - Frank Konietschke (2008)
# Aplicar al presente ejercicio.


library(nparcomp)

multi.comp<-nparcomp(Ab ~ cuenca, data=datos, asy.method = "probit",
 type = "Tukey", alternative = "two.sided", 
 plot.simci = TRUE, info = FALSE,correlation=TRUE)

names(multi.comp)
multi.comp$Analysis
multi.comp$Contrast
summary(multi.comp)





#==========================================
#-------------------------------------
# LIBRERÍAS REQUERIDAS
library(lattice)
library(ellipse)
library(plotrix)
require(SciViews)
require(stats)

#------------------------------------
# Figura 1
# Figuras de tortas (Función "pie")
datos<-read.csv2("Datos1.csv")
str(datos)

# Suma de las biomasas
datos1  <- colSums(datos[,5:9])
datos1

#Piechart
par(mfrow = c(2,2), mar = c(3, 3, 2, 1))
#
pie(datos1 , main = "Figura Circular Ordinaria")
#
pie(datos1 , col = gray(seq(0.4,1.0,length=6)),
clockwise=TRUE, main = "Escala de Grises", angle=45)
#
pie(datos1 , col = rainbow(6),clockwise = TRUE, main="Colores de Arcoiris")

# 3D
pie3D(datos1 , labels = names(datos1), explode = 0.1,  main = "Figura Circular en 3D", labelcex=0.8)
#
par(bentos)


#-----------------------------------------
# Grupos Funcionales de Invertebrados Acuáticos "gfun"

datos<-read.csv2("Datos1.csv")

# Datos del tramo A
datosA=datos[1:10,]
datosA

# Datos de Biomasa total del tramo A
tramoA=datosA$BIOM.TOT
tramoA
names(tramoA) <- datosA[,3]
tramoA

# Datos del tramo A
datosB=datos[11:20,]
datosB

# Datos de Biomasa total del tramo B
tramoB=datosB$BIOM.TOT
tramoB
names(tramoB) <- datosB[,3]
tramoB

# Tabla de biomasa total para los dos tramos
tramos <- cbind(tramoA, tramoB)
tramos

# Panel con 4 figuras de barras
par(mfrow = c(2,2), mar = c(3, 3, 2, 1))
barplot(tramoB , main = "Biomasas")
barplot(tramos)
barplot(t(tramos), col = gray(c(0.5,1)))
barplot(t(tramos), beside = TRUE)
par(mfrow = c(1,1)



#-------------------------------
#Example 2
# Figuras de columnas y desviaciones, para grupos funcionales 
# Por tramo y periodos climáticos.
datos<-read.csv2("Datos2.csv")
head(datos)

# Promedios y desviaciones por cada GF
datos.m <- tapply(datos$Ab, INDEX=bentos$GF, FUN=mean)
datos.de <- tapply(bentos$Ab, INDEX=bentos$GF, FUN=sd)

# Tabla de de medias y desviaciones por cada GFA
datos1<- cbind(Bent.m, Bent.de)
datos1

# Figura de barras con líneas acotadas

par(mfrow = c(2,1), mar = c(3, 5, 2, 1))
barplot(Bent.m, xlab = "GFA", ylab = "Abuandancia de GFA", ylim=c(0,400))
arrows(bp, Bent.m, bp, Bent.m + Bent.de, lwd = 1.5,angle=90,length=0.1)

barplot(datos.m, xlab = "GFA",ylab = "Abundancias (Indv)", col=rainbow(9), ylim=c(0,700))
arrows(bp, datos.m, bp, datos.m + datos.de, lwd = 1.5, angle=90,length=0.1)
box()
par(mfrow = c(1,1))


#--------------------------------
# Figuras de tiras
# Figura de líneas acotadas con errores estandar (.es)
# Error estandar = desviación estandar/raiz de tamaño de la muestra (.le)

# Tamaño de la muestra (tm) para los datos de abundancias por cada GF = 8 muestreos
datos.tm <- tapply(datos$Ab, INDEX=datos$GF, FUN=length)

# 
datos.ee <- Bent.de / sqrt(datos.tm) 


# Operación "random jittering (variación)" 
# 
# 
stripchart(datos$Biom ~ datos$GF, vert = TRUE,
   pch=1, method = "jitter", jit = 0.05, 
   xlab = "Grupos Funcionales",ylab = "Abundancias (Indv)")
points (1:5,datos.m, pch = 16, cex = 1.5)


# Líneas acotadas simbolizan los errores estandar (ee) 
arrows (1:8, datos.m,1:8, datos.m + datos.de, lwd = 1.5,
angle=90, length=0.1)
arrows (1:8, datos.m,1:8, datos.m - datos.de, lwd = 1.5,
angle=90, length=0.1)





#---------------------------------------
#Figuras de Cajas y Bigotes
datos<-read.csv2("Datos2.csv") 
library(latti)
#
par(mfrow = c(2,2), mar = c(3, 5, 2, 1))
#
boxplot(Ab~GF, data = datos, ylab ="Abundancia (Indv)",
cex.lab=1.3)
#
boxplot(Ab~GF,  notch=T, data = datos, ylab="")
#
boxplot(Ab~GF * Lluvia, data = datos, ylab="Biomasa (Biom)",
cex.lab=1.3)
#
boxplot(Ab~GF*Lluvia, names= c("P1/C-F","P2/C-F","P1/C-R","P2/C-R",
   "P1/D","P2/D","P1/R","P2/R","P1/T","P2/T"),data = datos,
ylab ="")

par(mfrow = c(1,1))



