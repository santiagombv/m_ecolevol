## Cargando paquetes y set de datos
library(lavaan)
datura <- read.table("datura1.txt", header = T)

## estandarización de todas las variables
datos <- data.frame(scale(datura))

## modelado de la estructura cusal mediante regresiones

# sólo daño afecta a crecimiento
m1<- lm(crec ~ dano, data = datos)
summary(m1)

# daño y crecimiento afectan a floración
m2 <- lm(flor.dias ~ dano + crec, data = datos)
summary(m2)

# finalmente daño, crecimiento y floración afectan a la adecuación
m3 <- lm(frutos ~ dano + crec + flor.dias, data = datos)
summary(m3)

#########################################################
## modelado en lavaan

# modelo saturado
sat <- '
frutos ~ crec + flor.dias + dano
flor.dias ~ crec + dano
crec ~ dano'
fit.sat <- sem(model = sat, data = datos, missing = "listwise",
               fixed.x = F)
summary(fit.sat, standardized = T, rsq = T, fit.measures = T)

# modelo con efectos indirectos del daño,
# y sin efecto del crecimiento en la floración
exdano <- '
frutos ~ crec + flor.dias
flor.dias ~ dano
crec ~ dano'
fit.exdano <- sem(model=exdano, data=datos, missing="listwise",
                  fixed.x=FALSE)
summary(fit.exdano, standardized = TRUE, rsq = T, fit.measures = T)

# modelo con efecto del crecimiento en la floración
# y sin efecto del daño en la floración
exdano2<-'
frutos ~ crec + flor.dias
flor.dias ~ crec
crec ~ dano'
fit.exdano2 <- sem(model=exdano2, data=datos, missing="listwise",
                   fixed.x=FALSE)
summary(fit.exdano2, standardized=TRUE, rsq=T, fit.measures=T)

