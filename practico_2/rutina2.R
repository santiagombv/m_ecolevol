#### Cargando paquetes y el set de datos
library(mgcv)
library(knitr)

data <- read.table("cyclop.txt", header = TRUE)
data <- na.omit(data) # evitar los datos faltantes

### Caso 1. Datos de conteos.
#### Ajustar el modelo.   

fit1 <- gam(exported.pollinaria ~ s(nectary.depth, bs = "cr") +
              s(flower.number, bs = "cr"),
            data = data, family = quasipoisson)
summary(fit1) ## revisar el parámetro Scale est (phi).

# Al ser Scale est < 2 podemos reajustar por Poisson 
fit2 <- gam(exported.pollinaria ~ s(nectary.depth, bs = "cr") +
              s(flower.number, bs = "cr"),
            data = data, family = poisson)
summary(fit2)

#### Gráficos univariados (para selección lineal y cuadrática)
## mediante el paquete visreg
library(visreg)

layout(matrix(1:2,1,2))
visreg(fit2, xvar = "nectary.depth", scale = "response", partial = FALSE,
       ylim = range(data$exported.pollinaria))
points(exported.pollinaria ~ nectary.depth, data = data, pch = 19)#puntos observados

visreg(fit2, xvar = "flower.number", scale = "response", partial = FALSE,
       ylim = range(data$exported.pollinaria))
points(exported.pollinaria ~ flower.number, data = data, pch = 19)#puntos observados
layout(1)

## mediante los paquetes tidymv + ggplot2
library(patchwork)
library(tidymv)
library(ggplot2)
g1 <- plot_smooths(fit2, series = nectary.depth, transform = exp) +
  geom_point(data=data, aes(x = nectary.depth, y = exported.pollinaria)) + 
  theme_linedraw()
g2 <- plot_smooths(fit2, series = nectary.depth, transform = exp) +
  geom_point(data=data, aes(x = nectary.depth, y = exported.pollinaria)) +
  theme_linedraw()
g1 + g2

#### Gráfico para superficies (selección correlacional)  
## mediante la función vis.gam de mgcv
layout(matrix(1:2,1,2))
vis.gam(fit2, view = c("flower.number", "nectary.depth"), type = "response",
        plot.type = "contour", color = "cm")
points(data$flower.number, data$nectary.depth, pch = 16, col ="grey")

vis.gam(fit2, view = c("flower.number", "nectary.depth"), type = "response",
        plot.type = "persp", color = "cm", phi = 10, theta = -45, 
        ticktype = "detailed") -> pp
layout(1)


## mediante la función predict_gam de tidymv
library(tidymv)
library(viridis)
library(ggplot2)

## exp es usada para invertir el enlace típico de los modelos Poisson
ggplot(predict_gam(fit2), aes(flower.number, nectary.depth, z = exp(fit))) + 
  geom_raster(aes(fill = exp(fit))) +
  geom_contour(colour = "white") +
  scale_fill_viridis(name = "pol", option = "viridis") + 
  theme_minimal() + geom_point(data = data, aes(flower.number, nectary.depth, z=NULL), color = "grey")         

### Caso 2. Datos de proporciones.
#### Ajustar el modelo   

resp <- cbind(data$exported.pollinaria, data$flower.number - data$exported.pollinaria)
fit3 <- gam(resp ~ s(nectary.depth, bs = "cr") +
              s(flower.number, bs = "cr"),
            data = data, family = quasibinomial)
summary(fit3) ## revisar el parámetro Scale est (phi).

# Al ser Scale est < 2 nos quedamos con el modelo que controla
# la sobredispersión

#### Gráficos univariados (para selección lineal y cuadrática)
## mediante el paquete visreg
library(visreg)

layout(matrix(1:2,1,2))
visreg(fit3, xvar = "nectary.depth", scale = "response", partial = FALSE,
       ylim = range(data$prop.pol))
points(prop.pol ~ nectary.depth, data = data, pch = 19)#puntos observados

visreg(fit3, xvar = "flower.number", scale = "response", partial = FALSE,
       ylim = range(data$prop.pol))
points(prop.pol ~ flower.number, data = data, pch = 19)#puntos observados
layout(1)

## paquete tidymv + ggplot2
library(patchwork)
library(tidymv)
library(ggplot2)
library(boot) #necesaria para inv.logit

g1 <- plot_smooths(fit3, series = nectary.depth, transform = inv.logit) +
  geom_point(data=data, aes(x = nectary.depth, y = prop.pol)) + 
  theme_linedraw()
g2 <- plot_smooths(fit3, series = flower.number, transform = inv.logit) +
  geom_point(data=data, aes(x = flower.number, y = prop.pol)) +
  theme_linedraw()
g1 + g2

#### Gráfico para superficies (selección correlacional)  
## mediante la función vis.gam de mgcv
layout(matrix(1:2,1,2))
vis.gam(fit3, view = c("flower.number", "nectary.depth"), type = "response",
        plot.type = "contour", color = "cm")
points(data$flower.number, data$nectary.depth, pch = 16, col ="grey")

vis.gam(fit3, view = c("flower.number", "nectary.depth"), type = "response",
        plot.type = "persp", color = "cm", phi = 10, theta = -45, 
        ticktype = "detailed") -> pp
layout(1)


## mediante la función predict_gam de tidymv
library(tidymv)
library(viridis)
library(ggplot2)
library(boot)

## inv.logit es usada para invertir el enlace típico de los modelos Binomiales
ggplot(predict_gam(fit3), aes(flower.number, nectary.depth, z = inv.logit(fit))) + 
  geom_raster(aes(fill = exp(fit))) +
  geom_contour(colour = "white") +
  scale_fill_viridis(name = "pol", option = "viridis") + 
  theme_minimal() + geom_point(data = data, aes(flower.number, nectary.depth, z=NULL), color = "grey")         


### END ###



