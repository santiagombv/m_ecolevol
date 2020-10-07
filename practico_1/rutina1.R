# Práctico 1. Estimación de la selección fenotípica.

# Modificado de Palacio, F.X., Ordano, M. & Benitez-Vieyra, S. 2019. 
# Measuring natural selection on multivariate phenotypic traits: 
# a protocol for verifiable and reproducible analyses of natural selection. 
# Israel Journal of Ecology and Evolution 65, 130–136.



#### Cargando paquetes y el set de datos
library(car)
library(boot)
library(visreg)
library(ggplot2)
library(knitr)

data <- read.table("cyclop.txt", header = TRUE)
data <- na.omit(data)


# Paso 1. Obtener el éxito reproductivo relativo y 
# estandarizar los rasgos fenotípicos (media = 0, varianza = 1)
W <- data$exported.pollinaria
wrel <- W/mean(W)
x1 <- data$flower.number
x2 <- data$nectary.depth
z1 <- (x1 - mean(x1))/sd(x1)
z2 <- (x2 - mean(x2))/sd(x2)

# Paso 2. Colinearidad de los rasgos fenotípicos   
# Factores de inflación de la varianza en el modelo (lineal) de Lande & Arnold.     
lin.grad <- lm(wrel ~ z1 + z2)
vif(lin.grad)

# Paso 3. Construir y chequear el modelo (lineal y completo) de Lande & Arnold
lin.grad <- lm(wrel ~ z1 + z2)
summary(lin.grad)

nonlin.grad <- lm(wrel ~ z1 + z2 + I(0.5*z1^2) + I(0.5*z2^2) + z1:z2)
summary(nonlin.grad)

layout(matrix(1:4, 2, 2))
plot(nonlin.grad)
layout(1)

# Paso 4. Intervalos de confianza
# Para la función *grad*, los datos deben tener el éxito reproductivo 
# relativo en la primer columna y las variabes estandarizadas en las 
# subsecuentes columnas. La función retorna un vector con los gradientes 
# lineales, cuadráticos y correlacionales, y representa la entrada de la función boot.   

grad <- function(data, original = c(1:nrow(data))){
  data <- data[original, ]
  vars  <- colnames(data)[-1]
  colnames(data)[1] <- "Wrel"
  model.lin <- as.formula(paste("Wrel", paste(vars, collapse=" + "), sep=" ~ "))
  m1 <- lm(formula = model.lin, data = data)
  part1 <- paste("(", paste(vars, collapse=" + "), ")^2", sep = "")
  part2 <- paste("I(0.5*(", vars, "^2))", sep = "", collapse = " + ")
  model.qua <- as.formula <- paste("Wrel", paste(part1, part2, sep = " + "), sep = " ~ ")
  m2 <- lm(formula = model.qua, data = data)
  sel.grad<-c(m1$coefficients[-1], m2$coefficients[-c(1:ncol(data))])
  return(sel.grad)
}

newdata <- data.frame(wrel, z1, z2)
selection.gradients <- grad(data = newdata)
boot.grad <- boot(data = newdata, statistic = grad, R = 999)

# Intervalos de confianza BCA para cada gradiente.   
CI <- list()
for(i in 1:length(boot.grad$t0)){
  CI[[i]] <- boot.ci(boot.grad, conf = 0.95, type = "bca", index = i)$bca[4:5]
}

names(CI) <- names(boot.grad$t0)
CI <- cbind(selection.gradients, do.call(rbind, CI))
colnames(CI) <-c("sel. gradients", "lower.ci", "upper.ci")
CI

# Paso 5. Graficando los resultados del modelo de Lande & Arnold.

# Selección lineal
new.z1 <- seq(-4, 4, length = 500)
plot(z1, wrel, pch = 19, cex = 1.5, col = "gray70", ylab = "Relative fitness")
pred.z1 <- predict(lin.grad, newdata = data.frame(z1 = new.z1, z2 = mean(z2)), 
                   se.fit = TRUE)
lines(new.z1, pred.z1$fit, lwd = 2, col = "blue")
lines(new.z1, pred.z1$fit + 2*pred.z1$se.fit, lty = 3, col = "blue")
lines(new.z1, pred.z1$fit - 2*pred.z1$se.fit, lty = 3, col = "blue")

new.z2 <- seq(-4, 4, length = 500)
plot(z2, wrel, pch = 19, cex = 1.5, col = "gray70", ylab = "Relative fitness")
pred.z2 <- predict(lin.grad, newdata = data.frame(z1 = mean(z1), z2 = new.z2),
                   se.fit = TRUE)
lines(new.z2, pred.z2$fit, lwd = 2, col = "blue")
lines(new.z2, pred.z2$fit + 2*pred.z2$se.fit, lty = 3, col = "blue")
lines(new.z2, pred.z2$fit - 2*pred.z2$se.fit, lty = 3, col = "blue")

# 5.2. Selección no-lineal   
new.z1 <- seq(-4, 4, length = 500)
plot(z1, wrel, pch = 19, cex = 1.5, col = "gray70", ylab = "Relative fitness")
pred.z1 <- predict(nonlin.grad, newdata = data.frame(z1 = new.z1, z2 = mean(z2)), 
                   se.fit = TRUE)
lines(new.z1, pred.z1$fit, lwd = 2, col = "blue")
lines(new.z1, pred.z1$fit + 2*pred.z1$se.fit, lty = 3, col = "blue")
lines(new.z1, pred.z1$fit - 2*pred.z1$se.fit, lty = 3, col = "blue")

new.z2 <- seq(-4, 4, length = 500)
plot(z2, wrel, pch = 19, cex = 1.5, col = "gray70", ylab = "Relative fitness")
pred.z2 <- predict(nonlin.grad, newdata = data.frame(z1 = mean(z1), z2 = new.z2), 
                   se.fit = TRUE)
lines(new.z2, pred.z2$fit, lwd = 2, col = "blue")
lines(new.z2, pred.z2$fit + 2*pred.z2$se.fit, lty = 3, col = "blue")
lines(new.z2, pred.z2$fit - 2*pred.z2$se.fit, lty = 3, col = "blue")

# 5.3. Selección correlacional   
nonlin.grad <- lm(wrel ~ z1 + z2 + I(0.5*z1^2) + I(0.5*z2^2) + z1:z2)
visreg2d(nonlin.grad, xvar = "z1", yvar = "z2", scale = "response", plot.type = "image")
visreg2d(nonlin.grad, xvar = "z1", yvar = "z2", scale = "response", plot.type = "persp")

# 6. Otras medidas comunes en estudios de selección fenotípica
# Oportunidad para la selección I
op <- var(wrel)
op

# Diferenciales de selección lineal s_i
s1 <- lm(wrel ~ z1)
summary(s1)

s2 <- lm(wrel ~ z2)
summary(s2)