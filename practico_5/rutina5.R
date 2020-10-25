## Práctico 5

## (SOLO RANDOM SKEWERS)

## ingreso de datos
tapi<-read.table("TAPI.txt", header = TRUE)
horq<-read.table("HORQ.txt", header = TRUE)

Gtapi <- cov(tapi, use = "complete.obs")
Ghorq <- cov(horq, use = "complete.obs")

## la función random skewers

# Argumentos
# G1 y G2 son matrices de varianza-covarianza.
# rep es el número de random skewers a aplicar.

# Output.
# La salida es un data.frame de dos columnas. La columna ratio es el 
# cociente entre la longitud de los vectores respuesta correspondientes 
# a G1 y G2. La columna correlation contiene medidas de la co-linealidad
# de los vecores respuest en el espacio morfométrico multivariado, que 
# equivale al coseno del ángulo entre los vectores respuesta.

random.skewers<-function(G1, G2, rep){
  N1 <- ncol(G1)
  N2 <- ncol(G2)
  if(N1 != N2) stop("G1 and G2 have different dimensions") else N <- N1
  rd.sk <- function(G1, G2){
    RS1 <- runif(n=N)		
    sign <- c(-1,1)
    RS2 <- vector(length=N)
    for (i in 1:N) RS2[i] <- sample(sign, size=1)*RS1[i]
    RS3 <- RS2/as.vector(sqrt((t(RS2)%*%RS2)))
    dz1 <- G1%*%RS3
    dz2 <- G2%*%RS3
    res <- c(sum(dz1^2)/sum(dz2^2), 
             (t(dz1)%*%dz2)/sqrt((t(dz1)%*%dz1)*(t(dz2)%*%dz2)))
    return(res)
  }
  RES <- data.frame(t(replicate(rep, rd.sk(G1=G1, G2=G2))))
  colnames(RES) <- c("ratio", "correlation")
  return(RES)
}

## la función sim.vec.cor

# Argumentos.
# k es el número de rasgos en el análisis de random skewers.
# rep es el número de pares de vectores simulados.

# Output.
# Un vector de longitud igual a rep conteniendo la correlación entre los
# vectores creados al azar. 

sim.vec.cor <- function(k, rep){
  vec.cor <- function(k){
    RS1 <- runif(n=k)*sample(c(-1,1), size = k, replace = TRUE)
    RS1 <- RS1/as.vector(sqrt((t(RS1)%*%RS1))) 
    RS2 <- runif(n=k)*sample(c(-1,1), size = k, replace = TRUE)
    RS2 <- RS2/as.vector(sqrt((t(RS2)%*%RS2)))
    ang <- acos(t(RS1)%*%RS2/as.vector(sqrt((t(RS1)%*%RS1)*(t(RS2)%*%RS2))))
    res <- ifelse(ang > pi/2, cos(pi-ang), cos(ang))
    return(res)
  }
  replicate(rep, vec.cor(k))
}

# Crear los random skewers.
RS <- random.skewers(G1 = Gtapi, G2 = Ghorq, rep = 1000)

# Obtener la respuesta media para la proporción entre la 
# longitud de los vectores de respuesta
mean(RS$ratio)

# Obtener la correlación media entre los vectores de respuesta
mean(RS$correlation)

# Estimar el percentil 95% de la correlación esperada
# entre vectores aleatorios
SC <- sim.vec.cor(k = 6, rep = 1000)
quantile(SC, probs = 0.95)