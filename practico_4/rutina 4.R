## ingreso de datos 

Data <- read.table("daturaHER2.txt", header=TRUE)
head(Data) #vistazo de los primeros valores
Data$animal <- as.factor(Data$animal)

## ingreso del pedigree
Ped <- read.table("daturaPED.txt", header=T)
head(Ped) # inspección de los primeros individuos

## reordenar el pedigree y asignar factores
library(MasterBayes)
Ped<-orderPed(Ped)
Ped$ID<-as.factor(Ped$ID)
Ped$MOTHER<-as.factor(Ped$MOTHER)
Ped$FATHER<-as.factor(Ped$FATHER)

## PRIMER MODELO BAYESIANO
library(MCMCglmm)
prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R =
                 list(V = 1,nu = 0.002))
model1 <- MCMCglmm(arcsinR2 ~ 1, random = ~animal, pedigree =
                     Ped, data = Data, prior = prior1, verbose = FALSE)

# diagnósticos numéricos. Idealmente, las correlaciones
# para lag>0 deben ser lo más bajas posibles
autocorr(model1$VCV)

#diagnósticos gráficos
#examinar variación alrededor de la media e histograma
plot(model1$Sol)

# examinar variación alrededor de la media e histograma
# y que no haya tendencias o serpenteos.
plot(model1$VCV)

## SEGUNDO MODELO, AUMENTANDO ITERACIONES
model2 <- MCMCglmm(arcsinR2 ~ 1, random = ~animal, pedigree = Ped,
                   data = Data, nitt = 130000, thin = 100, burnin = 30000,
                   prior = prior1, verbose = FALSE)

autocorr(model2$VCV)

## estimar la varianza posterior y el intervalo de credibilidad
posterior.mode(model2$VCV)
HPDinterval(model2$VCV)

## estimación de la heredabilidad
post.her<-model2$VCV[,"animal"] /
  (model2$VCV[,"animal"]+ model2$VCV[,"units"])# cadena
posterior.mode(post.her) # HEREDABILIDAD
plot(post.her) # examen gráfico


## MODELO CON EFECTOS MATERNOS

# prior
prior3 <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V =
                                                                  1,nu = 0.002)), R = list(V = 1,nu = 0.002))  
# modelo
model3 <- MCMCglmm(arcsinR2 ~ 1, random = ~animal + Dam, pedigree =
                     Ped, data = Data, nitt = 130000, thin = 100, burnin = 30000, prior =
                     prior3, verbose = FALSE)  

# cadena
post.her3 <- model3$VCV[, "animal"]/
  (model3$VCV[,"animal"] + model3$VCV[, "Dam"] + model3$VCV[, "units"])

HPDinterval(post.her3, 0.95)# intervalo de credibilidad

posterior.mode(post.her3) # heredabilidad

## COMPARACIÓN DE MODELOS

model2$DIC # sin efectos maternos
model3$DIC # con efectos maternos
