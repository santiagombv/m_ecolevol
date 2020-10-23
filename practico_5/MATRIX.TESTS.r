G.MATRICES <- function(DATA,Selected.Pop,Nos.of.Traits,STAT,FORMULA, Method)
{
# Function to calculate the G matrices
  Nos.of.Populations <- length(Selected.POP)
  N                  <- matrix(0,Nos.of.Populations)
  G.pop              <- array(0,dim=c(Nos.of.Traits,Nos.of.Traits,Nos.of.Populations))
  for ( Ith.Pop in 1:Nos.of.Populations)
{
   GROUP 			       <- DATA[DATA$POP==Selected.POP[Ith.Pop],]						# Set to GROUP
   N[Ith.Pop]        <- length(unique(GROUP$SIRE))
   G.pop[,,Ith.Pop]  <- MATRIX.CALCULATOR(GROUP, Nos.of.Traits, STAT, FORMULA, Method)
}
 Out <- list(G.pop,N)
 return(Out)
} # End of function
##########################################################################
COMMON.MATRIX <- function(G.pop, N, Nos.of.Traits, Nos.of.Populations)
 {
# Function that calculates the matrix of equality
  G.Equal   <- matrix(0,Nos.of.Traits,Nos.of.Traits)
  for ( Ith.Pop in 1:Nos.of.Populations){ G.Equal <- G.Equal + G.pop[,,Ith.Pop]*N[Ith.Pop]} # End of Ith.Pop
  Ntotal   <- sum(N)
  G.Equal  <- G.Equal/Ntotal
  X        <- array(0,dim=c(Nos.of.Traits,Nos.of.Traits,Nos.of.Populations))
  for(i in 1:Nos.of.Populations) { X[,,i] <- G.Equal} # Set MLE for all populations
  G.Equal  <- X
  return(G.Equal)
} # End of function
################################################################################
PROPORTIONAL.MATRICES <- function(Smatrices, n, Nos.of.Traits, N.matrices, Tolerance)
{
# Coding for Model 1 from
#Manly Bryan, F. J. & Rayner, J. C. W. 1987. 
#The comparison of sample covariance matrices using likelihood ratio tests. 
#Biometrika 74: 841-847.
# Smatrices: array containing sample matrices
# N = total sample size
# n = vector of sample size for each matrix
# Nos.of.Traits = number of traits = rank of matrix
# N.matrices = number of matrices
# c.j  =proportionality constants  
# Tolerance = minimum accuracy 
	N     <- sum(n)
  c.j  <- rep(1, times= N.matrices) 
	Diff <- 1
	LL   <- -1e-5
	Iteration.Counter <- 0
#	while (Diff > Tolerance)
for( rep in 1:10)
{
	Old.LL 	<- LL
	Omega.prime 	<- matrix(0,Nos.of.Traits,Nos.of.Traits)		 # Initialize matrix in eqn 3 of Manly and Rayner
# Calculate Equation 3
	for (j in 1:N.matrices)	{Omega.prime <- Omega.prime+n[j]*Smatrices[,,j]/(N*c.j[j]^2)}  
# Calculate Equation 4
	d <- det(Omega.prime); if(d=="NaN"){ print(c("Determinant=0", rep))}
if (d!="NaN" && d>0)
{
for (j in 2:N.matrices)
	{
	    c.jtemp <- ginv(Omega.prime)%*%Smatrices[,,j]
			c.j[j]  <- sqrt(sum(diag(c.jtemp))/Nos.of.Traits)
	}
# LL Equation 5
   #K 			<- -0.5*N*Nos.of.Traits*log(2*pi)-0.5*N*Nos.of.Traits
   #n.log.c.j 	<- n*c.j ; print(c.j)
   #LL			<-  K - 0.5*N*log(det(Omega.prime))-Nos.of.Traits*sum(n.log.c.j[2:N.matrices])	
	 #D#iff 		<- (Old.LL-LL)^2
   #p#rint(c(Old.LL,LL,Diff))
}else
{
  G.pop <- array(0,dim=c(Nos.of.Traits,Nos.of.Traits,N.matrices))
  for (Ith.Pop in 1:N.matrices) {G.pop[,,Ith.Pop] <- c.j[Ith.Pop]^2*Omega.prime}
  return(G.pop) # return proportionality estimates.
}
}# End of while
# Get MLE of matrices
  G.pop <- array(0,dim=c(Nos.of.Traits,Nos.of.Traits,N.matrices))
  for (Ith.Pop in 1:N.matrices) {G.pop[,,Ith.Pop] <- c.j[Ith.Pop]^2*Omega.prime}
  return(G.pop) # return proportionality estimates.  Note that c^2 is the value by which matrices multiplied
} # End of function
#################################################################################
CPC.MODEL <- function(Smatrices, n.g, Nos.of.Traits, Nos.of.Populations, iter)
{
# The coding for this function is from
# Trendafilov, N. T. 2010. Stepwise estimation of common principal components.
# Computational Statistics & Data Analysis 54: 3446-3457.
# Original coding in matlab kindly supplied by Dr Trendafilov
# Smatrices: Sample matrices
# n.g : vector with sample sizes
# Nos.of.Traits: Number of traits
# N.matrices: Number of matrices
# iter: Number of iterations Start with five
  n 			<- n.g/sum(n.g)
	D			  <- matrix(0,Nos.of.Traits,Nos.of.Populations)
	CPC			<- matrix(0,Nos.of.Traits,Nos.of.Traits)
	Qw			<- diag(Nos.of.Traits)
	s			  <- matrix(0,Nos.of.Traits,Nos.of.Traits)
	for (m in 1:Nos.of.Populations){s <- s + n[m]*Smatrices[,,m]}
	Ematrix <- eigen(s)
	q0 			<- Ematrix$vectors
	d0 			<- diag(Ematrix$values)
	if (d0[1,1]< d0[Nos.of.Traits,Nos.of.Traits]) {q0 <- q0[,seq(from=Nos.of.Traits, to=1, by=-1)]}
	for (ncomp in 1:Nos.of.Traits )	# Iterate over traits to get components
{
  	q 		<- q0[,ncomp]
  	d 		<- matrix(0,1,Nos.of.Populations)
  	for (m in 1:Nos.of.Populations){ d[m] <- t(q)%*%Smatrices[,,m]%*%q }
# Find optimum number of iterations based on first pass
  for ( i in 1:iter)  # Iterate until converged
{
	s <- matrix(0,Nos.of.Traits,Nos.of.Traits)
   for ( m in 1:Nos.of.Populations){s <- s + n.g[m]*Smatrices[,,m]/d[m]}
	w <- s%*%q
	if (ncomp!=1) { w <- Qw%*%w }
	q <- w%*%solve(((t(w)%*%w)^.5))
 	for (m in 1:Nos.of.Populations){ d[m] <- t(q)%*%Smatrices[,,m]%*%q	}
} # end of i
	D[ncomp,] 		<- d # Eigenvalues
	CPC[,ncomp] 	<- q #Eigenvectors
  	Qw 				<- Qw - q%*%t(q)
} # end of ncomp
#	LL <- 0; for (i in 1:N.matrices){LL <- LL +t(CPC[,i])%*%Smatrices[,,i]%*%CPC[,i]} 
# Get MLE of matrices
  G.pop <- array(0,dim=c(Nos.of.Traits,Nos.of.Traits,Nos.of.Populations))
  for (Ith.Pop in 1:Nos.of.Populations) {G.pop[,,Ith.Pop] <- t(CPC)%*%diag(D[,Ith.Pop])%*%CPC}
  return(G.pop) # return
}
########################################################################################
CHI2.VALUE <- function(G.pop, Nos.of.Populations,N, G.est)
{
# Function to get Chi2 value 
# G.pop = Observed G matrices
# G.est = Estimated matrices (Equal, proportional, CPC)
# Note that in some cases Chi2 cannot be calculates. In this case Nans returned
# and are eliminated from the calculations.  Warnings are given but can be ignored
  Chi2   <- 0
  for(Ith.Pop in 1:Nos.of.Populations){ Chi2 <- Chi2 + 
                       sum(N*log(det(G.est[,,Ith.Pop])/det(G.pop[,,Ith.Pop])))}
  return(Chi2)
}
#################################################################################
MATRIX.TESTS <- function(DATA, Nos.of.Rands, Selected.POP, Nos.of.Traits, STAT, FORMULA, Method)
{
# Function to do hiarchical tests
# Probabilities are calculated using the "jump up" approach"
# Thus probability calculated relative to unrelated structure
# Cases in which the determinant could not be calculated were dropped from the calculations
# The warnings refer to these cases and can be ignored
# The total number of randomizations and the number actually used (N.actual) are reported
  Nos.of.Populations    <- length(Selected.POP)
# Get the unique SIRE POP pairs
  Pop <- DATA$POP;  SIRE <- DATA$SIRE
	Unique.Combinations 	<- data.frame(Pop,SIRE)
  Unique.Combinations 	<- unique.data.frame(Unique.Combinations)  # Get set of unique combinations using SPLUS/R function
# Get the Chi2 for the original populations
   Out      <- G.MATRICES(DATA,Selected.Pop,Nos.of.Traits,STAT,FORMULA,Method)
   G.pop    <- Out[[1]]
   N        <- Out[[2]]
# Equal matrices
  G.equal        <- COMMON.MATRIX(G.pop,N,Nos.of.Traits,Nos.of.Populations)
  Chi2.equal.obs <- CHI2.VALUE(G.pop, Nos.of.Populations,N, G.equal)
  #DF <- 0.5*Nos.of.Traits*(Nos.of.Traits+1)-0.5*(Nos.of.Traits*(Nos.of.Traits+1))
  #P.equal <- 1-pchisq(Chi2.equal.obs,df=DF) # p-value for stat 

# Proportional matrices
  G.prop         <- PROPORTIONAL.MATRICES(G.pop, N, Nos.of.Traits, Nos.of.Populations, Tolerance=1e-5)
Chi2.prop.obs    <- CHI2.VALUE(G.pop, Nos.of.Populations, N, G.prop)
 # DF <- 0.5*Nos.of.Traits*(Nos.of.Traits+1)-(0.5*(Nos.of.Traits*(Nos.of.Traits+1))+Nos.of.Populations-1)
 # P.prop <- 1-pchisq(Chi2.prop.obs,df=DF) # p-value for stat 

# CPC Model
  G.CPC          <- CPC.MODEL(G.pop, N, Nos.of.Traits, Nos.of.Populations, iter=10) 
  Chi2.CPC.obs   <- CHI2.VALUE(G.pop, Nos.of.Populations, N, G.CPC)
  #DF <- 0.5*Nos.of.Traits*(Nos.of.Traits+1)-(0.5*(Nos.of.Traits*(Nos.of.Traits-1))+Nos.of.Populations*Nos.of.Traits)
  #P.CPC <- 1-pchisq(Chi2.CPC.obs,df=DF) # p-value for stat
#print(c(P.equal,P.prop,P.CPC))
# Iterate over combinations of randomized populations
  Chi2.rand		   <- matrix(0,Nos.of.Rands,3)		# Initialize matrix to hlod pairs of vector correlations
  Chi2.rand[1,1] <- Chi2.equal.obs
  Chi2.rand[1,2] <- Chi2.prop.obs
  Chi2.rand[1,3] <- Chi2.CPC.obs
# Print Observed Chi2 values
#  print("********Chi^2 for the observed data********")
#  print(Chi2.rand[1,])
#  print(" ")
  for ( i in 2:Nos.of.Rands)								# Iterate over replicates
{
  NEW.DATA <- RANDOMIZE(DATA, Unique.Combinations) # Randomize data set
  Out      <- G.MATRICES(NEW.DATA,Selected.Pop,Nos.of.Traits,STAT,FORMULA,Method)
  G.pop    <- Out[[1]]
  N        <- Out[[2]]
  G.equal  <- COMMON.MATRIX(G.pop,N,Nos.of.Traits,Nos.of.Populations)
  Chi2.rand[i]   <- CHI2.VALUE(G.pop, Nos.of.Populations,N, G.equal)
# Equal matrices
  G.equal        <- COMMON.MATRIX(G.pop,N,Nos.of.Traits,Nos.of.Populations)
  Chi2.rand[i,1] <- CHI2.VALUE(G.pop, Nos.of.Populations,N, G.equal)
# Proportional matrices
  G.prop         <- PROPORTIONAL.MATRICES(G.pop, N, Nos.of.Traits, Nos.of.Populations, Tolerance=1e-5)
  Chi2.rand[i,2] <- CHI2.VALUE(G.pop, Nos.of.Populations,N, G.prop)
# CPC Model
  G.CPC          <- CPC.MODEL(G.pop, N, Nos.of.Traits, Nos.of.Populations, iter=10) 
  Chi2.rand[i,3] <- CHI2.VALUE(G.pop, Nos.of.Populations,N, G.CPC)
} # End of i loop
# Now eliminate runs with NaN
  Out      <- matrix(0,5,1)
  Out[1]   <- Nos.of.Rands  # Maximum number of replicates
  k <- 1
for ( i in 1:3)
{                             
  X        <- Chi2.rand[,i]; X <- X[X!="NaN"]      # Exclude NaN
  N.reps   <- length(X)  
  k <- k+1
  Out[k]   <- N.reps                               # Number of replicates used
  k <-k+1
  Out[k] <- length(X[X>=Chi2.rand[1,i]])/N.reps	# Calculate probability
} # End of i loop

  print("*************************************************************************")
  print("Probabilities are calculated using the jump up approach")
  print("Thus the probability is calculated relative to unrelated structure")
  print("Cases in which the determinant could not be calculated are dropped from the calculations")
  print("The warnings refer to these cases and can be ignored")
  print("The total number of randomizations and the number actually used (N.actual) are reported")
  print("*************************************************************************")
  print("Results of randomization")
  print("Number of randomizations (includes observed)")
  print(Out[1])
  print("Comparing Equal with unrelated")
  print("N.actual   Prob "); print(Out[2:3])
  print("Comparing Proportional with unrelated")
  print("N.actual   Prob "); print(Out[4:5])
  print("Comparing CPC with unrelated")
  print("N.actual   Prob "); print(Out[6:7])
    return(Out)  #Return 7 items
}