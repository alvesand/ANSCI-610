#----------------------- ANSCI 610 L13 - Practice 3 ----------------------------#
#---------------------- Mixed model analysis with R ----------------------------#
#-------------------------------------------------------------------------------#
#Packages required for this class 
#install.packages("rrBLUP")
#install.packages("EMMREML")
#install.packages("pedigree")

#A simple function for creating incidence matrices
incidence = function(x){
  temp = data.frame(table(x))
  temp = subset(temp, Freq >0)
  n_class = dim(temp)[1]
  inc = matrix(0, length(x), n_class)
  colnames(inc) = temp[,1]
  for (i in 1:length(x)){
    j = which(x[i]==temp[,1])
    inc[i,j] = 1
  }
  return(inc)
}

#---------------------------------------------------------------------------#
#-------------- Create the A inverse ---------------------------------------#
#---------------------------------------------------------------------------#
library(pedigree)
Ai = makeAinv(ped=ysim[,1:3])
Ai <- read.table("Ainv.txt")
nInd <- nrow(ysim[,1:3])
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd
dim(Ainv)
colnames(Ainv) = rownames(Ainv) <- ysim$animal

ysim2 = subset(ysim, y>0) #Maintain only animals with phenotype

#Getting an idea for starting points within the parameter space
summary(lm(y~as.factor(sex)+feed,data=ysim[which(ysim$y>0),]))
guess1 = (28.22^2)/2

#Building the X matrix
X=cbind(incidence(ysim2$sex),ysim2$feed)
head(X)

#Building the Z matrix
Z=incidence(ysim$animal)
Z = Z[-which(ysim$y==0),] #delete the animals without information for Z
head(Z[1:10,1:10])
head(Z[1:10,121:131])
dim(Z)
#check if A and Z are conformable
dim(Z) == dim(Ainv)

#---------------------------------------------------------------------------------#
# Using our function for VCE using the REML method and the EM algorithm
#---------------------------------------------------------------------------------#
mme = emreml610(1000,guess = c(guess1,guess1),criter=10^-5,Ainv=Ainv,X=X,Z=Z,y=ysim2$y)

#Extract the solutions
blue = mme$blue
blup = mme$blup
mme$G
mme$R
 
mme$G/(mme$G+mme$R) #Compute the estimated heritability

#Check the correlation with true breeding values
par(mfrow =c(1,1))
plot(blup$ebv, ysim$tbv, xlab = "EBV", ylab= "TBV", main = "", pch = 16, col = rgb(0.5,0.5,0.9,0.6))
cor(blup$ebv, ysim$tbv)
mean(blup$acc)

#---------------------------------------------------------------------------------#
# Using EMMREML for VCE using the REML method and the EMMA algorithm
#---------------------------------------------------------------------------------#
library(EMMREML)
mme2 = emmreml(y=ysim2$y, X=X, Z=Z, K=solve(Ainv),varbetahat=T,varuhat=T, PEVuhat=T, test=FALSE)
mme2$Vu
mme2$Ve
blue2 = mme2$betahat
blup2 = mme2$uhat
head(ysim)

#---------------------------------------------------------------------------------
#Using rrBlup for VCE with the Maximum Likelihood (ML) method
#---------------------------------------------------------------------------------
library(rrBLUP)
mme3 = mixed.solve(y=ysim2$y, Z=Z, K=solve(Ainv), X=X, method="ML",
                   bounds=c(1e-09, 1e+09), SE=T, return.Hinv=FALSE)

mme3$Vu
mme3$Ve
blup3 = mme3$u
blue3 = mme3$beta
  
#---------------------------------------------------------------------------------#
# Comparing all solutions
#---------------------------------------------------------------------------------#
results = data.frame(blup1 = blup$ebv, blup2 = blup2, blup3 = blup3)

cor(results)

