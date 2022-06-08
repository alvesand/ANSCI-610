#-------------------------------------------------------------------------------------------#
setwd("C:\\Users\\Badger\\Desktop\\Class GS")
#-------------------------------------------------------------------------------------------#
# Loading packages
#-------------------------------------------------------------------------------------------#
#install.packages("gaston")
#install.packages("BGLR")
#install.packages("nadiv")
#install.packages("MCMCglmm")
#install.packages("doBy")
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("snpStats")

library("gaston")
library("BGLR")
library("snpStats")
library("nadiv")
library("MCMCglmm")

#-------------------------------------------------------------------------------------------#
# Importing data
#-------------------------------------------------------------------------------------------#
PHE <- read.table("pheno.txt",  header = T, sep = ",")
PHE$Animal <- as.factor(PHE$Animal)
PHE$Sex <- as.factor(PHE$Sex)
PHE$Phe <- as.numeric(PHE$Phe)
rownames(PHE) = PHE$Animal
boxplot(PHE$Phe~PHE$Sex)

#-------------------------------------------------------------------------------------------#
Pedig <- read.table("PED.txt",  header = T, sep = ",")
Pedig$Animal <- as.factor(Pedig$Animal)
Pedig$Sire <- as.factor(Pedig$Sire)
Pedig$Dam <- as.factor(Pedig$Dam)
#-------------------------------------------------------------------------------------------#
M <- as.matrix(read.table("GEN.txt", header = T, sep = ","))
rownames(M) <- PHE$Animal[as.numeric(levels(PHE$Animal))>110000]
dim(M)
#-------------------------------------------------------------------------------------------#
# Quality Control of SNPs
#-------------------------------------------------------------------------------------------#
markers = new("SnpMatrix", (M+1))
snpstat<-col.summary(markers)

maf = snpstat$MAF #Minnor allele frequency for the markers
call.rate = snpstat$Call.rate #Call rate  for the markers
hwe.p = 1-pnorm(abs(snpstat$z.HWE)) #HWE test p-values
p = snpstat$RAF #Frequency of the second allele

snps.ok<-which(maf>0.05&call.rate>0.9&hwe.p>1e-5)
M = M[,snps.ok]
p = p[snps.ok]
rm(markers)

#-------------------------------------------------------------------------------------------#
# Getting the Pseudo Phenotype
#-------------------------------------------------------------------------------------------#
mod = lm(Phe~Sex, data = PHE)
summary(mod)
PHE$Adj_Y = mod$residuals

Ped = prepPed(Pedig)
colnames(Ped) = c("animal","sire","dam")
colnames(PHE) = c("animal","Sex", "Phe", "Adj_Y")

flat_prior <- list(R = list(V = 1, nu = 0), G = list(G1 = list(V = 1, nu = 0)))
mod2 = MCMCglmm(Phe~Sex, random =~ animal, data = PHE, pedigree = Ped, verbose=T,
                nitt=10000, thin=2, burnin=2000, pr=T, family = "gaussian", prior =flat_prior)

h2 = mod2$VCV[,1]/(mod2$VCV[,1]+mod2$VCV[,2])

par(mfrow=c(1,2))
plot(as.numeric(h2), type = "l", ylab = "h2 values", main = "Trace plot")
plot(density(h2), main = "h2 marginal distribution")
polygon(density(h2), col = rgb(0.2,0.4,0.8,0.6))
abline(v = mean(h2), col = "red", lty = 2)

head(mod2$Sol[,1:10])
blup = data.frame("animal" = rep(0,nrow(Ped)),
                  "ebv" = rep(0,nrow(Ped)), "sd" = rep(0,nrow(Ped)))

for(i in 1:(ncol(mod2$Sol)-2)){
  blup[i,1] = Ped$animal[i]
  blup[i,2] = mean(mod2$Sol[,i+2])
  blup[i,3] = sd(mod2$Sol[,i+2])
}

PHE = merge(PHE, blup, by = "animal")
head(PHE)
cor(PHE$Phe, PHE$ebv)
cor(PHE$Phe, PHE$Adj_Y)

ind = intersect(PHE$animal, rownames(M))
rownames(PHE) = PHE$animal
head(ind)

PHE_G <- PHE[ind,]
M = M[ind,]

head(PHE_G)
head(rownames(M))

#-------------------------------------------------------------------------------------------#
# Computation of Genomic Relationship Matrix (Van Raden 2008)
#-------------------------------------------------------------------------------------------#
P = matrix(2*p, nrow=1000, ncol=length(p), byrow=TRUE)
Z = M - P
scal = 2*sum(p*(1-p))
G <- tcrossprod(Z)/scal
rm(Z, P)
par(mfrow=c(1,1))
hist(G[upper.tri(G, diag = F)], breaks = 30, col = "palegreen", main = "")

#-------------------------------------------------------------------------------------------#
# Numerator Relationship Matrix for the genotyped animals
#-------------------------------------------------------------------------------------------#
rownames(Pedig) = Pedig$Animal
ind = intersect(Pedig$Animal, rownames(M)) 
Ped = Pedig[ind,]
Ped = prepPed(Ped)
A = makeA(ped=Ped)
A = A[ind,ind]
A = as.matrix(A)
table(A[upper.tri(A, diag = F)])

#-------------------------------------------------------------------------------------------#
# Design Matrix for Fixed Effects
#-------------------------------------------------------------------------------------------#
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

X <- incidence(PHE_G$Sex)
head(X)

#-------------------------------------------------------------------------------------------#
# Pedigree-Based BLUP using only genotyped animals - (gaston, AI-REML)
#-------------------------------------------------------------------------------------------#
A_BLUP <- lmm.aireml(Y = PHE_G$Phe,
                     X = X,
                     K = A, 
                     verbose = FALSE)
#-------------------------------------------------------------------------------------------#
# Genomic BLUP - (gaston, AI-REML)
#-------------------------------------------------------------------------------------------#
G_BLUP <- lmm.aireml(Y = PHE_G$Phe,
                     X = X,
                     K = G, 
                     verbose = FALSE)
#-------------------------------------------------------------------------------------------#
# Bayesian GBLUP - RKHS (BGLR)
#-------------------------------------------------------------------------------------------#
Bayes_G <- BGLR(y=PHE_G$Phe, # y is vector with the phenotype of interest
                ETA=list(list( ~ Sex, data=PHE_G, model="FIXED"),
                         list(K=G, model='RKHS')), # ETA is a list containing the prior distributions
                nIter=2000, # The chain size
                burnIn=1000, # The number of iterations that will be discard
                thin=2, # Sampling every 2 iterations
                verbose=TRUE)
#-------------------------------------------------------------------------------------------#
# Bayesian Ridge Regression (gaussian prior)
#-------------------------------------------------------------------------------------------#
Bayes_R <- BGLR(y=PHE_G$Phe,
                ETA=list(list( ~ Sex, data=PHE_G, model="FIXED"),
                         list(X=M, model='BRR')),
                nIter=2000,
                burnIn=1000,
                thin=2,
                verbose=TRUE)

#-------------------------------------------------------------------------------------------#
# Bayes A (scaled-t)
#-------------------------------------------------------------------------------------------#
Bayes_A <- BGLR(y=PHE_G$Phe,
                ETA=list(list( ~ Sex, data=PHE_G, model="FIXED"),
                         list(X=M, model='BayesA')),
                nIter=2000,
                burnIn=1000,
                thin=2,
                verbose=TRUE)

#-------------------------------------------------------------------------------------------#
# Bayes B (Point of mass at zero + scaled-t slab)
#-------------------------------------------------------------------------------------------#
Bayes_B <- BGLR(y=PHE_G$Phe,
                ETA=list(list( ~ Sex, data=PHE_G, model="FIXED"),
                         list(X=M, model='BayesB')),
                nIter=2000,
                burnIn=1000,
                thin=2,
                verbose=TRUE)

#-------------------------------------------------------------------------------------------#
# Bayes C pi (Gaussian mixture)
#-------------------------------------------------------------------------------------------#
Bayes_C <- BGLR(y=PHE_G$Phe,
                ETA=list(list( ~ Sex, data=PHE_G, model="FIXED"),
                         list(X=M, model='BayesC')),
                nIter=2000,
                burnIn=1000,
                thin=2,
                verbose=TRUE)

#-------------------------------------------------------------------------------------------#
# Bayesian Lasso (double exponential)
#-------------------------------------------------------------------------------------------#
Bayes_L <- BGLR(y=PHE_G$Phe,
                ETA=list(list( ~ Sex, data=PHE_G, model="FIXED"),
                         list(X=M, model='BL')),
                nIter=2000,
                burnIn=1000,
                thin=2,
                verbose=TRUE)

#-------------------------------------------------------------------------------------------#
# Comparing the marker effects for all models
#-------------------------------------------------------------------------------------------#
betaRR = Bayes_R$ETA[[2]]$b
betaBA = Bayes_A$ETA[[2]]$b
betaBB = Bayes_B$ETA[[2]]$b
betaBC = Bayes_C$ETA[[2]]$b
betaBL = Bayes_L$ETA[[2]]$b

par(mfrow=c(3,2))
plot(betaRR^2, type = "o", cex = 0.7, col = rgb(0.5,0.5,0.5,0.6), ylab = "squared marker effects",  main = "Bayesian Ridge Regression")
plot(betaBA^2, type = "o", cex = 0.7, col = rgb(0.1,0.3,0.7,0.6), ylab = "squared marker effects",  main = "Bayes A - Scaled t")
plot(betaBB^2, type = "o", cex = 0.7, col = rgb(0.5,0.7,0.1,0.6), ylab = "squared marker effects",  main = "Bayes B - Scaled t mixture")
plot(betaBC^2, type = "o", cex = 0.7, col = rgb(0.7,0.1,0.5,0.6), ylab = "squared marker effects",  main = "Bayes C - Gaussian mixture")
plot(betaBL^2, type = "o", cex = 0.7, col = rgb(0.4,0.1,0.4,0.6), ylab = "squared marker effects",  main = "Bayesian LASSO")

#-------------------------------------------------------------------------------------------#
# Correlation between GEBVs and Phenotype
#-------------------------------------------------------------------------------------------#
GEBVs <- data.frame(Phen=PHE_G$Phe,
                    A_BLUP  = A_BLUP$BLUP_omega,
                    G_BLUP  = G_BLUP$BLUP_omega,
                    Bayes_G = Bayes_G$yHat,
                    Bayes_R = Bayes_R$yHat,
                    Bayes_L = Bayes_L$yHat,
                    Bayes_A = Bayes_A$yHat,
                    Bayes_B = Bayes_B$yHat,
                    Bayes_C = Bayes_C$yHat)
                    
round(cor(GEBVs),2)

#-------------------------------------------------------------------------------------------#
# Running Genomic Prediction using k-fold (k = 5) cross-validation method
#-------------------------------------------------------------------------------------------#

#Bayesian Ridge Regression (Adjusted Phenotype)
folds   <- 5
y       <- PHE_G$Adj_Y
sets    <- sample(5,nrow(PHE_G),replace = T)
mse_BRR <- numeric()
cor_BRR  <- numeric()

for(fold in 1:folds){
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  fm<-BGLR(y=yNa,
           ETA=list(list(X=M, model='BRR')),
           nIter=2000,
           burnIn=1000,
           thin=2,
           verbose=FALSE)
  cor_BRR[fold] <- cor(fm$yHat[fm$whichNa], y[whichNa])
  mse_BRR[fold] <- mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  cat("Fold: ", fold, "\n")
}

#Bayes BLUP A (Adjusted Phenotype)
mse_blupA <- numeric()
cor_blupA  <- numeric()
for(fold in 1:folds){
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  fm<-BGLR(y=yNa,
           ETA=list(list(K=A, model='RKHS')),
           nIter=2000,
           burnIn=1000,
           thin=2,
           verbose=FALSE)
  cor_blupA[fold] <- cor(fm$yHat[fm$whichNa], y[whichNa])
  mse_blupA[fold] <- mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  cat("Fold: ", fold, "\n")
}

#Bayesian GBLUP (Adjusted Phenotype)
mse_blupG <- numeric()
cor_blupG  <- numeric()
for(fold in 1:folds){
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  fm<-BGLR(y=yNa,
           ETA=list(list(K=G, model='RKHS')),
           nIter=2000,
           burnIn=1000,
           thin=2,
           verbose=FALSE)
  cor_blupG[fold] <- cor(fm$yHat[fm$whichNa], y[whichNa])
  mse_blupG[fold] <- mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  cat("Fold: ", fold, "\n")
}

#Bayes A (Adjusted Phenotype)
mse_BA <- numeric()
cor_BA  <- numeric()
for(fold in 1:folds){
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  fm<-BGLR(y=yNa,
           ETA=list(list(X=M, model='BayesA')),
           nIter=2000,
           burnIn=1000,
           thin=2,
           verbose=FALSE)
  cor_BA[fold] <- cor(fm$yHat[fm$whichNa], y[whichNa])
  mse_BA[fold] <- mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  cat("Fold: ", fold, "\n")
}

#Bayes B (Adjusted Phenotype)
mse_BB <- numeric()
cor_BB  <- numeric()
for(fold in 1:folds){
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  fm<-BGLR(y=yNa,
           ETA=list(list(X=M, model='BayesB')),
           nIter=2000,
           burnIn=1000,
           thin=2,
           verbose=FALSE)
  cor_BB[fold] <- cor(fm$yHat[fm$whichNa], y[whichNa])
  mse_BB[fold] <- mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  cat("Fold: ", fold, "\n")
}


#Bayes LASSO (Adjusted Phenotype)
mse_BL <- numeric()
cor_BL  <- numeric()
for(fold in 1:folds){
  yNa<-y
  whichNa<-which(sets==fold)
  yNa[whichNa]<-NA
  fm<-BGLR(y=yNa,
           ETA=list(list(X=M, model='BL')),
           nIter=2000,
           burnIn=1000,
           thin=2,
           verbose=FALSE)
  cor_BL[fold] <- cor(fm$yHat[fm$whichNa], y[whichNa])
  mse_BL[fold] <- mean((fm$yHat[fm$whichNa]-y[whichNa])^2)
  cat("Fold: ", fold, "\n")
}

#-------------------------------------------------------------------------------------------#
# Basic descriptive statistics Genomic Prediction Accuracy and Mean Square Error
#-------------------------------------------------------------------------------------------#
result = data.frame(cor = c(cor_blupA, cor_blupG, cor_BRR, cor_BA, cor_BB,
                           cor_BL), mse = c(mse_blupA, mse_blupG, mse_BRR,
                                            mse_BA, mse_BB, mse_BL),
                   method = c(rep("BLUP_A", 5),rep("BLUP_G", 5),rep("BRR", 5),rep("Bayes_A", 5),
                              rep("Bayes_B", 5),rep("BLASSO", 5)),
                   fold = rep(seq(1,5,1),6), 
                   names = c(rep("A", 5),rep("B", 5),rep("C", 5),rep("D", 5),
                   rep("E", 5),rep("F", 5))) 
head(result)

library(doBy)
summaryBy(cor~method, data = result, FUN = c(mean,sd))
summaryBy(mse~method, data = result, FUN = c(mean,sd))

par(mfrow=(c(2,1)))
boxplot(cor~names, data=result,
        col=(c("gold","darkgreen", "steelblue", "darkgray", "red", "violet")),
        main="Prediction Correlation", xlab="Method", names=c("BLUP_A", "BLUP_G", "BRR",
                                                              "Bayes_A", "Bayes_B", "BLASSO"))

boxplot(mse~names, data=result,
        col=(c("gold","darkgreen", "steelblue", "darkgray", "red", "violet")),
        main="Mean squared error", xlab="Method", names=c("BLUP_A", "BLUP_G", "BRR",
                                                          "Bayes_A", "Bayes_B", "BLASSO"))

#-------------------------------------------------------------------------------------------#




