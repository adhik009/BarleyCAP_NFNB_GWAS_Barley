rm(list =ls())

library(rrBLUP)
geno <- read.table("WSU_genotype.hmp.txt",header=TRUE,as.is=TRUE,check.names=FALSE)
M <- t(geno[,-c(1:4)])  #makes a marker file
hist(M) # Show that Genotypes were Downloaded correctly 
dim(M) 
m <- ncol(M)  #number of markers
n <- nrow(M)  #number of lines
W <- scale(M,scale=FALSE)  #centers each marker by its mean
W[is.na(W)] <- 0  	   #imputing missing values with pop mean
cov.W <- cov(t(W))  

###reads in your pheno
trait_data<-read.table("WSU_pheno.txt",header=TRUE,as.is=TRUE,check.names=FALSE)

A1 <- A.mat(M,shrink=TRUE) #Making the relationship matrix A1

#Create a dataframe with traitdata and genotype identifiers
env <- trait_data$Env
data2 <- data.frame(y=trait_data$Mean,env=trait_data$Env,gid=trait_data$Line) 


##### CHECK IF THE CALCULATED MARKER EFFECTS ARE THE SAME USING MIXED.SOLVE and KIN.BLUP
#__________________________________________________________________________________


#calculate the marker effects from kin.blup function
ans <- kin.blup(data2,K=A1,geno="gid",pheno="y",fixed="env")
write.csv(ans,file="WSU_BLUP.txt")

#Calculate the marker effects using results from mixed.solve
mkr_effects<-mixed.solve(ans$g,Z=W,SE=F,return.Hinv=F)
write.csv(cbind(geno[,1:4], mkr_effects), file="WSU_mkr_effects.csv")

#-------------- GUIDE FOR LOOKING AT RESULTS FROM THE TWO MIXED MODEL SOLUTIONS ------------

##Results from kin.blup###
#$Vg REML estimate of the genetic variance
#$Ve REML estimate of the error variance
#$g BLUP solution for the genetic values
#$resid residuals
#$pred predicted genetic values, averaged over the fixed effects

### Results from mixed.solve###
#$Vu estimator for σ2u
#$Ve estimator for σ2e
#$beta BLUE(β)
#$u BLUP(u)
#$LL maximized log-likelihood (full or restricted, depending on method)

