setwd("/Users/aniladhikari/Google Drive/TCAP/2-row")  #Set WD

library(rrBLUP)
geno <- read.table("2row_genotype.hmp.txt",header=TRUE,as.is=TRUE,check.names=FALSE)
str(geno) #Loaded the Gneotypes
head(geno)
M <- t(geno[,-c(1:4)])  #makes a marker file
hist(M) # Show that Genotypes were Dowlonaded correctly 
dim(M) #shows that we got all of the 
m <- ncol(M)  #number of markers
n <- nrow(M)  #number of lines

pheno <- read.table("Pheno_2row_new.txt",header=TRUE,as.is=TRUE,check.names=FALSE)
pheno[order(pheno$Mean),] #sorting by phenotypic values
hist(pheno$Mean)

#look at population structure
W <- scale(M,scale=FALSE)  #centers each marker by its mean
W[is.na(W)] <- 0  	   #imputing missing values with pop mean
cov.W <- cov(t(W))         #n x n covariance matrix
eig.result <- eigen(cov.W)
eig.vec <- eig.result$vectors
lambda <- eig.result$values
plot(lambda/sum(lambda),ylab="Fraction of Total Variance")
plot(eig.vec[,1],eig.vec[,2],xlab="PC1",ylab="PC2")


#P+K or K only model
#look at inflation
dev.new()
dev2=dev.cur()
dev.new()
dev.set(dev2)
P.K<-GWAS(pheno, geno[,-2], fixed="Env", K=NULL, n.PC=0, min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)

dev.new()
#look at GWAS hit
hit <- which.max(P.K[,4])
hit.region <- P.K[(hit-40):(hit+5),]
plot(hit.region[,3]/200,hit.region[,4],xlab="6H position (cM)",ylab="-log(p)")
P.K[hit,]
write.table(P.K,"emma.result_2row.txt", sep="\t")
