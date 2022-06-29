################################################################################
# STEP 1 LOADING RAW DATA

# Clean workspace
rm(list=ls())

# get current working directory
getwd()

# Set desired working directory
setwd("C:/Users/edoar/Desktop/DNA_RNA_Project/working_dir")

# Load minfi package
library(minfi)

# Load Illumina manifest
load("C:/Users/edoar/Desktop/DNA_RNA_Project/working_dir/Illumina450Manifest_clean.RData")

# Set input data directory
baseDir <- ("C:/Users/edoar/Desktop/DNA_RNA_Project/input_data_report")

# Reading sheet containing pheno-data
targets <- read.metharray.sheet(baseDir)

# Reading methylation array data
RGset <- read.metharray.exp(targets = targets)

# Checking output
RGset


################################################################################
# STEP 2 STORING RED AND GREEN FLUORESCENCE

# Red fluorescence
Red <- data.frame(getRed(RGset))
# Check df
head(Red)

# Green fluorescence
Green <- data.frame(getGreen(RGset))
# Check df
head(Green)


################################################################################
# STEP 3 INSPECTING ADDRESS AND DESIGN TYPE

# Red fluorescence dataframe inspection at "42639338"
Red[rownames(Red)=="42639338",]

# Red fluorescence dataframe inspection at "42639338"
Green[rownames(Green)=="42639338",]

# Store the address annotation
Addr42639338 <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="42639338",]

# Check the Infinium design type columns
Addr42639338$Infinium_Design_Type


################################################################################
# STEP 4 CONVERTING FLUORESCENCE INTENSITY INTO METHYLATION SIGNAL

# Create MSet.raw object from the RGset object
MSet.raw <- preprocessRaw(RGset)

# Check output
MSet.raw


################################################################################
# STEP 5 QUALITY CHECKS

# QCPLOT
qc <- getQC(MSet.raw)
plotQC(qc)

# NEGATIVE CONTROLS
controlStripPlot(RGset, controls="NEGATIVE")

# PVALUE DETECTION
detP <- detectionP(RGset)

failed <- detP>0.05
summary(failed)


################################################################################
# STEP 6 COMPUTING BETA AND M VALUES

# Store separately data from WT samples and MUT samples
WT <- targets[targets$Group=='WT',"Array"]
MUT <- targets[targets$Group=='MUT',"Array"]

## BETA VALUES
# Compute beta values
beta <- getBeta(MSet.raw)
WT_beta <- beta[, targets$Array %in% WT]
MUT_beta <- beta[, targets$Array %in% MUT]

# Compute beta values mean and density
WT_beta_mean <- apply(WT_beta,1,mean)
MUT_beta_mean <- apply(MUT_beta,1,mean)
WT_beta_d <- density(WT_beta_mean,na.rm=T)
MUT_beta_d <- density(MUT_beta_mean,na.rm=T)

# Plot beta values density
plot(WT_beta_d, col='red',main='Mean beta values density')
lines(MUT_beta_d,col='blue')
legend('topright', pch=c(15,15), col=c('red','blue'), legend=c('WT','MUT'))

## M VALUES
# Compute M values 
M <- getM(MSet.raw)
WT_M <- M[, targets$Array %in% WT]
MUT_M <- M[, targets$Array %in% MUT]

# Compute M values mean and density 
WT_M_mean <- apply(WT_M,1,mean)
MUT_M_mean <- apply(MUT_M,1,mean)
WT_M_d <- density(WT_M_mean,na.rm=T)
MUT_M_d <- density(MUT_M_mean,na.rm=T)

# Plot M values density
plot(WT_M_d, col='red',main='Mean M values density')
lines(MUT_M_d,col='blue')
legend('topright', pch=c(15,15), col=c('red','blue'), legend=c('WT','MUT'))


################################################################################
# STEP 7. NORMALIZATION

# Raw values (already computed in the previous step)
beta
# Normalization
preprocessQuantile_results <- preprocessQuantile(RGset)
# Normalized beta
normalized_beta <- getBeta(preprocessQuantile_results)

### RAW
## INFINIUM I DESIGN RAW
# Extract Infinium I annotated IDs
InfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
InfI_IDs <- droplevels(InfI)

# Cross-reference between InfI_IDs and raw beta
beta_I <- beta[rownames(beta) %in% InfI_IDs$IlmnID,]

# Computing mean and standard deviation
beta_I_mean <- apply(beta_I,1, mean)
beta_I_sd <- apply(beta_I, 1, sd)

# Computing densities for plots
beta_I_mean_d <- density(beta_I_mean, na.rm=T)
beta_I_sd_d <- density(beta_I_sd, na.rm=T)

## INFINIUM II DESIGN RAW
# Extract Infinium I annotated IDs
InfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
InfII_IDs <- droplevels(InfII)

# Cross-reference between InfII_IDs and raw beta
beta_II <- beta[rownames(beta) %in% InfII_IDs$IlmnID,]

# Computing mean and standard deviation
beta_II_mean <- apply(beta_II,1, mean)
beta_II_sd <- apply(beta_II, 1, sd)

# Computing densities for plots
beta_II_mean_d <- density(beta_II_mean, na.rm=T)
beta_II_sd_d <- density(beta_II_sd, na.rm=T)

### NORMALIZED
## INFINIUM I DESIGN NORMALIZED
# Cross-reference between InfI_IDs and raw beta
normal_beta_I <- normalized_beta[rownames(normalized_beta) %in% InfI_IDs$IlmnID,]

# Computing mean and standard deviation
normal_beta_I_mean <- apply(normal_beta_I,1, mean)
normal_beta_I_sd <- apply(normal_beta_I, 1, sd)

# Computing densities for plots
normal_beta_I_mean_d <- density(normal_beta_I_mean, na.rm=T)
normal_beta_I_sd_d <- density(normal_beta_I_sd, na.rm=T)

## INFINIUM II DESIGN NORMALIZED
# Cross-reference between InfII_IDs and raw beta
normal_beta_II <- normalized_beta[rownames(normalized_beta) %in% InfII_IDs$IlmnID,]

# Computing mean and standard deviation
normal_beta_II_mean <- apply(normal_beta_II,1, mean)
normal_beta_II_sd <- apply(normal_beta_II, 1, sd)

# Computing densities for plots
normal_beta_II_mean_d <- density(normal_beta_II_mean, na.rm=T)
normal_beta_II_sd_d <- density(normal_beta_II_sd, na.rm=T)

### PLOTS

# Raw beta mean plot
par(mfrow=c(3,2), pin = c(2,2))
plot(beta_I_mean_d, col = "green", main = "Mean Beta (raw)", xlim=c(-0.1,1.1),ylim = c(-0.5,6))
lines(beta_II_mean_d, col = "orange")
legend("topright", legend=c("InfI", "InfII"),
       col=c("green", "orange"), lty=1:1, cex=0.9)

# Normalized beta mean plot
plot(normal_beta_I_mean_d,col="blue",main=" Mean Beta (normalized)", xlim = c(-0.1,1.1), ylim = c(-0.5
                                                                                                  ,6))
lines(normal_beta_II_mean_d,col="red")
legend("topright", legend=c("InfI", "InfII"),
       col=c("blue", "red"), lty=1:1, cex=0.9)

# Raw beta sd plot
plot(beta_I_sd_d, col = "green", main = "St.Dev Beta (raw)", xlim = c(-0.05,0.45), ylim = c(-5,90))
lines(beta_II_sd_d, col = "orange")
legend("topright", legend=c("InfI", "InfII"),
       col=c("green", "orange"), lty=1:1, cex=0.9)

# Normalized beta sd plot
plot(normal_beta_I_sd_d,col="blue",main="St.Dev Beta (normalized)", xlim = c(-0.05,0.45), ylim = c(-5,90))
lines(normal_beta_II_sd_d,col="red")
legend("topright", legend=c("InfI", "InfII"),
       col=c("blue", "red"), lty=1:1, cex=0.9)

# Raw beta boxplot
boxplot(beta,main="Raw beta")
# Normalized beta boxplot
boxplot(normalized_beta, main = "preprocessQuantile normalized beta")

## OPTIONAL:

# WT slide arrays
WT_group <- colnames(WT_beta)

# MUT slide arrays
MUT_group <- colnames(MUT_beta)

# Vector color-mapping beta-values object columns to WT arrays and MUT arrays
colors = c("pink", "cadetblue1", "pink", "pink", "cadetblue1", "cadetblue1", "pink", "cadetblue1")

# Raw data Plot
boxplot(beta, col=colors, main="WT and MUT raw beta values")
legend("topright", inset=c(-0.4,0), legend=c("MUT", "WT"),
       col=c("cadetblue1", "pink"), lty=1:1, cex=0.9, xpd=TRUE)

# Normalized data Plot
boxplot(normalized_beta, col=colors, main="WT and MUT normalized beta values")
legend("topright", inset=c(-0.4,0), legend=c("MUT", "WT"),
       col=c("cadetblue1", "pink"), lty=1:1, cex=0.9, xpd=TRUE)


################################################################################
# STEP 8 PCA - PRINCIPAL COMPONENT ANALYSIS

# Perform PCA and check results
pca_results <- prcomp(t(normalized_beta),scale=T)
print(summary(pca_results))

# Extract and plot variances from PCA results
library(factoextra)
fviz_eig(pca_results, addlabels = TRUE, barfill = "grey", barcolor = "black")

# Plot first two components
group <- factor(targets$Group)
palette(c("red", "blue"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=19, col=c(group),xlab="PC1",ylab="PC2", xlim = c(-700, 700), ylim = c(-700, 700))
text(pca_results$x[,1], pca_results$x[,2],labels=rownames(pca_results$x),cex=0.6,pos=2)
legend("topright",legend=c("DS","WT"),col=c("red","blue"),pch=19)


################################################################################
# STEP 9 IDENTIFYING DIFFERENTIALLY METHYLATED PROBEDS

# Define a function to apply t.test() on each row
t_test_function <- function(x) {
  t_test <- t.test(x ~ targets$Group)
  return(t_test$p.value)
}

# Compute p-values
pValues_t <- apply(normalized_beta, 1, t_test_function)
beta_pValues <- data.frame(normalized_beta, pValues_t)


################################################################################
# STEP 10 MULTIPLE TEST CORRECTION

# Significant results without correction
raw_pValues <- beta_pValues[beta_pValues$pValues_t<=0.05,]

# Applying Bonferroni
Bonferroni_pValues <- p.adjust(beta_pValues$pValues_t, "bonferroni")

# Applying BH
BH_pValues <- p.adjust(beta_pValues$pValues_t, "BH")

# Compare
ttest_final <- data.frame(pValues_t, BH_pValues, Bonferroni_pValues)
dim(raw_pValues)
dim(ttest_final[ttest_final$Bonferroni_pValues<=0.05,])
dim(ttest_final[ttest_final$BH_pValues<=0.05,])

# Plot t-test final results
par(mfrow=c(1,1))
boxplot(ttest_final, ylim = c(-0.1, 1.1), col = c("red3", "seagreen3", "seashell2"))
legend("topright", inset=c(-0.7,0), legend=c("raw", "BH", "Bonferroni"),
       col=c("red3", "seagreen3", "seashell2"), lty=1:1, cex=0.9, xpd=TRUE)


################################################################################
# STEP 11 VOLCANO PLOT AND MANHATTAN PLOT

## VOLCANO PLOT
# WT group
WT_group <- beta_pValues[,targets$Group=="WT"]
WT_group_mean <- apply(WT_group, 1, mean)

# MUT group
MUT_group <- beta_pValues[,targets$Group=="MUT"]
MUT_group_mean <- apply(MUT_group, 1, mean)

# Compute delta
delta <- WT_group_mean - MUT_group_mean

# t-test p-values dataframe and highlight delta > 0.1
VolcPlot <- data.frame(delta, -log10(ttest_final$pValues_t))
Highlight <- VolcPlot[abs(VolcPlot[,1])>0.1 & VolcPlot[, 2]>(-log10(0.05)),]

# Plot
plot(VolcPlot[,1], VolcPlot[,2], pch=16, cex=0.4, ylab="p-Value (-log)", xlab="Delta",ylim = c(0, 8)) -log10(0.05)
abline(a=-log10(0.05),b=0,col="yellow")
points(Highlight[,1], Highlight[,2],pch=16,cex=0.4,col="cyan")

## MANHATTAN PLOT
# Load qqman package
library(qqman)

# Annotate p-value to probe genomic location
ttest_final_manh <- data.frame(rownames(ttest_final), ttest_final)
colnames(ttest_final_manh)
colnames(ttest_final_manh)[1] <- "IlmnID"
colnames(ttest_final_manh)
ttest_final_annotated <- merge(ttest_final_manh, Illumina450Manifest_clean,by="IlmnID")

# Curate format of input
input_Manhattan <- ttest_final_annotated[colnames(ttest_final_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_t")]
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)

# Plot
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_t",col=rainbow(24))


################################################################################
# STEP 12. HEATMAP

# Load gplots
library(gplots)

# Create input matrix for the heatmap
beta_pValues_ordered <- beta_pValues[order(beta_pValues$pValues_t),]
input_heatmap=as.matrix(beta_pValues_ordered[1:100,1:8])

# Vector color-mapping to WT arrays and MUT arrays according to group membership
colorbar <- c("blue","orange","blue","blue","orange","orange","blue","orange")
palette <- colorRampPalette(c("green","black","red"))(100)

# Complete linkage
heatmap.2(input_heatmap,col=palette,Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")

# Single linkage
heatmap.2(input_heatmap,col=palette,Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Single linkage")

# Average linkage
heatmap.2(heatmap,col=palette,Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Average linkage")

