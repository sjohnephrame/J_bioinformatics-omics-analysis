#==========================================================================================#
#CTRL VS OGT KO MICE BRAINS #WT-OGT_KO-MouseBrain-Mecano-10OutOf50-Fr1-24-Combined
#Dr Dias_OGT KO MICE BRAIN - Total proteomics analysis
#==========================================================================================#

# Clear environment
#rm(list = ls())

# Attach libraries
library(readxl)
library(tidyverse)
library(limma)
library(EnvStats)   # used to get geometric means
library(missForest) # Imputation
library(biomaRt) # library for mapping between annotations
library(qvalue)

# Load helper functions
path.to.functions = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\functions2.R"
source(path.to.functions)
#source(path/to/my_functions.R): In R, the source() function executes the script 
#at the specified path, making its functions and variables available in the 
#current environment.
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

#=========================#
# Proteomics Preprocessing ----
#=========================#

dat1k = read_xlsx("C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_DR DIAS OGT KO MICE\\WT-OGT_KO-MouseBrain-Mecano-10OutOf50-Fr1-24-Combined (1) - HIGH AND MED FDR.xlsx")
colnames(dat1k)
data2k<-dat1k

colnames(dat1k)
head(dat1k)
#dat1<-dat1[,-(20:23)]
colnames(dat1k)
names(dat1k)[names(dat1k) == 'Gene Symbol'] <- 'Symbol'
colnames(dat1k)
head(dat1k)
data2newk<-dat1k
library(dplyr)
dat1k <- dat1k %>%
  relocate(Symbol, .after = Accession)
colnames(dat1k)
head(dat1k)
data2new1k<-dat1k


dim(dat1k)
#[1] 8476   11 before removing duplicates 

new_dat1<-dat1k
head(new_dat1)
new_dat1<-new_dat1[,-(1)] #Remove FDR
head(new_dat1)
new_dat1<-new_dat1[,-(3)] #Remove Coverage
head(new_dat1)
new_dat1<-new_dat1[,-(9)] #Remove modfn
head(new_dat1) 

dim(new_dat1) #[1] 8476    8#without accession


#add_symbols - adds symbols to refseq id and removes duplicates
head(new_dat1)
new_dat2<-new_dat1
new_dat3<-new_dat1
head(new_dat3)
new_dat3<-new_dat3[,-(2)] #Remove symbol because not all symbols present; will re-add symbols
head(new_dat3)
new_dat2 = add_symbols(new_dat2) # This returns unique gene symbols by numbering any duplicates. Be careful!
new_dat3 = add_symbols(new_dat3) # This returns unique gene symbols by numbering any duplicates. Be careful!
# Duplicates are made unique by appending '.2', '.3', etc. 
# sum(grepl("\\.", dat$Symbol))
head(new_dat1)
head(new_dat2)
head(new_dat3)
new_dat4<-new_dat3

#but we need like matrix with no symbols

head(new_dat4)
new_dat4<-as.data.frame(new_dat4)
head(new_dat4)
rownames(new_dat4) <- new_dat4$Accession #convert column value to row names
head(new_dat4)
new_dat4.1<-new_dat4
head(new_dat4.1)
new_dat4<-new_dat4[,-(1)] #Remove Accession
head(new_dat4)

# This might look weird, but there will be NaN for some values that had no data
# and this will change them to NA, which is easier to handle.
new_dat4[is.na(new_dat4)] <- NA
dim(new_dat4) #[1] 7794    7 #without accession
new_dat4.1[is.na(new_dat4.1)] <- NA
dim(new_dat4.1) #[1] 7794   8 #with accession

# Assess missingness of proteins

head(new_dat4)
dim(new_dat4) #[1] 7794 7 #without accession, with symbol
new_dat5<-new_dat4[,-1] #without accession, without symbol
head(new_dat4) #with symbols
head(new_dat5) #without symbols
dim(new_dat4) #with symbols [1] 7794   7
dim(new_dat5) #without symbols [1] 7794   6

head(new_dat4) #with symbols
table(apply(new_dat4, 1, function(x) sum(is.na(x)))) #with symbols
# 0    6 
#7104  690 

head(new_dat5) #without symbols
table(apply(new_dat5, 1, function(x) sum(is.na(x)))) #with symbols
#0    6 
#7104  690 

#-------------------
#Code explanation:
#-------------------

#apply(new_dat1, 1, function(x) sum(is.na(x))):

#new_dat1: Your data frame.
#1: Indicates that the function should be applied to rows. 2 for cols
#function(x) sum(is.na(x)): A function that calculates the number of NA values in each row.

#----------------------------------------------------------------------------------------------


#thus, there is no missing value with symbol
#so ideally good to map symbols here before proceeding further

#For now let us omit NA with symbols

# Reduce to proteins with no missing data.

dim(new_dat4) #7868 rows 7 cols in new_dat4 #with symbols
head(new_dat4)
dat_comp1 = na.omit(new_dat4) #with symbols and no missing values
dim(dat_comp1) #with symbols and no missing values
#[1] 7104    7 #19 becauseof Symbol
head(new_dat4)
dim(new_dat4) #7794 rows 7 cols in new_dat1 #with symbols, without omit
dim(dat_comp1) #[1] 7104   7 #7 because of symbol, with omit
mapping1<-dat_comp1 
head(mapping1)
dim(mapping1)  #[1] 7171 7 #7 because of symbol, with omit
#mapping1<-mapping1[,-(2:19)]

#without symbols
head(new_dat5)
dim(new_dat5)
dat_comp2 = na.omit(new_dat5)
dim(dat_comp2)
#[1] 7104 6 #without symbols

#It is same with/without symbols

# log2 transform - to reduce skewness of a measurement variable 
#to make data more symmetrical, which helps it meet the assumptions of statistical models

head(dat_comp2)
dim(dat_comp2) #[1] 7104    6 after removing symbol, without NA


head(new_dat5)
dim(new_dat5) #[1] 7794    6 after remving symbols # with NAs

head(dat_comp2)
dim(dat_comp2) #[1] 7104    6 after removing symbol # without NAs
dat_comp_mat1 = as.matrix(log2(dat_comp2)) # without NAs #log2 transformation
head(dat_comp_mat1)
dim(dat_comp_mat1) #[1] 7104   6 after removing symbol # without NAs #log2 transformation

head(new_dat5)
dim(new_dat5)# with NAs
dat_mat1 = as.matrix(log2(new_dat5)) # with NAs
head(dat_mat1)
dim(dat_mat1) #[1] 7868    6 with NAs_log
dim(dat_comp_mat1) #[1] 7104   6  without NAs_log
head(dat_comp_mat1)
head(dat_mat1)
dim(dat_mat1) #[1] 7794    6 with NAs_log
dim(dat_comp_mat1) #[1] 7104   6  without NAs_log

#density plot - plot is applied to log10 transformed data in code 

head(dat_comp2)#not log transformed, not normalized, omitted NA
dim(dat_comp2)#[1] 7104 6 without NAs_ not log, not normalized
for (i in 1:ncol(dat_comp2)){
  if (i==1){
    plot(density(log10(dat_comp2[,1])),main="Density plot - OGT KO Mice brain \n Dr. Dias Collaboration Project",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp2[,i]))
    lines(den$x,den$y,col="blue")}
}

#uniform normal distribution - suitable for limma

#------------------------------------------------------------------------------------------------#
# Quantile normalize - to remove any technical variation before differential expression analysis
#------------------------------------------------------------------------------------------------#

head(dat_comp_mat1)
dim(dat_comp_mat1) #[1] 7104   6  without NAs_log
dat_comp_mat_quant1 = limma::normalizeQuantiles(dat_comp_mat1) # without NAs, log
head(dat_comp_mat_quant1)
dim(dat_comp_mat_quant1) #7104    6

head(dat_mat1)
dim(dat_mat1) #7794    6
dat_mat_quant1 = limma::normalizeQuantiles(dat_mat1) # with NAs, log
head(dat_mat_quant1)
dim(dat_mat_quant1)# 7794    6 with NAs, log, normalize #9206   18

head(dat_comp_mat_quant1) # without NAs, log, normalize 
dim(dat_comp_mat_quant1) #7104    6 # without NAs, log, normalize #8334   18

head(dat_comp2)
dim(dat_comp2) # 7104    6 without NAs, not log,
dat_comp1_norm = limma::normalizeQuantiles(dat_comp2) # without NAs, not log, normalized
dim(dat_comp1_norm) #7104   6 # without NAs, not log, normalized
head(dat_comp1_norm)


#-------------------------------------------------------------------------------------------------
# Median normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------
head(dat_comp_mat1)
head(dat_mat1)
dim(dat_comp_mat1) #[1] 7104    6 without NAs_log_matrix; not normalized
dim(dat_mat1) #[1] 7794    6 with NAs_log_matrix
dat_comp_mat_median1 = limma::normalizeMedianValues(dat_comp_mat1) # without NAs, log
dat_mat_median1 = limma::normalizeMedianValues(dat_mat1) # with NAs
head(dat_mat_median1) #[1] 7794    6 # with NAs, log, median norm
dim(dat_mat_median1) #[1] 7794    6 # with NAs, log, median norm
head(dat_comp_mat_median1) #[1] 7104    6 # without NAs, log, median norm
dim(dat_comp_mat_median1) #[1] 7104    6  # without NAs, log, median norm


#https://www.youtube.com/watch?v=ecjN6Xpv6SE&t=279s


## Batch Effect Correction ##
# Here is where one should run diagnostics for batch effects (i.e., systematic differences in measurements due to technical factors rather than biological signal), if the experimental design introduces batch effects

# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation
head(dat_comp_mat_quant1)
dim(dat_comp_mat_quant1) # without NAs, log, normalize #7171 6
#do a box plot
# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation

col <- c("blue","blue","blue","orange","orange","orange")
col
fill <- c("blue","orange")
fill

#par(mar = c(3, 3, 3, 10), mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2, xpd = TRUE)
boxplot(dat_comp_mat1, main="Boxplot before normalization - OGT KO Mice brain \n Dr. Dias Collaboration Project",col=col)
boxplot(dat_comp_mat_quant1, main="Boxplot Quantile normalized Intensities - OGT KO Mice brain \n Dr. Dias Collaboration Project",col=col)
boxplot(dat_comp_mat_median1, main="Boxplot Median normalized Intensities - OGT KO Mice brain \n Dr. Dias Collaboration Project",col=col)

#inset depends on margin of plot window

#=================================#
# Differential Expression Analysis ----
#     using limma-trend - # This assumes uniform and normal distribution of data (seen in density plot)
#=================================#

dat1 = dat_comp_mat_quant1
head(dat1)
dim(dat1) #[1] 7171   6 #removed NA, without Symbol, norm, log

# Change column names to not contain "-", since this will conflict with naming syntax in contrast matrix

colnames(dat1)
head(dat1)

# Getting a data.frame with 'Accession' and 'Symbol' columns

head(mapping1)
head(dat1)
dim(dat1) #[1] 7104 6 Without symbol
dim(mapping1) #[1] 7104 7 With Symbol
head(mapping1)
#both dat1 and mapping1 match
#so, try to separate symbol and accession for mapping1

mapping1.1<-mapping1
head(mapping1.1)
dim(mapping1.1) #[1] 7104   7
head(mapping1) #[1] 7104   7
dim(mapping1)
#make dim(mapping1) 8334 2
mapping2<-as.data.frame(cbind(rownames(mapping1),mapping1$Symbol))
dim(mapping2) #7171 2
head(mapping2)
colnames(mapping2) <- c('Accession','Symbol')
colnames(mapping2)
head(mapping2)
dim(mapping2) #[1] 7104    2

#----------------------------------------------------------------------------------------
#Code to add entrez gene id - need for AMEND
#----------------------------------------------------------------------------------------

#library(org.Hs.eg.db)
#dat4<-dat1
#Accession<-rownames(dat4)
#dim(dat4)

#idfound <-  Accession %in% mappedRkeys(org.Hs.egUNIPROT)
#y <- dat4[idfound,]
#table(idfound)
#y

#We add Entrez Gene IDs to the annotation
#egUNIPROT <- toTable(org.Hs.egUNIPROT)
#head(egREFSEQ)

#m <- match(y$genes$RefSeqID, egREFSEQ$accession)
#y$genes$EntrezGene <- egREFSEQ$gene_id[m]

## Now use Entrez Gene IDs to find gene symbol
#egSYMBOL <- toTable(org.Hs.egSYMBOL)
#head(egSYMBOL)
#m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
#y$genes$Symbol <- egSYMBOL$symbol[m]
#head(y$genes)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# Design matrix
head(dat1)
colnames(dat1) <- c('CTRL_1','CTRL_2','CTRL_3','OGTKO_1','OGTKO_2','OGTKO_3')
head(dat1)
txa = extract_string(colnames(dat1), "_", 1)
txa
Xa = model.matrix(~0 + txa)
Xa
colnames(Xa) = sort(unique(txa))
Xa #design matrix

# Contrast matrix
colnames(Xa)
#contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
contrastsa = c("OGTKO-CTRL")
#contrastsb = c("CTRL-OGTKO") Note: This is wrong. Must do trtmnt/ctrl so trtmnt 1st
#contrastsb
contrastsa
Xa
cm_data = makeContrasts(contrasts=contrastsa, levels=colnames(Xa))
cm_data #contrast matrix  #Note: Check with Dr. Thompson if it is correct!!!!!!
#cm_datab = makeContrasts(contrasts=contrastsb, levels=colnames(Xa)) #Did to check
#cm_datab

#Check out this recent paper by the authors which is precisely focused on explaining design matrices.

#Basically, limma (and other software) are focused on using linear models to analyse your data, 
#because these tools let you represent/model your study, taking into account any variables which are
#of interest to you.

#For example, you may want to compare expression of gene X between young and old samples. 
#A t-test will simply compare the means of the groups. However, sex may influence 
#the expression of gene X. A linear model lets you take into account sex when comparing
#the expression between young and old groups.

#Having said this, design matrices are the way that you indicate the study 
#variables to model (e.g., age and sex), and contrast matrices are the way that you 
#indicate which variables you want to test (e.g. age) between which groups.
#https://www.biostars.org/p/9554210/

#Linear models describe a continuous response variable as a function of one or more predictor variables. 
#eg linear regression


# Fitting models. The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)

head(mapping2)
head(dat1)
Xa
cm_data
dea = top_table_results(mapping2, dat1, Xa, cm_data, annotate = TRUE, trend = TRUE)
head(dea)
head(dea$tt) #8 columns #log FC 
gene.symbola = extract_string(dea$tt$Symbol, "\\.", 1)
head(gene.symbola)
head(dea$tt$Symbol)
# both match
tmpa = aggregate(dea$tt[,-c(7,8)], by = list(Symbol = gene.symbola), FUN = median) #remove last 2 columns
head(dea$tt)
head(tmpa)
tmpa$Accession = dea$tt$Accession[match(tmpa$Symbol, gene.symbola)] # Get the first listed Accession ID for each gene symbol
head(tmpa)
head(dea$tt)
dea$tt = tmpa
head(dea$tt)
head(tmpa)

##################################################################################################################
######################################################################################################################
#logFCprot - if you need an excel sheet of the differentially expressed values - For heatmap
head(dea$tt)
logFCprotKO<-dea$tt
head(logFCprotKO)
dim(logFCprotKO) #[1] 6845    8
logFCprotKO<-as.data.frame(logFCprotKO)
head(logFCprotKO)
type(logFCprotKO)
library("writexl") 
#install.packages("writexl")
#write_xlsx(logFCprotKO,"C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\logFCprotKOrep.xlsx")
dim(logFCprotKO)

# Fitting models - Simplified code - The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)
#de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)
#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
#tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median)
#tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
#de$tt = tmp

#head(de$tt)

##################################################################################################################
######################################################################################################################

#MDS - Multi-dimensional scaling plot - similar to PCA
head(dea$tt) #values after comparisons with contrast matrix
head(dat1) #normailzed and filtered values - use this for mds plot
dim(dat1)
labels<-colnames(dat1)
labels
library(ggplot2)
col <- c("blue","blue","blue","orange","orange","orange")
col
#https://support.bioconductor.org/p/100859/
plotMDS(dat1, main="Multipledimensional scaling plot (MDS) - OGT KO Mice brain \n Dr. Dias Collaboration Project",col=col)


#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html

######################################################################################

# PCA - Principal Component Analysis

labels<-c("CTRL","CTRL","CTRL","OGTKO","OGTKO","OGTKO")
labels
get_pca
get_pca1 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA) - OGT KO Mice Brain - Afetr Quantile Normalization \n Dr. Dias Collaboration Project" )+
    scale_color_manual(values = c("blue","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

#dat1 = dat_comp_mat_quant1
head(dat1)  
PCA<-get_pca1(dat1,labels)
PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/

get_pca2 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA) - OGT KO Mice brain - Afetr Median Normalization \n Dr. Dias Collaboration Project" )+
    scale_color_manual(values = c("blue","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

PCA<-get_pca2(dat_comp_mat_median1,labels)
PCA$p

##########################################################################################

#do a box plot
# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation

col <- c("blue","blue","blue","orange","orange","orange")
col
fill <- c("blue","orange")
fill
#par(mar = c(3, 3, 3, 10), mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2, xpd = TRUE)
boxplot(dat_comp_mat_quant1, main="Boxplot quantile normalized Intensities",col=col)
#legend("right", title="Treatment",legend=c("GFP-0","GFP-10","OGTKD-605-0","OGTKD-605-10","OGTKD-606-0","OGTKD-606-10"), fill=fill, inset = c(-0.155,0))

#inset depends on margin of plot window

##################################################################################################################################

# Volcano plots for p-value and adjusted p-values

colnames(dea$tt) #logFC and pvalues
head(dea$tt)

#CTRL-OGTKO

rx <- c(-1, 1)*max(abs(dea$tt$logFC))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTKO.CTRL), -log10(dea$tt$adj.P.Val.OGTKO.CTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")


#adjusted p value
plot(dea$tt$logFC, -log10(dea$tt$adj.P.Val.OGTKO.CTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KO Mice Brain_Dr. Dias Collaboration Project \n Without logFC cutoff")

# only adj p-value cutoff
lfc <- 0
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$logFC)> lfc & -log(dea$tt$adj.P.Val.OGTKO.CTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$logFC)< (-lfc) & -log(dea$tt$adj.P.Val.OGTKO.CTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$logFC,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTKO.CTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$logFC,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTKO.CTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
#abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)


###################################################################################

colnames(dea$tt) #logFC and pvalues
head(dea$tt)

rx <- c(-1, 1)*max(abs(dea$tt$logFC))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTKO.CTRL), -log10(dea$tt$adj.P.Val.OGTKO.CTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")


#p-value
plot(dea$tt$logFC, -log10(dea$tt$P.Value.OGTKO.CTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10 p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of p-values - OGT KO Mice Brain_Dr. Dias Collaboration Project")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$logFC)> lfc & -log(dea$tt$P.Value.OGTKO.CTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$logFC)< (-lfc) & -log(dea$tt$P.Value.OGTKO.CTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$logFC,-log(dea$tt[sigGenes1,]$P.Value.OGTKO.CTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$logFC,-log(dea$tt[sigGenes2,]$P.Value.OGTKO.CTRL,10),pch=20,col="blue",cex=2)

#points(dea$tt[sigGenes1,]$logFC,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTKO.CTRL,10),pch=20,col="orange",cex=2)
#points(dea$tt[sigGenes2,]$logFC,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTKO.CTRL,10),pch=20,col="blue",cex=2)

abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)

colnames(dea$tt) #logFC and pvalues
head(dea$tt)

#adjusted p value
plot(dea$tt$logFC, -log10(dea$tt$adj.P.Val.OGTKO.CTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KO Mice Brain_Dr. Dias Collaboration Project")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$logFC)> lfc & -log(dea$tt$adj.P.Val.OGTKO.CTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$logFC)< (-lfc) & -log(dea$tt$adj.P.Val.OGTKO.CTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$logFC,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTKO.CTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$logFC,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTKO.CTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

# Combine significant up and down genes
sig_indices <- which(sigGenes1 | sigGenes2)
sig_data <- dea$tt[sig_indices, ]

# Sort by adjusted p-value and take the top 20
top_20 <- sig_data[order(sig_data$adj.P.Val.OGTKO.CTRL), ][1:5, ]

# Coordinates for the points
x_points <- top_20$logFC
y_points <- -log10(top_20$adj.P.Val.OGTKO.CTRL)

# Offset coordinates for the labels (adjust these values if labels overlap too much)
x_labels <- x_points + 0.15 
y_labels <- y_points + 0.2

# 1. Draw the arrows
# length = 0.05 controls the size of the arrowhead
arrows(x0 = x_labels, y0 = y_labels, 
       x1 = x_points, y1 = y_points, 
       length = 0.05, col = "grey30", lwd = 1)

# 2. Draw the text labels
text(x_labels, y_labels, 
     labels = top_20$Symbol, 
     pos = 4,        # 4 = right of the coordinate
     cex = 0.7,      # font size
     font = 2,       # bold
     offset = 0.2)   # space between arrow end and text




# 1. Identify all significant genes first
sig_indices <- which(sigGenes1 | sigGenes2)
sig_data <- dea$tt[sig_indices, ]

# 2. Sort by p-value and take the top 20
top_20 <- sig_data[order(sig_data$adj.P.Val.OGTKO.CTRL), ][1:20, ]

# 3. Add labels only for these 20
text(top_20$logFC, 
     -log10(top_20$adj.P.Val.OGTKO.CTRL), 
     labels = top_20$Symbol, 
     cex = 0.7,      # Smaller text
     pos = 4,        # Position to the right of the dot
     font = 2)       # Bold font for better visibility

# Define which genes to label (union of up and down regulated)
all_sig <- sigGenes1 | sigGenes2

# Add the labels using the 'Symbol' column
text(dea$tt[all_sig, ]$logFC, 
     -log10(dea$tt[all_sig, ]$adj.P.Val.OGTKO.CTRL), 
     labels = dea$tt[all_sig, ]$Symbol, 
     cex = 0.6,      # Adjust text size
     pos = 3,        # Position text above the point (1=below, 2=left, 3=above, 4=right)
     offset = 0.5)   # Distance from the point

install.packages("calibrate")
library(calibrate)
textxy(dea$tt[allSig,]$logFC, -log10(dea$tt[allSig,]$adj.P.Val.OGTKO.CTRL), 
       labs=rownames(dea$tt)[allSig], cex=0.6)

#############################################################################################################################################

#Prep for GSEA

dea
dea2<-dea
dea$tt
head(dea2)
dea2$tt
head(dea2$tt)

#CTRL VS OGT KO
updownOGTTENGFPTEN<-cbind(dea2$tt$logFC,dea2$tt$adj.P.Val.OGTKO.CTRL,dea2$tt$Symbol)
head(updownOGTTENGFPTEN)
head(dea2$tt)
updownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN)
colnames(updownOGTTENGFPTEN)<-c('logFC','adj.P.Val.OGTKO.CTRL','Symbol')
head(updownOGTTENGFPTEN)
rownames(updownOGTTENGFPTEN)<-updownOGTTENGFPTEN$Symbol
head(updownOGTTENGFPTEN)
updownOGTTENGFPTEN<-updownOGTTENGFPTEN[,-3]
head(updownOGTTENGFPTEN)
updownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN)
head(updownOGTTENGFPTEN)
#convert e to nmeric so can sort
updownOGTTENGFPTEN$adj.P.Val.OGTKO.CTRL <- as.numeric(formatC(updownOGTTENGFPTEN$adj.P.Val.OGTKO.CTRL, format = "e", digits = 2))
head(updownOGTTENGFPTEN)
sortupdownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN[order(updownOGTTENGFPTEN$adj.P.Val.OGTKO.CTRL,decreasing = TRUE),]) 
head(sortupdownOGTTENGFPTEN)
tail(sortupdownOGTTENGFPTEN)
dim(sortupdownOGTTENGFPTEN)
###################################################################################################################
###################################################################################################################

#NOT NEEDED FOR GSEA BUT CAN DO FOR KEGG AND GO
#remove genes not significant p>0.05 
remsortupdownOGTTENGFPTEN <- sortupdownOGTTENGFPTEN[!(sortupdownOGTTENGFPTEN$adj.P.Val.OGTKO.CTRL>0.05), ]
head(remsortupdownOGTTENGFPTEN)
remsortupdownOGTTENGFPTEN<-as.data.frame(remsortupdownOGTTENGFPTEN)
head(remsortupdownOGTTENGFPTEN)
dim(remsortupdownOGTTENGFPTEN) #1] 452   2
#463 proteins < p=0.05 and FDR Filtration 

# Filter log FC >1.5 and <-1.5 separately
head(remsortupdownOGTTENGFPTEN)

logFCless<-remsortupdownOGTTENGFPTEN
head(logFCless)
logFCless<-as.data.frame(logFCless)
head(logFCless)
type(logFCless)
#logFCless<-as.numeric(logFCless$logFC.OGTTEN.GFPTEN)
#head(logFCless)
logFCless <- transform(logFCless, logFC = as.numeric(logFC))
head(logFCless)
type(logFCless)
#logFCless<-subset(logFCless,(((logFCless$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCless$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCless<-subset(logFCless,(((logFCless$logFC)<(-0.58))))
dim(logFCless)
#0,2

logFCgreat<-remsortupdownOGTTENGFPTEN
head(logFCgreat)
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-as.numeric(logFCgreat$logFC.OGTTEN.GFPTEN)
#head(logFCgreat)
logFCgreat <- transform(logFCgreat, logFC = as.numeric(logFC))
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC)>(0.58))))


head(logFCgreat)
dim(logFCgreat) #61 2

#>0.5 and <-0.5 combined proteins

logFCcombo<-remsortupdownOGTTENGFPTEN
head(logFCcombo)
logFCcombo<-as.data.frame(logFCcombo)
head(logFCcombo)
type(logFCcombo)
#logFCcombo<-as.numeric(logFCcombo$logFC.OGTTEN.GFPTEN)
#head(logFCcombo)
logFCcombo <- transform(logFCcombo, logFC = as.numeric(logFC))
head(logFCcombo)
type(logFCcombo)
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCcombo<-as.data.frame(logFCcombo)
head(logFCcombo)
type(logFCcombo)
logFCcombo<-subset(logFCcombo,logFCcombo$logFC>0.58 | logFCcombo$logFC<(-0.58))
head(logFCcombo)
#logFCcombo <- logFCcombo[logFCcombo$logFC.OGTTEN.GFPTEN > 1.5 | logFCcombo$logFC.OGTTEN.GFPTEN < -1.5, ]
head(logFCcombo)
dim(logFCcombo) #61  2 combined >0.5 and <-0.5

###################################################################################################################
###################################################################################################################

#Prepare protlist for waterfall plot - rank-ordered protein list

#GSEA 606 PROT LIST
head(sortupdownOGTTENGFPTEN) #all proteins without removing based on significance
dim(sortupdownOGTTENGFPTEN)
prot.list <- sortupdownOGTTENGFPTEN$logFC             # rank-ordered protein list
head(prot.list)
prot.list<-as.numeric(prot.list)
head(prot.list)
head(sortupdownOGTTENGFPTEN)
names(prot.list) <- rownames(sortupdownOGTTENGFPTEN)
head(prot.list)

# Waterfall plot - rank-ordered protein list

# Barplot of ranked fold changes (waterfall plot)
#barplot(sort(gene.list, decreasing = T),axisnames=FALSE,main="Plot of ranked gene Fold changes")

head(prot.list)
col <- ifelse(prot.list >0, "orange", "blue")
barplot(sort(prot.list, decreasing = T),axisnames=FALSE,main="Plot of ranked Fold changes of total proteome - OGT KO Mice Brain Samples \n DR. Dias Collaboration Project",  col=col, border=col, ylim=c(-4,4))

###############################################################################################################################################

###GO analysis

#goana: Gene Ontology or KEGG Pathway Analysis
#In limma: Linear Models for Microarray Data

#topGO package provides tools for testing GO terms while accounting for the topology of the GO graph. Different test statistics and different methods for eliminating local similarities and dependencies between GO terms can be implemented and applied.

#In R, topKEGG is a function that extracts the most significant KEGG pathways from kegga output. 
#kegga function to gather pathway enrichment for my dataset.

library(org.Mm.eg.db)
#library(edgeR)
library(GO.db)

head(remsortupdownOGTTENGFPTEN)
dim(remsortupdownOGTTENGFPTEN) #[1] 452   2 Filtered for adj. p. value
d.go1 <- remsortupdownOGTTENGFPTEN
head(d.go1) #subsetted already for p<0.05
#d.go.DE <- subset(d.go,PValue<0.05) 
d.go.DE1<-d.go1
head(d.go.DE1)
d.entrez.id1 <- mapIds(org.Mm.eg.db, keys=rownames(d.go.DE1),column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id1)
head(d.entrez.id1)
all(rownames(d.go.DE1)==names(d.entrez.id1)) 
go.test1 <- goana(d.entrez.id1,species="Mm")
go.results1 <- topGO(go.test1, sort = "DE", number = Inf)
head(go.results1)
sum(go.results1$P.DE<10^(-5)) #461
sum(go.results1$P.DE<0.05) #2896

head(remsortupdownOGTTENGFPTEN)
dim(remsortupdownOGTTENGFPTEN) #[1] 452   2 Filtered for adj. p. value
head(logFCgreat)
dim(logFCgreat) #61  2
d.go2 <- logFCgreat
head(d.go2) #subsetted already for p<0.05, FDR, log FC >1.5 and <-1.5
#d.go.DE <- subset(d.go,PValue<0.05) 
d.go.DE2<-d.go2
head(d.go.DE2)
dim(d.go.DE2)
d.entrez.id2 <- mapIds(org.Mm.eg.db, keys=rownames(d.go.DE2),column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id2)
head(d.entrez.id2)
all(rownames(d.go.DE2)==names(d.entrez.id2)) 
go.test2 <- goana(d.entrez.id2,species="Mm")
go.results2 <- topGO(go.test2, sort = "DE", number = Inf)
head(go.results2)
sum(go.results2$P.DE<10^(-5)) #248
sum(go.results2$P.DE<0.05) #1574


#topBP <- topGO(go.test1, ontology = "BP", number = 20)
#print(topBP)

library(ggplot2)

# Get top 20 pathways
df_plot <- topGO(go.test1, sort = "DE", number = 20)
head(go.test1)
# Add -log10 P-value for the x-axis
df_plot$logP <- -log10(df_plot$P.DE)
head(df_plot)
# Reorder Term by significance so the plot is ranked
df_plot$Term <- reorder(df_plot$Term, df_plot$logP)
head(df_plot)

#Dot plot - Top 20 Enriched GO Terms - OGT KD 605 10 min vs GFP - 10 min
#Without filtering
ggplot(df_plot, aes(x = logP, y = Term, size = DE, color = logP)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 20 Enriched GO Terms - Before Filtering \n CTRL vs. OGT KO Mice Brain - Dr. Dias Collaboration Project",
       x = "-log10(P-value)",
       y = "GO Term",
       size = "Gene Count",
       color = "Significance") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # This centers the title
    axis.title = element_text(face = "bold"), # This bolds the Y axis label
    axis.text.y = element_text(face = "bold")   # Bolds the individual terms
  ) 

#Terms like "intracellular anatomical structure" and "protein binding" are so broad that they include almost half the genome. While the 
#-values are astronomically low (10^-202), they aren't telling you much about the specific biology of your experiment.


#Focus on Biological Process (BP)
#Cellular Component (CC) and Molecular Function (MF) often give 
#generic structural results. Biological Process usually contains
#the pathways researchers actually care about (e.g., "glucose metabolic process").

# Filter for just Biological Process
head(go.test1)
go.bp <- topGO(go.test1, ontology = "BP", sort = "DE", number = 50)
head(go.bp, 20)

#N: The total number of genes in the entire background "universe" 
#that are annotated to that specific GO term.
#DE: The number of genes from your Differentially Expressed (DE) list 
#that are annotated to that GO term.

#Let's grab terms that have fewer than 1,000 total genes (N < 1000). This removes 
#the "cellular process" noise and finds the specific pathways.

# Filter for Biological Process and smaller, more specific terms
go.specific <- go.results1[go.results1$Ont == "BP" & go.results1$N < 1000, ]

# Take the top 20 of THESE results
plot_data <- head(go.specific, 20)
head(plot_data)
# Create -log10 P-value for the plot
plot_data$logP <- -log10(plot_data$P.DE)
plot_data$Term <- reorder(plot_data$Term, plot_data$logP)
head(plot_data)
print(plot_data)

library(ggplot2)

#if (!require(ggtext)) install.packages("ggtext")
library(ggtext)

#Filtered for BP and N < 1000 genes
ggplot(plot_data, aes(x = logP, y = Term)) +
  geom_point(aes(size = DE, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Gene Ontology (GO) - Top 20 Specific Biological Processes",
    # Use HTML tags: <b> for bold, <br> for new line
    subtitle = "<b>CTRL vs. OGT KO Mice Brain Samples</b><br>Filtered for N < 1000 genes",
    x = "-log10(P-value)",
    y = NULL,
    size = "Genes in List",
    color = "Significance"
  ) +
  theme(
    # This line tells ggplot to interpret the HTML tags in the subtitle
    plot.subtitle = element_markdown(hjust = 0.5, lineheight = 1.2), 
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(face = "bold")
  )

#Filetr N<500
# Filter for Biological Process and smaller, more specific terms
go.specific2 <- go.results1[go.results1$Ont == "BP" & go.results1$N < 500, ]

# Take the top 20 of THESE results
plot_data2 <- head(go.specific2, 20)

# Create -log10 P-value for the plot
plot_data2$logP <- -log10(plot_data2$P.DE)
plot_data2$Term <- reorder(plot_data2$Term, plot_data2$logP)

print(plot_data2)

library(ggtext)

#Update the Plot Code for mixed title

ggplot(plot_data2, aes(x = logP, y = Term)) +
  geom_point(aes(size = DE, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Gene Ontology (GO) - Top 20 Specific Biological Processes",
    # Use HTML tags: <b> for bold, <br> for new line
    subtitle = "<b>CTRL vs. OGT KO Mice Brain Samples</b><br>Filtered for N < 500 genes",
    x = "-log10(P-value)",
    y = NULL,
    size = "Genes in List",
    color = "Significance"
  ) +
  theme(
    # This line tells ggplot to interpret the HTML tags in the subtitle
    plot.subtitle = element_markdown(hjust = 0.5, lineheight = 1.2), 
    
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(face = "bold")
  )

#if (!require(gt)) install.packages("gt")

library(gt)
library(dplyr)
#install.packages("webshot2")
library(webshot2)
# Prepare and style the table

plot_data %>%
  # Target only the relevant columns
  dplyr::select(Term, Ont, N, DE, P.DE) %>% 
  gt() %>%
  # 1. Bold the Title and Subtitle
  tab_header(
    title = md("**Gene Ontology (GO) Top 20 Specific Biological Processes**"),
    # Use <br> for the line break and move the ** markers
    subtitle = md("**OGT KD 605 10 min vs GFP 10 min** <br> Filtered for N < 1000 genes")
  ) %>%
  # 2. Bold ALL Column Headings
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  # 3. Bold ONLY the 'Term' column body
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Term)
  ) %>%
  # 4. Final Formatting (Scientific notation and alignment)
  fmt_scientific(columns = P.DE, decimals = 2) %>%
  opt_row_striping() %>%
  cols_align(align = "center", columns = c(Ont, N, DE, P.DE))%>%
  cols_width(P.DE ~ px(150)) %>% # Adjust 150 higher if you want even more space
  gtsave("C:\\Users\\sophi\\OneDrive\\Desktop\\J_TOTAL PROTEOME OGT KO mice brain PLOTS\\J_Updated figures_04052026\\GO_TABLE1_CTRL VS. OGT KO MICE BRAIN_N LESS 1000.png") # This saves it to your working directory

plot_data2 %>%
  # Target only the relevant columns
  dplyr::select(Term, Ont, N, DE, P.DE) %>% 
  gt() %>%
  # 1. Bold the Title and Subtitle
  tab_header(
    title = md("**Gene Ontology (GO) Top 20 Specific Biological Processes**"),
    # Use <br> for the line break and move the ** markers
    subtitle = md("**OGT KD 605 10 min vs GFP 10 min** <br> Filtered for N < 500 genes")
  ) %>%
  # 2. Bold ALL Column Headings
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  # 3. Bold ONLY the 'Term' column body
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Term)
  ) %>%
  # 4. Final Formatting (Scientific notation and alignment)
  fmt_scientific(columns = P.DE, decimals = 2) %>%
  opt_row_striping() %>%
  cols_align(align = "center", columns = c(Ont, N, DE, P.DE))%>%
  cols_width(P.DE ~ px(150)) %>% # Adjust 150 higher if you want even more space
  gtsave("C:\\Users\\sophi\\OneDrive\\Desktop\\J_TOTAL PROTEOME OGT KO mice brain PLOTS\\J_Updated figures_04052026\\GO_TABLE1_CTRL VS. OGT KO MICE BRAIN_N LESS 500.png") # This saves it to your working directory



###KEGG analysis

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("KEGGREST")
library(KEGGREST)
#kegg.test <- kegga(d.entrez.id,species="Mm")
#kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
#head(kegg.results)
#sum(kegg.results$P.DE<10^(-5))

d.kegg1 <- remsortupdownOGTTENGFPTEN
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
#d.kegg.DE <- subset(d.kegg,PValue<0.05)
d.kegg.DE1<-d.kegg1
head(d.kegg.DE1)
dim(d.kegg.DE1)
head(d.entrez.id1)
all(rownames(d.kegg.DE1)==names(d.entrez.id1)) 

kegg.test1 <- kegga(d.entrez.id1,species="Mm")
kegg.results1 <- topKEGG(kegg.test1, sort = "DE", number = Inf)
head(kegg.results1)
sum(kegg.results1$P.DE<10^(-5))#10
sum(kegg.results1$P.DE<0.05) #74

d.kegg2 <- logFCgreat
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
#d.kegg.DE <- subset(d.kegg,PValue<0.05)
d.kegg.DE2<-d.kegg2
head(d.kegg.DE2)
head(d.entrez.id2)
all(rownames(d.kegg.DE2)==names(d.entrez.id2)) 

kegg.test2 <- kegga(d.entrez.id2,species="Mm")
kegg.results2 <- topKEGG(kegg.test2, sort = "DE", number = Inf)
head(kegg.results2)
sum(kegg.results2$P.DE<10^(-5)) #4
sum(kegg.results2$P.DE<0.05) #39

######################################################################################################################
######################################################################################################################

### J_GSEA analysis old code with prefixes 

# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
library(fgsea)

#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v2023.2.Hs.symbols.gmt") #got lot of warnings
#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
#all.gene.sets <- gmtPathways("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\m2.all.v2023.1.Mm.symbols.gmt")        #m2curated gene sets
#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\m2.all.v2023.1.Mm.symbols.gmt")        #m2curated gene sets
all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\m2.all.v2026.1.Mm.symbols.gmt")        #m2curated gene sets

class(all.gene.sets)
length(all.gene.sets)
all.gene.sets[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
library(tidyverse)
all.gene.sets %>% head() %>% lapply(head)

head(prot.list)
#head(d.entrez.id1)
#d.entrez.id2<-
### Now run fgsea 
fgseaRes <- fgsea(pathways = all.gene.sets, stats = prot.list, minSize=15, maxSize=500,eps=0)
#warnings()
head(fgseaRes)
head(fgseaRes[order(pval), ])
sum(fgseaRes[, padj < 0.05])#279
fgseaRes1<-fgseaRes[order(padj), ]
head(fgseaRes1)


# Make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#dev.off() #run if plots get messed up
library(ggplot2)
p<-plotGseaTable(all.gene.sets[topPathways], prot.list, fgseaRes,gseaParam = 0.5, pathwayLabelStyle=list(size=7),
                 valueStyle = list(size=7),
                 axisLabelStyle = list(size=4))
p

#ggsave(p, width=10, height=8, file="table.png")

plotGseaTable(all.gene.sets[topPathwaysUp], prot.list, fgseaRes,gseaParam = 0.5, pathwayLabelStyle=list(size=7),
              valueStyle = list(size=7),
              axisLabelStyle = list(size=6))

plotGseaTable(all.gene.sets[rev(topPathwaysDown)], prot.list, fgseaRes,gseaParam = 0.5, pathwayLabelStyle=list(size=7),
              valueStyle = list(size=7),
              axisLabelStyle = list(size=6))

head(topPathways)

library(ggplot2)
topPathways
# Make a few Enrichment Plots

#Top 4 up

#Plot1
plotEnrichment(all.gene.sets[["MARKEY_RB1_ACUTE_LOF_DN"]],prot.list) + labs(title="MARKEY_RB1_ACUTE_LOF_DN")
#Plot2
plotEnrichment(all.gene.sets[["CHEN_METABOLIC_SYNDROM_NETWORK"]],prot.list) + labs(title="CHEN_METABOLIC_SYNDROM_NETWORK")
#Plot3
plotEnrichment(all.gene.sets[["LEE_AGING_CEREBELLUM_UP"]],prot.list) + labs(title="LEE_AGING_CEREBELLUM_UP")
#Plot4
#plotEnrichment(all.gene.sets[["WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]],prot.list) + labs(title="WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA")

#Top 4 down

#Plot1
plotEnrichment(all.gene.sets[["MIKKELSEN_NPC_HCP_WITH_H3K27ME3"]],prot.list) + labs(title="MIKKELSEN_NPC_HCP_WITH_H3K27ME3")
#gseaplot(all.gene.sets[["REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"]],prot.list) + labs(title="REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")

#Plot2
plotEnrichment(all.gene.sets[["WP_OXIDATIVE_PHOSPHORYLATION"]],prot.list) + labs(title="WP_OXIDATIVE_PHOSPHORYLATION")
#Plot3
plotEnrichment(all.gene.sets[["REACTOME_RAB_GERANYLGERANYLATION"]],prot.list) + labs(title="REACTOME_RAB_GERANYLGERANYLATION")
#Plot4
#plotEnrichment(all.gene.sets[["LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN"]],prot.list) + labs(title="LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN")
#Plot5
#plotEnrichment(all.gene.sets[["ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP"]],prot.list) + labs(title="ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP")

######################################################################################################################
######################################################################################################################

#==================================================================#
### J_GSEA analysis - Final code - OGT KO MICE BRAIN
#==================================================================#

# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb

library(fgsea)
library(msigdbr)

# Pull all MSigDB collections (caution: this will be a huge list)
all_msigdb_df <- msigdbr(species = "Mus musculus")

#If you want the broadest possible coverage, 
#including the  Hallmark (H) and Reactome (C2) sets that were originally defined in humans.

# Retrieve mouse-native gene sets directly
# (No ortholog mapping required; uses native mouse identifiers)
all_msigdb_df <- msigdbr(species = "Mus musculus", db_species = "MM")

#all_gene_sets2 <- split(x = all_msigdb_df$gene_symbol, f = all_msigdb_df$gs_name)
all_gene_sets2 <- split(x = all_msigdb_df$gene_symbol, f = all_msigdb_df$gs_name)
length(all_gene_sets2) # Should be ~17,000+ vs 33000 for broad coverage
#[1] 17068

library(tidyverse)

#Quick "peek" at the first 6 pathways and the first 6 genes within each.

#all.gene.sets %>% head() %>% lapply(head) #2023 list
all_gene_sets2 %>% head() %>% lapply(head) #updated 2025 list

head(prot.list) 

fgseaRes2 <- fgsea(pathways = all_gene_sets2, stats = prot.list, minSize=15, maxSize=500,eps=0)
head(fgseaRes2) # updated 2025 list
head(fgseaRes2[order(pval), ])
sum(fgseaRes2[, padj < 0.05]) #1316 # updated 2025 list - so better to use

fgseaRes4<-fgseaRes2[order(padj), ]
head(fgseaRes4)

# Make a table plot for a bunch of selected pathways:
#if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)
library(grid)

topPathwaysUp <- fgseaRes2[ES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
#topPathwaysUp <- fgseaRes2[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes2[ES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#with prefixes 
p<-plotGseaTable(all_gene_sets2[topPathways], prot.list, fgseaRes2,gseaParam = 0.5, pathwayLabelStyle=list(size=7, fontface = "bold"), 
                 valueStyle = list(size=7),
                 axisLabelStyle = list(size=4))

# 2. Add the title
grid.arrange(p, top = textGrob("GSEA - CTRL vs. OGT KO MICE BRAIN SAMPLES", 
                               x = 0.62, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))

#Removing prefixes

# 1. Clean the names (remove prefixes)
fgseaRes2[, pathway := gsub("^[^_]*_", "", pathway)]

# 2. Sort by significance (padj) so the best version of a pathway is on top
fgseaRes2 <- fgseaRes2[order(padj)]

# 3. Remove duplicates based on the pathway name
# This keeps ONLY the most significant version of "OXIDATIVE_PHOSPHORYLATION"
#fgseaRes_unique <- fgseaRes2[!duplicated(pathway)]
fgseaRes_unique <- fgseaRes2

# 4. Now selecting Top 10 from the UNIQUE list
topPathwaysUp <- fgseaRes_unique[ES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes_unique[ES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
head(topPathways)
# 5. Making sure the pathway list names match the cleaned names
names(all_gene_sets2) <- gsub("^[^_]*_", "", names(all_gene_sets2))

# 6. Re-generating the plot
p <- plotGseaTable(all_gene_sets2[topPathways], 
                   prot.list, 
                   fgseaRes_unique, 
                   gseaParam = 0.5, 
                   pathwayLabelStyle = list(size=7, fontface = "bold"), 
                   valueStyle = list(size=7),
                   axisLabelStyle = list(size=4))

# 7. Drawing with the title
grid.arrange(p, top = textGrob("GSEA - CTRL vs. OGT KO Mice Brain Total Proteomics", 
                               x = 0.65, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))

head(topPathways)

library(ggplot2)
topPathways
# Make a few Enrichment Plots

#Top 4 up
head(prot.list)
#Plot1
plotEnrichment(all_gene_sets2[["EXTRACELLULAR_SPACE"]],prot.list) + labs(title="EXTRACELLULAR_SPACE")
plotEnrichment(all_gene_sets2[["EXTRACELLULAR_SPACE"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n EXTRACELLULAR_SPACE",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot2
plotEnrichment(all_gene_sets2[["CDC2_IL15_RESPONSE_UP"]],prot.list) + labs(title="CDC2_IL15_RESPONSE_UP")
plotEnrichment(all_gene_sets2[["CDC2_IL15_RESPONSE_UP"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n CDC2_IL15_RESPONSE_UP",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot3
plotEnrichment(all_gene_sets2[["RB1_ACUTE_LOF_DN"]],prot.list) + labs(title="RB1_ACUTE_LOF_DN")
plotEnrichment(all_gene_sets2[["RB1_ACUTE_LOF_DN"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n RB1_ACUTE_LOF_DN",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot4
plotEnrichment(all_gene_sets2[["CDC2_IFNA1_RESPONSE_UP"]],prot.list) + labs(title="CDC2_IFNA1_RESPONSE_UP")
plotEnrichment(all_gene_sets2[["CDC2_IFNA1_RESPONSE_UP"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n CDC2_IFNA1_RESPONSE_UP",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot4
plotEnrichment(all_gene_sets2[["CDC2_IFNA1_RESPONSE_UP"]],prot.list) + labs(title="CDC2_IFNA1_RESPONSE_UP")
plotEnrichment(all_gene_sets2[["CDC2_IFNA1_RESPONSE_UP"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n CDC2_IFNA1_RESPONSE_UP",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

topPathways

#topPathways Down plot

#Plot1
plotEnrichment(all_gene_sets2[["TRANSPORTER_COMPLEX"]],prot.list) + labs(title="TRANSPORTER_COMPLEX")
plotEnrichment(all_gene_sets2[["TRANSPORTER_COMPLEX"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n TRANSPORTER_COMPLEX",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot2
plotEnrichment(all_gene_sets2[["CATION_CHANNEL_COMPLEX"]],prot.list) + labs(title="CATION_CHANNEL_COMPLEX")
plotEnrichment(all_gene_sets2[["CATION_CHANNEL_COMPLEX"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n CATION_CHANNEL_COMPLEX",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot3
plotEnrichment(all_gene_sets2[["VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"]],prot.list) + labs(title="VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE")
plotEnrichment(all_gene_sets2[["VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot3
plotEnrichment(all_gene_sets2[["RESPIRATORY_CHAIN_COMPLEX"]],prot.list) + labs(title="RESPIRATORY_CHAIN_COMPLEX")
plotEnrichment(all_gene_sets2[["RESPIRATORY_CHAIN_COMPLEX"]], prot.list) + 
  labs(
    title = "GSEA - CTRL vs. OGT KO MICE BRAIN_TOTAL PROTEOMICS \n RESPIRATORY_CHAIN_COMPLEX",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )


#----------------------------------------------------------------------------------------------------

#Pretty gsea visualizations    

#library(org.Hs.eg.db)
library(org.Mm.eg.db)
#install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
head(prot.list)
type(prot.list)
library(enrichplot)

head(fgseaRes2)

topPathwaysUp1 <- fgseaRes_unique[ES > 0][head(order(padj), n=10)]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-2]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-(4:5)]
topPathwaysUp1
topPathwaysUp1<-as.data.frame(topPathwaysUp1)
topPathwaysUp1

topPathwaysDown1 <- fgseaRes_unique[ES < 0][head(order(padj), n=10)]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-2]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-(4:5)]
topPathwaysDown1
topPathwaysDown1<-as.data.frame(topPathwaysDown1)
topPathwaysDown1

topPathways1 <- rbind(topPathwaysUp1, rev(topPathwaysDown1)) #merge rows
topPathways1
topPathways1<-as.data.frame(topPathways1)
topPathways1

#https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package

#Pathways with underscore
ggplot(topPathways1, aes(reorder(pathway, NES), NES)) +
  # Removed color="black" to take away the borders around the bars
  geom_col(aes(fill = NES > 0)) + 
  scale_fill_manual(
    values = c("TRUE" = "#FF8C00", "FALSE" = "#4682B4"), 
    labels = c("TRUE" = "Upregulated", "FALSE" = "Downregulated"),
    name = "Direction"
  ) +
  coord_flip() +
  labs(
    x = "Pathway", 
    y = "Normalized Enrichment Score (NES)",
    title = "GSEA: Top 10 Up and Down Regulated Pathways \n CTRL vs. OGT KO Mice Brain Total Proteomics"
  ) + 
  theme_bw() + 
  theme(
    # lineheight = 1.5 creates the space between the two title lines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, lineheight = 1.2),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    # This keeps the border around the whole plot, but not the bars
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

#Pathways without underscore
ggplot(topPathways1, aes(x = reorder(gsub("_", " ", pathway), NES), y = NES)) +
  geom_col(aes(fill = NES > 0)) + 
  scale_fill_manual(
    values = c("TRUE" = "#FF8C00", "FALSE" = "#4682B4"), 
    labels = c("TRUE" = "Upregulated", "FALSE" = "Downregulated"),
    name = "Direction"
  ) +
  coord_flip() +
  labs(
    x = NULL, # Keep the axis label clean
    y = "Normalized Enrichment Score (NES)",
    title = "GSEA: Top 10 Up and Down Regulated Pathways \n CTRL vs. OGT KO Mice Brain Total Proteomics"
  ) + 
  theme_bw() + 
  theme(
    # Centered bold title with that 1.5 line spacing
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, lineheight = 1.2),
    # Bold text for everything else
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

#========================#
#J_HEATMAP
#========================#

head(prot.list)
head(sortupdownOGTTENGFPTEN) #all proteins without removing based on significance
dim(sortupdownOGTTENGFPTEN)

# 1. Filter for adjusted p-value < 0.05
sig_only <- sortupdownOGTTENGFPTEN[sortupdownOGTTENGFPTEN$adj.P.Val.OGTKO.CTRL < 0.045, ]
sig_only

# 2. Top 20 Up of sorted list
top_20_up <- sig_only[sig_only$logFC > 0, ]
head(top_20_up)
top_20_up$logFC <- as.numeric(top_20_up$logFC)
top_20_up <- top_20_up[order(-top_20_up$logFC), ]
head(top_20_up)
top_20_up <- head(top_20_up, 20)
top_20_up
#- sorts from lare to samll (descending order)

# 3. Top 20 Down (Bottom of the list)
top_20_down <- sig_only[sig_only$logFC < 0, ]
head(top_20_down)
top_20_down$logFC <- as.numeric(top_20_down$logFC)
top_20_down <- top_20_down[order(top_20_down$logFC), ]
head(top_20_down)
top_20_down <- head(top_20_down, 20)
top_20_down

# 4. Combining them
heatmap_df <- rbind(top_20_up, top_20_down)
heatmap_df
plot_mat <- cbind(
  LogFC = heatmap_df$logFC,
  Significance = heatmap_df$adj.P.Val.OGTKO.CTRL
)
rownames(plot_mat) <- rownames(heatmap_df)
plot_mat

# 1. Create an empty matrix with the correct dimensions
display_mat <- matrix("", nrow = nrow(plot_mat), ncol = ncol(plot_mat))

# 2. Fill the first column (LogFC)
display_mat[, 1] <- as.character(round(plot_mat[, 1], 2))

# 3. Fill the second column (Significance) 
display_mat[, 2] <- formatC(plot_mat[, 2], format = "e", digits = 2)

display_mat

# 1. Create a matrix for colors where both columns ARE the LogFC
color_mat <- cbind(plot_mat[, 1], plot_mat[, 1]) 

# 2. Creating the Heatmap
library(pheatmap)
pheatmap(color_mat, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         display_numbers = display_mat,  # Shows the text you want
         number_color = "black",
         color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(256),
         # Center the scale so 0 is white
         breaks = seq(-max(abs(plot_mat[,1])), max(abs(plot_mat[,1])), length.out = 257),
         labels_col = c(expression(bold("LogFC")), expression(bold("Adj. P-Value"))),
         angle_col = 0,  # 0 degrees = perfectly horizontal
         main = ("Heat Map - CTRL vs. OGT KO Mice Brain Total Proteomics"))


pheatmap(color_mat, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         # This flips the dendrogram and rows upside down
         clustering_callback = function(hc, mat){
           as.hclust(rev(as.dendrogram(hc)))
         },
         display_numbers = display_mat,
         number_color = "black",
         color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(256),
         breaks = seq(-max(abs(color_mat)), max(abs(color_mat)), length.out = 257),
         labels_col = c(expression(bold("LogFC")), expression(bold("Adj. P-Value"))),
         angle_col = 0,
         main = "Heat Map - CTRL vs. OGT KO Mice Brain Total Proteomics")

#=================================================================================#

#==============================#
#OGT KO Mice - Netwrok Analysis
#==============================#

head(sortupdownOGTTENGFPTEN)
dim(sortupdownOGTTENGFPTEN)
sortupdownOGTTENGFPTEN2<-sortupdownOGTTENGFPTEN
head(sortupdownOGTTENGFPTEN2)
dim(sortupdownOGTTENGFPTEN2)

# 1. Filter for adjusted p-value < 0.05 only
sig_only <- sortupdownOGTTENGFPTEN2[sortupdownOGTTENGFPTEN2$adj.P.Val.OGTKO.CTRL < 0.045, ]
head(sig_only)
dim(sig_only)

# 1. Force the logFC column to be numeric 
sortupdownOGTTENGFPTEN2$logFC.numeric <- as.numeric(as.character(sortupdownOGTTENGFPTEN2$logFC))
head(sortupdownOGTTENGFPTEN2)
dim(sortupdownOGTTENGFPTEN2)
# 2. Apply both filters (Adjusted P < 0.044 AND absolute LogFC > 0.58)
sig_only <- sortupdownOGTTENGFPTEN2[
  !is.na(sortupdownOGTTENGFPTEN2$adj.P.Val.OGTKO.CTRL) & # Remove NAs
    sortupdownOGTTENGFPTEN2$adj.P.Val.OGTKO.CTRL < 0.045 & 
    abs(sortupdownOGTTENGFPTEN2$logFC.numeric) > 0.58, 
]

# 3. Check how many proteins are left
dim(sig_only)
dim(sig_only2)
head(sig_only)
tail(sig_only)

# Get the list of significant gene names
genes_for_ppi <- rownames(sig_only)
head(genes_for_ppi)

library(STRINGdb)
# Initialize for Human (9606) with a medium confidence score (400)
string_db <- STRINGdb$new(version="12.0", species=10090, score_threshold=700)
# Species 10090 is the NCBI taxonomy ID for Mus musculus
#head(string_db)

# 9606 = Taxonomy ID for Homo sapiens (Humans).
#400= Medium confidence
#700= High confidence
#900= highest confidence
#150= low confidence


# Map symbols to STRING IDs
#prot_mapped <- string_db$map(sig_only, "gene_symbol", removeUnmappedRows = TRUE)

# Get the interactions
#interactions <- string_db$get_interactions(prot_mapped$STRING_id)

# 1. Move rownames into a column so STRINGdb can "see" them
sig_only$Symbol <- rownames(sig_only)
head(sig_only)

# 2. Map symbols to STRING IDs (now it will find the 'gene_symbol' column)
prot_mapped <- string_db$map(sig_only, "Symbol", removeUnmappedRows = TRUE)
nrow(prot_mapped)
head(prot_mapped)
# 3. Get the list of mapped IDs for the plot
hits <- prot_mapped$STRING_id
head(hits)
# 4. Draw the network
#string_db$plot_network(hits)

#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("STRINGdb", update = TRUE, ask = FALSE)

# 1. Get the graph object locally
graph <- string_db$get_graph()
library(igraph)
# 2. Extract the subgraph for your 'hits'
# Ensure 'hits' is a character vector of STRING IDs
subgraph <- igraph::induced_subgraph(graph, vids = V(graph)[name %in% hits])

# 3. Plot it using standard R graphics
#plot(subgraph, vertex.label = V(subgraph)$name, vertex.size = 5) #very cluttered, not working

# 1. Map the Gene Symbols to the graph nodes
# This creates a named vector: names = STRING IDs, values = Symbols
symbol_map <- setNames(prot_mapped$Symbol, prot_mapped$STRING_id)

# 2. Assign these symbols as a new attribute to your subgraph
V(subgraph)$label <- symbol_map[V(subgraph)$name]

# 3. Clean up the plot settings for a more professional look

# 1. Use the 'Fruchterman-Reingold' layout with more iterations for better spacing
# niter = 1000 gives the algorithm more time to push nodes apart
l <- layout_with_fr(subgraph, niter = 1000)

# 1. Create a color mapping based on LogFC
# Everything > 0.58 is Red, everything < -0.58 is Blue
node_colors <- ifelse(sig_only$logFC.numeric > 0, "#d7191c", "#2c7bb6")

# 2. Assign these colors to the graph nodes
# Match the order of your 'hits' with the 'sig_only' dataframe
V(subgraph)$color <- node_colors[match(V(subgraph)$name, prot_mapped$STRING_id)]


# 1. Remove "outliers" by selecting only the largest connected component
# This removes all the "floating" dots and focuses on the main interaction hub
components <- clusters(subgraph)
giant_component <- induced_subgraph(subgraph, V(subgraph)[components$membership == which.max(components$csize)])

library(igraph)

# 1. Identify the connected components (the "islands" of proteins)
# 'components()' replaces the deprecated 'clusters()'
comp <- components(subgraph)

# 2. Extract only the "Giant Component" (the largest connected hub)
# This removes all the "floating" outlier dots
giant_component <- induced_subgraph(subgraph, V(subgraph)[comp$membership == which.max(comp$csize)])

# 3. Map Gene Symbols to the labels for readability
V(giant_component)$label <- prot_mapped$Symbol[match(V(giant_component)$name, prot_mapped$STRING_id)]


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("RCy3")
library(RCy3)

cytoscapePing() # Confirms R can "see" the running Cytoscape app
cytoscapePing() # Should return "You are connected to Cytoscape!"
#layoutNetwork('force-directed')
#fitContent()

createNetworkFromIgraph(giant_component, 
                        title = "OGT_KO Mice Brain 2", 
                        collection = "STRING_Analysis")
# 1. Set the background color

setVisualPropertyDefault(list(visualProperty = "NETWORK_BACKGROUND_PAINT", value = "#FFFFFF"))


# 1. Update the table using the STRING IDs as the common key
loadTableData(prot_mapped, 
              data.key.column = "STRING_id", 
              table.key.column = "name")

# 2. Now the column 'logFC.numeric' will exist in Cytoscape. 
# You can run the size mapping:
#setNodeSizeMapping('logFC.numeric', c(0.58, 2.0), c(20, 80))
head(prot_mapped)
tail(prot_mapped)
# 1. Scale Node size by the actual LogFC (Higher LogFC = Bigger Node)
# This highlights the most significant biological changes visually
#setNodeSizeMapping('logFC.numeric', c(0.58, 2.0), c(20, 80))

# 2. Add a border to nodes to make them pop against the white background
setNodeBorderColorDefault('#333333')
setNodeBorderWidthDefault(2)

# 3. Handle Label Overlap (CRITICAL for this many nodes)
# This clears up the "mess" of text so labels don't sit on top of each other
setNodeLabelMapping('label')
setNodeLabelColorDefault('#000000')
setNodeFontSizeDefault(18)
setNodeFontSizeDefault(18)
setNodeShapeDefault('ELLIPSE')

# 2. Make the labels bold to make them stand out
setNodeFontFaceDefault("Arial-Bold")

# 1. Change the global edge color to a light-medium gray
setEdgeColorDefault('#000000')

# This manually sets the 'NODE_LABEL_FONT_FACE' property to Bold
setVisualPropertyDefault(list(visualProperty="NODE_LABEL_FONT_FACE", value="Arial,bold,plain"))

# 2. Make the edges slightly thinner to let the red/blue nodes "pop"
setEdgeLineWidthDefault(1.0)

# -4 (Deep Blue), 0 (Grey/White), +4 (Deep Red)
setNodeColorMapping('logFC.numeric', 
                    table.column.values = c(-4, 0, 4), 
                    colors = c('blue', '#D3D3D3', 'red'), 
                    mapping.type = 'c')

# -4 (Large), 0 (Small), +4 (Large)
setNodeSizeMapping('logFC.numeric', 
                   table.column.values = c(-4, 0, 4), 
                   sizes = c(80, 20, 80), 
                   mapping.type = 'c')

# Fit everything back into the window view
fitContent()


layoutNetwork('force-directed')

# 2. Reset the scale if things are still too far apart
# This fits the whole network perfectly into the window
fitContent()



