#========================================================#
#J_ImmVar data analysis using varianceparition in R
#By: Sophiya J. Hanigan
#========================================================#

#--------------------#
#Notes
#--------------------#

#Data download from GEO -  NCBI Gene Expression Omnibus

#GEO: https://www.ncbi.nlm.nih.gov/gds/?term=Immvar
#ImmVar data = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56035
#From the Paper: https://pubmed.ncbi.nlm.nih.gov/24786080/ : GEO/GSE56035
#CD4 T Cells
#CD14 Monocytes

#Series GSE56035: Superseries
#Includes the following subseries:
  #GSE56033	Immune Variation Project (ImmVar) [CD4]
  #GSE56034	Immune Variation Project (ImmVar) [CD14]

#https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.R

#================================================================================================================#

#================================================================================#
#Download the dataset
#=================================================================================#
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")

#if (!require("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
  #BiocManager::install("variancePartition")

  #if (!require("BiocManager", quietly = TRUE))
   # install.packages("BiocManager")
  #BiocManager::install("BiocParallel")

#------------------#
# Attach libraries
#------------------#
library(GEOquery)
library(variancePartition) 
library(BiocParallel) #For parallel computing
  #https://www.bioconductor.org/packages//release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#introduction
#Configuring Paralellel computation threads:
#BiocParallel Objective is to reduce the complexity faced when developing and using software that performs parallel computations. 

#parallel::detectCores() #My computer has 10 cores; I use windows
#param<-SnowParam(workers=4, progressbar = TRUE)
#register(param)

#snow= simplenetwork of workstation; used to speed up analysis time by giving multiple background R Workstations 
#to process genes in chunks instead of a single workstation processing 1 by 1
#workers=10, uses 10 cores instead of 2
#progressbar=TRUE shows time left to finish analysis
#register - activates the variable param
library(BiocManager)
library(readxl)
library(tidyverse)
library(limma)
library(EnvStats)   # used to get geometric means
library(missForest) # Imputation
library(biomaRt) # library for mapping between annotations
library(qvalue)

#-----------------------------------------------------------#
#Downloading GSE superseries: GSE56035
#-----------------------------------------------------------#

GSEsuper<-getGEO("GSE56035",GSEMatrix=TRUE)
#getGEO - A character string representing a GEO object for download and parsing. (eg., 'GDS505','GSE2','GSM2','GPL96')
#https://www.rdocumentation.org/packages/GEOquery/versions/2.38.4/topics/getGEO
head(GSEsuper)
print(names(GSEsuper))
print((GSEsuper))
#Obtained a list
#Extract very first item of list
expr_element<-GSEsuper[[1]] 
head(expr_element)
dim(expr_element)
#Features  Samples 
#21727      984

#-----------------------------------------------------------#
#Extracting Expression Matrix (Matrix of gene expression)
#-----------------------------------------------------------#

expr_matrix<-exprs(expr_element)
#exprs returns a (usually large!) matrix of expression values; se.exprs returns the corresponding matrix of standard errors, when available.
head(expr_matrix)
dim(expr_matrix)
#[1] 21727   984

#---------------------------------------------#
#Assess missing values in Expression Matrix
#---------------------------------------------#

table(is.na(expr_matrix))
# FALSE 
#21379368 
table(is.nan(expr_matrix))
#  FALSE 
#21379368
type(expr_matrix)
#double
#To see if there are any rows with all columns 0
head(expr_matrix)
ncol(expr_matrix) #[1] 984 #GEO Sample Accesion ID (no of patients)
nrow(expr_matrix) #21727 Microarray probe ID or genes
head(rowSums(expr_matrix==0))
#rowSums counts 0s in each row separately
#https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/rowsum
zerorow<-sum(rowSums(expr_matrix==0)==ncol(expr_matrix))

#where sum counts total no. of TRUEs
zerorow

#No gene with 0 across all columns

#-----------------------------------------------------------#
#Exptracting metadata
#-----------------------------------------------------------#

metadata<-pData(expr_element)
#where, pData used to extract phenotyping data
head(metadata)
colnames(metadata)
head(metadata$source_name_ch1)
table(metadata$source_name_ch1)
#healthy individual_CD4 T-Cells from PBMCs   healthy individual_monocytes from PBMCs 
#499                                         485
head(metadata)
table(metadata$title)
table(metadata$"cell type")
#monocytes from peripheral blood mononuclear cells (PBMCs) 
#485 
#T4 Naive cells from human peripheral blood mononuclear cells (PBMCs) 
#499 

##Thus we have 499 TCells and 485 Monocytes all from healthy individual

table(metadata$"age (yrs)")
table(metadata$"batch")
table(metadata$gender)
#female   male 
#269    184
#Had to combine gender and sex (did it later downstream)
#And other metadat like age, batch, gender etc.

head(expr_matrix)
#sampleNames: GSM1350856 GSM1350857 ... GSM1355987 (984 total)
#featureNames: 7896740 7896742 ... 7896756 (6 total) #featurename = probeID for microarray data
#https://support.bioconductor.org/p/117804/

head(expr_matrix) #Has the gene expression
head(metadata) #Has the metadata

#================================================================================================================#

#===============================================#
#Log transformation of expr_matrix
#===============================================#

#expr_matrix_log<-log2(expr_matrix+1)
#Warning message:
#NaNs produced

summary(as.vector(expr_matrix))
#summary(as.vector(expr_matrix))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-1359.05    57.93   132.75   353.14   312.44 25496.44 

#Values are raw data, have negative values
#Microarrays different than rna seq nad proteomics, so has -ve values because of background noise
#Signal = Intebsity-background; If background high can have a -ve no.

expr_matrix_filter<-expr_matrix
expr_matrix_filter[expr_matrix_filter<1]<-1 #Substitute values less than 1 with 1
#https://www.researchgate.net/post/Negative-number-in-gene-expression-microarray-any-thoughts
expr_matrix_log<-log2(expr_matrix_filter)
head(expr_matrix_log)
table(is.na(expr_matrix_log))
#FALSE 
#21379368

#So log transformed and no missing values 

summary(as.vector(expr_matrix_log))
#Now min = 0

summary((expr_matrix_log[,1:4]))
summary((expr_matrix[,1:4]))

colnames(metadata) 
cat(metadata$data_processing[1])

#The datasets was pre-filtered to keep only those probesets for which a 
#gene symbol could be found in the Affymetrix annotation. 
#CEL files were normalized using Affymetrix Power Tools on the predefined 
#probeset ID list mentioned above, and using the standard RMA workflow 
#(background adjustment, quantile normalization, median polish probeset 
#summarization). Output on 2^ scale.

#So the downloaded data is already quantile normalized! 
#Output on 2^ scale. means the log transformation was reversed

#================================================================================================================#

#======================#
#QC and Normalization
#======================#

#-------------------------------------------------------------------------------------------------#
#Boxplot prior to Normalization
#-------------------------------------------------------------------------------------------------#

col <- c("red","blue","darkgreen","purple","orange","pink","green","maroon","violet","deeppink")
col
boxplot(expr_matrix_log[,1:10], main="Boxplot to check normalization of data",col=col)
#boxplot(expr_matrix[,1:10], main="ImmVar_Boxplot to check Normalization",col=col)

#But data doesn't look very normalized in box plot, could be as the papers said this data was processed batch to batch

#I am not sure if I should run another quantile normalization
#As metadata data processing mentions data is quantile normalized
#Maybe since the data was changed back to 2^ scale normalization needed again?

#================#
#Normalization
#================#

#-------------------------------------------------------------------------------------------------#
# Quantile normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------#

expr_matrix_log_quant = limma::normalizeQuantiles(expr_matrix_log) # without NAs, log

#Quantile normalization is preferred for microarray data
#Because of dense fluorescence background intensities used
# Quantile normalization assumes most changes is due to technical variation rather than biological

apply(expr_matrix_log[,1:5],2,median)
apply(expr_matrix_log_quant[,1:5],2,median)

#>apply(expr_matrix_log[,1:5],2,median)
#GSM1350856 GSM1350857 GSM1350858 GSM1350859 GSM1350860 
#6.897158   6.872997   6.905648   6.825854   6.905595 
#> apply(expr_matrix_log_quant[,1:5],2,median)
#GSM1350856 GSM1350857 GSM1350858 GSM1350859 GSM1350860 
#7.054449   7.054449   7.054449   7.054449   7.054449 

#-------------------------------------------------------------------------------------------------#
# Median normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------#

expr_matrix_log_median = limma::normalizeMedianValues(expr_matrix_log) # without NAs, log

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#Boxplot After Normalization
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

boxplot(expr_matrix_log_quant[,1:10], main="Boxplot to check data_post quantile normalization",col=col)

boxplot(expr_matrix_log_median[,1:10], main="Boxplot to check data_post median normalization",col=col)

#Data looks much better and normalized post quantile and median normalization.

#We will use quantile normalized data for further analysis
#Quantile normalization is preferred for microarray data
#Because of dense fluorescence background intensities used
# Quantile normalization assumes most changes is due to technical variation rather than biological

#----------------------------------------------------------------------#
# PCA - Principal Component Analysis (dynamic titles instead of fixed)
# Plotted for Batch and Cell separately
#----------------------------------------------------------------------#

get_pca <- function(dat, labels, legend_title = 'Group', plot_title = 'PCA Plot') {
  require(ggplot2)
  
  #Running PCA on transposed matrix
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  #Extracting percentage variance for axes labels
  var_explained<-summary(pca_res)$importance[2,1:2]*100
  # R generates a summary table when we run PCA called 
  #Importance of components that has
  #Row 1 equal to standard deviation, 
  #row 2 equal to Proportion of variance, say .456
  #row 3 equal to cumulative proportion
  #[2,1:2]=row 2, column 1-2 (PC1 and PC2)
  # *100 converts to % say proport .456*100 = 45.6%
  
  lbl_pc1<-paste0("PC1(",round(var_explained[1],1),"%)")
  #var_explained[1] = PC1 say 45.6
  #prints as PC1 45.6%
  lbl_pc2<-paste0("PC2(",round(var_explained[2],1),"%)")
  
  
  #Plot PCA
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    #where, U <- data.frame(pca_res$x)
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(x=lbl_pc1, y=lbl_pc2, color = legend_title, title=plot_title)+
    #title=plot_title - sets main heading at top of chart
    #color = legend_title - sets header for colors as whatever the legend title we give cell type/batch etc.
    
    scale_color_brewer(palette="Set1")+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

#ImmVar metadata
colnames(metadata)
Cell_labels<-as.character(metadata$source_name_ch1) #Separates CD4 vs Monocytes
Batch_labels<-as.character(metadata$`batch:ch1`) #Separates CD4 vs Monocytes

#head(dat1)
#get_pca2 <- function(dat, labels, legend_title = 'Group', plot_title = 'PCA Plot') 
PCA_raw_cell<-get_pca(expr_matrix,Cell_labels,"Cell Type", "Principal Component Analysis (PCA) Based on Cell Type\n Before log2 transformation and normalization")
PCA_raw_cell$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/
PCA_raw_batch<-get_pca(expr_matrix,Batch_labels,"Batch Effect", "Principal Component Analysis (PCA) Based on Batch\n Before log2 transformation and normalization")
PCA_raw_batch$p

PCA_log_cell<-get_pca(expr_matrix_log,Cell_labels,"Cell Type", "Principal Component Analysis (PCA) Based on Cell Type \n After log2 transformation and before normalization")
PCA_log_cell$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/
PCA_log_batch<-get_pca(expr_matrix_log,Batch_labels,"Batch Effect", "Principal Component Analysis (PCA) Based on Batch \n After log2 transformation and Before normalization")
PCA_log_batch$p

PCA_quant_cell<-get_pca(expr_matrix_log_quant,Cell_labels,"Cell Type", "Principal Component Analysis (PCA) Based on Cell Type \n After log2 transformation and Quantile normalization")
PCA_quant_cell$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/
PCA_quant_batch<-get_pca(expr_matrix_log_quant,Batch_labels,"Batch Effect", "Principal Component Analysis (PCA) Based on Batch \n After log2 transformation and Quantile normalization")
PCA_quant_batch$p

#PCA_median_cell<-get_pca(expr_matrix_log_quant,Cell_labels,"Cell Type", "Principal Component Analysis (PCA) Based on Cell Type \n After log2 transformation and Median normalization")
#PCA_median_cell$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/
#PCA_median_batch<-get_pca(expr_matrix_log_quant,Batch_labels,"Batch Effect", "Principal Component Analysis (PCA) Based on Batch \n After log2 transformation and Median normalization")
#PCA_median_batch$p

#-------------------------------#
#Batch and Cell combo PCA plot
#-------------------------------#

get_pca_combo <- function(dat, color_labels, shape_labels, Color_title = 'Batch', Shape_title = 'Cell_Type', plot_title = 'PCA Plot') {
  require(ggplot2)
  
  #Running PCA on transposed matrix
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  #Extracting percentage variance for axes labels
  var_explained<-summary(pca_res)$importance[2,1:2]*100
  # R generates a summary table when we run PCA called 
  #Importance of components that has
  #Row 1 equal to standard deviation, 
  #row 2 equal to Proportion of variance, say .456
  #row 3 equal to cumulative proportion
  #[2,1:2]=row 2, column 1-2 (PC1 and PC2)
  # *100 converts to % say proport .456*100 = 45.6%
  
  lbl_pc1<-paste0("PC1(",round(var_explained[1],1),"%)")
  #var_explained[1] = PC1 say 45.6
  #prints as PC1 45.6%
  lbl_pc2<-paste0("PC2(",round(var_explained[2],1),"%)")
  
  
  #Plot PCA
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = color_labels, shape=shape_labels)) +
    #where, U <- data.frame(pca_res$x)
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(x=lbl_pc1, y=lbl_pc2, color = Color_title, shape=Shape_title, title=plot_title)+
    #title=plot_title - sets main heading at top of chart
    #color = legend_title - sets header for colors as whatever the legend title we give cell type/batch etc.
    
    scale_color_brewer(palette="Set1")+ #Can only use for upto 9
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca


PCA_raw_cell_batch<-get_pca_combo(dat=expr_matrix, 
                                  color_labels=Batch_labels, 
                                  shape_labels=Cell_labels, 
                                  Color_title = 'Batch Effect', 
                                  Shape_title = 'Cell Type', 
                                  plot_title = 'Principal Component Analysis (PCA) Based on Batch and Cell Type \n Before log2 transformation and normalization')

PCA_raw_cell_batch$p

# scale_shape_manual(values=c(16,17)) #16=circle, 17=triangle, was default here
#https://www.datanovia.com/en/blog/ggplot-point-shapes-best-tips/

PCA_log_cell_batch<-get_pca_combo(dat=expr_matrix_log, 
                                  color_labels=Batch_labels, 
                                  shape_labels=Cell_labels, 
                                  Color_title = 'Batch Effect', 
                                  Shape_title = 'Cell Type', 
                                  plot_title = 'Principal Component Analysis (PCA) Based on Batch and Cell Type \n After log2 transformation and before normalization')
PCA_log_cell_batch$p

PCA_quant_cell_batch<-get_pca_combo(dat=expr_matrix_log_quant, 
                                  color_labels=Batch_labels, 
                                  shape_labels=Cell_labels, 
                                  Color_title = 'Batch Effect', 
                                  Shape_title = 'Cell Type', 
                                  plot_title = 'Principal Component Analysis (PCA) Based on Batch and Cell Type \n After log2 transformation and quantile normalization')
PCA_quant_cell_batch$p

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#Density plots
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

head(expr_matrix_filter)#No missing values and w/o duplicate SiteID (median summarized); is a dataframe; Not matrix and not log transformed; not normalized; need for pca
dim(expr_matrix_filter)#[1] 12298   18 #No missing values and w/o duplicate SiteID (median summarized); is a dataframe; Not matrix and not log transformed; not normalized; need for pca

#pdat22.2 = not log transformed, not normalized, omitted NA
for (i in 1:ncol(expr_matrix_filter)){
  if (i==1){
    plot(density(log10(expr_matrix_filter[,1])),main="Density plot - ImmVar \n Not log transformed, not normalized, omitted NA", xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(expr_matrix_filter[,i]))
    lines(den$x,den$y,col="blue")}
} 

#Already saved image

#pdat23_mat_norm = Median normalization, matrix, log2 transformed, removed missing values and duplicates
for (i in 1:ncol(expr_matrix_log_quant)){
  if (i==1){
    plot(density((expr_matrix_log_quant[,1])),main="Density plot \n Quantile normalized, log2 transformed, Removed missing values and duplicates",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((expr_matrix_log_quant[,i]))
    lines(den$x,den$y,col="blue")}
}

for (i in 1:ncol(expr_matrix_log_quant)){
  if (i==1){
    plot(density((expr_matrix_log_quant[,1])),main="Density plot \n Quantile normalized, Palette colored, log2 transformed, Removed missing values and duplicates",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((expr_matrix_log_quant[,i]))
    lines(den$x,den$y,col= rainbow(ncol(expr_matrix_log_quant))[i])}
}

#pdat24_mat_norm = Quantile normalization, matrix, log2 transformed, removed missing values and duplicates

for (i in 1:ncol(expr_matrix_log_median)){
  if (i==1){
    plot(density((expr_matrix_log_median[,1])),main="Density plot  \n Median normalized, log2 transformed, Removed missing values and duplicates",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((expr_matrix_log_median[,i]))
    lines(den$x,den$y,col="blue")}
}

head(expr_matrix_log)
for (i in 1:ncol(expr_matrix_log)){
  if (i==1){
    plot(density((expr_matrix_log[,1])),main="Density plot  \n  log2 transformed, Removed missing values and duplicates, Not normalized",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((expr_matrix_log[,i]))
    lines(den$x,den$y,col="blue")}
}

#Roughly normal distribution
#Can use parametric linear models
#However here since we are comparing multiple groups, here the cell types T cells and monocytes
#from same person we will use lmm
#Because different cell types but with shared genetic background (from same person)

#======================================================================================#
#-----------------------#
#Preparing metadata
#-----------------------#
head(metadata)
colnames(metadata)
head(metadata$'cell type:ch1')
table(metadata$'cell type:ch1')
head(metadata$'Sex:ch1')
head(metadata$'gender:ch1')
head(metadata$'age (yrs):ch1')
head(metadata$'characteristics_ch1')
head(metadata$'characteristics_ch1.1')
head(metadata$'batch:ch1')
head(metadata$'phenotype markers:ch1')
head(metadata$'title')
head(expr_matrix_log)
#Extracting individual IDs from Title
Individual_id=as.factor(gsub("\\..*","",metadata$'title'))
head(Individual_id)
clean_metadata<-data.frame(
  CellType=as.factor(metadata$'cell type:ch1'),
  #If metadata$title has T4 Naive , gives Cell type = T_cell else gives Monocyte
  #grepl = Regular expression logical, searches through a list of text a particular pattern and gives true if found, False if not found
  Gender=as.factor(metadata$'gender:ch1'),
  Batch=as.factor(metadata$'batch:ch1'),
  Age=as.numeric(as.character((metadata$'age (yrs):ch1'))),
  Pheno=as.factor(metadata$'phenotype markers:ch1'),
  Individual=as.factor(Individual_id),
  
  #Check papaer for other like individual
  row.names = colnames(expr_matrix_log)
)

table(table(clean_metadata$Individual))

head(clean_metadata)

#Linear Mixed Model Formula

#continuous fixed effects - plain; categorical random effects (1|Variable) 
#lmm<-~(1|Individual)+(1|CellType)+(1|Gender)+(1|Batch)+Age

#No need Pheno as Cell Type sorted based on pheno
#Error


#Example:
# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
#form <- ~ Age + (1 | Individual) + (1 | Tissue) + (1 | Batch)
#https://bioconductor.posit.co/packages/3.23/bioc/vignettes/variancePartition/inst/doc/variancePartition.html

#--------------------#                                                   
#VariancePartition
#--------------------#
#varpartlog<-fitExtractVarPartModel(expr_matrix_log,lmm,clean_metadata,BPPARAM=param)
#varpartquant<-fitExtractVarPartModel(expr_matrix_log_quant,lmm,clean_metadata,BPPARAM=param)
#Gave error
#Variables contain NA's: Gender, Age 
#Model failed

#So need to filter samples to kepp those with no NA

colnames(metadata)
head(metadata$'gender:ch1')
head(metadata$'age (yrs):ch1')
Gender=as.factor(metadata$'gender:ch1')
table(is.na(Gender))
#FALSE  TRUE 
#453   531  # A lot missing
Sex=as.factor(metadata$'Sex:ch1')
table(is.na(Sex))
#table(is.na(Sex))
#FALSE  TRUE 
#  531   453 

#But gender False 453 and Sex False = 531 which is true for na in the other set
#Maybe problem in data curation

Combined_gender<-ifelse(is.na(metadata$'gender:ch1'),
                        as.character(metadata$'Sex:ch1'),
                        as.character(metadata$'gender:ch1'))
head(Combined_gender)
table(is.na(Combined_gender))
#FALSE 
#984

Age=as.numeric(as.character((metadata$'age (yrs):ch1')))
table(is.na(Age))
#FALSE  TRUE 
#938    46 

#Need to filter Age as well

# Age=as.numeric(as.character((metadata$'age (yrs):ch1'))),
table(is.na(Age)) #46 missing Age
#FALSE  TRUE 
#938    46 
938+46 #984
Clean_Age<-!(is.na(metadata$`age (yrs):ch1`))
head(Clean_Age)
length(Clean_Age) #984

clean_metadata_2<-data.frame(
  CellType=as.factor(metadata$'cell type:ch1'),
  #If metadata$title has T4 Naive , gives Cell type = T_cell else gives Monocyte
  #grepl = Regular expression logical, searches through a list of text a particular pattern and gives true if found, False if not found
  Gender=as.factor(Combined_gender), #---Used combined gender
  Batch=as.factor(metadata$'batch:ch1'),
  Age=as.numeric(as.character((metadata$'age (yrs):ch1'))),
  Pheno=as.factor(metadata$'phenotype markers:ch1'),
  Individual=as.factor(Individual_id),
  
  #Check papaer for other like individual
  row.names = colnames(expr_matrix_log)
)


head(clean_metadata_2)
dim(clean_metadata_2) #984   6
clean_metadata_filtered<-clean_metadata_2[Clean_Age,]
head(clean_metadata_filtered)
dim(clean_metadata_filtered)
#938  6
colnames(clean_metadata_filtered)

#-----------------------------#
#Quantile Normalized data
#-----------------------------#

head(expr_matrix_log_quant)
dim(expr_matrix_log_quant)
#[1] 21727   984

expr_matrix_log_quant_filter<-expr_matrix_log_quant[,Clean_Age]
head(expr_matrix_log_quant_filter)  
dim(expr_matrix_log_quant_filter)  
#[1] 21727   938

#--------------------------------#
#Data prior to Normalization
#--------------------------------#

head(expr_matrix_log)
dim(expr_matrix_log)
#[1] 21727   984

expr_matrix_log_filter<-expr_matrix_log[,Clean_Age]
head(expr_matrix_log_filter)  
dim(expr_matrix_log_filter)  
#[1] 21727   938

#--------------------#                                                   
#lmm
#--------------------#
colnames(clean_metadata_filtered)
lmm<-~(1|Individual)+(1|CellType)+(1|Gender)+(1|Batch)+Age

#====================#                                                   
#VariancePartition
#====================#

#-------------------#
#Parallel connction
#-------------------#

#closeAllConnections()
#parallel::detectCores() #My computer has 10 cores; I use windows
#param<-SnowParam(workers=6, progressbar = TRUE)
#register(param)
#options(warn=0)
#varpartquant<-fitExtractVarPartModel(expr_matrix_log_quant_filter,lmm,clean_metadata_filtered,BPPARAM=param)

#Error So I tried Serial Param, but must fix parallel param for increased efficiency

#-------------------#
#Serial Connection
#-------------------#

options(warn=0)
closeAllConnections()
parallel::detectCores() #My computer has 10 cores; I use windows
serial_param<-SerialParam(progressbar = TRUE)
#register(param)

#------------------------------------------------------#
#VariantPartition Analysis - Quantile Normalized Data
#------------------------------------------------------#

varpartquant<-fitExtractVarPartModel(expr_matrix_log_quant_filter,lmm,clean_metadata_filtered,BPPARAM=serial_param)

#Sort and Plot

varpartsorted_quant<-sortCols(varpartquant)
plotVarPart(varpartsorted_quant)

library(ggplot2)
plotVarPart(varpartsorted_quant)+
 labs(title = ('VariantPartition Analysis - Quantile Normalized Data'))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) # This centers the title
    
#------------------------------------------------------------------#
#VariantPartition Analysis - Data prior to Quantile Normalization
#------------------------------------------------------------------#

options(warn=0)
head(expr_matrix_log_filter)
dim(expr_matrix_log_filter)
varpartquant_prior<-fitExtractVarPartModel(expr_matrix_log_filter,lmm,clean_metadata_filtered,BPPARAM=serial_param)

#Sort and Plot

varpartsorted_prior<-sortCols(varpartquant_prior)
plotVarPart(varpartsorted_prior)

library(ggplot2)
plotVarPart(varpartsorted_prior)+
  labs(title = ('VariantPartition Analysis - Data Prior to Quantile Normalization'))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) # This centers the title

#================================================================================================================================#

#-----------------------------------------------------------------------#
#Mapping individual variables to gene names - Quantile Normalized data
#-----------------------------------------------------------------------#

head(expr_matrix_log_quant_filter)
dim(expr_matrix_log_quant_filter)
#[1] 21727   938
colnames(metadata)
head(metadata$title)
head(metadata$geo_accession)
head(metadata$data_processing)
#Affymetrix annotation
#The datasets was pre-filtered to keep only those probesets for which a gene
#symbol could be found in the Affymetrix annotation
head(metadata$type)

#RNA
head(colnames(expr_matrix_log_quant_filter))
head(metadata$geo_accession)


#Match
#geo_accession=patient ID
#GEO=Gene Expression Omnibus
#GSM=GEO Sample

head(rownames(expr_matrix_log_quant_filter))
# "7896740" "7896742" "7896744" "7896750" "7896754" "7896756"

#Micro array Probe_ID
#Platforms (1)	
#GPL6244	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56035

#Affymetrix Human Gene 1.0-ST Array Transcriptcluster Revision 8 annotation data (chip hugene10sttranscriptcluster)
#https://www.bioconductor.org/packages//2.10/data/annotation/html/hugene10sttranscriptcluster.db.html

#if (!require("BiocManager", quietly = TRUE))
 #install.packages("BiocManager")
#BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

#--------------------#
#Extracting probe ID
#--------------------#

head(expr_matrix_log_quant_filter)
dim(expr_matrix_log_quant_filter)
#[1] 21727   938
dim(expr_matrix_log_quant)
#[1] 21727   984

probe_id<-rownames(expr_matrix_log_quant_filter)
head(probe_id)
length(probe_id)
#[1] 21727


#--------------------#
#Map to gene symbol
#--------------------#

columns(hugene10sttranscriptcluster.db)
#mapIds: https://support.bioconductor.org/p/9138389/
gene_symbol<-mapIds(hugene10sttranscriptcluster.db,
                    keys = probe_id,
                    column = "SYMBOL",
                    keytype = "PROBEID",
                    multiVals = "first") #gives 1st symbol for IDs mapping to multiple symbols
head(gene_symbol)                    
# 7896740      7896742      7896744      7896750      7896754      7896756 
#"OR4F4"     "PCMTD2"     "OR4F29"           NA "SEPTIN7P13"     "FAM87B"
length(gene_symbol)
#21727

#Need to filter missing values

Symbol_not_NA<-!(is.na(gene_symbol))
head(Symbol_not_NA)

gene_symbol2<-gene_symbol[Symbol_not_NA]
head(gene_symbol2)
length(gene_symbol2)
#[1] 21219
gene_symbol2<-as.data.frame(gene_symbol2)
head(gene_symbol2)
gene_symbol2$Probe_ID<-rownames(gene_symbol2)
head(gene_symbol2)
colnames(gene_symbol2)<-c("Gene_Symbol","Probe_ID")
head(gene_symbol2)
gene_symbol2 <- gene_symbol2 %>%
  relocate(Gene_Symbol, .after = Probe_ID)
head(gene_symbol2)
rownames(gene_symbol2)<-NULL
head(gene_symbol2)

table(duplicated(gene_symbol2$Probe_ID))
#FALSE 
#21219 

table(duplicated(gene_symbol2$Gene_Symbol))
#FALSE  TRUE 
#19837  1382 

#Thus there are duplicated Symbols


#-------------------------------#
#Replace probeID with Symbol
#-------------------------------#

head(expr_matrix_log_quant_filter)
dim(expr_matrix_log_quant_filter)
#[1] 21727   938
head(rownames(expr_matrix_log_quant_filter))
expr_matrix_log_quant_map<-expr_matrix_log_quant_filter
head(expr_matrix_log_quant_map)
dim(expr_matrix_log_quant_map)
#[1] 21727   938
head(gene_symbol2) 
dim(gene_symbol2) 
dim(expr_matrix_log_quant_map)
#[1] 21727   938

#Filter expr_matrix_log_quant_map
expr_matrix_log_quant_map<-expr_matrix_log_quant_map[gene_symbol2$Probe_ID,] #only rows that are in gene_symbol2 kept
head(expr_matrix_log_quant_map)
dim(expr_matrix_log_quant_map)
#[1] 21219   938
dim(gene_symbol2)
#[1] 21219     2
#Both match
head(expr_matrix_log_quant_map)
library(limma)
head(gene_symbol2)
#Average Duplicate gene symbol for microarray data
expr_matrix_mapped_final<-avereps(expr_matrix_log_quant_map, ID=gene_symbol2$Gene_Symbol) 
head(expr_matrix_mapped_final)
dim(expr_matrix_mapped_final)
#[1] 19837   938
head(rownames(expr_matrix_mapped_final))

#-------------------------------------------------------------------#
#VariantPartition Analysis - Quantile Normalized and Mapped Data
#-------------------------------------------------------------------#

head(expr_matrix_mapped_final)
dim(expr_matrix_mapped_final)
#[1] 19837   938
head(clean_metadata_filtered)
dim(clean_metadata_filtered)
#[1] 938   6
dim(expr_matrix_mapped_final)
#[1] 19837   938
lmm
varpartquant_mapped<-fitExtractVarPartModel(expr_matrix_mapped_final,lmm,clean_metadata_filtered,BPPARAM=serial_param)

#Sort and Plot

varpartsorted_mapped<-sortCols(varpartquant_mapped)
plotVarPart(varpartsorted_mapped)

library(ggplot2)
plotVarPart(varpartsorted_mapped)+
  labs(title = ('VariantPartition Analysis - Quantile Normalized Data After Mapping to Gene Names'))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) # This centers the title

#==============================================================================================================#
#======================================================#
#Variant Partition Analysis for Individual genes
#======================================================#

head(expr_matrix_mapped_final)
dim(expr_matrix_mapped_final)
head(rownames(expr_matrix_mapped_final))
Random6<-head(rownames(expr_matrix_mapped_final))
Random6

#---------------------------------------------------#
#Individual variant analysis for First 6 genes
#---------------------------------------------------#

plotPercentBars(varpartquant_mapped[Random6,])
plotPercentBars(varpartquant_mapped[Random6,])+
  labs(title = ('VariantPartition Analysis - For The First 6 Individual Genes, after quantile normalization'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold",
    margin = margin(t=10)))

#---------------------------------------------------#
#Top 2 genes from each metadata category
#---------------------------------------------------#

head(varpartquant_mapped)
varpartquant_mapped_df<-as.data.frame(varpartquant_mapped)
head(varpartquant_mapped_df)
type(varpartquant_mapped_df$CellType)
colnames(clean_metadata)
lmm

top2_all<-c(rownames(varpartquant_mapped_df[order(varpartquant_mapped_df$CellType, decreasing = TRUE),])[1:2],
         rownames(varpartquant_mapped_df[order(varpartquant_mapped_df$Individual, decreasing = TRUE),])[1:2],
         rownames(varpartquant_mapped_df[order(varpartquant_mapped_df$Age, decreasing = TRUE),])[1:2],
         rownames(varpartquant_mapped_df[order(varpartquant_mapped_df$Gender, decreasing = TRUE),])[1:2],
         rownames(varpartquant_mapped_df[order(varpartquant_mapped_df$Batch, decreasing = TRUE),])[1:2]
         )  
top2_all

plotPercentBars(varpartquant_mapped[top2_all,])+
  labs(title = ('VariantPartition Analysis - For Top2 Genes from Each Metadata Group, after quantile normalization'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold",
    margin = margin(t=10)))+
  guides(fill=guide_legend(title = "Metadata Groups"))

#---------------------------------------------------#
#Genes from paper
#---------------------------------------------------#

Paper<-c('AGPAT1','KARS','TLR4','DUSP6','GSTM1','ZFP57','RPS4Y1','UTY','SULT1E1') 
Paper
plotPercentBars(varpartquant_mapped[Paper,])+
  labs(title = ('VariantPartition Analysis - For Genes from VariancePartition Paper, after quantile normalization'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold",
    margin = margin(t=10)))+
  guides(fill=guide_legend(title = "Metadata Groups"))

#Batch effect completely removed (different from seen in papaer)

#=====================================================================================================#

#-------------------------------------#
#PlotStratifyby
#-------------------------------------#

#------------#
#1. CellType
#------------#

#1. GLT1D1

GLT1D1<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["GLT1D1",]),
                   CellType = clean_metadata_filtered$CellType)
head(GLT1D1)
dim(GLT1D1)
#938  2
plotStratifyBy(GLT1D1,"CellType","Expression") + #(dataframe,x,y)+
labs(title = ('GLT1D1 Stratfication Plot_Cell Type'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#2. CPVL

CPVL<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["CPVL",]),
                   CellType = clean_metadata_filtered$CellType)
head(CPVL)
dim(CPVL)
#938  2
plotStratifyBy(CPVL,"CellType","Expression") + #(dataframe,x,y)+
  labs(title = ('CVPL Stratfication Plot_Cell Type'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#------------#
#2. Gender
#------------#

#1.PRKY

PRKY<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["PRKY",]),
                 Gender = clean_metadata_filtered$Gender)
head(PRKY$Gender)
PRKY$Gender<-toupper(PRKY$Gender)
PRKY$Gender<-as.factor(PRKY$Gender)
head(PRKY)
dim(PRKY)
#938  2
plotStratifyBy(PRKY,"Gender","Expression") + #(dataframe,x,y)+
  labs(title = ('PRKY Stratfication Plot_Gender'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#2.KDM5D

KDM5D<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["KDM5D",]),
                 Gender = clean_metadata_filtered$Gender)
head(KDM5D$Gender)
KDM5D$Gender<-toupper(KDM5D$Gender)
KDM5D$Gender<-as.factor(KDM5D$Gender)
head(KDM5D)
dim(KDM5D)
#938  2
plotStratifyBy(KDM5D,"Gender","Expression") + #(dataframe,x,y)+
  labs(title = ('KDM5D Stratfication Plot_Gender'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#----------#
#3. Batch
#----------#

#1.JCHAIN

JCHAIN<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["JCHAIN",]),
                   Batch = as.factor(clean_metadata_filtered$Batch))
head(JCHAIN)
dim(JCHAIN)
#938  2
plotStratifyBy(JCHAIN,"Batch","Expression") + #(dataframe,x,y)+
  labs(title = ('JCHAIN Stratfication Plot_Batch'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))


#2.DDX3Y

DDX3Y<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["DDX3Y",]),
                  Batch = as.factor(clean_metadata_filtered$Batch))
head(DDX3Y)
dim(DDX3Y)
#938  2
plotStratifyBy(DDX3Y,"Batch","Expression") + #(dataframe,x,y)+
  labs(title = ('DDX3Y Stratfication Plot_Batch'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#-------------------------------------#
#4. Age - Scatter plot using ggplot2
#-------------------------------------#

#1. AFF3

AFF3<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["AFF3",]),
                 Age = as.numeric(clean_metadata_filtered$Age))
head(AFF3)
ggplot(AFF3, aes(x=Age, y=Expression))+
  geom_point()+
  geom_smooth(method = "lm")+#lm=linear model
  theme_bw()+
  labs(title = ('AFF3 Scatter Plot_Age'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#Note: 

#head(AFF3$Age)
#AFF3$Age<-as.factor(AFF3$Age)
#head(AFF3)
#dim(AFF3)
#938  2
#plotStratify(Expression~Age, data=AFF3) + #(dataframe,x,y)+
 # labs(title = ('AFF3 Stratfication Plot_Age'))+
#  theme(plot.title = element_text(
#    hjust = 0.5, 
#    face = "bold" ))

#plotStratify Gives boxplot plot for all ages as 1 category

#So must use ggplot2

#2. IGFBP3

IGFBP3<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["IGFBP3",]),
                 Age = as.numeric(clean_metadata_filtered$Age))
head(IGFBP3)
ggplot(IGFBP3, aes(x=Age, y=Expression))+
  geom_point()+
  geom_smooth(method = "lm")+#lm=linear model
  theme_bw()+
  labs(title = ('IGFBP3 Scatter Plot_Age'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))


#----------------------------------------------------------------------------#
#5. Individual - Need to figure out Fig. 4f from VariancePartition Paper
#----------------------------------------------------------------------------#

#1. GSTM1

#Need to figure out the scatter plot for Individual like in ImmVar paper

head(clean_metadata_filtered$Individual)
GSTM1<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["GSTM1",]),
                  Individual = as.factor(clean_metadata_filtered$Individual)) #categorical factor
head(GSTM1)

ggplot(GSTM1, aes(x=Individual, y=Expression))+
  geom_point()+
  geom_smooth(method = "lm")+#lm=linear model
  theme_bw()+
  labs(title = ('GSTM1 Scatter Plot_Age'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#Standard scatter plot does not work.


head(clean_metadata_filtered$Individual)
GSTM1<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["GSTM1",]),
                  Individual = as.factor(clean_metadata_filtered$Individual)) #categorical factor
head(GSTM1)

#To replicated Figure 4f from Variance Partition papaer
#GSTM1 Individual Stratification Plot

#Find median to get the smooth line and prevent points jumping up and down

GSTM1_median<-GSTM1
head(GSTM1_median)
type(GSTM1_median$Expression)
type(GSTM1_median$Individual)

head(GSTM1_median)

GSTM1_median<-GSTM1_median%>%
  group_by(Individual)%>% #Groups Expression of each individual together
  mutate(Median_Expression=median(Expression))%>% #Finds median in a new column Median_Exprsn
  ungroup()%>% #ungroup for future calulation
  mutate(Individual_Sorted=reorder(Individual,Median_Expression)) #creates new column with sorted values of Individual based on median exprsn

head(GSTM1_median)
#Looks as ordered in Asccending order

ggplot(GSTM1_median,aes(x=Individual_Sorted,y=Expression))+ #y=raw expression values
  stat_summary(fun=median,geom="line")+
  geom_point(alpha=0.5)+
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black",fill=NA))+
  labs(title = ('GSTM1 Individual Stratification Plot'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))


#2. ZFP57

head(clean_metadata_filtered$Individual)
ZFP57<-data.frame(Expression=as.numeric(expr_matrix_mapped_final["ZFP57",]),
                  Individual = as.factor(clean_metadata_filtered$Individual)) #categorical factor

#To replicated Figure 4f from Variance Partition papaer
#ZFP57 Individual Stratification Plot

#Find median to get the smooth line and prevent points jumping up and down

ZFP57_median<-ZFP57
head(ZFP57_median)
type(ZFP57_median$Expression)
type(ZFP57_median$Individual)

head(ZFP57_median)

ZFP57_median<-ZFP57_median%>%
  group_by(Individual)%>% #Groups Expression of each individual together
  mutate(Median_Expression=median(Expression))%>% #Finds median in a new column Median_Exprsn
  ungroup()%>% #ungroup for future calulation
  mutate(Individual_Sorted=reorder(Individual,Median_Expression)) #creates new column with sorted values of Individual based on median exprsn

head(ZFP57_median)
#Looks as ordered in Asccending order

ggplot(ZFP57_median,aes(x=Individual_Sorted,y=Expression))+ #y=raw expression values
  stat_summary(fun=median,geom="line")+
  geom_point(alpha=0.5)+
  theme(panel.background = element_blank(),panel.border = element_rect(color = "black",fill=NA))+
  labs(title = ('ZFP57 Individual Stratification Plot'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold" ))

#=================================================================================================================#
#=================================================================================================================#


#=================================================================================================================#
#Mapping individual variables to gene names _ Data not Quantile Normalized
#=================================================================================================================#

#=================================================================================================================#
#Note:

#I am not sure if the data has already been quantile normalized
#While the author notes says the data is quantile normalized,
#Boxplot showed quantile normalization may be necessary
#However both PCA and Variant Partition Analysis didnt show batch effects

#So I am checking how the data changes for Individual genes with non-quantile normalized data1
#But this is something I need to verify

#=================================================================================================================#

head(expr_matrix_log_quant_filter)
dim(expr_matrix_log_quant_filter)
#[1] 21727   938
apply(expr_matrix_log_quant_filter[,1:5],2,median)

head(expr_matrix_log_filter)  
dim(expr_matrix_log_filter)  #Not quantile normalized
apply(expr_matrix_log_filter[,1:5],2,median)

#[1] 21727   938
colnames(metadata)
head(metadata$title)
head(metadata$geo_accession)
head(metadata$data_processing)
#Affymetrix annotation
#The datasets was pre-filtered to keep only those probesets for which a gene
#symbol could be found in the Affymetrix annotation
head(metadata$type)

#RNA
head(colnames(expr_matrix_log_filter))
head(metadata$geo_accession)


#Match
#geo_accession=patient ID
#GEO=Gene Expression Omnibus
#GSM=GEO Sample

head(rownames(expr_matrix_log_filter))
# "7896740" "7896742" "7896744" "7896750" "7896754" "7896756"

#Micro array Probe_ID
#Platforms (1)	
#GPL6244	[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56035

#Affymetrix Human Gene 1.0-ST Array Transcriptcluster Revision 8 annotation data (chip hugene10sttranscriptcluster)
#https://www.bioconductor.org/packages//2.10/data/annotation/html/hugene10sttranscriptcluster.db.html

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

#--------------------#
#Extracting probe ID
#--------------------#

head(expr_matrix_log_filter)
dim(expr_matrix_log_filter)
#[1] 21727   938
dim(expr_matrix_log)
#[1] 21727   984

probe_id<-rownames(expr_matrix_log_filter)
head(probe_id)
length(probe_id)
#[1] 21727


#--------------------#
#Map to gene symbol
#--------------------#

columns(hugene10sttranscriptcluster.db)
#mapIds: https://support.bioconductor.org/p/9138389/
gene_symbol<-mapIds(hugene10sttranscriptcluster.db,
                    keys = probe_id,
                    column = "SYMBOL",
                    keytype = "PROBEID",
                    multiVals = "first") #gives 1st symbol for IDs mapping to multiple symbols
head(gene_symbol)                    
# 7896740      7896742      7896744      7896750      7896754      7896756 
#"OR4F4"     "PCMTD2"     "OR4F29"           NA "SEPTIN7P13"     "FAM87B"
length(gene_symbol)
#21727

#Need to filter missing values

Symbol_not_NA<-!(is.na(gene_symbol))
head(Symbol_not_NA)

gene_symbol2<-gene_symbol[Symbol_not_NA]
head(gene_symbol2)
length(gene_symbol2)
#[1] 21219
gene_symbol2<-as.data.frame(gene_symbol2)
head(gene_symbol2)
gene_symbol2$Probe_ID<-rownames(gene_symbol2)
head(gene_symbol2)
colnames(gene_symbol2)<-c("Gene_Symbol","Probe_ID")
head(gene_symbol2)
gene_symbol2 <- gene_symbol2 %>%
  relocate(Gene_Symbol, .after = Probe_ID)
head(gene_symbol2)
rownames(gene_symbol2)<-NULL
head(gene_symbol2)

table(duplicated(gene_symbol2$Probe_ID))
#FALSE 
#21219 

table(duplicated(gene_symbol2$Gene_Symbol))
#FALSE  TRUE 
#19837  1382 

#Thus there are duplicated Symbols


#-------------------------------#
#Replace probeID with Symbol
#-------------------------------#

head(expr_matrix_log_filter)
dim(expr_matrix_log_filter)
#[1] 21727   938
head(rownames(expr_matrix_log_filter))
expr_matrix_log_map<-expr_matrix_log_filter
head(expr_matrix_log_map)
dim(expr_matrix_log_map)
#[1] 21727   938
head(gene_symbol2) 
dim(gene_symbol2) 
dim(expr_matrix_log_map)
#[1] 21727   938

#Filter expr_matrix_log_map
expr_matrix_log_map<-expr_matrix_log_map[gene_symbol2$Probe_ID,] #only rows that are in gene_symbol2 kept
head(expr_matrix_log_map)
dim(expr_matrix_log_map)
#[1] 21219   938
dim(gene_symbol2)
#[1] 21219     2
#Both match
head(expr_matrix_log_map)
library(limma)
head(gene_symbol2)
#Average Duplicate gene symbol for microarray data
expr_matrix_log_mapped_final<-avereps(expr_matrix_log_map, ID=gene_symbol2$Gene_Symbol) 
head(expr_matrix_log_mapped_final)
dim(expr_matrix_log_mapped_final)
#[1] 19837   938
head(rownames(expr_matrix_log_mapped_final))

#---------------------------------------------------------------------#
#VariantPartition Analysis - Not Quantile Normalized and Mapped Data
#---------------------------------------------------------------------#

head(expr_matrix_log_mapped_final)
dim(expr_matrix_log_mapped_final)
#[1] 19837   938
head(clean_metadata_filtered)
dim(clean_metadata_filtered)
#[1] 938   6
dim(expr_matrix_log_mapped_final)
#[1] 19837   938
lmm
varpart_log_mapped<-fitExtractVarPartModel(expr_matrix_log_mapped_final,lmm,clean_metadata_filtered,BPPARAM=serial_param)

#Sort and Plot

varpartsorted_log_mapped<-sortCols(varpart_log_mapped)
plotVarPart(varpartsorted_log_mapped)

library(ggplot2)
plotVarPart(varpartsorted_log_mapped)+
  labs(title = ('VariantPartition Analysis - After Mapping to Gene Names, Before Quantile Normalization'))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) # This centers the title

#==============================================================================================================#
#==========================================================================================================#
#Variant Partition Analysis for Individual genes (prior to Quantile Normalization for cross verification)
#==========================================================================================================#

head(expr_matrix_log_mapped_final)
dim(expr_matrix_log_mapped_final)
head(rownames(expr_matrix_log_mapped_final))
Random6<-head(rownames(expr_matrix_log_mapped_final))
Random6

#---------------------------------------------------#
#Individual variant analysis for First 6 genes
#---------------------------------------------------#

plotPercentBars(varpart_log_mapped[Random6,])
plotPercentBars(varpart_log_mapped[Random6,])+
  labs(title = ('VariantPartition Analysis - For The First 6 Individual Genes, prior to quantile normalization'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold",
    margin = margin(t=10)))

#---------------------------------------------------#
#Top 2 genes from each metadata category
#---------------------------------------------------#

head(varpart_log_mapped)
varpart_log_mapped_df<-as.data.frame(varpart_log_mapped)
head(varpart_log_mapped_df)
type(varpart_log_mapped_df$CellType)
colnames(clean_metadata)
lmm

top2_all<-c(rownames(varpart_log_mapped_df[order(varpart_log_mapped_df$CellType, decreasing = TRUE),])[1:2],
            rownames(varpart_log_mapped_df[order(varpart_log_mapped_df$Individual, decreasing = TRUE),])[1:2],
            rownames(varpart_log_mapped_df[order(varpart_log_mapped_df$Age, decreasing = TRUE),])[1:2],
            rownames(varpart_log_mapped_df[order(varpart_log_mapped_df$Gender, decreasing = TRUE),])[1:2],
            rownames(varpart_log_mapped_df[order(varpart_log_mapped_df$Batch, decreasing = TRUE),])[1:2]
)  
top2_all

plotPercentBars(varpart_log_mapped[top2_all,])+
  labs(title = ('VariantPartition Analysis - For Top2 Genes from Each Metadata Group, prior to quantile Normalization'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold",
    margin = margin(t=10)))+
  guides(fill=guide_legend(title = "Metadata Groups"))

#---------------------------------------------------#
#Genes from paper
#---------------------------------------------------#

Paper<-c('AGPAT1','KARS','TLR4','DUSP6','GSTM1','ZFP57','RPS4Y1','UTY','SULT1E1') 
Paper
plotPercentBars(varpart_log_mapped[Paper,])+
  labs(title = ('VariantPartition Analysis - For Genes from VariancePartition Paper, prior to quantile Normalization'))+
  theme(plot.title = element_text(
    hjust = 0.5, 
    face = "bold",
    margin = margin(t=10)))+
  guides(fill=guide_legend(title = "Metadata Groups"))

#Batch effect completely removed (different from seen in papaer)

#=====================================================================================================#

#--------------------------------------------------------------------------------------------------------------------#
#PlotStratifyby - Can do if needed to analyze influence of metatadata groups 
#on variance of individual genes
#---------------------------------------------------------------------------------------------------------------------#

#=====================================================================================================#

#----------------------------------------------------------------------------------------#
#Future: Downstream Differential Expression Analysis using Dream form VariancePartition
#----------------------------------------------------------------------------------------#

#DEA - using Dream in Variance Partition published in 2021
#https://academic.oup.com/bioinformatics/article/37/2/192/5878955

#=====================================================================================================#

