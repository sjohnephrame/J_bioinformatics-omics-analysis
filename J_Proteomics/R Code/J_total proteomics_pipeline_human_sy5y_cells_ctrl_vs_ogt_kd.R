#===================================================================================
#J_18-Plex TMT-Total Proteomics analysis_SY5Y OGT Knockdown (KD) Neuroblatoma cells
#===================================================================================

# Attach libraries
library(BiocManager)
library(readxl)
library(tidyverse)
library(limma)
library(EnvStats)   # used to get geometric means
library(missForest) # Imputation
library(biomaRt) # library for mapping between annotations
library(qvalue)

# Load helper functions
path.to.functions = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\functions2.R"
source(path.to.functions) #functions and variables made available in the current environment.
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

#=========================#
# Proteomics Preprocessing ----
#=========================#

dat1 = read_xlsx("C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_Total Proteome_R_High and Medium FDR.xlsx")
colnames(dat1)
data2<-dat1
colnames(dat1)
#dat1<-dat1[,-(20:23)]
colnames(dat1)
names(dat1)[names(dat1) == 'Gene Symbol'] <- 'Symbol'
colnames(dat1)
data2new<-dat1
library(dplyr)
dat1 <- dat1 %>%
  relocate(Symbol, .after = Accession)
colnames(dat1)
data2new1<-dat1

dim(dat1)
#[1] 10150    20 before removing duplicates #24 with description

#add_symbols - adds symbols to refseq id and removes duplicates
#We already have symbols
#Remove duplicates
d <- duplicated(dat1$Symbol)
head(d)
dat1 <- dat1[!d,]
nrow(dat1)
dim(dat1)

#[1] 9206   20 after removing duplicates ##24 with description, removed duplicates

head(dat1) #Has symbols already
dat1new2<-dat1
head(dat1)
head(dat1new2)
#dat3<-dat1
#dat3<-as.data.frame(dat3)
#dat3<-as.matrix(dat3)
#rownames(dat3) <- dat3$Accession
new_dat1 <- dat1
head(new_dat1)
new_dat1<-as.data.frame(new_dat1)
head(new_dat1)
rownames(new_dat1) <- new_dat1$Accession #convert column value to row names
head(new_dat1)
new_dat1.1<-new_dat1
head(new_dat1.1)
new_dat1<-new_dat1[,-1]
head(new_dat1) #Without Accession because accession is now rowname
head(new_dat1.1) #with Accession
dim(new_dat1) #[1] 9206   23 #without accession
dim(new_dat1.1) #[1] 9206   24 #with accession

# This might look weird, but there will be NaN for some values that had no data
# and this will change them to NA, which is easier to handle.

new_dat1[is.na(new_dat1)] <- NA
dim(new_dat1) #[1] 9206   23 #without accession
new_dat1.1[is.na(new_dat1.1)] <- NA
dim(new_dat1.1) #[1] 9206   24 #with accession

#--------------------------------
# Assess missingness of proteins
#--------------------------------

head(new_dat1)
dim(new_dat1) #[1] 9206   23 #without accession
new_dat2<-new_dat1[,-1]
head(new_dat1) #with symbols
head(new_dat2) #without symbols
dim(new_dat1) #with symbols [1] 9206   19/23 with Description
head(new_dat1)
#new_dat1<-new_dat1[,-1]
head(new_dat2)
dim(new_dat2) #without symbols [1] 9206   18/22 with Description
head(new_dat1) #with symbols
new_dat1<-new_dat1[,-(21:23)]
head(new_dat1) #with symbols
dim(new_dat1)
new_dat1<-new_dat1[,-(20)]
head(new_dat1) #with symbols

table(apply(new_dat1, 1, function(x) sum(is.na(x)))) #with symbols
#   0    1   18 
#8334    1  871 

#apply(new_dat1, 1, function(x) sum(is.na(x))):

#new_dat1: Your data frame.
#1: Indicates that the function should be applied to rows. 2 for cols
#function(x) sum(is.na(x)): A function that calculates the number of NA values in each row.

head(new_dat2) #without symbols 
new_dat2<-new_dat2[,-(19:22)]
head(new_dat2) #without symbols 

table(apply(new_dat2, 1, function(x) sum(is.na(x))))#without symbols
#0   18 
#8335  871

dim(new_dat1) #[1] 9206   19 in new_dat1 #with symbols
head(new_dat1)
symbol<-new_dat1$Symbol
symbol<-as.data.frame(symbol)
table(apply(symbol, 1, function(x) sum(is.na(x))))#check symbols for missing values
#   0    1 
#9205    1 
#thus, there is one missing value with symbol
#so ideally good to map symbols here before proceeding further

#For now let us omit NA with symbols

#before omitting missing values 
#new_dat1 has 9206 obs.

# Imputation doesn't seem to be warranted, since a protein with missing values is missing across all samples.

#----------------------------------------------
# Reduce to proteins with no missing data
#----------------------------------------------
dim(new_dat1) #9206 rows 19 cols in new_dat1 #with symbols
head(new_dat1)
dat_comp1 = na.omit(new_dat1) #with symbols
dim(dat_comp1) #with symbols
#[1] 8334   19 #19 becauseof Symbol

#Matches: table(apply(new_dat1, 1, function(x) sum(is.na(x)))) #with symbols
#   0    1   18 
#8334    1  871 

head(new_dat1)
dim(new_dat1) #9206 rows 19 cols in new_dat1 #with symbols, without omit
dim(dat_comp1) #[1] 8334   19 #19 because of symbol, with omit
mapping1<-dat_comp1 
head(mapping1)
dim(mapping1)  #[1] 8334   19 #19 because of symbol
#mapping1<-mapping1[,-(2:19)]

#without symbols
head(new_dat2)
dim(new_dat2)
dat_comp2 = na.omit(new_dat2)
dim(dat_comp2)
#[1] 8335   18

#-----------------------------------------------------------------------------------------
# log2 transform - to reduce skewness of a measurement variable 
#to make data more symmetrical, which helps it meet the assumptions of statistical models
#-----------------------------------------------------------------------------------------

dim(dat_comp1) #[1] 8334   19 before removing symbol
head(dat_comp1) 
dat_comp1.1<-dat_comp1
dat_comp1<-dat_comp1[,-1]
dim(dat_comp1) #[1] 8334   18 after removing symbol # without NAs
head(dat_comp1)
head(new_dat1)
dim(new_dat1) #[1] 9206   19 with NA, with symbol
new_dat1.1<-new_dat1
new_dat1<-new_dat1[,-1]
head(new_dat1)
dim(new_dat1) #[1] 9206   18 after remving symbols # with NAs

head(dat_comp1)
dim(dat_comp1) #[1] 8334   18 after removing symbol # without NAs
dat_comp_mat1 = as.matrix(log2(dat_comp1)) # without NAs
head(dat_comp_mat1)
head(new_dat1)
dim(new_dat1)# with NAs
dat_mat1 = as.matrix(log2(new_dat1)) # with NAs

head(dat_comp_mat1)
head(dat_mat1)
dim(dat_comp_mat1) #[1] 8334   18 without NAs_log_matrix
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix

dim(dat_comp1) #[1] 8334   18 without NAs_ not log
head(dat_comp1)

#density plot - plot is applied to log10 transformed data in code 
#dat_comp1 = na.omit(new_dat1) #with symbols
head(dat_comp1)#not log transformed, not normalized, omitted NA
dim(dat_comp1)#[1] 8334   18 without NAs_ not log, not normalized
for (i in 1:ncol(dat_comp1)){
  if (i==1){
    plot(density(log10(dat_comp1[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp1[,i]))
    lines(den$x,den$y,col="blue")}
}

#-------------------------------------------------------------------------------------------------
# Quantile normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------

head(dat_comp_mat1)
head(dat_mat1)
dim(dat_comp_mat1) #[1] 8334   18 without NAs_log_matrix
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
dat_comp_mat_quant1 = limma::normalizeQuantiles(dat_comp_mat1) # without NAs, log
dat_mat_quant1 = limma::normalizeQuantiles(dat_mat1) # with NAs
head(dat_mat_quant1) #[1] 9206  18 # with NAs, log, quantile norm
dim(dat_mat_quant1) #[1] 9206   18 # with NAs, log, quantile norm
head(dat_comp_mat_quant1) #[1] 8334   18 # without NAs, log, quantile norm
dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quantile norm

dim(dat_comp1) # [1] 8334   18 without NAs, not log,
head(dat_comp1)
dat_comp1_norm = limma::normalizeQuantiles(dat_comp1) # without NAs, not log, normalized
dim(dat_comp1_norm) ## [1] 8334   18 without NAs, not log, quantile normalized
head(dat_comp1_norm)

#https://www.youtube.com/watch?v=ecjN6Xpv6SE&t=279s

#-------------------------------------------------------------------------------------------------
# Median normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------
head(dat_comp_mat1)
head(dat_mat1)
dim(dat_comp_mat1) #[1] 8334   18 without NAs_log_matrix; not normalized
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
dat_comp_mat_median1 = limma::normalizeMedianValues(dat_comp_mat1) # without NAs, log
dat_mat_median1 = limma::normalizeMedianValues(dat_mat1) # with NAs
head(dat_mat_median1) #[1] 9206  18 # with NAs, log, median norm
dim(dat_mat_median1) #[1] 9206   18 # with NAs, log, median norm
head(dat_comp_mat_median1) #[1] 8334   18 # without NAs, log, median norm
dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm

##########################################################################################

#for erk data, do a box plot

## Batch Effect Correction ##
# Here is where one should run diagnostics for batch effects (i.e., systematic differences in measurements due to technical factors rather than biological signal), if the experimental design introduces batch effects

# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat_comp_mat_quant, main="Boxplot normalized Intensities")

col <- c("red","red","red","blue","blue","blue","darkgreen","darkgreen","darkgreen","green","green","green","purple","purple","purple","orange","orange","orange")
col
fill <- c("red","blue","darkgreen","green","purple","orange")
fill
#par(mar = c(3, 3, 3, 10), mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2, xpd = TRUE)
boxplot(dat_comp_mat1, main="Boxplot before normalization",col=col)
boxplot(dat_comp_mat_quant1, main="Boxplot Quantile normalized Intensities",col=col)
boxplot(dat_comp_mat_median1, main="Boxplot Median normalized Intensities",col=col)
#dim(dat_comp_mat1) #[1] 8334   18 without NAs_log_matrix
#dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
#legend("right", title="Treatment",legend=c("GFP-0","GFP-10","OGTKD-605-0","OGTKD-605-10","OGTKD-606-0","OGTKD-606-10"), fill=fill, inset = c(-0.155,0))

#inset depends on margin of plot window

#############################################################################################

################################################################################################################
#density plot - after normalization - not needed just for reference

head(dat_comp1_norm)# without NAs, not log, normalized
for (i in 1:ncol(dat_comp1_norm)){
  if (i==1){
    plot(density(log10(dat_comp1_norm[,1])),main="Density plot across study samples",xlab="Subjects",col="red",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp1_norm[,i]))
    lines(den$x,den$y,col="blue")}
}


#Normalization makes density plot a single line

#============================================================================================

# Imputation (if warranted) #In data science, imputation is the process of replacing 
#missing or unavailable values in a dataset with substituted values. 

any(is.na(dat_mat_quant))      # Checks for NAs or NaNs
#True
any(is.na(dat_comp_mat_quant))      # Checks for NAs or NaNs
#False

#IMPUTE = FALSE
#if(IMPUTE){
# dat_mat_quant_t = t(dat_mat_quant)
  # missForest takes rows as samples, columns as variables
# dat_mat_quant_imp_t = missForest::missForest(dat_mat_quant_t)
  # If a protein has missing values across all samples, imputation won't work.
# dat_mat_quant_imp = t(dat_mat_quant_imp_t$ximp)
#}

any(is.na(dat_comp_mat_quant1))      # Checks for NAs or NaNs
any(is.infinite(dat_comp_mat_quant1)) # Checks for Inf (often caused by log2(0))
any(is.na(dat_comp_mat_median1))      # Checks for NAs or NaNs
any(is.infinite(dat_comp_mat_median1)) # Checks for Inf (often caused by log2(0))
any(is.na(dat_mat_median1))      # Checks for NAs or NaNs
any(is.infinite(dat_mat_median1)) # Checks for Inf (often caused by log2(0))
any(is.na(dat_mat_quant1))      # Checks for NAs or NaNs
#True
any(is.infinite(dat_mat_quant1)) # Checks for Inf (often caused by log2(0))
dim(dat_mat_quant1) #[1] 9206   18 with NAs_log_matrix_quantile norm
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
#So no need imputataion for all (False) but dat_mat_quant1 (True)

# Imputation (if warranted) - can skip As we are going to use the NA omitted data dat_comp_mat_quant1

#IMPUTE = FALSE
#if(IMPUTE){
#  dat_mat_quant_t1 = t(dat_mat_quant1)
  # missForest takes rows as samples, columns as variables
# dat_mat_quant_imp_t1 = missForest::missForest(dat_mat_quant_t1)
  # If a protein has missing values across all samples, imputation won't work.
# dat_mat_quant_imp1 = t(dat_mat_quant_imp_t1$ximp)
#}

#dim(dat_mat_quant_imp)

#===============================================================================================

#==========================================================================================================#
# Differential Expression Analysis 
# using limma-trend - # This assumes uniform and normal distribution of data (seen in density plot)
#==========================================================================================================#

#When the library sizes are quite variable between samples, then the voom approach 
#is theoretically more powerful than limma-trend. 

#Library size could mean one of two things: the total number of reads that were sequenced 
#in the run or the total number of mapped reads

#Read = sequenced fragment of cDNA (obtained from RNA) here ptn

head(dat_mat_quant1) #[1] 9206  18 # with NAs, log, quantile norm
dim(dat_mat_quant1) #[1] 9206   18 # with NAs, log, quantile norm
head(dat_comp_mat_quant1) #[1] 8334   18 # without NAs, log, quantile norm
dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quantile norm

datold1<-dat1
dat1 = dat_comp_mat_quant1 #[1] 8334   18  # without NAs, log, quantile norm
head(dat1)
dim(dat1) #[1] 8334   18 #removed NA, without Symbol, norm, log

# Change column names to not contain "-", since this will conflict with naming syntax in contrast matrix
colnames(dat1)
colnames(dat1) = gsub(pattern = "-", replacement = "", x = colnames(dat1))
colnames(dat1)
colnames(dat1) <- c('GFP01','GFP02','GFP03','GFP101','GFP102','GFP103','T60501','T60502','T60503','T605101','T605102','T605103','T60601','T60602','T60603','T606101','T606102','T606103')
colnames(dat1)
head(dat1)

# Getting a data.frame with 'Accession' and 'Symbol' columns
#mapping = add_symbols(data.frame(Accession = rownames(dat)))
#head(mapping)

head(mapping1)
head(dat1)
dim(dat1) #[1] 8334   18 Without symbol
dim(mapping1) #[1] 8334   19 With Symbol
#both dat1 and mapping1 match
#so, try to separate symbol and accession for mapping1

head(mapping)
mapping1.1<-mapping1
head(mapping1.1)
dim(mapping1.1) #[1] 8334   19
head(mapping1) #[1] 8334   19
dim(mapping1)
#make dim(mapping1) 8334 2
mapping2<-as.data.frame(cbind(rownames(mapping1),mapping1$Symbol))
dim(mapping2)
head(mapping2)
head(mapping)
colnames(mapping2) <- c('Accession','Symbol')
colnames(mapping2)
head(mapping2)
head(mapping)
dim(mapping2) #[1] 8334    2

# Design matrix

head(dat1)
colnames(dat1)<-c('GFPCTRL_1','GFPCTRL_2','GFPCTRL_3','GFPTEN_1','GFPTEN_2','GFPTEN_3','OGTCTRL_1','OGTCTRL_2','OGTCTRL_3','OGTTEN_1','OGTTEN_2','OGTTEN_3','ROGTCTRL_1','ROGTCTRL_2','ROGTCTRL_3','ROGTTEN_1','ROGTTEN_2','ROGTTEN_3')
colnames(dat1)
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
contrastsa = c("GFPTEN-GFPCTRL", "OGTTEN-OGTCTRL", "ROGTTEN-ROGTCTRL")
contrastsa
Xa
cm_dataa = makeContrasts(contrasts=contrastsa, levels=colnames(Xa))
cm_dataa #contrast matrix

#Try to make another contrast matrix for OGT VS GFP
contrastsb = c("GFPTEN-GFPCTRL", "OGTTEN-OGTCTRL", "ROGTTEN-ROGTCTRL","OGTCTRL-GFPCTRL","OGTTEN-GFPTEN","ROGTCTRL-GFPCTRL","ROGTTEN-GFPTEN")
contrastsb
Xa
cm_data = makeContrasts(contrasts=contrastsb, levels=colnames(Xa))
cm_data #contrast matrix

#----------------------------------------------------------------------------------------------------------------
#Function to do DE

##' Builds a limma model and gets the topTable results, with gene symbols added
#' 
#' @param dat (data.frame): a data.frame with an Accession and a Symbol column
#' @param dat_mat (matrix): a data.frame with an genes as rows
#' @param design (matrix): a design matrix for limma
#' @param cnt_mat (matrix): a contrast matrix for limma
#' @param annotate (logical): flag indicating if data should be annotated
#' @param trend (logical): flag indicating if mean-variance trend should be accounted for in eBayes (limma-trend method)
#' 
#' @return list with the topTable results and the model
#'

#top_table_results <- function(dat, dat_mat, design, cnt_mat, annotate = F, robust = F, trend = F) {
  ###
# if(0){ # For debugging
#   dat = mapping; dat_mat = dat[[1]]; design = DM[[1]]; cnt_mat = cnt_mat[[1]]; annotate = FALSE; trend = TRUE
#  }
  ###
#  require(limma)
  
# mod <- get_model(dat_mat, design, cnt_mat, robust = robust, trend = trend)
# tt <- topTable(mod, number = 100000, adjust = 'BH') # topTable - Extract a table of the top-ranked genes from a linear model fit.
#  tt = tt[,!colnames(tt) %in% c("F", "P.Value", "adj.P.Val")]
  
#  cnt_mat_names = colnames(cnt_mat)
#  p.vals = vector(mode = "list", length = length(cnt_mat_names))
#  for(i in 1:ncol(cnt_mat)){
#    tt.temp = topTable(mod, coef = cnt_mat_names[i], number = 100000, adjust = 'BH')
#    p.vals[[i]] = tt.temp[match(row.names(tt), row.names(tt.temp)), c("P.Value", "adj.P.Val")]
#    colnames(p.vals[[i]]) = str_replace(paste(colnames(p.vals[[i]]), cnt_mat_names[i], sep = "."), "-", ".")
#  }
#  tt = cbind(tt, do.call(cbind, p.vals))
  
#  if(annotate) {
#    tt <- add_annotation(tt, dat)	
#  } else {
#    tt$Symbol = rownames(tt)
#    tt$Accession = rownames(tt)
# }
#  return(list(tt = tt, mod = mod))
#} # end top_table_results

#-------------------------------------------------------------------------------------------------------------------

# Fitting models. The trend=TRUE argument indicates that the mean-variance trend 
#will be accounted for in eBayes (limma-trend method)
de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)

#In proteomics, the limma-trend method is often preferred because it accounts for 
#the fact that low-intensity proteins usually have more "noise" than high-intensity ones.

################################

head(mapping2)
head(dat1)
Xa
cm_data
dea = top_table_results(mapping2, dat1, Xa, cm_data, annotate = TRUE, trend = TRUE)
head(dea)
head(dea$tt) #24 columns #log FC 
head(dea$tt$Symbol[grepl("\\.", dea$tt$Symbol)])
dim(dea$tt) #[1] 8334   24
table(grepl("\\.", dea$tt$Symbol))
#grepl(): Short for "grep logical." It returns a vector of TRUE or FALSE values 
#rather than indices, which is perfect for filtering data frames.
#FALSE  TRUE 
# 8333     1
#1 gene symbols have "."; eg: "Uchl5"   "Uchl5.2
which(grepl("\\.", dea$tt$Symbol))
#[1] 910
dea$tt$Symbol[910]
dea$tt$Symbol[911]
#[1] "ap1s2; Ap1s2; AP1S2; ap1s2.L; LOC101341861; LOC111157096"
dea$tt[910, ]
#GFPTEN.GFPCTRL OGTTEN.OGTCTRL ROGTTEN.ROGTCTRL OGTCTRL.GFPCTRL OGTTEN.GFPTEN ROGTCTRL.GFPCTRL ROGTTEN.GFPTEN
#910      0.1035608         1.0138      -0.02977226      0.02104583     0.9312853         1.722384       1.589051
#AveExpr P.Value.GFPTEN.GFPCTRL adj.P.Val.GFPTEN.GFPCTRL P.Value.OGTTEN.OGTCTRL adj.P.Val.OGTTEN.OGTCTRL
#910 9.462007              0.7116239                0.8991319            0.002167602               0.02122121
#P.Value.ROGTTEN.ROGTCTRL adj.P.Val.ROGTTEN.ROGTCTRL P.Value.OGTCTRL.GFPCTRL adj.P.Val.OGTCTRL.GFPCTRL
#910                 0.915183                          1               0.9399809                 0.9727805
#P.Value.OGTTEN.GFPTEN adj.P.Val.OGTTEN.GFPTEN P.Value.ROGTCTRL.GFPCTRL adj.P.Val.ROGTCTRL.GFPCTRL
#910           0.004020953              0.01541427             1.460305e-05               0.0001330074
#P.Value.ROGTTEN.GFPTEN adj.P.Val.ROGTTEN.GFPTEN Accession
#910           3.518695e-05             0.0002457465    F6SFB5
#Symbol
#910 ap1s2; Ap1s2; AP1S2; ap1s2.L; LOC101341861; LOC111157096

#dea$tt$Symbol[910] <- "Ap1s2"

grep("^ap1s2", dea$tt$Symbol, value = TRUE, ignore.case = TRUE)
#The ^ ensures it starts with Acly (so you don't accidentally get other genes containing those letters).
#ignore.case = TRUE handles "ACLY", "Acly", or "acly"
dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]
#There are 2 rows with apls2; They are Case sensitive

#gene.symbola = extract_string(dea$tt$Symbol, "\\.", 1) # This doesn't work for case sensitive rows
#head(gene.symbola)
#length(gene.symbola)
head(dea$tt$Symbol)
length(dea$tt$Symbol)
#gene.symbola[grep("ap1s2", gene.symbola, ignore.case = TRUE)]
#grep("^ap1s2", gene.symbola, value = TRUE, ignore.case = TRUE)
#2 rows "ap1s2; Ap1s2; AP1S2; ap1s2" and "AP1S2"  
# 1. Take ONLY the first symbol before any semicolon or dot
clean_symbols = extract_string(dea$tt$Symbol, "[;.]", 1)
# 1. Clean the symbols - Take ONLY the first symbol before any semicolon or dot
#AND force them all to Uppercase
# This makes "ap1s2" and "AP1S2" both become "AP1S2"
clean_symbols_up = toupper(extract_string(dea$tt$Symbol, "[;.]", 1))

# 2. Check the count now—it should be LOWER than 8334
length(unique(clean_symbols))
length(unique(clean_symbols_up))
length(clean_symbols)
length(clean_symbols_up)
grep("^ap1s2", clean_symbols, value = TRUE, ignore.case = TRUE)
grep("^ap1s2", clean_symbols_up, value = TRUE, ignore.case = TRUE)
# both match
colnames(dea$tt)
tmpa = aggregate(dea$tt[,-c(23,24)], by = list(Symbol = clean_symbols_up), FUN = median)
#removing old symbols, accession
#adding the new symbols (gene.symbola) with say ABC.1 , ABC.2 extracted to ABC, ABC
#Finding median of ABC
#Be Careful using median as it can lose P value significance But in this case we had only one row with dot
dea$tt$Symbol[910]
tmpa$Symbol[910]
#Your new tmpa is now sorted alphabetically by Symbol (A-Z).
head(dea$tt)
length(dea$tt$Symbol)
head(tmpa)
length(tmpa$Symbol)
length(dea$tt$Symbol)

tmpa$Accession = dea$tt$Accession[match(tmpa$Symbol, clean_symbols_up)] 

#WE have two tables that are out of order:
#  tmpa: Sorted alphabetically.
#dea$tt: Sorted by P-value (Significant genes first).
#When you run: match(tmpa$Symbol, gene.symbola)
#R is saying:
#  "Hey, tmpa$Symbol says the first gene is A2M. 
#Go look through the original gene.symbola list and tell me which row A2M was originally in.
#If A2M was originally in row 500, the match function returns the number 500.
#Then, tmpa$Accession = dea$tt$Accession[500]
#"Go to the original table dea$tt$Accession, grab the Accession ID from row 500, 
#and paste it into the first row of my new table tmpa"
#match(A,B) takes a list of things you have (Table A) and finds 
#where they are located in a master list (Table B).

head(tmpa)
head(dea$tt)
nrow(tmpa) #[1] 8330
dim(dea$tt) #[1] 8334   24
dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]
tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
dea$tt = tmpa
head(dea$tt)
head(tmpa)

##################################################################################################################
######################################################################################################################
#logFCprot
head(dea$tt)
logFCprot<-dea$tt
head(logFCprot)
logFCprot<-as.data.frame(logFCprot)
head(logFCprot)
type(logFCprot)
library("writexl") 
#install.packages("writexl")
#write_xlsx(logFCprot,"C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\logFCprot.xlsx")
dim(logFCprot) #[1] 8330   24

#================================================================================================
#==============#
#MDS Plot
#==============#
head(dea)
head(dea$tt) #values after comparisons with contrast matrix
head(dat1) #normailzed and filtered values - use this for mds plot
#dat1 = dat_comp_mat_quant1 #[1] 8334   18  # without NAs, is log, quantile norm
#boxplot(dat_comp_mat1, main="Boxplot before normalization",col=col)
#boxplot(dat_comp_mat_quant1, main="Boxplot Quantile normalized Intensities",col=col)
#boxplot(dat_comp_mat_median1, main="Boxplot Median normalized Intensities",col=col)
dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
head(dat_comp_mat1)
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm
dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quntile norm

dim(dat1)
labels<-colnames(dat1)
labels
library(ggplot2)
col <- c("red","red","red","blue","blue","blue","darkgreen","darkgreen","darkgreen","green","green","green","purple","purple","purple","orange","orange","orange")
col
#https://support.bioconductor.org/p/100859/
#plotMDS(dat1, main="Multipledimensional scaling plot (MDS) \n After quantile normalization",col=col)
plotMDS(dat_comp_mat1, main="Multipledimensional scaling plot (MDS) \n Before normalization",col=col)
plotMDS(dat_comp_mat_quant1, main="Multipledimensional scaling plot (MDS) \n After quantile normalization",col=col)
plotMDS(dat_comp_mat_median1, main="Multipledimensional scaling plot (MDS) \n After median normalization",col=col)

#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html

######################################################################################

#Get PCA

get_pca1 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA) \n After Quantile normalization" )+
    scale_color_manual(values = c("red","blue","darkgreen","green","purple","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

head(dat1) 
#dat1 = dat_comp_mat_quant1 #[1] 8334   18  # without NAs, is log, quantile norm
labels<-c("GFP-0","GFP-0","GFP-0","GFP-10","GFP-10","GFP-10","OGTKD-605-0","OGTKD-605-0","OGTKD-605-0","OGTKD-605-10","OGTKD-605-10","OGTKD-605-10","OGTKD-606-0","OGTKD-606-0","OGTKD-606-0","OGTKD-606-10","OGTKD-606-10","OGTKD-606-10")
labels
#head(dat1)
PCA<-get_pca1(dat1,labels)
PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/

#===========================#
#PCA - Before Normalization
#===========================#
dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
head(dat_comp_mat1)

get_pca1 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA) \n Before normalization" )+
    scale_color_manual(values = c("red","blue","darkgreen","green","purple","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca


labels<-c("GFP-0","GFP-0","GFP-0","GFP-10","GFP-10","GFP-10","OGTKD-605-0","OGTKD-605-0","OGTKD-605-0","OGTKD-605-10","OGTKD-605-10","OGTKD-605-10","OGTKD-606-0","OGTKD-606-0","OGTKD-606-0","OGTKD-606-10","OGTKD-606-10","OGTKD-606-10")
labels
#head(dat1)
PCA<-get_pca1(dat_comp_mat1,labels)
PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/

#===================================#
#PCA - After Quantile Normalization
#===================================#

get_pca1 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA) \n After Quantile normalization" )+
    scale_color_manual(values = c("red","blue","darkgreen","green","purple","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

head(dat1)  
dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
head(dat_comp_mat1)
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm
dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quntile norm

labels<-c("GFP-0","GFP-0","GFP-0","GFP-10","GFP-10","GFP-10","OGTKD-605-0","OGTKD-605-0","OGTKD-605-0","OGTKD-605-10","OGTKD-605-10","OGTKD-605-10","OGTKD-606-0","OGTKD-606-0","OGTKD-606-0","OGTKD-606-10","OGTKD-606-10","OGTKD-606-10")
labels
#head(dat1)
PCA<-get_pca1(dat_comp_mat_quant1,labels)
PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/


#===================================#
#PCA - After Median Normalization
#===================================#

get_pca1 <- function(dat, labels, legend_title = 'Treatment') {
  require(ggplot2)
  
  pca_res <- prcomp(t(dat), center = T, scale. = F)
  
  U <- data.frame(pca_res$x)
  
  p <- ggplot(data = U, aes(x = PC1, y = PC2, color = labels )) +
    geom_point(size = 3, alpha = 0.5) + 
    theme_bw() +
    labs(color = legend_title, title="Principal Component Analysis (PCA) \n After Median normalization" )+
    scale_color_manual(values = c("red","blue","darkgreen","green","purple","orange"))+
    theme(plot.title = element_text(hjust=0.5))
  
  return(list(p = p , summ = summary(pca_res)))
} # end get_pca

head(dat1)  
labels<-c("GFP-0","GFP-0","GFP-0","GFP-10","GFP-10","GFP-10","OGTKD-605-0","OGTKD-605-0","OGTKD-605-0","OGTKD-605-10","OGTKD-605-10","OGTKD-605-10","OGTKD-606-0","OGTKD-606-0","OGTKD-606-0","OGTKD-606-10","OGTKD-606-10","OGTKD-606-10")
labels
#head(dat1)
PCA<-get_pca1(dat_comp_mat_median1,labels)
PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/


##########################################################################################

#for erk data, do a box plot
# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat_comp_mat_quant, main="Boxplot normalized Intensities")

col <- c("red","red","red","blue","blue","blue","darkgreen","darkgreen","darkgreen","green","green","green","purple","purple","purple","orange","orange","orange")
col
fill <- c("red","blue","darkgreen","green","purple","orange")
fill
#par(mar = c(3, 3, 3, 10), mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2, xpd = TRUE)
boxplot(dat_comp_mat1, main="Boxplot before normalization",col=col)
boxplot(dat_comp_mat_quant1, main="Boxplot Quantile normalized Intensities",col=col)
boxplot(dat_comp_mat_median1, main="Boxplot Median normalized Intensities",col=col)

#legend("right", title="Treatment",legend=c("GFP-0","GFP-10","OGTKD-605-0","OGTKD-605-10","OGTKD-606-0","OGTKD-606-10"), fill=fill, inset = c(-0.155,0))

#inset depends on margin of plot window

#############################################################################################

### Density plots by sample
###-----------------------------------------------

head(dat1)
dataa<-dat1 #quantile normalized
head(dataa)
head(dat)
for (i in 1:ncol(dat)){
  if (i==1){
    plot(density(log10(dat[,1])),main="Density plot across study samples",xlab="Subjects",col="red",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat[,i]))
    lines(den$x,den$y,col="blue")}
}

for (i in 1:ncol(dataa)){
  if (i==1){
    plot(density(dataa[,1]),main="Density plot across study samples \n quantile normalized",xlab="Subjects",col="red")}
  else {
    den <- density(dataa[,i])
    lines(den$x,den$y,col="blue")}
}

head(dat_comp1_norm)#not log transformed, normalized
for (i in 1:ncol(dat_comp1_norm)){
  if (i==1){
    plot(density(log10(dat_comp1_norm[,1])),main="Density plot across study samples",xlab="Subjects",col="red",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp1_norm[,i]))
    lines(den$x,den$y,col="blue")}
}

head(dat_comp1)#not log transformed, not normalized, omitted NA
dim(dat_comp1)
for (i in 1:ncol(dat_comp1)){
  if (i==1){
    plot(density(log10(dat_comp1[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp1[,i]))
    lines(den$x,den$y,col="blue")}
}

head(dat1)  
dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
head(dat_comp_mat1)
dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm
dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quntile norm

for (i in 1:ncol(dat_comp_mat1)){
  if (i==1){
    plot(density((dat_comp_mat1[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((dat_comp_mat1[,i]))
    lines(den$x,den$y,col="blue")}
}

for (i in 1:ncol(dat_comp_mat_quant1)){
  if (i==1){
    plot(density((dat_comp_mat_quant1[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((dat_comp_mat_quant1[,i]))
    lines(den$x,den$y,col="blue")}
}

for (i in 1:ncol(dat_comp_mat_median1)){
  if (i==1){
    plot(density((dat_comp_mat_median1[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((dat_comp_mat_median1[,i]))
    lines(den$x,den$y,col="blue")}
}

#####################################################################################################

#==================================================#
# Volcano plots for adjusted p-values
#==================================================#

#Final Code
#==================================================#
colnames(dea$tt) #logFC and pvalues

#####################################################################################################################################################################

#GFPTEN.GFPCTRL

rx <- c(-1, 1)*max(abs(dea$tt$GFPTEN.GFPCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.GFPTEN.GFPCTRL), -log10(dea$tt$adj.P.Val.GFPTEN.GFPCTRL))))

#https://biostatsquid.com/volcano-plot/

#adjusted p value
plot(dea$tt$GFPTEN.GFPCTRL, -log10(dea$tt$adj.P.Val.GFPTEN.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - GFP 10 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$GFPTEN.GFPCTRL)> lfc & -log(dea$tt$adj.P.Val.GFPTEN.GFPCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$GFPTEN.GFPCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.GFPTEN.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.GFPTEN.GFPCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.GFPTEN.GFPCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#############################################################################################################################################
#OGTTEN.OGTCTRL

rx <- c(-1, 1)*max(abs(dea$tt$OGTTEN.OGTCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTTEN.OGTCTRL), -log10(dea$tt$adj.P.Val.OGTTEN.OGTCTRL))))

#adjusted p value
plot(dea$tt$OGTTEN.OGTCTRL, -log10(dea$tt$adj.P.Val.OGTTEN.OGTCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KD 605-10 min vs. OGT KD 605-0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$OGTTEN.OGTCTRL)> lfc & -log(dea$tt$adj.P.Val.OGTTEN.OGTCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$OGTTEN.OGTCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.OGTTEN.OGTCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$OGTTEN.OGTCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTTEN.OGTCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$OGTTEN.OGTCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTTEN.OGTCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)


#############################################################################################################################################
colnames(dea$tt)
#ROGTTEN.ROGTCTRL

rx <- c(-1, 1)*max(abs(dea$tt$ROGTTEN.ROGTCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.ROGTTEN.ROGTCTRL), -log10(dea$tt$adj.P.Val.ROGTTEN.ROGTCTRL))))

#adjusted p value
plot(dea$tt$ROGTTEN.ROGTCTRL, -log10(dea$tt$adj.P.Val.ROGTTEN.ROGTCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KD 606-10 min vs. OGT KD 606-0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$ROGTTEN.ROGTCTRL)> lfc & -log(dea$tt$adj.P.Val.ROGTTEN.ROGTCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$ROGTTEN.ROGTCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.ROGTTEN.ROGTCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$ROGTTEN.ROGTCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.ROGTTEN.ROGTCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$ROGTTEN.ROGTCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.ROGTTEN.ROGTCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#############################################################################################################################################
colnames(dea$tt)
#OGTCTRL.GFPCTRL

rx <- c(-1, 1)*max(abs(dea$tt$OGTCTRL.GFPCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTCTRL.GFPCTRL), -log10(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL))))

#adjusted p value
plot(dea$tt$OGTCTRL.GFPCTRL, -log10(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KD 605-0 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$OGTCTRL.GFPCTRL)> lfc & -log(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$OGTCTRL.GFPCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$OGTCTRL.GFPCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTCTRL.GFPCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$OGTCTRL.GFPCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTCTRL.GFPCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#############################################################################################################################################
colnames(dea$tt)
#OGTTEN.GFPTEN

rx <- c(-1, 1)*max(abs(dea$tt$OGTTEN.GFPTEN))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTTEN.GFPTEN), -log10(dea$tt$adj.P.Val.OGTTEN.GFPTEN))))

#adjusted p value
plot(dea$tt$OGTTEN.GFPTEN, -log10(dea$tt$adj.P.Val.OGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KD 605-10 min vs. GFP 10 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$OGTTEN.GFPTEN)> lfc & -log(dea$tt$adj.P.Val.OGTTEN.GFPTEN,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$OGTTEN.GFPTEN)< (-lfc) & -log(dea$tt$adj.P.Val.OGTTEN.GFPTEN,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$OGTTEN.GFPTEN,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTTEN.GFPTEN,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$OGTTEN.GFPTEN,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTTEN.GFPTEN,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#############################################################################################################################################
colnames(dea$tt)
head(dea$tt)
#ROGTCTRL.GFPCTRL

#If your input data was already log2 transformed, 
#then the "Fold Change" you calculated is actually the log2 Fold Change

rx <- c(-1, 1)*max(abs(dea$tt$ROGTCTRL.GFPCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.ROGTCTRL.GFPCTRL), -log10(dea$tt$adj.P.Val.ROGTCTRL.GFPCTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")


#p-value
#plot(dea$tt$ROGTCTRL.GFPCTRL, -log10(dea$tt$P.Value.ROGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
 #    xlim=rx, ylim=ry, xaxt="n",
  #   xlab="Fold change", ylab="-log10  p-value")
#abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
#abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
#axis(1, seq(-6,6,1))
#title("Volcano plot of p-values - OGT KD 606-0 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
#lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
#pval <- 0.05 
#-log10(pval)
#-log10(0.01)
#>1.3 is good
# Selecting interesting genes
#sigGenes1 <- ((dea$tt$ROGTCTRL.GFPCTRL)> lfc & -log(dea$tt$P.Value.ROGTCTRL.GFPCTRL,10) > -log10(pval))   
#sigGenes2 <- ((dea$tt$ROGTCTRL.GFPCTRL)< (-lfc) & -log(dea$tt$P.Value.ROGTCTRL.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
#points(dea$tt[sigGenes1,]$ROGTCTRL.GFPCTRL,-log(dea$tt[sigGenes1,]$P.Value.ROGTCTRL.GFPCTRL,10),pch=20,col="orange",cex=2)
#points(dea$tt[sigGenes2,]$ROGTCTRL.GFPCTRL,-log(dea$tt[sigGenes2,]$P.Value.ROGTCTRL.GFPCTRL,10),pch=20,col="blue",cex=2)
#abline(h=-log10(pval),col="brown4",lty=2)
#abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
#mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)

#adjusted p value
plot(dea$tt$ROGTCTRL.GFPCTRL, -log10(dea$tt$adj.P.Val.ROGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values - OGT KD 606-0 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$ROGTCTRL.GFPCTRL)> lfc & -log(dea$tt$adj.P.Val.ROGTCTRL.GFPCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$ROGTCTRL.GFPCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.ROGTCTRL.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$ROGTCTRL.GFPCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.ROGTCTRL.GFPCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$ROGTCTRL.GFPCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.ROGTCTRL.GFPCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
# Keep lfc at 0.58 for the position, but use "1.5" for the text
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#The antilog of 0.58 (in base 2) is 1.5.

#############################################################################################################################################
colnames(dea$tt)
#ROGTTEN.GFPTEN

rx <- c(-1, 1)*max(abs(dea$tt$ROGTTEN.GFPTEN))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.ROGTTEN.GFPTEN), -log10(dea$tt$adj.P.Val.ROGTTEN.GFPTEN))))

#adjusted p value
plot(dea$tt$ROGTTEN.GFPTEN, -log10(dea$tt$adj.P.Val.ROGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of adjusted p-values -  OGT KD 606-10 min vs. GFP 10 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
# Selecting interesting genes
sigGenes1 <- ((dea$tt$ROGTTEN.GFPTEN)> lfc & -log(dea$tt$adj.P.Val.ROGTTEN.GFPTEN,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$ROGTTEN.GFPTEN)< (-lfc) & -log(dea$tt$adj.P.Val.ROGTTEN.GFPTEN,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$ROGTTEN.GFPTEN,-log(dea$tt[sigGenes1,]$adj.P.Val.ROGTTEN.GFPTEN,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$ROGTTEN.GFPTEN,-log(dea$tt[sigGenes2,]$adj.P.Val.ROGTTEN.GFPTEN,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#############################################################################################################################################

#PRep for GSEA

#OGT KD 605
dea
dea2<-dea
#dea2 <- estimateDisp(dea2)-didnt work	  ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)
#dat.filtered$samples$lib.size <- colSums(dat.filtered$counts)
#dat.filtered
### PLOTS 
#summary(de <- decideTestsDGE(dea2$tt, p=0.05, adjust="BH"))
dea$tt
head(dea2)
dea2$tt
updownOGTTENGFPTEN<-cbind(dea2$tt$OGTTEN.GFPTEN,dea2$tt$adj.P.Val.OGTTEN.GFPTEN,dea2$tt$Symbol)
head(updownOGTTENGFPTEN)
head(dea2$tt)
updownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN)
colnames(updownOGTTENGFPTEN)<-c('logFC.OGTTEN.GFPTEN','adj.P.Val.OGTTEN.GFPTEN','Symbol')
head(updownOGTTENGFPTEN)
rownames(updownOGTTENGFPTEN)<-updownOGTTENGFPTEN$Symbol
head(updownOGTTENGFPTEN)
updownOGTTENGFPTEN<-updownOGTTENGFPTEN[,-3]
head(updownOGTTENGFPTEN)
updownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN)
head(updownOGTTENGFPTEN)
#convert e to nmeric so can sort
updownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN <- as.numeric(formatC(updownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN, format = "e", digits = 2))
head(updownOGTTENGFPTEN)
sortupdownOGTTENGFPTEN<-as.data.frame(updownOGTTENGFPTEN[order(updownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN,decreasing = TRUE),]) 
head(sortupdownOGTTENGFPTEN)

#OGTKD 606 FOR GSEA
dea
dea2<-dea
#dea2 <- estimateDisp(dea2)-didnt work	  ## Calculates genewise dispersion parameter adjusted using bayesian empirical method(Recommended, edgeR)
#dat.filtered$samples$lib.size <- colSums(dat.filtered$counts)
#dat.filtered
### PLOTS 
#summary(de <- decideTestsDGE(dea2$tt, p=0.05, adjust="BH"))
dea$tt
head(dea2)
dea2$tt
updownROGTTENGFPTEN<-cbind(dea2$tt$ROGTTEN.GFPTEN,dea2$tt$adj.P.Val.ROGTTEN.GFPTEN,dea2$tt$Symbol)
head(updownROGTTENGFPTEN)
head(dea2$tt)
head(updownROGTTENGFPTEN)
updownROGTTENGFPTEN<-as.data.frame(updownROGTTENGFPTEN)
head(updownROGTTENGFPTEN)
colnames(updownROGTTENGFPTEN)<-c('logFC.ROGTTEN.GFPTEN','adj.P.Val.ROGTTEN.GFPTEN','Symbol')
head(updownROGTTENGFPTEN)
rownames(updownROGTTENGFPTEN)<-updownROGTTENGFPTEN$Symbol
head(updownROGTTENGFPTEN)
updownROGTTENGFPTEN<-updownROGTTENGFPTEN[,-3]
head(updownROGTTENGFPTEN)
updownROGTTENGFPTEN<-as.data.frame(updownROGTTENGFPTEN)
head(updownROGTTENGFPTEN)
#convert e to nmeric so can sort
updownROGTTENGFPTEN$adj.P.Val.ROGTTEN.GFPTEN <- as.numeric(formatC(updownROGTTENGFPTEN$adj.P.Val.ROGTTEN.GFPTEN, format = "e", digits = 2))
head(updownROGTTENGFPTEN)
type(updownROGTTENGFPTEN)
updownROGTTENGFPTEN$logFC.ROGTTEN.GFPTEN <- as.numeric(formatC(updownROGTTENGFPTEN$logFC.ROGTTEN.GFPTEN, format = "e", digits = 2))
head(updownROGTTENGFPTEN)
type(updownROGTTENGFPTEN)
sortupdownROGTTENGFPTEN<-as.data.frame(updownROGTTENGFPTEN[order(updownROGTTENGFPTEN$logFC.ROGTTEN.GFPTEN,decreasing = TRUE),]) 
head(sortupdownROGTTENGFPTEN)

###################################################################################################################
###################################################################################################################

#NOT NEEDED FOR GSEA BUT CAN DO FOR KEGG AND GO
#remove genes not significant p>0.05 
#OGTKD 605 10 MIN VS GFP 10 MIN
remsortupdownOGTTENGFPTEN <- sortupdownOGTTENGFPTEN[!(sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN>0.05), ]
head(remsortupdownOGTTENGFPTEN)
remsortupdownOGTTENGFPTEN<-as.data.frame(remsortupdownOGTTENGFPTEN)
head(remsortupdownOGTTENGFPTEN)
dim(remsortupdownOGTTENGFPTEN)
#3348 proteins < p=0.05 and FDR Filtration 

#OGTKD 606 10 MIN VS GFP 10 MIN
remsortupdownROGTTENGFPTEN <- sortupdownROGTTENGFPTEN[!(sortupdownROGTTENGFPTEN$adj.P.Val.ROGTTEN.GFPTEN>0.05), ]
head(remsortupdownROGTTENGFPTEN)
remsortupdownROGTTENGFPTEN<-as.data.frame(remsortupdownROGTTENGFPTEN)
head(remsortupdownROGTTENGFPTEN)
dim(remsortupdownROGTTENGFPTEN)
#[1] 4652    2 < p=0.05 and FDR Filtration 


# Filter log FC >1.5
head(remsortupdownOGTTENGFPTEN)
logFCgreat<-remsortupdownOGTTENGFPTEN
head(logFCgreat)
logFCgreat<-as.data.frame(logFCgreat)
#logFCgreat<-logFCgreat[,-2]
head(logFCgreat)
logFCgreat<-subset(logFCgreat,(logFCgreat$logFC.OGTTEN.GFPTEN)>1.5)
head(logFCgreat)
dim(logFCgreat)

logFCless<-remsortupdownOGTTENGFPTEN
head(logFCless)
logFCless<-as.data.frame(logFCless)
head(logFCless)
type(logFCless)
#logFCless<-as.numeric(logFCless$logFC.OGTTEN.GFPTEN)
#head(logFCless)
logFCless <- transform(logFCless, logFC.OGTTEN.GFPTEN = as.numeric(logFC.OGTTEN.GFPTEN))
head(logFCless)
type(logFCless)
#logFCless<-subset(logFCless,(((logFCless$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCless$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCless<-subset(logFCless,(((logFCless$logFC.OGTTEN.GFPTEN)<(-1.5))))

logFCgreat<-remsortupdownOGTTENGFPTEN
head(logFCgreat)
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-as.numeric(logFCgreat$logFC.OGTTEN.GFPTEN)
#head(logFCgreat)
logFCgreat <- transform(logFCgreat, logFC.OGTTEN.GFPTEN = as.numeric(logFC.OGTTEN.GFPTEN))
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))


head(logFCgreat)
dim(logFCgreat)

logFCgreat<-remsortupdownOGTTENGFPTEN
head(logFCgreat)
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-as.numeric(logFCgreat$logFC.OGTTEN.GFPTEN)
#head(logFCgreat)
logFCgreat <- transform(logFCgreat, logFC.OGTTEN.GFPTEN = as.numeric(logFC.OGTTEN.GFPTEN))
head(logFCgreat)
type(logFCgreat)
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
#logFCgreat<-subset(logFCgreat,(((logFCgreat$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCgreat<-as.data.frame(logFCgreat)
head(logFCgreat)
type(logFCgreat)
logFCgreat<-subset(logFCgreat,logFCgreat$logFC.OGTTEN.GFPTEN>1.5 | logFCgreat$logFC.OGTTEN.GFPTEN<(-1.5))
head(logFCgreat)
#logFCgreat <- logFCgreat[logFCgreat$logFC.OGTTEN.GFPTEN > 1.5 | logFCgreat$logFC.OGTTEN.GFPTEN < -1.5, ]
head(logFCgreat)
dim(logFCgreat)


#p-value
plot(logFCgreat$logFC.OGTTEN.GFPTEN, -log10(logFCgreat$adj.P.Val.OGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10 p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Volcano plot of p-values - GFP 10 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 1.5
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((logFCgreat$logFC.OGTTEN.GFPTEN)> lfc & -log(logFCgreat$adj.P.Val.OGTTEN.GFPTEN,10) > -log10(pval))   
sigGenes2 <- ((logFCgreat$logFC.OGTTEN.GFPTEN)< (-lfc) & -log(logFCgreat$adj.P.Val.OGTTEN.GFPTEN,10) > -log10(pval))   
# Identifying the selected genes
points(logFCgreat[sigGenes1,]$logFC.OGTTEN.GFPTEN,-log(logFCgreat[sigGenes1,]$adj.P.Val.OGTTEN.GFPTEN,10),pch=20,col="orange",cex=2)
points(logFCgreat[sigGenes2,]$logFC.OGTTEN.GFPTEN,-log(logFCgreat[sigGenes2,]$adj.P.Val.OGTTEN.GFPTEN,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)

###################################################################################################################
###################################################################################################################

#Prepare protlist for waterfall plot

#GSEA 606 PROT LIST
head(sortupdownROGTTENGFPTEN)
prot.list <- sortupdownROGTTENGFPTEN$logFC.ROGTTEN.GFPTEN             # rank-ordered gene list
head(prot.list)
prot.list<-as.numeric(prot.list)
head(prot.list)
head(sortupdownROGTTENGFPTEN)
names(prot.list) <- rownames(sortupdownROGTTENGFPTEN)
head(prot.list)

#GSEA 605 PROT LIST
head(sortupdownOGTTENGFPTEN)
prot.list2 <- sortupdownOGTTENGFPTEN$logFC.OGTTEN.GFPTEN             # rank-ordered gene list
head(prot.list2)
prot.list2<-as.numeric(prot.list2)
head(prot.list2)
head(sortupdownOGTTENGFPTEN)
names(prot.list2) <- rownames(sortupdownOGTTENGFPTEN)
head(prot.list2)

# Waterfall plot

#OGT KD 606

# Barplot of ranked fold changes (waterfall plot)
#barplot(sort(gene.list, decreasing = T),axisnames=FALSE,main="Plot of ranked gene Fold changes")

head(prot.list)
col <- ifelse(prot.list >0, "orange", "blue")
barplot(sort(prot.list, decreasing = T),axisnames=FALSE,main="Plot of ranked Fold changes of total proteome - OGT KD - 606 10 min vs. GFP - 10 min",  col=col, border=col, ylim=c(-4,4))

barplot(sort(prot.list2, decreasing = T),axisnames=FALSE,main="Plot of ranked Fold changes of total proteome - OGT KD - 605 10 min vs. GFP - 10 min", col=col, border=col, ylim=c(-4,4))


#barplot(lib.size,xaxt="n",xlab="Study samples",ylab="Library size values",main="Barplot of library size",
 #       col=c(rep("red",3),rep("blue",3),rep("darkgreen",3),rep("green",3),rep("purple",3),rep("orange",3)),sub="(Red line represents mean)")
#abline(h=mean(lib.size),lwd=2,col="red")
#legend("bottomright", title="Treatment",legend=c("GFP-0","GFP-10","OGTKD-605-0","OGTKD-605-10","OGTKD-606-0","OGTKD-606-10"), fill=fill)


###############################################################################################################################################

###GO analysis
###====================================###

#goana: Gene Ontology or KEGG Pathway Analysis
#In limma: Linear Models for Microarray Data

#topGO package provides tools for testing GO terms while accounting for the topology of the GO graph. Different test statistics and different methods for eliminating local similarities and dependencies between GO terms can be implemented and applied.

#In R, topKEGG is a function that extracts the most significant KEGG pathways from kegga output. 
#kegga function to gather pathway enrichment for my dataset.

library(org.Hs.eg.db)
#library(edgeR)

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("GO.db")

library(GO.db)

d.go <- d.tmp
head(d.go)
#d.go.DE <- subset(d.go,FDR<0.05) 
d.go.DE <- subset(d.go,PValue<0.05) #subsetted based on p value
head(d.go.DE)
d.entrez.id <- mapIds(org.Hs.eg.db, keys=d.go.DE$hgnc_symbol,column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id)
head(d.entrez.id)
all(d.go.DE$hgnc_symbol==names(d.entrez.id)) 
go.test <- goana(d.entrez.id,species="Hs")
go.results <- topGO(go.test, sort = "DE", number = Inf)
head(go.results)
sum(go.results$P.DE<10^(-5))

#OGT KD 605 10 MIN VS GFP 10 MIN

#d.go1 <- dea$tt
#head(d.go1)
#d.go.DE <- subset(d.go,FDR<0.05) 
head(remsortupdownOGTTENGFPTEN)
d.go1 <- remsortupdownOGTTENGFPTEN
head(d.go1) #subsetted already for p<0.05
#d.go.DE <- subset(d.go,PValue<0.05) 
d.go.DE1<-d.go1
head(d.go.DE1)
d.entrez.id1 <- mapIds(org.Hs.eg.db, keys=rownames(d.go.DE1),column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id1)
head(d.entrez.id1)
all(rownames(d.go.DE1)==names(d.entrez.id1)) 
go.test1 <- goana(d.entrez.id1,species="Hs")
go.results1 <- topGO(go.test1, sort = "DE", number = Inf)
head(go.results1, 20)
sum(go.results1$P.DE<10^(-5)) #487
sum(go.results1$P.DE<0.05) #2532
#In the goana output from the limma package, the columns represent the following:
#Term: The full descriptive name of the Gene Ontology (GO) term (e.g., "cell cycle" or "metabolic process").
#Ont: The specific GO ontology or "aspect" the term belongs to:
  #BP: Biological Process (e.g., DNA replication).
  #CC: Cellular Component (e.g., mitochondrion).
  #MF: Molecular Function (e.g., enzyme activity).
#N: The total number of genes in the entire background "universe" 
    #that are annotated to that specific GO term.
#DE: The number of genes from your Differentially Expressed (DE) list 
    #that are annotated to that GO term.
#P.DE: The p-value for the over-representation of that GO term in your 
#DE gene set, typically calculated using a hypergeometric test 
#(equivalent to Fisher's exact test). 
#www.rdocumentation.org
#www.rdocumentation.org
 

#topBP <- topGO(go.test1, ontology = "BP", number = 20)
#print(topBP)

library(ggplot2)

# Get top 20 pathways
df_plot <- topGO(go.test1, sort = "DE", number = 20)

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
  labs(title = "Top 20 Enriched GO Terms - Before Filtering \n OGT KD 605 10 min vs GFP - 10 min",
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

ggplot(plot_data, aes(x = logP, y = Term)) +
  geom_point(aes(size = DE, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Gene Ontology - Top 20 Specific Biological Processes",
    subtitle = "OGT KD 605 10 min vs GFP 10 min \n Filtered for N < 1000 genes",
    x = "-log10(P-value)",
    y = NULL,
    size = "Genes in List",
    color = "Significance"
  ) +
  theme(
    # Bolds the pathway names (Y-axis)
    axis.text.y = element_text(size = 10, face = "bold", color = "black"), 
    # Bolds the X-axis numbers (-log10 P-values)
    axis.text.x = element_text(face = "bold"),
    # Centers the titles
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold")
  )


#if (!require(ggtext)) install.packages("ggtext")
library(ggtext)

#Update the Plot Code for mixed title

ggplot(plot_data, aes(x = logP, y = Term)) +
  geom_point(aes(size = DE, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Gene Ontology (GO) - Top 20 Specific Biological Processes",
    # Use HTML tags: <b> for bold, <br> for new line
    subtitle = "<b>OGT KD 605 10 min vs GFP 10 min</b><br>Filtered for N < 1000 genes",
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
    subtitle = "<b>OGT KD 605 10 min vs GFP 10 min</b><br>Filtered for N < 500 genes",
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
install.packages("webshot2")
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
  gtsave("C:\\Users\\sophi\\OneDrive\\Desktop\\J_Total proteome SY5Y plots\\J_QUANTILE NORMALIZATION\\GO_TABLE1_OGT KD 605 10 MIN VS GFP 10MIN_N LESS 1000.png") # This saves it to your working directory

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
  gtsave("C:\\Users\\sophi\\OneDrive\\Desktop\\J_Total proteome SY5Y plots\\J_QUANTILE NORMALIZATION\\GO_TABLE1_OGT KD 605 10 MIN VS GFP 10MIN_N LESS 500.png") # This saves it to your working directory

#OGT KD 606 10 MIN VS GFP 10 MIN

#d.go1 <- dea$tt
#head(d.go1)
#d.go.DE <- subset(d.go,FDR<0.05) 
head(remsortupdownROGTTENGFPTEN)
d.go2 <- remsortupdownROGTTENGFPTEN
head(d.go2) #subsetted already for p<0.05
#d.go.DE <- subset(d.go,PValue<0.05) 
d.go.DE2<-d.go2
head(d.go.DE2)
head(d.go.DE1)
d.entrez.id2 <- mapIds(org.Hs.eg.db, keys=rownames(d.go.DE2),column="ENTREZID",keytype="SYMBOL")
length(d.entrez.id2)
head(d.entrez.id2)
all(rownames(d.go.DE2)==names(d.entrez.id2)) 
go.test2 <- goana(d.entrez.id2,species="Hs")
go.results2 <- topGO(go.test2, sort = "DE", number = Inf)
head(go.results2, 20)
sum(go.results2$P.DE<10^(-5)) #761
sum(go.results2$P.DE<0.05) #3410
#In the goana output from the limma package, the columns represent the following:
#Term: The full descriptive name of the Gene Ontology (GO) term (e.g., "cell cycle" or "metabolic process").
#Ont: The specific GO ontology or "aspect" the term belongs to:
#BP: Biological Process (e.g., DNA replication).
#CC: Cellular Component (e.g., mitochondrion).
#MF: Molecular Function (e.g., enzyme activity).
#N: The total number of genes in the entire background "universe" 
#that are annotated to that specific GO term.
#DE: The number of genes from your Differentially Expressed (DE) list 
#that are annotated to that GO term.
#P.DE: The p-value for the over-representation of that GO term in your 
#DE gene set, typically calculated using a hypergeometric test 
#(equivalent to Fisher's exact test). 
#www.rdocumentation.org
#www.rdocumentation.org


#topBP <- topGO(go.test1, ontology = "BP", number = 20)
#print(topBP)

library(ggplot2)

# Get top 20 pathways
df_plot2 <- topGO(go.test2, sort = "DE", number = 20)

# Add -log10 P-value for the x-axis
df_plot2$logP <- -log10(df_plot2$P.DE)
head(df_plot2)
# Reorder Term by significance so the plot is ranked
df_plot2$Term <- reorder(df_plot2$Term, df_plot2$logP)
head(df_plot2)

#Dot plot - Top 20 Enriched GO Terms - OGT KD 605 10 min vs GFP - 10 min
#Without filtering
ggplot(df_plot2, aes(x = logP, y = Term, size = DE, color = logP)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top 20 Enriched GO Terms - Before Filtering \n OGT KD 606 10 min vs GFP - 10 min",
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
go.bp2 <- topGO(go.test2, ontology = "BP", sort = "DE", number = 50)
head(go.bp2, 20)

#N: The total number of genes in the entire background "universe" 
#that are annotated to that specific GO term.
#DE: The number of genes from your Differentially Expressed (DE) list 
#that are annotated to that GO term.

#Let's grab terms that have fewer than 1,000 total genes (N < 1000). This removes 
#the "cellular process" noise and finds the specific pathways.

# Filter for Biological Process and smaller, more specific terms
go.specific3 <- go.results2[go.results2$Ont == "BP" & go.results2$N < 1000, ]

# Take the top 20 of THESE results
plot_data3 <- head(go.specific3, 20)
head(plot_data3)
# Create -log10 P-value for the plot
plot_data3$logP <- -log10(plot_data3$P.DE)
plot_data3$Term <- reorder(plot_data3$Term, plot_data3$logP)
head(plot_data3)
print(plot_data3)

library(ggplot2)

#if (!require(ggtext)) install.packages("ggtext")
library(ggtext)

#Update the Plot Code for mixed title

ggplot(plot_data3, aes(x = logP, y = Term)) +
  geom_point(aes(size = DE, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Gene Ontology (GO) - Top 20 Specific Biological Processes",
    # Use HTML tags: <b> for bold, <br> for new line
    subtitle = "<b>OGT KD 606 10 min vs GFP 10 min</b><br>Filtered for N < 1000 genes",
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
go.specific4 <- go.results2[go.results2$Ont == "BP" & go.results2$N < 500, ]

# Take the top 20 of THESE results
plot_data4 <- head(go.specific4, 20)

# Create -log10 P-value for the plot
plot_data4$logP <- -log10(plot_data4$P.DE)
plot_data4$Term <- reorder(plot_data4$Term, plot_data4$logP)

print(plot_data4)

library(ggtext)

#Update the Plot Code for mixed title

ggplot(plot_data4, aes(x = logP, y = Term)) +
  geom_point(aes(size = DE, color = logP)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(
    title = "Gene Ontology (GO) - Top 20 Specific Biological Processes",
    # Use HTML tags: <b> for bold, <br> for new line
    subtitle = "<b>OGT KD 606 10 min vs GFP 10 min</b><br>Filtered for N < 500 genes",
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

plot_data3 %>%
  # Target only the relevant columns
  dplyr::select(Term, Ont, N, DE, P.DE) %>% 
  gt() %>%
  # 1. Bold the Title and Subtitle
  tab_header(
    title = md("**Gene Ontology (GO) Top 20 Specific Biological Processes**"),
    # Use <br> for the line break and move the ** markers
    subtitle = md("**OGT KD 606 10 min vs GFP 10 min** <br> Filtered for N < 1000 genes")
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
  gtsave("C:\\Users\\sophi\\OneDrive\\Desktop\\J_Total proteome SY5Y plots\\J_QUANTILE NORMALIZATION\\GO_TABLE1_OGT KD 606 10 MIN VS GFP 10MIN_N LESS 1000.png") # This saves it to your working directory

plot_data4 %>%
  # Target only the relevant columns
  dplyr::select(Term, Ont, N, DE, P.DE) %>% 
  gt() %>%
  # 1. Bold the Title and Subtitle
  tab_header(
    title = md("**Gene Ontology (GO) Top 20 Specific Biological Processes**"),
    # Use <br> for the line break and move the ** markers
    subtitle = md("**OGT KD 606 10 min vs GFP 10 min** <br> Filtered for N < 500 genes")
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
  gtsave("C:\\Users\\sophi\\OneDrive\\Desktop\\J_Total proteome SY5Y plots\\J_QUANTILE NORMALIZATION\\GO_TABLE1_OGT KD 606 10 MIN VS GFP 10MIN_N LESS 500.png") # This saves it to your working directory

#=======================#            
###KEGG analysis
#=======================#
            
d.kegg1 <- remsortupdownOGTTENGFPTEN
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
#d.kegg.DE <- subset(d.kegg,PValue<0.05)
d.kegg.DE1<-d.kegg1
head(d.kegg.DE1)
head(d.entrez.id1)
all(rownames(d.kegg.DE1)==names(d.entrez.id1)) 

kegg.test1 <- kegga(d.entrez.id1,species="Hs")
kegg.results1 <- topKEGG(kegg.test1, sort = "DE", number = Inf)
head(kegg.results1)
sum(kegg.results1$P.DE<10^(-5))#13
sum(kegg.results1$P.DE<0.05) #67

d.kegg2 <- logFCgreat
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
#d.kegg.DE <- subset(d.kegg,PValue<0.05)
d.kegg.DE2<-d.kegg2
head(d.kegg.DE2)
head(d.entrez.id2)
all(rownames(d.kegg.DE2)==names(d.entrez.id2)) 

kegg.test2 <- kegga(d.entrez.id2,species="Hs")
kegg.results2 <- topKEGG(kegg.test2, sort = "DE", number = Inf)
head(kegg.results2)
sum(kegg.results2$P.DE<10^(-5)) #0
sum(kegg.results2$P.DE<0.05) #6

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("KEGGREST")
library(KEGGREST)
#kegg.test <- kegga(d.entrez.id,species="Mm")
#kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
#head(kegg.results)
#sum(kegg.results$P.DE<10^(-5))

######################################################################################################################
######################################################################################################################

### J_GSEA analysis for OGT 606 10 MIN VS GFP 10 MIN

# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
library(fgsea)

all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v2023.2.Hs.symbols.gmt") #got lot of warnings
all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
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
warnings()
head(fgseaRes)
head(fgseaRes[order(pval), ])
sum(fgseaRes[, padj < 0.05])#289
fgseaRes1<-fgseaRes[order(padj), ]
head(fgseaRes1)


# Make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
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
plotEnrichment(all.gene.sets[["GOMF_TRANSPORTER_ACTIVITY"]],prot.list) + labs(title="GOMF_TRANSPORTER_ACTIVITY")
#Plot2
plotEnrichment(all.gene.sets[["GOCC_ORGANELLE_INNER_MEMBRANE"]],prot.list) + labs(title="GOCC_ORGANELLE_INNER_MEMBRANE")
#Plot3
plotEnrichment(all.gene.sets[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],prot.list) + labs(title="KEGG_OXIDATIVE_PHOSPHORYLATION")
#Plot4
plotEnrichment(all.gene.sets[["WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]],prot.list) + labs(title="WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA")

#Top 4 down

#Plot1
plotEnrichment(all.gene.sets[["REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"]],prot.list) + labs(title="REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")
#gseaplot(all.gene.sets[["REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"]],prot.list) + labs(title="REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")

#Plot2
plotEnrichment(all.gene.sets[["REACTOME_SELENOAMINO_ACID_METABOLISM"]],prot.list) + labs(title="REACTOME_SELENOAMINO_ACID_METABOLISM")
#Plot3
plotEnrichment(all.gene.sets[["REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"]],prot.list) + labs(title="REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")
#Plot4
plotEnrichment(all.gene.sets[["REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY"]],prot.list) + labs(title="REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY")

#----------------------------------------------------------------------------------------------------

#Pretty gsea visualizations    #07/15/2024 

library(org.Hs.eg.db)
#install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
head(prot.list)

head(prot.list)
type(prot.list)


#GSEA 606 PROT LIST


library(enrichplot)

head(fgseaRes)
#listpathway<-fgseaRes[,c("pathway","padj","ES")]
#head(listpathway)
#topPathwaysUp1 <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
#listpathway <- listpathway[order(listpathway$padj),] #to sort based on padj
#head(listpathway)
#listpathway<-as.matrix(listpathway)

topPathwaysUp1 <- fgseaRes[ES > 0][head(order(padj), n=10)]
head(topPathwaysUp1)
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

topPathwaysDown1 <- fgseaRes[ES < 0][head(order(padj), n=10)]
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
ggplot(topPathways1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05&NES>0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()+
  theme(text = element_text(face = "bold"))

#Change legend





