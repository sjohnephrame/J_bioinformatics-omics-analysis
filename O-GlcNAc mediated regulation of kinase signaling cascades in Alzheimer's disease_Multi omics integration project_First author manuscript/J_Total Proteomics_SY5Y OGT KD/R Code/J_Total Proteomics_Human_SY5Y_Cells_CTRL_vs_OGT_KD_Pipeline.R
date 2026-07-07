#===================================================================================
#J_18-Plex TMT-Total Proteomics analysis_SY5Y OGT Knockdown (KD) Neuroblatoma cells
#===================================================================================

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# 2. Install the STRINGdb package
#BiocManager::install("STRINGdb")

# Clear environment
#rm(list = ls())

# Attach libraries
library(BiocManager)
library(readxl)
library(tidyverse)
library(limma)
library(EnvStats)   # used to get geometric means
library(missForest) # Imputation
library(biomaRt) # library for mapping between annotations
library(qvalue)
library(ggplot2)
# Load helper functions
path.to.functions = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\functions2.R"
source(path.to.functions)
#source(path/to/my_functions.R): In R, the source() function executes the script 
#at the specified path, making its functions and variables available in the 
#current environment.
extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))
#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
#strsplit(x, "\\."): Breaks a string like "ENSG000001.2" into two parts: [1] "ENSG000001" and [2] "2".

#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
#code explanation:


#In programming, a function is a named, self-contained block of code 
#designed to perform a specific task.


#extract_string = function(x, k, pos) unlist(lapply(strsplit(x, k), function(y) y[pos]))

#The overall result is a list containing the elements extracted from each string at the specified position. 

#Usage: DEA => txa = extract_string(colnames(dat1), "_", 1) 

#eg: GFP_2: O/P: GFP

#function(x, k, pos) = text (x)=colnames(dat1), a separator (k=_),
#and the position we want to grab the piecce from (pos=1).
#This sets up the "controls." You need to tell R what to split (x), where to cut it (k), and which piece to keep (pos).
#unlist= Make the list flat and easier to use = Otherwise it'll be like folders
#[1] 
#apple etc - Not easier to use
#lapply= list apply 
#l=The result is always a list.
#apply: It applies a specific function to each element Of the list.
#Format: lapply(X, FUNCTION, ...) 
#strsplit(x, k) = Split string X at separator k 
#"messy." It returns a List (folders)
#lapply(strsplit(x, k), function(y) y[pos]): Since strsplit created folders, 
#lapply reaches into each folder (y) and grabs the specific piece from the pos we want
#Here, strsplit(x="colnames(dat1)", k="_") = Split string X at separator k => GFP_2 will be: "GFP" "2"
#Here pos=1 so lapply(strsplit(x="colnames(dat1)", k="_"), function(y) y[pos=1]): O/P: GFP


#unlist(...): This converts the list of extracted elements back into a single vector.


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

#=========================#
# Proteomics Preprocessing ----
#=========================#
path.to.file = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_updated_11924\\J_R ANALYSIS\\J_SAM\\J_PROTEOMICS R CODE\\proteomic_example_data - Copy.xlsx"
sheet.name = NULL
sheet.range = "A5:AW5801"
# Proteomic data 
dat = read_xlsx(path = path.to.file, sheet = sheet.name, range = sheet.range)

# Column names has all of the sample information we need
colnames(dat)
data1<-dat
# e.g., F142: Fractionation number... a specific 'run' of the Mass Spec machine
# e.g., UA4427: Unique specimen ID
# OGA-KO, OGT-KO, TMG: treatments

head(dat)
colnames(dat) #Has no symbols
dat = add_symbols(dat) # This returns unique gene symbols by numbering any duplicates. Be careful!
colnames(dat) #Has Symbol

#=======================================================================================#

#-----------------------------------#
#add_symbols - Function explanation
#-----------------------------------#

# Duplicates are made unique by appending '.2', '.3', etc. 
# sum(grepl("\\.", dat$Symbol))
#' Add HUGO gene symbols to a data.frame with genes and RefSeq gene ids.
#' 
#' Assumes that data.frame passed has an Acc column with RefSeq gene ids. 
#' 
#' @author Jeffrey A. Thompson
#' 
#' @param df (data.frame): data.frame with genes as rows
#' @param acc_col_name (character): name of column with accession number
#' 
#' @return data.frame with new Symbol column
#'

#add_symbols <- function(df, acc_col_name = 'Accession') {
 # accs <- as.matrix(df[,acc_col_name])[,1]
  
  ## Many of the accessions will have a version, indicated by a dot and
  ## a number at the end. We'll strip these off.
  #dot_locs <- unlist(gregexpr(pattern = '\\.', accs))
  #df[,acc_col_name] <- ifelse(dot_locs == -1, accs, substr(accs, 1, dot_locs - 1))
  
  #accs <- as.matrix(df[,acc_col_name])[,1]
  
  ## translate RefSeq accessions to gene symbols
  #symbols <- unlist(translate(accs, org.Mm.egREFSEQ2EG, org.Mm.egSYMBOL))
  
  ## Make symbols unique by numbering any duplicates, be aware that this
  ## will break matching these symbols with anything else though. 
  #dup_symbols <- list()
  #for(i in which(duplicated(symbols))) {
  #  if(is.null(dup_symbols[[symbols[i]]])) {
  #    dup_symbols[[symbols[i]]] <- 1
  #    symbols[i] <- paste0(symbols[i], '.', dup_symbols[[symbols[i]]] + 1)
  #  } else {
  #    dup_symbols[[symbols[i]]] = dup_symbols[[symbols[i]]] + 1
  #    symbols[i] <- paste0(symbols[i], '.', dup_symbols[[symbols[i]]] + 1)
  #  }
  #}
  
  ## Create data.frame with symbols and accession and then join back to
  ## original data.frame.
  #symbols <- data.frame(Acc = names(symbols), Symbol = symbols)
  #colnames(symbols)[1] <- acc_col_name
  #df <- symbols %>% left_join(df)
  
  #return(df)
#}

#Only works for refseq ID
#=======================================================================================#

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
#Already have symbols
#remove duplicates
d <- duplicated(dat1$Symbol)
head(d)
sum(d)
#dup<-dat1[d, ]

d2 <- duplicated(dat1$Accession)
head(d2)
sum(d2)
#So no duplicated Accession. this means duplicated symbols can be isoform etc..
#dup2<-dat1[d, ]

#dat1[grep("SPTAN1", dat1$Symbol), ]
#Accession  Symbol `GFP0-1` `GFP0-2` `GFP0-3` `GFP10-1` `GFP10-2` `GFP10-3`
#<chr>      <chr>     <dbl>    <dbl>    <dbl>     <dbl>     <dbl>     <dbl>
#  1 Q13813     SPTAN1   71631.   69746.   71508.    74220.    71882.    71402.
#  2 A0A0D9SF54 SPTAN1   70934.   69253.   70904.    73485.    71514.    71046.

#dat1[grep("DST", dat1$Symbol), ]

#Accession  Symbol `GFP0-1` `GFP0-2` `GFP0-3` `GFP10-1` `GFP10-2` `GFP10-3`
#<chr>      <chr>     <dbl>    <dbl>    <dbl>     <dbl>     <dbl>     <dbl>
 # 1 Q03001     DST      21948.   22882.   22319.    22669.    21692.    21490.
 # 2 F6QMI7     DST      20261.   21897.   21345.    21351.    20854     20625.
 # 3 A0A0U1RQJ2 DST       9024.   10140.    9660.     9437.     9394.     9113.

#Looks like we do not need median summarization of duplicates as accession = no dups
#and symbol dups can be isoforms!

#Same Site, Multiple Rows (Phospho): SUM intensities -> Log2 -> DEA.
#Same Symbol, Multiple Accessions (Total): DEA -> MEDIAN LogFCs

#In total proteomics, if the same Accession (e.g., O43719) appears twice, 
#it is usually because the software found it in two different "Protein Groups."
#Action: Median summarize the intensities.
#Why not Sum? In total proteomics, summing can lead to "double-counting" peptides 
#that are shared between the two rows. 
#The median is the safest way to find the "true" abundance of that protein without 
#artificially inflating the numbers.

#===========================================================================================#
#Aggressive Duplicate filtering by just removing all duplicates - not needed as aggresssive 
#===========================================================================================#

table(duplicated(dat1$Symbol))
#FALSE  TRUE 
#9206   944 
table(duplicated(dat1$Accession))
#FALSE 
#10150 

#No need to remove duplicates for now

#dat1.2 <- dat1[!d,]
#nrow(dat1.2)
#dim(dat1.2)

#[1] 9206   20 after removing duplicates ##24 with description, removed duplicates

#But this is very aggressiveBecause keeps only the 1st row and removes the rest of the rows

#head(dat1) #Has symbols already
#head(dat)

#======================================================================================================#
#So Summarize duplicates by Median intensity per site (most common)- not needed as accession no dups
#======================================================================================================#

table(duplicated(dat1$Symbol))
#FALSE  TRUE 
#9206   944 
table(duplicated(dat1$Accession))
#FALSE 
#10150 

#Grabbing DST protein symbol to check

dat1 %>% 
  dplyr::filter(Symbol == "DST") %>% 
  dplyr::select(Symbol, Accession, where(is.numeric))

# So Just isoforms; Leave as is

#dat1 <- dat1 %>% 
# group_by(Symbol) %>%
#  summarize(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) 
dat1<-as.data.frame(dat1)
dim(dat1) #[1] 10150    24 #After Summarizing duplicates by Median intensity per site 
head(dat1)

#pdat16.2 <- pdat15_filtered_check %>% 
 # group_by(SiteID) %>%
  #summarize(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) 
#pdat16.2<-as.data.frame(pdat16.2)
#dim(pdat16.2) #[1] 12298    19 #After Summarizing duplicates by Median intensity per site 
#head(pdat16.2)

#========================#
#Duplicate Sites Check
#========================#
#dat1[grep("DST", dat1$Symbol), ]
#dat1[grep("SPTAN1", dat1$Symbol), ]
#dat1.2[grep("DST", dat1.2$Symbol), ]
#dat1.2[grep("SPTAN1", dat1.2$Symbol), ]


# Unique ID for each biological replicate

colnames(dat)
colnames(dat) = gsub(pattern = "Abundances \\(Normalized\\): ", replacement = "", x = colnames(dat)) #remove Abundances Normalized from column name
colnames(dat)
UA.id = unlist(gregexpr(pattern = "UA", text = colnames(dat)[-c(1,2)]))
# Biological replicate labels
brep_labs = substr(colnames(dat)[-c(1,2)], UA.id, UA.id + 5)
# Each biological replicate has 3 technical replicates
table(brep_labs) #check why he has technical replicates... may need for phospho-proteomics

#We have 3 bioligcal replicates
#GFP0 - 3; GFP10 - 3, etc...
#No technical replicates 

colnames(dat)
# Treatment IDs
tmp = colnames(dat)[-c(1,2)]
tmp
tx_labs = unlist(lapply(strsplit(tmp, ", "), function(x) x[4]))
tx_labs #For getting wt vs ogt/oga

# Denote technical replicates 
frac_labs = substr(colnames(dat)[-c(1,2)], 1, 4)
frac_labs
sample_labs = paste(brep_labs, tx_labs, sep = "_")
sample_labs
colnames(dat) <- c(colnames(dat)[1:2], sample_labs)
colnames(dat)
table(sample_labs)

#Sofar just renaming 
head(dat)

# Reduce the technical replicates with geometric means. This also puts Accession IDs as rownames
new_dat = get_geo_means(dat, sample_labs)
colnames(new_dat)
colnames(dat)

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#' Get geometric means.
#' 
#' Note that if there are missing values for all technical replicates but one
#' then just that one will be returned.
#' 
#' @author Jeffrey A. Thompson
#' 
#' @param dat (data.frame): data.frame with genes as rows
#' @param sample_ids (character): which sample ids to find the mean for
#' @param cols (integer): which columns have the numeric data
#' @param label_col (character): name of column with gene ids
#' @param rowname_length (integer): max length of a gene name, used to make unique row names
#' 
#' @return data.frame with geometric means
#'
#get_geo_means <- function(dat, sample_ids, cols = 3:ncol(dat), label_col = 'Accession', rowname_length = 15) {
#  unique_ids <- unique(sample_ids)
  
#  new_dat <- matrix(0L, nrow = nrow(dat), ncol = length(unique_ids))
  
#  for(i in 1:nrow(dat)) {
#    for(j in 1:length(unique_ids)) {
#      new_dat[i,j] <- geoMean(as.numeric(dat[i,cols][which(sample_ids == unique_ids[j])]), na.rm=T)
#    }
#    #message(i, appendLF = T)
#  }
  
#  new_dat <- data.frame(new_dat)
#  colnames(new_dat) <- unique_ids
#  rownames(new_dat) <- make.unique(substr(as.character(as.matrix((dat[,label_col]))), 1, rowname_length))
  
#  return(new_dat)
#}

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------

#for loop
fruits <- list("apple", "banana", "cherry")

for (x in fruits) {
  print(x)
}
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------


#This section of the code is only needed when you have technical replicates.
#will need for phospho-proteomics

#but new dat now looks like matrix with no symbols

head(new_dat)
dim(new_dat) #[1] 5376   16
head(dat1)
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
dim(new_dat1) #[1] 10150   23 #without accession
dim(new_dat1.1) #[1] 10150   24 #with accession

# This might look weird, but there will be NaN for some values that had no data
# and this will change them to NA, which is easier to handle.
dim(new_dat) #[1] 5376   16
new_dat[is.na(new_dat)] <- NA
new_dat1[is.na(new_dat1)] <- NA
dim(new_dat) #[1] 5376   16
dim(new_dat1) #[1] 9206   23 #without accession
new_dat1.1[is.na(new_dat1.1)] <- NA
dim(new_dat1.1) #[1] 10150   24 #with accession

#--------------------------------
# Assess missingness of proteins
#--------------------------------

table(apply(new_dat, 1, function(x) sum(is.na(x))))
# 0   16 
#4881  544 

head(new_dat1)
dim(new_dat1) #[1] 10150   23 #without accession
new_dat2<-new_dat1[,-1] #Removing Symbols
head(new_dat1) #with symbols
head(new_dat2) #without symbols
dim(new_dat1) #with symbols [1] 10150   19/23 with Description
head(new_dat1)
#new_dat1<-new_dat1[,-1]
head(new_dat2)
dim(new_dat2) #without symbols [1] 10150   18/22 with Description
head(new_dat1) #with symbols
new_dat1<-new_dat1[,-(21:23)]
head(new_dat1) #with symbols
dim(new_dat1)
new_dat1<-new_dat1[,-(20)]
head(new_dat1) #with symbols

table(apply(new_dat1, 1, function(x) sum(is.na(x)))) #with symbols
#0    1   18   19 
#8923  276  914   37 


#Old
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
#9199  951

dim(new_dat1) #[1] 10150   19 in new_dat1 #with symbols
head(new_dat1)
symbol<-new_dat1$Symbol
symbol<-as.data.frame(symbol)
table(apply(symbol, 1, function(x) sum(is.na(x))))#check symbols for missing values
#   0    1 
#9837   313  
#thus, there are 313 missing values with symbol
#so ideally good to map symbols before proceeding further with DEA
head(dat1)
table(is.na(dat1$Accession))
table(duplicated(dat1$Accession))
#Both False


# Imputation doesn't seem to be warranted, since a protein with missing values is missing across all samples.

#----------------------------------------------
# Reduce to proteins with no missing data
#----------------------------------------------
dat_comp = na.omit(new_dat)
dim(dat_comp)
#[1] 4881   16

dim(new_dat1) #10150 rows 19 cols in new_dat1 #with symbols
head(new_dat1)
#dat_comp1 = na.omit(new_dat1) #with symbols - Too aggressive
##need to remove these all 0 rows before proceeding
type(new_dat1$`GFP0-1`)
#double
new_dat1 %>%
  summarise(across(where(is.double), ~ sum(is.na(.))))
#see blank rows across all numeric columns - Total = 951

new_dat1 %>%
  filter(if_any(where(is.double), is.na))
#Gives the rows with numeric values = NA (951)

new_dat1 <- new_dat1 %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))
#If the value is NA, give me a 0; otherwise, give me the value itself (.)
#Converted all NA rows to all 0 rows

new_dat1 %>%
  filter(if_any(where(is.double), is.na))
#Gives the rows with numeric values = NA (951)
#Now all NA rows turned to all 0 rows

table(is.na(new_dat1$Symbol))
table(is.na(new_dat1$`GFP0-1`))
table(is.na(new_dat1$`GFP0-2`))
table(is.na(new_dat1$`GFP0-3`))
table(is.na(new_dat1$`GFP10-1`))
table(is.na(new_dat1$`GFP10-2`))
table(is.na(new_dat1$`GFP10-3`))
table(is.na(new_dat1$`605-01`))
table(is.na(new_dat1$`605-02`))
table(is.na(new_dat1$`605-03`))
table(is.na(new_dat1$`605-10-1`))
table(is.na(new_dat1$`605-10-2`))
table(is.na(new_dat1$`605-10-3`))
table(is.na(new_dat1$`606-01`))
table(is.na(new_dat1$`606-02`))
table(is.na(new_dat1$`606-03`))
table(is.na(new_dat1$`606-10-1`))
table(is.na(new_dat1$`606-10-2`))
table(is.na(new_dat1$`606-10-3`))

dim(new_dat1)#[1] 10150    19 # Includes all 0 rows
allzero <- new_dat1 %>%
  filter(if_any(where(is.numeric), ~ .x == 0))
dim(allzero) #951

new_dat1 <- new_dat1 %>%
  filter(if_any(where(is.numeric), ~ .x > 0))

dim(new_dat1)
#9199
10150-951
#9199

table(is.na(new_dat1$Symbol)) #Ignore NA Symbol
#FALSE  TRUE 
#8923   276 
#Need to map the missing symbol later

dat_comp1<-new_dat1
dim(dat_comp1) #with symbols
#[1] 9199   19 #19 becauseof Symbol

#Matches: table(apply(new_dat1, 1, function(x) sum(is.na(x)))) #with symbols
#   0    1   18 
#8334    1  871 

head(new_dat1)
dim(new_dat1) #9199 rows 19 cols in new_dat1 #with symbols, without omit
dim(dat_comp1) #[1] 9199   19 #19 because of symbol, with omit
mapping1<-dat_comp1 
head(mapping1)
dim(mapping1)  #[1] 9199   19 #19 because of symbol
#mapping1<-mapping1[,-(2:19)]

#without symbols and before removing all 0 rows
head(new_dat2)
dim(new_dat2)
#[1] 10150    18

#dat_comp2 = na.omit(new_dat2)
dat_comp2<-dat_comp1
head(dat_comp1)
dim(dat_comp1)
#[1] 9199   19
#Removed all 0 rows

dat_comp2<-dat_comp2[,-c(1)]
head(dat_comp2)
dim(dat_comp2)
#[1] 9199   18 #After removing Symbol column
dim(dat_comp1)
#[1] 9199   19 #Before removing symbol column

#-----------------------------------------------------------------------------------------
# log2 transform - to reduce skewness of a measurement variable 
#to make data more symmetrical, which helps it meet the assumptions of statistical models
#-----------------------------------------------------------------------------------------
head(new_dat)
head(dat_comp)
dim(new_dat) #5376   16   # with NAs
dim(dat_comp) #4835   16  # without NAs
dat_comp_mat = as.matrix(log2(dat_comp)) # without NAs
dat_mat = as.matrix(log2(new_dat)) # with NAs

head(dat_comp_mat) #matrix and log transformed, without NA
dim(dat_comp_mat)#4835   16 #matrix and log transformed, without NA
head(dat_mat)
dim(dat_mat) #[1] 5376   16 #matrix and log transformed, with NA

dim(dat_comp1) #[1] 9199   19 before removing symbol
head(dat_comp1) 
dat_comp1.1<-dat_comp1
dat_comp1<-dat_comp1[,-1]
dim(dat_comp1) #[1] 9199   18 after removing symbol # without NAs
head(dat_comp1)
head(new_dat1)
dim(new_dat1) #[1] 9199   19 without all 18 cols NA, with symbol 
dim(dat_comp1) #[1] 9199   18 after removing symbol # without NAs
new_dat1.1<-new_dat1
new_dat1<-new_dat1[,-1]
head(new_dat1)
dim(new_dat1) #[1] 9199   18 after remving symbols # without NAs

head(dat_comp1)
dim(dat_comp1) #[1] 9199   18 after removing symbol # without NAs
dat_comp_mat1 = as.matrix(log2(dat_comp1+1)) # without NAs
head(dat_comp_mat1)
head(new_dat1)
dim(new_dat1)# with NAs
dat_mat1 = as.matrix(log2(new_dat1+1)) # with NAs

head(dat_comp_mat1)
head(dat_mat1)
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix
dim(dat_mat1) #[1] #9199   18 without NAs_log_matrix (Was previously with NAs, with filtering new way we only removed all 18 cols NAs)

dim(dat_comp1) #[1] #9199   18 without NAs_ not log
head(dat_comp1)

#density plot - plot is applied to log10 transformed data in code 
#dat_comp1 = na.omit(new_dat1) #with symbols
head(dat_comp1)#not log transformed, not normalized, omitted NA
dim(dat_comp1)#[1] 9199   18 without NAs_ not log, not normalized
for (i in 1:ncol(dat_comp1)){
  if (i==1){
    plot(density(log10(dat_comp1[,1])),main="Density plot across study samples - Total Proteomics_SY5Y",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp1[,i]))
    lines(den$x,den$y,col="blue")}
}

#Roughly normally distributed so  can use Limma trend

#dat_comp1.1<-dat_comp1 #with symbol
#dim(dat_comp1.1)#[1] 8334   19 #19 becauseof symbol
#new_dat1.1<-new_dat1
#dim(new_dat1.1)
#head(dat_comp1)
#dat_comp1<-dat_comp1[,-1]
#head(dat_comp1)
#dim(dat_comp1) #[1] 8334   18

#head(new_dat1.1)
#head(new_dat1)
#new_dat1<-new_dat1[,-1]
#head(new_dat1)
#head(new_dat1.1)
#dim(new_dat1)
#dim(new_dat1.1)

#dat_comp_mat1 = as.matrix(log2(dat_comp1)) # without NAs
#dat_mat1 = as.matrix(log2(new_dat1)) # with NAs
#dim(dat_comp_mat1)#[1] 8335   18

#dat_comp_mat1 = as.matrix(log2(dat_comp1)) # without NAs
#dat_mat1 = as.matrix(log2(new_dat1)) # with NAs
#dim(dat_comp_mat1)#[1] 8335   18

#-------------------------------------------------------------------------------------------------
# Quantile normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------

head(dat_comp_mat) #matrix and log transformed, without NA
dim(dat_comp_mat)#4835   16 #matrix and log transformed, without NA
head(dat_mat)
dim(dat_mat) #[1] 5376   16 #matrix and log transformed, with NA
dat_comp_mat_quant = limma::normalizeQuantiles(dat_comp_mat) # without NAs, log, norm
dat_mat_quant = limma::normalizeQuantiles(dat_mat) # with NAs, log, norm

head(dat_comp_mat1)  #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_mat1)
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_mat1) #[1] 9199   18 without NAs_log_matrix (Before used to be with NAs)
dat_comp_mat_quant1 = limma::normalizeQuantiles(dat_comp_mat1) # without NAs, log
dat_mat_quant1 = limma::normalizeQuantiles(dat_mat1) # with NAs
head(dat_mat_quant1) #[1] 9199  18 # with NAs, log, quantile norm
dim(dat_mat_quant1) #[1] 9199   18 # with NAs, log, quantile norm
head(dat_comp_mat_quant1) #[1] 9199   18 # without NAs, log, quantile norm
dim(dat_comp_mat_quant1) #[1] 9199   18  # without NAs, log, quantile norm

dim(dat_comp1) # [1] 9199   18 without NAs, not log,
head(dat_comp1)
dat_comp1_norm = limma::normalizeQuantiles(dat_comp1) # without NAs, not log, normalized
dim(dat_comp1_norm) ## [1] 9199   18 without NAs, not log, quantile normalized
head(dat_comp1_norm)

#https://www.youtube.com/watch?v=ecjN6Xpv6SE&t=279s

#-------------------------------------------------------------------------------------------------
# Median normalize - to remove any technical variation before differential expression analysis
#-------------------------------------------------------------------------------------------------
head(dat_comp_mat1)  #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_mat1)
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized #[1] 9206   18 with NAs_log_matrix (Before)
dat_comp_mat_median1 = limma::normalizeMedianValues(dat_comp_mat1) # without NAs, log
dat_mat_median1 = limma::normalizeMedianValues(dat_mat1) # without NA; log (used to be with NAs)
head(dat_mat_median1) #[1] 9199   18 # without NAs, log, median norm #[1] 9206  18 # with NAs, log, median norm (Before)
dim(dat_mat_median1) #[1] 9199   18 # without NAs, log, median norm #[1] 9206   18 # with NAs, log, median norm (Before)
head(dat_comp_mat_median1) #[1] 9199   18 # without NAs, log, median norm -- USED!!!
dim(dat_comp_mat_median1) #[1] 9199   18  # without NAs, log, median norm -- USED!!!

#=========================================================#
#Data to use for normalizing for phosphoproteomics:
#=========================================================#
Normalize<-dat_comp_mat_median1
head(Normalize)
dim(Normalize) #[1] 9199   18  # without NAs, log, median norm
Normalize<-as.data.frame(Normalize)
head(Normalize)
dim(Normalize) #[1] 9199   18  # without NAs, log, median norm

library(dplyr)
library(tibble)
Normalize2<-as.data.frame(Normalize)
head(Normalize2)
dim(Normalize2) #[1] 9199   18  # without NAs, log, median norm


Normalize2 <- Normalize2 %>%
  rownames_to_column(var = "Accession")
head(Normalize2)
dim(Normalize2)

identical(Normalize, Normalize2)
all.equal(Normalize, Normalize2)

head(Normalize2)
head(Normalize)

dim(Normalize2) #[1] 9199   19  # without NAs, log, median norm
dim(Normalize) #[1] 9199   18  # without NAs, log, median norm


#9199  19 - use this for normalizing phosphoproteomiccs
library(writexl)
###write_xlsx(Normalize2, "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\Normalize2v4.xlsx")
#write_xlsx(Normalize2, "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\Normalize2v4_copy.xlsx")
#write_xlsx(Normalize, "C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\Normalize2v3_copy.xlsx")
##########################################################################################

#===========#
#Box Plot
#===========#

# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation

par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat_comp_mat_quant, main="Boxplot normalized Intensities")

col <- c("red","red","red","blue","blue","blue","darkgreen","darkgreen","darkgreen","green","green","green","purple","purple","purple","orange","orange","orange")
col
fill <- c("red","blue","darkgreen","green","purple","orange")
fill
#par(mar = c(3, 3, 3, 10), mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2, xpd = TRUE)
head(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
dim(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
head(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
dim(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized

boxplot(dat_comp_mat1, main="Boxplot Before Normalization_Total Proteomics_SY5Y",col=col)
boxplot(dat_comp_mat_quant1, main="Boxplot Quantile Normalized Intensities_Total Proteomics_SY5Y",col=col)
boxplot(dat_comp_mat_median1, main="Boxplot Median Normalized Intensities_Total Proteomics_SY5Y",col=col)

#dim(dat_comp_mat1) #[1] 8334   18 without NAs_log_matrix
#dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
#legend("right", title="Treatment",legend=c("GFP-0","GFP-10","OGTKD-605-0","OGTKD-605-10","OGTKD-606-0","OGTKD-606-10"), fill=fill, inset = c(-0.155,0))

#inset depends on margin of plot window

#############################################################################################

################################################################################################################
#density plot - after normalization - not needed just for reference
dim(dat_comp1_norm) ## [1] 9199   18 without NAs, not log, quantile normalized
head(dat_comp1_norm)## [1] 9199   18 without NAs, not log, quantile normalized
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
head(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
dim(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
head(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
dim(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
head(dat_mat_median1) #[1] 9199   18 # without NAs, log, median norm #[1] 9206  18 # used to be with NAs, log, median norm (Before)
dim(dat_mat_median1) #[1] 9199   18 # without NAs, log, median norm #[1] 9206   18 # used to be with NAs, log, median norm (Before)
head(dat_mat_quant1) #[1] 9199  18 # without NAs, log, quantile norm #[1] 9206   18 # used to be with NAs, log, median norm (Before)
dim(dat_mat_quant1) #[1] 9199   18 # without NAs, log, quantile norm #[1] 9206   18 # used to be with NAs, log, median norm (Before)

head(dat_comp_mat1)  #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized #[1] 9206   18 with NAs_log_matrix (Before)
dim(dat_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized #[1] 9206   18 with NAs_log_matrix (Before)

any(is.na(dat_comp_mat_quant1))      # Checks for NAs or NaNs # without NAs_log_matrix; quantile normalized
any(is.infinite(dat_comp_mat_quant1)) # Checks for Inf (often caused by log2(0))
any(is.na(dat_comp_mat_median1))      # Checks for NAs or NaNs #without NAs_log_matrix; median normalized
any(is.infinite(dat_comp_mat_median1)) # Checks for Inf (often caused by log2(0))
any(is.na(dat_mat_median1))      # Checks for NAs or NaNs # without NAs, log, median norm #[1] 9206   18 # used to be with NAs, log, median norm (Before)
any(is.infinite(dat_mat_median1)) # Checks for Inf (often caused by log2(0))
any(is.na(dat_mat_quant1))      # Checks for NAs or NaNs # without NAs, log, quantile norm #[1] 9206   18 # used to be with NAs, log, median norm (Before)
any(is.infinite(dat_mat_quant1)) # Checks for Inf (often caused by log2(0)) # without NAs, log, quantile norm #[1] 9206   18 # used to be with NAs, log, median norm (Before)

#So no need imputataion for all (False) 
#-----------------------------------------#

# (If dat_mat_quant1 (True)

# Imputation (if warranted) - can skip As we are going to use the NA omitted data dat_comp_mat_quant1)

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
#================================================================================================
#y1-created dgelist and normalized...

#MDS AFTER REMOVING GENE COUNTS =0
#head(y1)
#plotMDS(y1, main="Multipledimensional scaling plot (MDS2)")
#head(dea)
#head(dat1)
#plotMDS(dat1, main="Multipledimensional scaling plot (MDS2)")
#=========================================================================================#
#==============#
#MDS Plot
#==============#
#head(dea)
#head(dea$tt) #values after comparisons with contrast matrix
#dim(dea$tt) # [1] 8477   24
#head(dat1) #normailzed and filtered values - use this for mds plot

#dat1 = dat_comp_mat_quant1 #[1] 8334   18  # without NAs, is log, quantile norm
#boxplot(dat_comp_mat1, main="Boxplot before normalization",col=col)
#boxplot(dat_comp_mat_quant1, main="Boxplot Quantile normalized Intensities",col=col)
#boxplot(dat_comp_mat_median1, main="Boxplot Median normalized Intensities",col=col)
#dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
#head(dat_comp_mat1)
#dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
#dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm
#dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quntile norm

head(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
dim(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
head(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
dim(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
head(dat1)
dim(dat1)
labels<-colnames(dat_comp_mat1)
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
#PCA<-get_pca1(dat1,labels)
#PCA$p
#https://stackoverflow.com/questions/33640492/change-point-colors-and-color-of-frame-ellipse-around-points
#https://www.geeksforgeeks.org/how-to-change-position-of-ggplot-title-in-r/

#===========================#
#PCA - Before Normalization
#===========================#
head(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
dim(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
head(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
dim(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized

#dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
#head(dat_comp_mat1)

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

#head(dat1)  
#dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
#head(dat_comp_mat1)
#dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
#dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm
#dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quntile norm

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

#=====================================================================================================================#

### Density plots by sample
###-----------------------------------------------
#head(dat1)
#dataa<-dat1 #quantile normalized
#head(dataa)
#head(dat)
#for (i in 1:ncol(dat)){
#  if (i==1){
#    plot(density(log10(dat[,1])),main="Density plot across study samples",xlab="Subjects",col="red",ylim=c(0,0.8))}
#  else {
#    den <- density(log10(dat[,i]))
#    lines(den$x,den$y,col="blue")}
#}

#for (i in 1:ncol(dataa)){
#  if (i==1){
#    plot(density(dataa[,1]),main="Density plot across study samples \n quantile normalized",xlab="Subjects",col="red")}
#  else {
#    den <- density(dataa[,i])
#    lines(den$x,den$y,col="blue")}
#}

#head(dat_comp1_norm)#not log transformed, normalized
#for (i in 1:ncol(dat_comp1_norm)){
#  if (i==1){
#    plot(density(log10(dat_comp1_norm[,1])),main="Density plot across study samples",xlab="Subjects",col="red",ylim=c(0,0.8))}
#  else {
#    den <- density(log10(dat_comp1_norm[,i]))
#    lines(den$x,den$y,col="blue")}
#}

#head(dat_comp1)#not log transformed, not normalized, omitted NA
#dim(dat_comp1)
#for (i in 1:ncol(dat_comp1)){
#  if (i==1){
#    plot(density(log10(dat_comp1[,1])),main="Density plot across study samples",xlab="Subjects",col="blue",ylim=c(0,0.8))}
#  else {
#    den <- density(log10(dat_comp1[,i]))
#    lines(den$x,den$y,col="blue")}
#}

#head(dat1)  
#dim(dat_comp_mat1) #[1] 8334   18 without NAs_is log_matrix; not normalized
#head(dat_comp_mat1)
#dim(dat_mat1) #[1] 9206   18 with NAs_log_matrix
#dim(dat_comp_mat_median1) #[1] 8334   18  # without NAs, log, median norm
#dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quntile norm
#density plot - plot is applied to log10 transformed data in code 
#dat_comp1 = na.omit(new_dat1) #with symbols

head(dat_comp1)#not log transformed, not normalized, omitted NA
dim(dat_comp1)#[1] 9199   18 without NAs_ not log, not normalized
for (i in 1:ncol(dat_comp1)){
  if (i==1){
    plot(density(log10(dat_comp1[,1])),main="Density plot across study samples - Total Proteomics_SY5Y",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density(log10(dat_comp1[,i]))
    lines(den$x,den$y,col="blue")}
}

#Roughly normally distributed so  can use Limma trend

head(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
dim(dat_comp_mat1) #[1] 9199   18 without NAs_log_matrix; not normalized
head(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
dim(dat_comp_mat_quant1) #[1] 9199   18 without NAs_log_matrix; quantile normalized
head(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
dim(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized

for (i in 1:ncol(dat_comp_mat1)){
  if (i==1){
    plot(density((dat_comp_mat1[,1])),main="Density plot across study samples after log2 transformation before normalization_for reference only",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((dat_comp_mat1[,i]))  #Already log tansformed so removed log
    lines(den$x,den$y,col="blue")}
}

for (i in 1:ncol(dat_comp_mat_quant1)){
  if (i==1){
    plot(density((dat_comp_mat_quant1[,1])),main="Density plot across study samples_After Quantile Normalization",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((dat_comp_mat_quant1[,i]))
    lines(den$x,den$y,col="blue")}
}

for (i in 1:ncol(dat_comp_mat_median1)){
  if (i==1){
    plot(density((dat_comp_mat_median1[,1])),main="Density plot across study samples_After Median Normalization",xlab="Subjects",col="blue",ylim=c(0,0.8))}
  else {
    den <- density((dat_comp_mat_median1[,i]))
    lines(den$x,den$y,col="blue")}
}

#------------------------------------------------
#Black and White box plot for sample Sam's data
#-------------------------------------------------
## Batch Effect Correction ##
# Here is where one should run diagnostics for batch effects (i.e., systematic differences in measurements due to technical factors rather than biological signal), if the experimental design introduces batch effects
# This dataset doesn't have any batch effects. 
# The measurements for each biological replicate were aggregated across technical replicates, therefore all samples incorporate the same variability due to technical factors

#for erk data, do a box plot
# boxplot: intensities of all 16 channels after data preprocessing and normalihttp://127.0.0.1:40129/graphics/17b44cfd-ded5-4327-9e49-56e1f2f75d97.pngzation
#par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#boxplot(dat_comp_mat_quant, main="Boxplot normalized Intensities")

#par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#boxplot(dat_comp_mat_quant1, main="Boxplot normalized Intensities")

# without NAs, not log, normalized - shows need to log transform for boxplot
#par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#boxplot(dat_comp1_norm, main="Boxplot normalized Intensities")

#==========================================================================================================#
# Differential Expression Analysis 
# using limma-trend - # This assumes uniform and normal distribution of data (seen in density plot)
#==========================================================================================================#

#When the library sizes are quite variable between samples, then the voom approach 
#is theoretically more powerful than limma-trend. 

#Library size could mean one of two things: the total number of reads that were sequenced 
#in the run or the total number of mapped reads

#Read = sequenced fragment of cDNA (obtained from RNA) here ptn

datold<-dat
dat = dat_comp_mat_quant 

head(dat_mat_quant1) #[1] 9206  18 # with NAs, log, quantile norm
dim(dat_mat_quant1) #[1] 9206   18 # with NAs, log, quantile norm
head(dat_comp_mat_quant1) #[1] 8334   18 # without NAs, log, quantile norm
dim(dat_comp_mat_quant1) #[1] 8334   18  # without NAs, log, quantile norm
head(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized
dim(dat_comp_mat_median1) #[1] 9199   18 without NAs_log_matrix; median normalized

datold1<-dat1
#dat1 = dat_comp_mat_quant1 #[1] 8334   18  # without NAs, log, quantile norm
dat1 = dat_comp_mat_median1 #[1] 8334   18  # without NAs, log, quantile norm
head(dat1)
dim(dat1) #[1] 9199   18 #removed NA, without Symbol, norm, log

# Change column names to not contain "-", since this will conflict with naming syntax in contrast matrix
colnames(dat)
colnames(dat) = gsub(pattern = "-", replacement = "", x = colnames(dat)) #removing "-dash" as R can misinterprt it as subtraction
colnames(dat)

colnames(dat1)
colnames(dat1) = gsub(pattern = "-", replacement = "", x = colnames(dat1))
colnames(dat1)
colnames(dat1) <- c('GFP01','GFP02','GFP03','GFP101','GFP102','GFP103','T60501','T60502','T60503','T605101','T605102','T605103','T60601','T60602','T60603','T606101','T606102','T606103')
colnames(dat1)
head(dat1)

#=======================#
#Mapping
#=======================#

# Getting a data.frame with 'Accession' and 'Symbol' columns
mapping = add_symbols(data.frame(Accession = rownames(dat)))
head(mapping)
table(is.na(mapping))

head(mapping1)
head(dat1)
dim(dat1) #[1] 9199   18 Without symbol
dim(mapping1) #[1] 9199   19 With Symbol
#both dat1 and mapping1 match
#so, try to separate symbol and accession for mapping1

head(mapping)
mapping1.1<-mapping1
head(mapping1.1)
dim(mapping1.1) #[1] 9199   19
head(mapping1) #[1] 9199   19
dim(mapping1) #[1] 9199   19
#make dim(mapping1) 9199 2
mapping2<-as.data.frame(cbind(rownames(mapping1),mapping1$Symbol))
dim(mapping2)
head(mapping2)
head(mapping)
colnames(mapping2) <- c('Accession','Symbol')
colnames(mapping2)
head(mapping2)
head(mapping)
dim(mapping2) #[1] 9199    2

#Checking for missing Symbols:
table(is.na(mapping2$Symbol))
#FALSE  TRUE 
#8923   276 

table(is.na(mapping2$Accession))
#FALSE 
#9199 

needs_mapping2 <- is.na(mapping2$Symbol)
#This creates a logical vector (TRUE/FALSE) where TRUE corresponds to rows where Symbol is NA.
head(needs_mapping2)
head(mapping2)
missing_ids2 <- unique(mapping2$Accession[needs_mapping2])
#This subsets the Accession column using that logical vector,
#keeping only the Accessions where needs_mapping2 is TRUE.
#This ensures that if an Accession ID has multiple missing symbols,
#it only appears once in the final missing_ids2 list.
head(missing_ids2)
length(missing_ids2)
#[1] 276 (old=322 before replacing blank)

#-----------------------------------------------------------------------#
#Used Biomart first but had more missing symbols (747) (322 UNIQUE NAs)
#-----------------------------------------------------------------------#

library(biomaRt)

#BiocManager::install("fastmap")
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#useMart: You are telling R to go online and connect to the Ensembl "market" (Mart).
#dataset = "hsapiens_gene_ensembl": You are specifically asking for the Human (Homo sapiens) gene database.
#mart
#bm_map2 <- getBM(attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol"), 
#               filters = "uniprotsptrembl", 
#              values = (pmapping2$Accession), 
#             mart = mart)

mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast") # Options: 'useast', 'uswest', 'asia'
mart


# Initialize the BioMart dataset with the US East mirror
#mart <- useEnsembl(biomart = "genes", 
 #                  dataset = "hsapiens_gene_ensembl", 
  #                 mirror = "useast")

# Print the object to verify connection
#mart

#library(biomaRt)

# Bypass useEnsembl and connect directly to a dedicated mirror host
#mart <- useMart(biomart = "genes", 
 #               dataset = "hsapiens_gene_ensembl", 
  #              host = "https://useast.ensembl.org") # Or "https://asia.ensembl.org"

#library(biomaRt)

# Establish a strict connection directly to the main server block
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
 #               dataset = "hsapiens_gene_ensembl", 
  #              host = "https://www.ensembl.org",
   #             )
#mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
 #                        dataset = "hsapiens_gene_ensembl",
  #                       host = "http://www.ensembl.org")

head(mapping2)
dim(mapping2)

bm_map <- getBM(attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol"), 
                filters = "uniprot_gn_id", 
                values = (mapping2$Accession), 
                mart = mart)
head(bm_map) #uniprotswissprot uniprotsptrembl hgnc_symbol

#bm_map2 <- getBM(
# attributes = c("uniprotswissprot", "uniprotsptrembl", "uniparc", "hgnc_symbol"),
#filters = "uniprot_gn_id", # Or use "uniparc" if your input list contains UPI IDs
#values = pmapping2$UniprotAcc,
#mart = mart
#)

# Combine both columns into one 'Accession' for matching
head(bm_map) #uniprotswissprot uniprotsptrembl hgnc_symbol
dim(bm_map) #[1] 54422     3
bm_clean <- bm_map %>%
  pivot_longer(cols = c(uniprotswissprot, uniprotsptrembl), 
               values_to = "Accession") %>%
  filter(Accession != "", !is.na(Accession)) %>%
  dplyr::select(Accession, hgnc_symbol) %>%  # Explicitly use dplyr
  distinct()

#bm_clean <- bm_map %>%
# pivot_longer(cols = c(uniprotswissprot, uniprotsptrembl), 
#            values_to = "UniprotAcc") %>%
#filter(UniprotAcc != "", !is.na(UniprotAcc)) %>%
#dplyr::select(UniprotAcc, hgnc_symbol)   # Explicitly use dplyr

#pivot_longer(...): This moves data from "wide" to "long." 
#It takes the two separate columns for UniProt IDs 
#(uniprotswissprot and uniprotsptrembl) and stacks them into a single new column
#called UniprotAcc.

#filter(UniprotAcc != "", !is.na(UniprotAcc)): Since many genes might only have 
#one type of ID (or none), the pivot creates a lot of empty or NA rows. 
#This line scrubs those out so you only keep valid IDs.
#dplyr::select(UniprotAcc, hgnc_symbol): This drops any extra metadata columns
#from the original biomaRt query, leaving you with just the two columns you 
#actually need for a mapping.
#dplyr::distinct(): This is the final cleanup. Because one gene symbol can be 
#associated with multiple protein isoforms (or the same ID might have appeared 
#in both original columns), this ensures that every UniProt-to-Symbol pair is 
#unique, preventing your data from "exploding" (duplicating rows) when you 
#later join this to another dataset.

head(bm_clean)
dim(bm_clean) # with distinct [1] 46110   2 
head(mapping2)
dim(mapping2)#[1] 9199     2
mapping4 <- mapping2 %>%
  left_join(bm_clean, by = "Accession")

head(mapping4)
dim(mapping4) #[1] 9232     3?? - without distincct #[1] 12315     3 - with distinct
#[1] 9232     3

dim(mapping2)
#[1] 9199     2

#pmapping 4 nrows increased due to dupliates

# Need to Deduplicate bm_clean first
head(bm_clean)
table(duplicated(bm_clean$Accession))
#FALSE  TRUE 
#46055  55

bm_clean_unique <- bm_clean %>%
  distinct(Accession, .keep_all = TRUE)
table(duplicated(bm_clean_unique$Accession))
#FALSE 
#46055 

# Joining
mapping4 <- mapping2 %>%
  left_join(bm_clean_unique, by = "Accession")

# Result changed to 12298 rows
dim(mapping4) #[1] 9199     3
dim(mapping2) #[1] 9199     2 
#Rows match
head(mapping4)
head(mapping2)

table(is.na(mapping2$Symbol))
table(is.na(mapping4$Symbol))
#match
#FALSE  TRUE 
#8923   276
table(is.na(mapping4$hgnc_symbol))
#FALSE  TRUE 
#8804   395 


mapping4 %>%
  group_by(Symbol) %>%
  filter(n() > 1) %>%
  arrange(Symbol)

#Different accesion same symbol could be isomers

head(mapping4)
table(is.na(mapping4$Symbol))
#FALSE  TRUE 
#8923   276 
table(is.na(mapping4$hgnc_symbol))
#FALSE  TRUE 
#8804   395 

# Replace literal dashes with NA
#hgnc_symbol blank, so replace with NA
mapping4$hgnc_symbol[mapping4$hgnc_symbol == "-"] <- NA
mapping4$hgnc_symbol[mapping4$hgnc_symbol == ""] <- NA
mapping4$Symbol[mapping4$Symbol == "-"] <- NA
mapping4$Symbol[mapping4$Symbol == ""] <- NA
table(is.na(mapping4$hgnc_symbol))
#FALSE  TRUE 
#8797   402 
table(is.na(mapping4$Symbol))
#FALSE  TRUE 
#8923   276 

library(dplyr)

mapping5 <- mapping4 %>%
  mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol), Symbol, hgnc_symbol))
#Symbol from original mass spec file; hgnc_symbol - biomart mapped and updated nas with symbol from mass spec file

head(mapping5) #Main to use for joining for mapping integration
dim(mapping5)
#[1] 9199    3

table(is.na(mapping4$Symbol))
#FALSE  TRUE 
#8923   276 
table(is.na(mapping5$Symbol))
#FALSE  TRUE 
#8923   276
table(is.na(mapping4$hgnc_symbol))
#FALSE  TRUE 
#8804   395 
table(is.na(mapping5$hgnc_symbol))
#FALSE  TRUE 
#9056   143

needs_mapping3 <- is.na(mapping5$hgnc_symbol)
#This creates a logical vector (TRUE/FALSE) where TRUE corresponds to rows where Symbol is NA.
head(needs_mapping3) # TRUE FALSE
head(mapping5)
missing_ids3 <- unique(mapping5$Accession[needs_mapping3])
#This subsets the Accession column using that logical vector,
#keeping only the Accessions where needs_mapping2 is TRUE.
#This ensures that if an Accession ID has multiple missing symbols,
#it only appears once in the final missing_ids2 list.
head(missing_ids2)
length(missing_ids2) #Before mapping missing IDs to Symbols from original MAss Spec File
#[1] 276 
head(missing_ids3) #After mapping some missing IDs to Symbols from original MAss Spec File
length(missing_ids3)
#[1] 143 
head(mapping2)
dim(mapping2) #[1] 9199    2 - original file for mapping but befroe any mapping done
table(is.na(mapping2$Symbol))
#FALSE  TRUE 
#8923   276

#---------------------------------#
#USe missing_ids3 for mapping
#---------------------------------#

bm_map <- getBM(attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol"), 
                filters = "uniprot_gn_id", 
                values = (missing_ids3), 
                mart = mart)
head(bm_map) #uniprotswissprot uniprotsptrembl hgnc_symbol
#0 rows

bm_map <- getBM(
  attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol"), 
  filters = "uniprotswissprot",  # Changed filter
  values = missing_ids3, 
  mart = mart
)

head(bm_map) #uniprotswissprot uniprotsptrembl hgnc_symbol
#0 rows

#So:

#-----------------------------------------------------------------------#
# 1. Query Org.HS.db for ONLY those specific IDs
#-----------------------------------------------------------------------#

# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
head(mapping5)
head(missing_ids3)
length(missing_ids3)
map_res <- mapIds(org.Hs.eg.db, 
                  keys = missing_ids3, 
                  column = "SYMBOL", 
                  keytype = "UNIPROT", 
                  multiVals = "first")
#Symbol: This is the information I want back." You are specifically asking for the HGNC Gene Symbol (like ELOA).
#Keytype: This tells R: "This is what my input IDs are."
#multiVals = "first"

head(map_res)
# Add it back to your table

head(missing_ids3) #List of accessions
length(missing_ids3) #143
missing_idsmap<-as.data.frame(missing_ids3)
head(missing_idsmap)

#pmapping6<-pmapping2
#head(pmapping6)
missing_idsmap$Symbol <- map_res
head(missing_idsmap)
any(is.na(missing_idsmap))
table(is.na(missing_idsmap$missing_ids3))
#FALSE 
#143
table(is.na(missing_idsmap$Symbol))
#FALSE  TRUE 
#13     130
head(missing_idsmap) #Has 13 more symbols mapped

#missing_idsmap like bm clean but org hs

head(missing_idsmap) #143  2 #Has 13 more symbols mapped
dim(missing_idsmap)  #143  2
head(missing_ids3) #was used for missing_idsmap
length(missing_ids3) #143

#Org.Hs standardized mapping table
org_hs_clean <- missing_idsmap %>%
  dplyr::select(Accession = missing_ids3, symbol_org_hs = Symbol) %>%
  filter(!is.na(symbol_org_hs))
head(org_hs_clean)
dim(org_hs_clean)#[1] 13   2 #Will use later for joining - only has 13 mapped symbols

head(missing_idsmap)
dim(missing_idsmap) #143  2
needs_mapping4 <- is.na(missing_idsmap$Symbol)
head(needs_mapping4) #True or False
length(needs_mapping4)

missing_ids4 <- unique(missing_idsmap$missing_ids3[needs_mapping4])
head(missing_ids4)
length(missing_ids4)
#[1] 130
143-13

#-----------------------------------------------------------------------#
# 2. Query UniProt for ONLY those specific IDs
#-----------------------------------------------------------------------#

library(UniProt.ws)
up <- UniProt.ws(taxId = 9606)
keytypes(up)
columns(up)
# Querying "GENE_NAMES" to get the symbols (like ELOA)

#extra_map <- select(up, 
#                   keys = missing_ids, 
#                  columns = c("uniparc_id","xref_ensembl","gene_names", "gene_synonym","gene_primary"),
#                 keytype = "UniProtKB")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("UniProt.ws")

head(mapping5)
dim(mapping5)
extra_map <- AnnotationDbi::select(up,
                                   keys = missing_ids4,
                                   columns = c("xref_refseq","uniparc_id","xref_ensembl","gene_names", "gene_synonym", "gene_primary"),
                                   keytype = "UniProtKB")

head(extra_map)
colnames(extra_map)
table(is.na(extra_map))
#FALSE  TRUE 
#348   692 

nrow(extra_map) #[1] 130

nrow(extra_map) == length(missing_ids4)
#[1] TRUE

#pmapping4 has bm_clean=biomart joined
#missing_idsmap like bm clean but org hs
#extra_map = uniprot.ws mapping (Gene.Names)

#Uniprot.ws standardized mapping table
head(extra_map)

uniprot_ws_clean <- extra_map %>%
  mutate(symbol_uni_ws = sub(" .*", "", Gene.Names)) %>%
  dplyr::select(UniprotAcc = From, symbol_uni_ws) %>%
  filter(!is.na(symbol_uni_ws))

head(uniprot_ws_clean)
dim(uniprot_ws_clean) #[1] 11  2

#Change column names
head(uniprot_ws_clean)
dim(uniprot_ws_clean) #[1] 11  2
names(uniprot_ws_clean)[names(uniprot_ws_clean) == 'UniprotAcc'] <- 'Accession'
head(uniprot_ws_clean)
dim(uniprot_ws_clean) #[1] 11  2

ids_to_fix <- extra_map$From[is.na(extra_map$Gene.Names)]
head(ids_to_fix) #List of accession to fix
length(ids_to_fix) 
#119
130-11

# Copy these from your console
cat(ids_to_fix, sep = ", ") 

head(ids_to_fix)
length(ids_to_fix)

#=====================================================================================#
#-------------------------#
#Uniprot.ws troubleshoot
#-------------------------#

#extra_map <- AnnotationDbi::select(up,
#                                   keys = mapping5$Accession,
#                                   columns = c("xref_refseq","uniparc_id","xref_ensembl","gene_names", "gene_synonym", "gene_primary"),
#                                   keytype = "UniProtKB")

#intersect(c("xref_refseq","uniparc_id","xref_ensembl","gene_names", "gene_synonym", "gene_primary"), columns(up))


#extra_map_test <- AnnotationDbi::select(
#  up,
#  keys = mapping5$Accession[1:5],  # Subset to just the first 5 elements
#  columns = c("gene_names"),
#  keytype = "UniProtKB"
#)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# Update BiocManager first
#BiocManager::install(version = "3.22") # Or your system's latest stable version

# Force a clean install of the core querying packages
#BiocManager::install(c("AnnotationDbi", "UniProt.ws", "httr2"), force = TRUE, ask = FALSE)

#extra_map <- mapIds(up, keys = missing_ids4, 
#                    column = "UNIPARC_ID", # Choose one column for mapIds
#                    keytype = "UNIPROTKB")

#extra_map <- AnnotationDbi::select(up,
#                                   keys = missing_ids3,
#                                   columns = c("xref_refseq","uniparc_id","xref_ensembl","gene_names", "gene_synonym", "gene_primary"),
#                                   keytype = "UniProtKB")

#Note: None of this worked May 17th but all were fine May 18th 2026
#=====================================================================================#

#-----------------------------------------------------------------------#
#Manual map
#-----------------------------------------------------------------------#

# Create the manual lookup table
# Create the manual lookup table for your exact list of 119 accessions
manual_map <- data.frame(
  UniprotAcc = c(
    "A0A087WY61", "A0A0A0MTS2", "A0A087X2D8", "A0A087WTS8", "E7EVH7", 
    "A0A087WVV2", "A0A0A0MTH9", "F8VV04", "A0A2R8Y611", "B1AKJ5", 
    "H3BNC9", "A0A0A0MT60", "A0A087WXU3", "A0A1W2PNV4", "E7EST9", 
    "F8W914", "F6TLX2", "A0A096LPI6", "C9J840", "J3KR44", 
    "E9PH62", "H7BXF4", "J3KNN5", "A0A1W2PNZ9", "A0A087X211", 
    "F5H5P2", "A0A087WXS7", "A0A1W2PS05", "H3BN98", "A0A494C169", 
    "F6TB26", "A0A087WTA5", "D6RB85", "A0A087WUC6", "F5GWE5", 
    "A0A0A0MR59", "A0A0B4J2E5", "A0A087WUM0", "K7EJH8", "Q5SQP8", 
    "M0R2C6", "H0Y7F4", "H7BYZ3", "A0A2R8Y422", "F8WAK2", 
    "E9PSI1", "F8WAN9", "A0A0A0MRB7", "A0A3B3ITR4", "S4R400", 
    "A0A0A0MR09", "A0A087WUQ6", "A6NML8", "H0YHG0", "A0A0A0MSJ9", 
    "A0A087WT12", "J3QSX6", "H0YIV9", "H0YGG7", "A0A087WV05", 
    "D6R9U4", "A0A087WVF8", "A0A0C4DGZ1", "J3KPL8", "A0A286YEW0", 
    "A0A0D9SFF8", "C9JLM9", "Q6ZSU2", "Q5JUV4", "A6NE70", 
    "A0A494C0C4", "V9GYY5", "J3KQS5", "F6Y5H0", "Q9BZD3", 
    "G3V4G9", "W4VSQ3", "A0A1W2PQJ7", "A0A087WVE9", "A0A087WU02", 
    "A0A087X060", "B5MD58", "A0A2R8Y600", "K7EP59", "K7ENM7", 
    "A0A087WVS8", "A0A0A0MRG8", "A0A1W2PP04", "A0A087WSZ7", "A0A087X1N7", 
    "A0A0U1RRB3", "F6WH68", "G3V3Q6", "A0A087WTQ6", "H7C3M1", 
    "F6QB42", "H0YEF3", "H7C5V1", "A0A0A0MS59", "E7EP41", 
    "A0A1W2PRW1", "K7N7A0", "A0A1B0GUF4", "A0A1W2PPT3", "A0A1W2PQ13", 
    "A0A087X1I5", "A0A0A6YYF3", "Q8IWY7", "A0A1W2PRV3", "U3KQK5", 
    "H0Y2W4", "A0A075B6Q3", "A6NJE6", "A0A2R8YEI5", "O00370", 
    "A0A2R8YCV2", "A9X3U0", "A0A087X0T7", "M0R3E8"
  ),
  ManualSymbol = c(
    "IPO4", "ZNF208", "OR4F5", "OR4F5", "MUC2", 
    "MUC5B", "ZNF208", "POTEKP", "POTEKP", "IPO4", 
    "CAMSAP1", "ZNF208", "NBPF14", "MUC2", "KMT2C", 
    "ABCC1", "C1orf159", "ZNF257", "ZNF208", "OR4F5", 
    "OR4F29", "EPHA1", "ZNF257", "MUC16", "OR4F29", 
    "TBC1D3", "MUC5B", "MUC16", "FAM234A", "TBC1D3", 
    "NBPF14", "MUC5B", "POTEF", "MUC5B", "NBPF20", 
    "ZNF705A", "LPA", "MUC5B", "SMN2", "VTI1A", 
    "AMY1A", "ZNF705A", "ZNF813", "POTEF", "LPA", 
    "OR4F16", "GOLGA6B", "ZNF705A", "FRG2C", "FRG2C", 
    "ZNF705A", "MUC5B", "MUC5B", "ZNF718", "ZNF705A", 
    "OR4F5", "SEC22B", "OR2T1", "ZNF727", "MUC5B", 
    "NUTM2B", "MUC5B", "SEC22B", "ZNF99", "NUTM2B", 
    "ZNF99", "ZNF780A", "PCDHB16", "NUP85", "PRAMEF1", 
    "TBC1D3", "ZNF543", "ZNF736", "PRR20A", "HCN3", 
    "ZNF705D", "CTAGE5", "MUC16", "MUC5B", "OR4F16", 
    "NBPF10", "HDAC5", "POTEF", "GAGE12F", "GAGE12I", 
    "MUC5B", "ZNF705A", "MUC2", "MMP24", "ZNF727", 
    "SRGAP2C", "NBPF10", "SRGAP2C", "NOTCH2NL", "NOTCH2NL", 
    "KIR2DL4", "ZNF735", "GTF2IRD2", "ZNF705A", "MUC16", 
    "MUC16", "FCGBP", "FCGBP", "MUC16", "MUC2", 
    "ZNF735", "KRTAP4-11", "TEX41", "MUC2", "KRTAP4-11", 
    "ZNF705E", "ZNF844", "ZNF844", "POTEF", "LINE1", 
    "POTEF", "SMN1", "SMN1", "AMY2A"
  ),
  stringsAsFactors = FALSE
)

head(manual_map)
dim(manual_map)
head(ids_to_fix)
length(ids_to_fix)
(manual_map$UniprotAcc) == ids_to_fix
#All true

table(is.na(manual_map$ManualSymbol))
#FALSE   
#119     

library(dplyr)

# Filter for missing symbols
na_symbols <- manual_map %>% 
  filter(is.na(ManualSymbol))

# Display the results
print(na_symbols)
#o rows

head(manual_map)

# Update the specific ID in the existing manual_map dataframe - Not needed here
#manual_map$ManualSymbol[manual_map$UniprotAcc == "A0A1P0AZG4"] <- "LCOR"
#manual_map$ManualSymbol[manual_map$UniprotAcc == "A0A2Q2TH77"] <- "GOLGA2"
##manual_map$ManualSymbol[manual_map$UniprotAcc == "A0A494C1P3"] <- "MOSC domain containing protein"
#manual_map$ManualSymbol[manual_map$UniprotAcc == "A0A494C1P3"] <- "MOSC2"

# Verify the change
#anual_map %>% filter(UniprotAcc == "A0A1P0AZG4")
#manual_map %>% filter(UniprotAcc == "A0A2Q2TH77")
#manual_map %>% filter(UniprotAcc == "A0A494C1P3")

head(manual_map)
table(is.na(manual_map$ManualSymbol))
#FALSE 
#119
dim(manual_map) #119   2

#Change Column name to Accession:
names(manual_map)[names(manual_map) == 'UniprotAcc'] <- 'Accession'
head(manual_map)
dim(manual_map) #[1] 119  2

#---------------------------------------#
#Summary of mapping files
#---------------------------------------#

#Thus,

#mapping4 <- mapping2 %>%
#  left_join(bm_clean_unique, by = "Accession")
head(bm_clean)
dim(bm_clean) #[1] 46110   2
head(bm_clean_unique)
dim(bm_clean_unique) #[1] 46055   2 #Biomart mapping
#mapping5 <- mapping4 %>%
#  mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol), Symbol, hgnc_symbol))
#Took mapping 4 and mapped Symbol from Mass Spec for mapping 5
head(mapping5)
table(is.na(mapping5$hgnc_symbol))
#FALSE  TRUE 
#9056   143

#needs_mapping3 <- is.na(mapping5$hgnc_symbol)
#This creates a logical vector (TRUE/FALSE) where TRUE corresponds to rows where Symbol is NA.
#head(needs_mapping3) # TRUE FALSE
#head(mapping5)
#missing_ids3 <- unique(mapping5$Accession[needs_mapping3])
head(missing_ids3)
length(missing_ids3) #143

#map_res <- mapIds(org.Hs.eg.db, 
#                  keys = missing_ids3, 
#                  column = "SYMBOL", 
#                  keytype = "UNIPROT", 
#                  multiVals = "first")
#missing_idsmap<-as.data.frame(missing_ids3)
#head(missing_idsmap)
#missing_idsmap$Symbol <- map_res #(from org.hs.db)

head(missing_idsmap) # Has mapped org.hs.db
table(is.na(missing_idsmap$Symbol))
dim(missing_idsmap) #[1] 143   2 includes 13 form org.hs.db and all NAs
head(org_hs_clean)
dim(org_hs_clean) #[1] 13   2 #(from org.hs.db)
143-13
head(extra_map)
dim(extra_map) #[1] 130   2
143-13
head(uniprot_ws_clean)
dim(uniprot_ws_clean) #[1] 11   2
130-11
head(manual_map)
dim(manual_map) #[1] 119   2
130-11

dim(org_hs_clean) + dim(uniprot_ws_clean) + dim(manual_map)
#[1] 13   2         #[1] 11   2             #[1] 119   2

(nrow(org_hs_clean) + nrow(uniprot_ws_clean) + nrow(manual_map)) == (length(missing_ids3))

#[1] 13   2         #[1] 11   2             #[1] 119   2

head(org_hs_clean)
head(uniprot_ws_clean)
head(manual_map)

dim(org_hs_clean) #13
dim(uniprot_ws_clean) #11
dim(manual_map) #119
13+11+119
#143

#mapping5 has bm_clean=biomart joined; also has missing symbols form biomart mapped to Symbol from MAss Spec
#missing_idsmap has org hs mapped with NAs; org_hs_clean = mapping table without NAs
#extra_map = uniprot.ws mapping (Gene.Names); uniprot_ws_clean = mapping table without NAs
#manual_map is  all manual map

#-----------------------------------------------------------------------#
#Integration of all mapping:
#-----------------------------------------------------------------------#

library(dplyr)

# 1. Standardize your intermediate mapping tables for joining
# From org.Hs.eg.db
head(org_hs_clean) 
dim(org_hs_clean)#[1] 13   2 

# From UniProt.ws (cleaning Gene.Names to get just the first symbol)
head(uniprot_ws_clean) 
dim(uniprot_ws_clean) #[1] 11  2

head(manual_map) 
dim(manual_map) #[1] 119  2

nrow(org_hs_clean)+nrow(uniprot_ws_clean)+nrow(manual_map)
#[1] 143

#table(is.na(pmapping4$UniprotAcc))
#FALSE 
#12298 
#table(is.na(pmapping4$hgnc_symbol))
#FALSE  TRUE 
#11568   747

#needs_mapping2 <- is.na(pmapping4$hgnc_symbol)
#head(needs_mapping2)

#missing_ids2 <- unique(pmapping4$UniprotAcc[needs_mapping2])
#missing_ids3 <- unique(mapping5$Accession[needs_mapping3])
head(missing_ids3)
length(missing_ids3)
#[1] 143


# 2. Sequential Join
head(mapping5)
mapping10 <- mapping5 %>%
  # Join automatic sources
  left_join(org_hs_clean, by = "Accession") %>%
  left_join(uniprot_ws_clean, by = "Accession") %>%
  # Join manual source
  left_join(manual_map, by = "Accession") %>%
  
  # 3. Create 'source' and final 'hgnc_symbol' columns
  mutate(
    Symbol = coalesce(hgnc_symbol, symbol_org_hs, symbol_uni_ws, ManualSymbol),
    
    source = case_when(
      !is.na(hgnc_symbol) ~ "biomart",
      !is.na(symbol_org_hs) ~ "org_hs_db",
      !is.na(symbol_uni_ws) ~ "uniprot_ws",
      !is.na(ManualSymbol) ~ "manual",
      TRUE ~ "accession_only"
    ))

#This case_when logic acts like a priority ladder. 
#It evaluates each row from top to bottom and assigns a label based on the
#first column that isn't empty (NA).
#Here is the step-by-step explanation:
# !is.na(hgnc_symbol) ~ "biomart":
# R looks at the original BioMart column first.
#If a symbol is found here, it tags it "biomart" and stops looking at 
#other sources for this row. This preserves your primary data source.
#!is.na(symbol_org_hs) ~ "org_hs_db":
#  If BioMart was empty, R checks the org.Hs.eg.db results.
#If a match is found here, it tags it "org_hs_db".
#!is.na(symbol_uni_ws) ~ "uniprot_ws":
#  If both previous sources were empty, it checks the UniProt.ws results.
#If found, it tags it "uniprot_ws".
#!is.na(ManualSymbol) ~ "manual":
#  If all automated tools failed, it looks at your manual map (the 110 IDs we fixed).
#If found, it tags it "manual".
#TRUE ~ "accession_only":
#  This is the "catch-all." If every single source above was NA,
#it assigns the label "accession_only".
#In your coalesce function, this corresponds to using the UniProt Accession itself as the symbol.


# Final fallback: use Accession ID if everything else is NA
#hgnc_symbol = coalesce(hgnc_symbol, UniprotAcc)
#) %>%

# Clean up temporary columns
#dplyr::select(-symbol_org_hs, -symbol_uni_ws, -ManualSymbol)

# 4. Verification
#table(pmapping_final$source)

head(mapping10) #has all columns including source
table(is.na(mapping10$Symbol))
#FALSE 
#9199
table(is.na(mapping10$Accession))
#FALSE 
#9199 
# Count how many are tagged as accession_only
table(mapping10$source == "accession_only")
#FALSE 
#9199 # So mapped all Symbols  

mapping11<-mapping10
head(mapping11)
dim(mapping11)#[1] 9199    7
mapping11<-mapping11[,-c(3:7)]
head(mapping11)
dim(mapping11)
#[1] 9199    2
table(is.na(mapping11$Symbol))
#FALSE 
#9199

table(is.na(mapping11$Accession))
#FALSE 
#9199


#-----------------------------------------------------------------------#
#Final files for DE analysis
#-----------------------------------------------------------------------#

head(mapping10) #J_Has all columns for mapping reference
dim(mapping10) #9199  7
head(mapping11) #Has only 2 columns for DEA
dim(mapping11) #9199  2

head(dat1) #Median normalized matrix for DEA
#where, dat1 = dat_comp_mat_median1 #[1] 9199   18  # without NAs, log, median norm #Median normalized matrix for DEA
dim(dat1) #[1] 9199   18 #removed NA, without Symbol, norm, log
head(mapping11) #9199  2
dim(mapping11) #9199  2

#---------------------#
# Design matrix
#---------------------#

colnames(dat)
tx = extract_string(colnames(dat), "_", 2)
tx
X = model.matrix(~0 + tx)
X
colnames(X) = sort(unique(tx))
X #design matrix - with treatment labels

# Contrast matrix
#make sure to not include any numbers as it would cause problems
contrasts = c("TMG-saline", "OGTKO-OGTWT", "OGAKO-OGAWT")
contrasts
cm_dat = makeContrasts(contrasts = contrasts, levels = colnames(X))
cm_dat #contrast matrix - the comparisons

# Design matrix
head(dat1)
colnames(dat1)<-c('GFPCTRL_1','GFPCTRL_2','GFPCTRL_3','GFPTEN_1','GFPTEN_2','GFPTEN_3','OGTCTRL_1','OGTCTRL_2','OGTCTRL_3','OGTTEN_1','OGTTEN_2','OGTTEN_3','ROGTCTRL_1','ROGTCTRL_2','ROGTCTRL_3','ROGTTEN_1','ROGTTEN_2','ROGTTEN_3')
colnames(dat1)
head(dat1)
head(dat)
#tx = extract_string(colnames(dat), "_", 2)
tx
txa = extract_string(colnames(dat1), "_", 1)
txa
#X = model.matrix(~0 + tx)
Xa = model.matrix(~0 + txa)
Xa
#colnames(X) = sort(unique(tx))
colnames(Xa) = sort(unique(txa))
Xa #design matrix
X

# Contrast matrix
head(dat)
colnames(X)
contrasts
colnames(Xa)
#contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
#contrasts = c("TMG-saline", "OGTKO-OGTWT", "OGAKO-OGAWT")
contrastsa = c("GFPTEN-GFPCTRL", "OGTTEN-OGTCTRL", "ROGTTEN-ROGTCTRL")
contrastsa
Xa
#cm_dat = makeContrasts(contrasts = contrasts, levels = colnames(X))
cm_dataa = makeContrasts(contrasts=contrastsa, levels=colnames(Xa))
cm_dataa #contrast matrix
cm_dat

#Try to make another contrast matrix for OGT VS GFP
#contrasts = c("TMG-saline", "OGTKO-OGTWT", "OGAKO-OGAWT")
contrastsb = c("GFPTEN-GFPCTRL", "OGTTEN-OGTCTRL", "ROGTTEN-ROGTCTRL","OGTCTRL-GFPCTRL","OGTTEN-GFPTEN","ROGTCTRL-GFPCTRL","ROGTTEN-GFPTEN")
contrastsb
Xa
#cm_dat = makeContrasts(contrasts = contrasts, levels = colnames(X))
cm_data = makeContrasts(contrasts=contrastsb, levels=colnames(Xa))
cm_data #contrast matrix
cm_dat

#Check with Sam the contrast matrix

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

#------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#General Function:
#julia <- function(x, y){ # function declaration - function name, parameters etc...
#return =  (x^2 + y^2) # function definition - body of the function
#}
#a=julia(1,2) #function call
#a #5

#square<-function(a){
#  square=a*a
#  return(square)
#}

#square(3)

#------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#Function to do DE

#' This functions does all the steps to use limma to fit an eBayes model.
#' 
#' @param mat (matrix): matrix of gene expression w/ genes as rows
#' @param design: design matrix from model.matrix
#' @param cnt_mat: contrast matrix from MakeContrasts
#' @param trend (logical): flag indicating if mean-variance trend should be accounted for in eBayes (limma-trend method)
#' 
#' @return differential expression data.frame
#'
#get_model = function(mat, design, cnt_mat, robust = F, trend = F) {
# require(limma)
  
# fit = lmFit(mat, design)
# fit_contrast = contrasts.fit(fit, cnt_mat)
# if(robust) {
#   fit_ebayes = eBayes(fit_contrast, trend = T, robust = T)
# } else {
#   fit_ebayes = eBayes(fit_contrast, trend = trend, robust = F)	
# }
  
# return(fit_ebayes)
#}


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


#=========================================================================================#

# Fitting models. The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)
#de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)
#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
#tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median)
#tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
#head(tmp$Accession)
#de$tt = tmp

#head(de$tt)

#=========================================================================================#

#Went deeply through sample code to understand DEA better
#--------------------------------------------------------------#

## Fitting models. The trend=TRUE argument indicates that the mean-variance trend 
##will be accounted for in eBayes (limma-trend method)
#de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)

##In proteomics, the limma-trend method is often preferred because it accounts for 
##the fact that low-intensity proteins usually have more "noise" than high-intensity ones.

#head(de)
#head(de$tt$Symbol)

#table(grepl("\\.", de$tt$Symbol))
##grepl(): Short for "grep logical." It returns a vector of TRUE or FALSE values 
##rather than indices, which is perfect for filtering data frames.
##FALSE  TRUE 
##4743    92 
##92 gene symbols have "."; eg: "Uchl5"   "Uchl5.2

#head(de$tt$Symbol[grepl("\\.", de$tt$Symbol)])

#grep("^Acly", de$tt$Symbol, value = TRUE, ignore.case = TRUE)
##The ^ ensures it starts with Acly (so you don't accidentally get other genes containing those letters).
##ignore.case = TRUE handles "ACLY", "Acly", or "acly"
#grep("^Pcid2", de$tt$Symbol, value = TRUE, ignore.case = TRUE)
#grep("^Uchl5", de$tt$Symbol, value = TRUE, ignore.case = TRUE)
##[1] "Uchl5"   "Uchl5.2"

#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)

##Uchl5.2": The function splits it at the dot into ["Uchl5", "2"]. It then grabs the 1st piece: "Uchl5".

##In proteomics datasets, the difference between Uchl5 and Uchl5.2 (or any gene with a .2, .3 suffix)
##typically indicates that multiple unique protein sequences (Accession IDs)
##mapped to the same single Gene Symbol. 

## Sometimes gene symbols come with extra version numbers or suffixes (like Gene.1). 
##gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
##uses your extract_string function to chop off everything after the first period (.), 
##leaving just the clean gene name (e.g., "OGT").

#head(gene.symbol)
#grep("^Uchl5", gene.symbol, value = TRUE, ignore.case = TRUE)
##[1] "Uchl5"   "Uchl5"; so '.' removed

#colnames(de$tt)
#any(duplicated(de$tt$Symbol))
#tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median) #removing "Accession" and "Symbol" 
#tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
#head(tmp$Accession)
#head(de$tt)
#head(tmp)
#de$tt = tmp

#head(de$tt)
#grep("^Uchl5", de$tt$Symbol, value = TRUE, ignore.case = TRUE)
#de$tt[grep("Uchl5", de$tt$Symbol, ignore.case = TRUE), ]

#=========================================================================================#

#=========================================================================================#

# Fitting models. The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)
de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)
gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median)
tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
head(tmp$Accession)
de$tt = tmp

head(de$tt)

#=========================================================================================#

head(mapping11)
dim(mapping11) #9199  2
head(dat1)
dim(dat1) #9199  18 (Median Normalized)
Xa
cm_data
#de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)
dea = top_table_results(mapping11, dat1, Xa, cm_data, annotate = TRUE, trend = TRUE)
head(dea)
head(dea$tt) #24 columns #log FC 
ncol(dea$tt) #[1] 24
nrow(dea$tt) #[1] 9199
nrow(dat1) #[1] 9199
head(dea$tt$Symbol[grepl("\\.", dea$tt$Symbol)])
dim(dea$tt) #[1] 9199   24
table(grepl("\\.", dea$tt$Symbol))
#FALSE 
#9199 

#For old way of 8330 with aggressive duplicate symbol Filter 
#----------------------------------------------------------------#

#grepl(): Short for "grep logical." It returns a vector of TRUE or FALSE values 
#rather than indices, which is perfect for filtering data frames.
#FALSE  TRUE 
# 8333     1
#1 gene symbols have "."; eg: "Uchl5"   "Uchl5.2
#which(grepl("\\.", dea$tt$Symbol))
#[1] 910
#dea$tt$Symbol[910]
#dea$tt$Symbol[911]
##[1] "ap1s2; Ap1s2; AP1S2; ap1s2.L; LOC101341861; LOC111157096"
#dea$tt[910, ]
##GFPTEN.GFPCTRL OGTTEN.OGTCTRL ROGTTEN.ROGTCTRL OGTCTRL.GFPCTRL OGTTEN.GFPTEN ROGTCTRL.GFPCTRL ROGTTEN.GFPTEN
##910      0.1035608         1.0138      -0.02977226      0.02104583     0.9312853         1.722384       1.589051
##AveExpr P.Value.GFPTEN.GFPCTRL adj.P.Val.GFPTEN.GFPCTRL P.Value.OGTTEN.OGTCTRL adj.P.Val.OGTTEN.OGTCTRL
##910 9.462007              0.7116239                0.8991319            0.002167602               0.02122121
##P.Value.ROGTTEN.ROGTCTRL adj.P.Val.ROGTTEN.ROGTCTRL P.Value.OGTCTRL.GFPCTRL adj.P.Val.OGTCTRL.GFPCTRL
##910                 0.915183                          1               0.9399809                 0.9727805
##P.Value.OGTTEN.GFPTEN adj.P.Val.OGTTEN.GFPTEN P.Value.ROGTCTRL.GFPCTRL adj.P.Val.ROGTCTRL.GFPCTRL
##910           0.004020953              0.01541427             1.460305e-05               0.0001330074
##P.Value.ROGTTEN.GFPTEN adj.P.Val.ROGTTEN.GFPTEN Accession
##910           3.518695e-05             0.0002457465    F6SFB5
##Symbol
##910 ap1s2; Ap1s2; AP1S2; ap1s2.L; LOC101341861; LOC111157096

##dea$tt$Symbol[910] <- "Ap1s2"

grep("^ap1s2", dea$tt$Symbol, value = TRUE, ignore.case = TRUE)
#The ^ ensures it starts with Acly (so you don't accidentally get other genes containing those letters).
#ignore.case = TRUE handles "ACLY", "Acly", or "acly"

dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]
#There are 2 rows with apls2; They were Case sensitive in old way
#But we need to aggregate duplicates Symbols after DEA

#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
##Uchl5.2": The function splits it at the dot into ["Uchl5", "2"]. It then grabs the 1st piece: "Uchl5".

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
length(clean_symbols)
# 1. Clean the symbols - Take ONLY the first symbol before any semicolon or dot
#AND force them all to Uppercase
# This makes "ap1s2" and "AP1S2" both become "AP1S2"
clean_symbols_up = toupper(extract_string(dea$tt$Symbol, "[;.]", 1))

# 2. Check the count now—it should be LOWER than 8334
length(unique(clean_symbols)) #8477
length(unique(clean_symbols_up)) #8477
length(clean_symbols) #9199
length(clean_symbols_up) #9199
head(clean_symbols_up)
#head(gene.symbol)
grep("^ap1s2", clean_symbols, value = TRUE, ignore.case = TRUE)
grep("^ap1s2", clean_symbols_up, value = TRUE, ignore.case = TRUE)
grep("^hand2", clean_symbols_up, value = TRUE, ignore.case = TRUE)
grep("^C10orf88", clean_symbols_up, value = TRUE, ignore.case = TRUE)
grep("^EFNA5", clean_symbols_up, value = TRUE, ignore.case = TRUE)

colnames(dea$tt)
#tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median)
tmpa = aggregate(dea$tt[,-c(23,24)], by = list(Symbol = clean_symbols_up), FUN = median) #But isoforms are median aggregated

#The code calculates the median value for each protein (represented by clean_symbols_up)
#across all numeric columns in the results table, effectively collapsing multiple rows
#belonging to the same protein into a single representative row. 

#where, clean_symbols_up = toupper(extract_string(dea$tt$Symbol, "[;.]", 1))

#ISoforms are median aggregated but we will fix this later


#removing old symbols, accession
#adding the new symbols (gene.symbola) with say ABC.1 , ABC.2 extracted to ABC, ABC
#Finding median of ABC
#Be Careful using median as it can lose P value significance But in this case we had only one row with dot

dea$tt$Symbol[910]
tmpa$Symbol[910]
#Your new tmpa is now sorted alphabetically by Symbol (A-Z).
head(dea$tt)
length(dea$tt$Symbol) #9199
head(tmpa) #Has no Accession
length(dea$tt$Symbol) #9199
length(tmpa$Symbol) #8477
9199-8477 #722 probably aggregated by median 
head(clean_symbols_up)
tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]

#tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
tmpa$Accession = dea$tt$Accession[match(tmpa$Symbol, clean_symbols_up)] 
head(tmpa)
tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]

#Explanation:
#--------------#
#clean_symbols_up = toupper(extract_string(dea$tt$Symbol, "[;.]", 1))
#Clean symbols up has order of dea$tt$symbol
#match(tmpa$Symbol, clean_symbols_up) = gives the position
#WE have two tables that are out of order:
#  tmpa: Sorted alphabetically.
#dea$tt: Sorted by P-value (Significant genes first).
#When you run: match(tmpa$Symbol, clean_symbols_up)
#R is saying:
#  "Hey, tmpa$Symbol says the first gene is A2M. 
#Go look through the original clean_symbols_up list and tell me which row A2M was originally in.
#If A2M was originally in row 500, the match function returns the number 500.
#Then, tmpa$Accession = dea$tt$Accession[500]
#"Go to the original table dea$tt$Accession, grab the Accession ID from row 500, 
#and paste it into the first row of my new table tmpa"
#match(A,B) takes a list of things you have (Table A) and finds 
#where they are located in a master list (Table B).

#head(tmp$Accession)
head(tmpa$Accession)
head(tmpa)
tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]
head(dea$tt)
head(tmpa)
nrow(tmpa) #[1] 8477; old=8330
dim(dea$tt) #[1] 9199; old=8334   24
dea$tt[grep("ap1s2", dea$tt$Symbol, ignore.case = TRUE), ]
tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
#de$tt = tmp
dea$tt2 = dea$tt
dea$tt = tmpa
#head(de$tt)
head(dea$tt2)
head(tmpa)
head(dea$tt)

dim(dea$tt2) #[1] 9199   24
dim(tmpa) #[1] 8477   24
dim(dea$tt) #[1] 8477   24

head(dea)


##Check

tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
tmpa[grep("Hand2; HAND2", tmpa$Symbol, ignore.case = TRUE), ]
tmpa[grep("HAND2", tmpa$Symbol, ignore.case = TRUE), ]
tmpa[grep("C10orf88", tmpa$Symbol, ignore.case = TRUE), ]
tmpa[grep("Efna5; EFNA5", tmpa$Symbol, ignore.case = TRUE), ]

# Extract rows where the symbol contains at least one lowercase letter
mixed_case_rows <- tmpa[grep("[a-z]", tmpa$Symbol, ignore.case = FALSE), ]

# View the matching symbols
unique(mixed_case_rows$Symbol)

#This is better
#Because ; removed form Symbols


#----------------------------------------------------------------------------------------------------------#
#Original 8334 table's apls2 - why we chose clean_symbols_up for the original 8334 way

#Symbol GFPTEN.GFPCTRL OGTTEN.OGTCTRL ROGTTEN.ROGTCTRL OGTCTRL.GFPCTRL OGTTEN.GFPTEN
#370  ap1s2      0.1035608         1.0138      -0.02977226      0.02104583     0.9312853
#371  AP1S2      0.1035608         1.0138      -0.02977226      0.02104583     0.9312853
#ROGTCTRL.GFPCTRL ROGTTEN.GFPTEN  AveExpr P.Value.GFPTEN.GFPCTRL adj.P.Val.GFPTEN.GFPCTRL
#370         1.722384       1.589051 9.462007              0.7116239                0.8991319
#371         1.722384       1.589051 9.462007              0.7116239                0.8991319
#P.Value.OGTTEN.OGTCTRL adj.P.Val.OGTTEN.OGTCTRL P.Value.ROGTTEN.ROGTCTRL
#370            0.002167602               0.02122121                 0.915183
#371            0.002167602               0.02122121                 0.915183
#adj.P.Val.ROGTTEN.ROGTCTRL P.Value.OGTCTRL.GFPCTRL adj.P.Val.OGTCTRL.GFPCTRL
#370                          1               0.9399809                 0.9727805
#371                          1               0.9399809                 0.9727805
#P.Value.OGTTEN.GFPTEN adj.P.Val.OGTTEN.GFPTEN P.Value.ROGTCTRL.GFPCTRL
#370           0.004020953              0.01541427             1.460305e-05
#371           0.004020953              0.01541427             1.460305e-05
#adj.P.Val.ROGTCTRL.GFPCTRL P.Value.ROGTTEN.GFPTEN adj.P.Val.ROGTTEN.GFPTEN Accession
#370               0.0001330074           3.518695e-05             0.0002457465    F6SFB5
#371               0.0001330074           3.518695e-05             0.0002457465    P56377

#----------------------------------------------------------------------------------------------------------#

#===========================================================================================================#

# Checking why clean_symbols_up could be better than gene.symbola
#--------------------------------------------------------------------#
# Fitting models. The trend=TRUE argument indicates that the mean-variance trend will be accounted for in eBayes (limma-trend method)
#de = top_table_results(mapping, dat, X, cm_dat, annotate = TRUE, trend = TRUE)
#gene.symbol = extract_string(de$tt$Symbol, "\\.", 1)
#tmp = aggregate(de$tt[,-c(11,12)], by = list(Symbol = gene.symbol), FUN = median)
#tmp$Accession = de$tt$Accession[match(tmp$Symbol, gene.symbol)] # Get the first listed Accession ID for each gene symbol
#head(tmp$Accession)
#de$tt = tmp

#head(de$tt)

#with gene.symbol directly instead of clean_symbol_Up
#--------------------------------------------------------#

#head(mapping11)
#head(dat1)
#Xa
#cm_data
#dea = top_table_results(mapping11, dat1, Xa, cm_data, annotate = TRUE, trend = TRUE)
#head(dea)
#head(dea$tt) #24 columns #log FC 
#ncol(dea$tt) #24
#gene.symbola = extract_string(dea$tt$Symbol, "\\.", 1)
#head(gene.symbola)
#head(dea$tt$Symbol)
# both match
#tmpa = aggregate(dea$tt[,-c(23,24)], by = list(Symbol = gene.symbola), FUN = median)
#head(dea$tt)
#head(tmpa)
#tmpa$Accession = dea$tt$Accession[match(tmpa$Symbol, gene.symbola)] # Get the first listed Accession ID for each gene symbol
#head(tmpa)
#head(dea$tt)
#dea$tt = tmpa
#head(dea$tt)
#head(tmpa)
#dim(tmpa)

#Check

#tmpa[grep("ap1s2", tmpa$Symbol, ignore.case = TRUE), ]
#tmpa[grep("Hand2; HAND2", tmpa$Symbol, ignore.case = TRUE), ]
#tmpa[grep("HAND2", tmpa$Symbol, ignore.case = TRUE), ]
#tmpa[grep("C10orf88", tmpa$Symbol, ignore.case = TRUE), ]
#tmpa[grep("Efna5; EFNA5", tmpa$Symbol, ignore.case = TRUE), ]

# Extract rows where the symbol contains at least one lowercase letter
#mixed_case_rows <- tmpa[grep("[a-z]", tmpa$Symbol, ignore.case = FALSE), ]

# View the matching symbols
#unique(mixed_case_rows$Symbol)

#===========================================================================================================#

#===========================================================================================================#

#logFCprot
head(dea$tt2)
head(tmpa)
head(dea$tt)
dim(dea$tt2) #[1] 9199   24
dim(tmpa) #[1] 8477   24
dim(dea$tt) #[1] 8477   24
logFCprot<-dea$tt
head(logFCprot) 
dim(logFCprot) #[1] 8477   24
logFCprot<-as.data.frame(logFCprot)
head(logFCprot)
type(logFCprot)
library("writexl") 
#install.packages("writexl")
#write_xlsx(logFCprot,"C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\logFCprot_new_05182026.2.xlsx")
dim(logFCprot) #[1] 8477   24; #old = [1] 8330   24 #Median summarized  isoforms

#=========================================================================================================================================================================================#
#============================================================================================================================================================================================#

#==================================================#
# Volcano plots for p-value and adjusted p-values
#==================================================#

#Final Code
#==================================================#
head(dea$tt) #logFC and pvalues
dim(dea$tt) #logFC and pvalues
colnames(dea$tt) #logFC and pvalues

#==================#
#GFPTEN.GFPCTRL
#==================#

rx <- c(-1, 1)*max(abs(dea$tt$GFPTEN.GFPCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.GFPTEN.GFPCTRL), -log10(dea$tt$adj.P.Val.GFPTEN.GFPCTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#-----------------#
#p-value
#-----------------#

plot(dea$tt$GFPTEN.GFPCTRL, -log10(dea$tt$P.Value.GFPTEN.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10 p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values - GFP 10 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$GFPTEN.GFPCTRL)> lfc & -log(dea$tt$P.Value.GFPTEN.GFPCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$GFPTEN.GFPCTRL)< (-lfc) & -log(dea$tt$P.Value.GFPTEN.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes1,]$P.Value.GFPTEN.GFPCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes2,]$P.Value.GFPTEN.GFPCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#--------------------#
#adjusted p value
#--------------------#

plot(dea$tt$GFPTEN.GFPCTRL, -log10(dea$tt$adj.P.Val.GFPTEN.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values - GFP 10 min vs. GFP 0 min")

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

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------------------------------------#
#============================================================#
#GFPTEN.GFPCTRL_Total proteomics (For Sample Calculation)
#============================================================#

#rx <- c(-1, 1)*max(abs(dea$tt$GFPTEN.GFPCTRL))*1.1 #log2FC
#ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.GFPTEN.GFPCTRL), -log10(dea$tt$adj.P.Val.GFPTEN.GFPCTRL))))

###https://biostatsquid.com/volcano-plot/

###par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
###par(las=1, xaxs="i", yaxs="i")

#----------#
#p-value
#----------#
#plot(dea$tt$GFPTEN.GFPCTRL, -log10(dea$tt$P.Value.GFPTEN.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
#    xlim=rx, ylim=ry, xaxt="n",
#   xlab="Fold change", ylab="-log10 p-value")
#abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
#abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
#axis(1, seq(-6,6,1))
#title("Volcano plot of p-values - GFP 10 min vs. GFP 0 min")

### Log2 fold change and p-value cutoff
#lfc <- 0.58
###lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
#pval <- 0.05 
#-log10(pval)
#-log10(0.01)
###>1.3 is good
### Selecting interesting genes
#sigGenes1 <- ((dea$tt$GFPTEN.GFPCTRL)> lfc & -log(dea$tt$P.Value.GFPTEN.GFPCTRL,10) > -log10(pval))   
#sigGenes2 <- ((dea$tt$GFPTEN.GFPCTRL)< (-lfc) & -log(dea$tt$P.Value.GFPTEN.GFPCTRL,10) > -log10(pval))   
### Identifying the selected genes
#points(dea$tt[sigGenes1,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes1,]$P.Value.GFPTEN.GFPCTRL,10),pch=20,col="orange",cex=2)
#points(dea$tt[sigGenes2,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes2,]$P.Value.GFPTEN.GFPCTRL,10),pch=20,col="blue",cex=2)
#abline(h=-log10(pval),col="brown4",lty=2)
#abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
#mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
###mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
#mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#--------------------#
#adjusted p value
#--------------------#
#plot(dea$tt$GFPTEN.GFPCTRL, -log10(dea$tt$adj.P.Val.GFPTEN.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
#    xlim=rx, ylim=ry, xaxt="n",
#   xlab="log2 fold change", ylab="-log10  p-value")
#abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
#abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
#axis(1, seq(-6,6,1))
#title("Volcano plot of adjusted p-values - GFP 10 min vs. GFP 0 min")

### Log2 fold change and p-value cutoff
#lfc <- 0.58
###lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
#pval <- 0.05 
### Selecting interesting genes
#sigGenes1 <- ((dea$tt$GFPTEN.GFPCTRL)> lfc & -log(dea$tt$adj.P.Val.GFPTEN.GFPCTRL,10) > -log10(pval))   
#sigGenes2 <- ((dea$tt$GFPTEN.GFPCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.GFPTEN.GFPCTRL,10) > -log10(pval))   
### Identifying the selected genes
#points(dea$tt[sigGenes1,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.GFPTEN.GFPCTRL,10),pch=20,col="orange",cex=2)
#points(dea$tt[sigGenes2,]$GFPTEN.GFPCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.GFPTEN.GFPCTRL,10),pch=20,col="blue",cex=2)
#abline(h=-log10(pval),col="brown4",lty=2)
#abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
#mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
###mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
### Keep lfc at 0.58 for the position, but use "1.5" for the text
#mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)
#-------------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#==================#
#OGTTEN.OGTCTRL
#==================#
rx <- c(-1, 1)*max(abs(dea$tt$OGTTEN.OGTCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTTEN.OGTCTRL), -log10(dea$tt$adj.P.Val.OGTTEN.OGTCTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#----------------------------#
#p-value
#----------------------------#
plot(dea$tt$OGTTEN.OGTCTRL, -log10(dea$tt$P.Value.OGTTEN.OGTCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values - OGT KD 605-10 min vs. OGT KD 605-0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$OGTTEN.OGTCTRL)> lfc & -log(dea$tt$P.Value.OGTTEN.OGTCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$OGTTEN.OGTCTRL)< (-lfc) & -log(dea$tt$P.Value.OGTTEN.OGTCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$OGTTEN.OGTCTRL,-log(dea$tt[sigGenes1,]$P.Value.OGTTEN.OGTCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$OGTTEN.OGTCTRL,-log(dea$tt[sigGenes2,]$P.Value.OGTTEN.OGTCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#----------------------------#
#adjusted p value
#----------------------------#
plot(dea$tt$OGTTEN.OGTCTRL, -log10(dea$tt$adj.P.Val.OGTTEN.OGTCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values - OGT KD 605-10 min vs. OGT KD 605-0 min")

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

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#======================#
#ROGTTEN.ROGTCTRL
#======================#

colnames(dea$tt)
rx <- c(-1, 1)*max(abs(dea$tt$ROGTTEN.ROGTCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.ROGTTEN.ROGTCTRL), -log10(dea$tt$adj.P.Val.ROGTTEN.ROGTCTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#-----------------------------#
#p-value
#-----------------------------#
plot(dea$tt$ROGTTEN.ROGTCTRL, -log10(dea$tt$P.Value.ROGTTEN.ROGTCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values - OGT KD 606-10 min vs. OGT KD 606-0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$ROGTTEN.ROGTCTRL)> lfc & -log(dea$tt$P.Value.ROGTTEN.ROGTCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$ROGTTEN.ROGTCTRL)< (-lfc) & -log(dea$tt$P.Value.ROGTTEN.ROGTCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$ROGTTEN.ROGTCTRL,-log(dea$tt[sigGenes1,]$P.Value.ROGTTEN.ROGTCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$ROGTTEN.ROGTCTRL,-log(dea$tt[sigGenes2,]$P.Value.ROGTTEN.ROGTCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#-----------------------------#
#adjusted p value
#-----------------------------#
plot(dea$tt$ROGTTEN.ROGTCTRL, -log10(dea$tt$adj.P.Val.ROGTTEN.ROGTCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values - OGT KD 606-10 min vs. OGT KD 606-0 min")

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

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#=======================#
#OGTCTRL.GFPCTRL
#=======================#
colnames(dea$tt)
rx <- c(-1, 1)*max(abs(dea$tt$OGTCTRL.GFPCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTCTRL.GFPCTRL), -log10(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#----------------------------#
#p-value
#----------------------------#

plot(dea$tt$OGTCTRL.GFPCTRL, -log10(dea$tt$P.Value.OGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 Fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values - OGT KD 605-0 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$OGTCTRL.GFPCTRL)> lfc & -log(dea$tt$P.Value.OGTCTRL.GFPCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$OGTCTRL.GFPCTRL)< (-lfc) & -log(dea$tt$P.Value.OGTCTRL.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$OGTCTRL.GFPCTRL,-log(dea$tt[sigGenes1,]$P.Value.OGTCTRL.GFPCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$OGTCTRL.GFPCTRL,-log(dea$tt[sigGenes2,]$P.Value.OGTCTRL.GFPCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#----------------------------#
#adjusted p value
#----------------------------#

plot(dea$tt$OGTCTRL.GFPCTRL, -log10(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values - OGT KD 605-0 min vs. GFP 0 min")

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

#-------------------------------------------------------------------------------------------------------------------------------------------#
###Only p value cutoff and no lfc cutoff
#lfc <- 0.58
###lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
#pval <- 0.05 

### Selecting interesting genes excluding logFC
#sigGenes1 <- ((dea$tt$OGTCTRL.GFPCTRL)> lfc & -log(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL,10) > -log10(pval))   
#sigGenes2 <- ((dea$tt$OGTCTRL.GFPCTRL)< (-lfc) & -log(dea$tt$adj.P.Val.OGTCTRL.GFPCTRL,10) > -log10(pval))   
### Identifying the selected genes
#points(dea$tt[sigGenes1,]$OGTCTRL.GFPCTRL,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTCTRL.GFPCTRL,10),pch=20,col="orange",cex=2)
#points(dea$tt[sigGenes2,]$OGTCTRL.GFPCTRL,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTCTRL.GFPCTRL,10),pch=20,col="blue",cex=2)
#abline(h=-log10(pval),col="brown4",lty=2)
###abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
#mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
###mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
### Keep lfc at 0.58 for the position, but use "1.5" for the text
#mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)
#-------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#=======================================#
#OGTTEN.GFPTEN
#=======================================#
colnames(dea$tt)
rx <- c(-1, 1)*max(abs(dea$tt$OGTTEN.GFPTEN))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.OGTTEN.GFPTEN), -log10(dea$tt$adj.P.Val.OGTTEN.GFPTEN))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#---------------#
#p-value
#---------------#
plot(dea$tt$OGTTEN.GFPTEN, -log10(dea$tt$P.Value.OGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values - OGT KD 605-10 min vs. GFP 10 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$OGTTEN.GFPTEN)> lfc & -log(dea$tt$P.Value.OGTTEN.GFPTEN,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$OGTTEN.GFPTEN)< (-lfc) & -log(dea$tt$P.Value.OGTTEN.GFPTEN,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$OGTTEN.GFPTEN,-log(dea$tt[sigGenes1,]$P.Value.OGTTEN.GFPTEN,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$OGTTEN.GFPTEN,-log(dea$tt[sigGenes2,]$P.Value.OGTTEN.GFPTEN,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#---------------------#
#adjusted p value
#---------------------#
plot(dea$tt$OGTTEN.GFPTEN, -log10(dea$tt$adj.P.Val.OGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values - OGT KD 605-10 min vs. GFP 10 min")

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

#--------------------------------------------------------------------------------------------------------------------------#
# only p-value cutoff and no lfc
#lfc <- 0.58
###lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
#pval <- 0.05 
### Selecting interesting genes
#sigGenes1 <- ((dea$tt$OGTTEN.GFPTEN)> lfc & -log(dea$tt$adj.P.Val.OGTTEN.GFPTEN,10) > -log10(pval))   
#sigGenes2 <- ((dea$tt$OGTTEN.GFPTEN)< (-lfc) & -log(dea$tt$adj.P.Val.OGTTEN.GFPTEN,10) > -log10(pval))   
### Identifying the selected genes
#points(dea$tt[sigGenes1,]$OGTTEN.GFPTEN,-log(dea$tt[sigGenes1,]$adj.P.Val.OGTTEN.GFPTEN,10),pch=20,col="orange",cex=2)
#points(dea$tt[sigGenes2,]$OGTTEN.GFPTEN,-log(dea$tt[sigGenes2,]$adj.P.Val.OGTTEN.GFPTEN,10),pch=20,col="blue",cex=2)
#abline(h=-log10(pval),col="brown4",lty=2)
###abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
#mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
###mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
#mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)
#--------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#=========================#
#ROGTCTRL.GFPCTRL
#=========================#
#If your input data was already log2 transformed, 
#then the "Fold Change" you calculated is actually the log2 Fold Change
colnames(dea$tt)
head(dea$tt)
rx <- c(-1, 1)*max(abs(dea$tt$ROGTCTRL.GFPCTRL))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.ROGTCTRL.GFPCTRL), -log10(dea$tt$adj.P.Val.ROGTCTRL.GFPCTRL))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#----------------------#
#p-value
#----------------------#

plot(dea$tt$ROGTCTRL.GFPCTRL, -log10(dea$tt$P.Value.ROGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values - OGT KD 606-0 min vs. GFP 0 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$ROGTCTRL.GFPCTRL)> lfc & -log(dea$tt$P.Value.ROGTCTRL.GFPCTRL,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$ROGTCTRL.GFPCTRL)< (-lfc) & -log(dea$tt$P.Value.ROGTCTRL.GFPCTRL,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$ROGTCTRL.GFPCTRL,-log(dea$tt[sigGenes1,]$P.Value.ROGTCTRL.GFPCTRL,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$ROGTCTRL.GFPCTRL,-log(dea$tt[sigGenes2,]$P.Value.ROGTCTRL.GFPCTRL,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#----------------------#
#adjusted p value
#----------------------#

plot(dea$tt$ROGTCTRL.GFPCTRL, -log10(dea$tt$adj.P.Val.ROGTCTRL.GFPCTRL), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values - OGT KD 606-0 min vs. GFP 0 min")

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

#-------------------------------------------------------------------------------------------------------------------------------------------------#

#================================#
#ROGTTEN.GFPTEN
#================================#
colnames(dea$tt)
rx <- c(-1, 1)*max(abs(dea$tt$ROGTTEN.GFPTEN))*1.1 #log2FC
ry <- c(0, ceiling(max(-log10(dea$tt$P.Value.ROGTTEN.GFPTEN), -log10(dea$tt$adj.P.Val.ROGTTEN.GFPTEN))))

#https://biostatsquid.com/volcano-plot/

#par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
#par(las=1, xaxs="i", yaxs="i")

#---------------------------#
#p-value
#---------------------------#

plot(dea$tt$ROGTTEN.GFPTEN, -log10(dea$tt$P.Value.ROGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="Fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of p-values -  OGT KD 606-10 min vs. GFP 10 min")

# Log2 fold change and p-value cutoff
lfc <- 0.58
#lfc=expression in 10 min/expression in 0 min; 5/5=1; 10/5=2
pval <- 0.05 
-log10(pval)
-log10(0.01)
#>1.3 is good
# Selecting interesting genes
sigGenes1 <- ((dea$tt$ROGTTEN.GFPTEN)> lfc & -log(dea$tt$P.Value.ROGTTEN.GFPTEN,10) > -log10(pval))   
sigGenes2 <- ((dea$tt$ROGTTEN.GFPTEN)< (-lfc) & -log(dea$tt$P.Value.ROGTTEN.GFPTEN,10) > -log10(pval))   
# Identifying the selected genes
points(dea$tt[sigGenes1,]$ROGTTEN.GFPTEN,-log(dea$tt[sigGenes1,]$P.Value.ROGTTEN.GFPTEN,10),pch=20,col="orange",cex=2)
points(dea$tt[sigGenes2,]$ROGTTEN.GFPTEN,-log(dea$tt[sigGenes2,]$P.Value.ROGTTEN.GFPTEN,10),pch=20,col="blue",cex=2)
abline(h=-log10(pval),col="brown4",lty=2)
abline(v=c(-lfc,lfc),col="deeppink2",lty=2)
mtext(paste("pval = ",round(pval,2)), side=2,at=-log10(pval),cex=0.8,line=0.5,las=1)
#mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex=1,line=0.2)
mtext(c("- 1.5 fold", "+ 1.5 fold"), side = 3, at = c(-0.58, 0.58), cex = 1, line = 0.2)

#---------------------------#
#adjusted p value
#---------------------------#

plot(dea$tt$ROGTTEN.GFPTEN, -log10(dea$tt$adj.P.Val.ROGTTEN.GFPTEN), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="log2 fold change", ylab="-log10  p-value")
abline(v=seq(-6,6,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-6,6,1))
title("Total Proteomics_SY5Y - Volcano plot of adjusted p-values -  OGT KD 606-10 min vs. GFP 10 min")

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

#========================================================================================================================================================================================#

#===================================================#
#Prep for GSEA and other Downstream Analyses
#===================================================#

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
dim(dea2$tt) #8477 24 
colnames(dea2$tt)

#-------------#
#For loop
#-------------#

head(dea2$tt)

# Define your 7 exact treatment comparisons
contrasts <- c("OGTTEN.GFPTEN", "OGTCTRL.GFPCTRL", "OGTTEN.OGTCTRL", 
               "ROGTTEN.GFPTEN", "ROGTCTRL.GFPCTRL", "ROGTTEN.ROGTCTRL", "GFPTEN.GFPCTRL")

for (comp in contrasts) {
  
  #comp: This is a temporary "placeholder" variable.
  #The first time the loop runs, comp equals "OGTTEN.GFPTEN". 
  #The second time, comp automatically changes to "OGTCTRL.GFPCTRL", 
  #and so on, until it hits the end of the list.
  

  # 1. Recreate your exact initial 'cbind' step
  # This matches your original: cbind(dea2$tt$COMP, dea2$tt$adj.P.Val.COMP, dea2$tt$Symbol)
  raw_matrix <- cbind(dea2$tt[[comp]], 
                      dea2$tt[[paste0("adj.P.Val.", comp)]], 
                      dea2$tt$Symbol)
  
  # 2. Convert to data frame and set your exact custom column names
  updown_df <- as.data.frame(raw_matrix, stringsAsFactors = FALSE)
  colnames(updown_df) <- c(paste0('logFC.', comp), paste0('adj.P.Val.', comp), 'Symbol')
  
  # 3. Assign your Symbol column to rownames and then drop the Symbol column (-3)
  rownames(updown_df) <- updown_df$Symbol
  updown_df <- updown_df[, -3]
  updown_df <- as.data.frame(updown_df)
  
  # 4. Apply your exact formatC and numeric formatting step for p-values
  p_col_name <- paste0('adj.P.Val.', comp)
  updown_df[[p_col_name]] <- as.numeric(formatC(as.numeric(updown_df[[p_col_name]]), format = "e", digits = 2))
  
  # 5. Apply your exact formatC and numeric formatting step for logFC
  fc_col_name <- paste0('logFC.', comp)
  updown_df[[fc_col_name]] <- as.numeric(formatC(as.numeric(updown_df[[fc_col_name]]), format = "e", digits = 2))
  
  # 6. Sort by decreasing LogFC exactly like your final line of code
  sorted_df <- as.data.frame(updown_df[order(updown_df[[fc_col_name]], decreasing = TRUE), ])
  
  # 7. Create the exact global variable name (e.g., sortupdown.OGTTEN.GFPTEN) - ordered by decreasing log FC
  var_name <- paste0("sortupdown.", comp)
  assign(var_name, sorted_df, envir = .GlobalEnv)
}

#print head and dim of all using for loop

head(sortupdown.GFPTEN.GFPCTRL)
dim(sortupdown.GFPTEN.GFPCTRL)
contrasts <- c("OGTTEN.GFPTEN", "OGTCTRL.GFPCTRL", "OGTTEN.OGTCTRL", 
               "ROGTTEN.GFPTEN", "ROGTCTRL.GFPCTRL", "ROGTTEN.ROGTCTRL", "GFPTEN.GFPCTRL")

#for (a in contrasts) {
#print(head(sortupdown.a))
#  }

#for (a in contrasts) {
#  print(dim(sortupdown.a))
#}

for (a in contrasts) {
  # 1. Construct the variable name string
  var_name <- paste0("sortupdown.", a)
  
  # 2. Print a clean separator header
  cat("\n#========================================#\n")
  cat("Dataset:", var_name, "\n")
  
  # 3. Fetch the data frame using get()
  df_to_print <- get(var_name)
  
  # 4. Print the dimensions (Rows, Columns)
  cat("Dimensions (Rows x Columns):", paste(dim(df_to_print), collapse = " x "), "\n")
  cat("#========================================#\n")
  
  # 5. Print the top 6 rows
  print(head(df_to_print))
}

#msg <- "Analyzing data"
#print(msg)
#cat(msg) # cat = Clean print (cat stands for concatenate and print)

#---------------------------#
#Traditional repeated code
#---------------------------#

#Different treatment comparisons 
updown.OGTTEN.GFPTEN<-cbind(dea2$tt$OGTTEN.GFPTEN,dea2$tt$adj.P.Val.OGTTEN.GFPTEN,dea2$tt$Symbol)           #OGTTEN.GFPTEN
updown.OGTCTRL.GFPCTRL<-cbind(dea2$tt$OGTCTRL.GFPCTRL,dea2$tt$adj.P.Val.OGTCTRL.GFPCTRL,dea2$tt$Symbol)     #OGTCTRL.GFPCTRL
updown.OGTTEN.OGTCTRL<-cbind(dea2$tt$OGTTEN.OGTCTRL,dea2$tt$adj.P.Val.OGTTEN.OGTCTRL,dea2$tt$Symbol)        #OGTTEN.OGTCTRL
updown.ROGTTEN.GFPTEN<-cbind(dea2$tt$ROGTTEN.GFPTEN,dea2$tt$adj.P.Val.ROGTTEN.GFPTEN,dea2$tt$Symbol)        #ROGTTEN.GFPTEN
updown.ROGTCTRL.GFPCTRL<-cbind(dea2$tt$ROGTCTRL.GFPCTRL,dea2$tt$adj.P.Val.ROGTCTRL.GFPCTRL,dea2$tt$Symbol)  #ROGTCTRL.GFPCTRL
updown.ROGTTEN.ROGTCTRL<-cbind(dea2$tt$ROGTTEN.ROGTCTRL,dea2$tt$adj.P.Val.ROGTTEN.ROGTCTRL,dea2$tt$Symbol)  #ROGTTEN.ROGTCTRL
updown.GFPTEN.GFPCTRL<-cbind(dea2$tt$GFPTEN.GFPCTRL,dea2$tt$adj.P.Val.GFPTEN.GFPCTRL,dea2$tt$Symbol)        #GFPTEN.GFPCTRL

#OGTTEN.GFPTEN
#OGTCTRL.GFPCTRL
#OGTTEN.OGTCTRL
#ROGTTEN.GFPTEN
#ROGTCTRL.GFPCTRL
#ROGTTEN.ROGTCTRL
#GFPTEN.GFPCTRL

head(updown.OGTTEN.GFPTEN)
head(updown.OGTCTRL.GFPCTRL)
head(updown.OGTTEN.OGTCTRL)
head(updown.ROGTTEN.GFPTEN)
head(updown.ROGTCTRL.GFPCTRL)
head(updown.ROGTTEN.ROGTCTRL)
head(updown.GFPTEN.GFPCTRL)

#---------------------#
#1. OGTTEN.GFPTEN
#---------------------#

#head(dea2$tt$OGTCTRL.GFPCTRL)
#head(dea2$tt$adj.P.Val.OGTCTRL.GFPCTRL)
head(dea2$tt)

head(updown.OGTTEN.GFPTEN)
updown.OGTTEN.GFPTEN<-as.data.frame(updown.OGTTEN.GFPTEN)
colnames(updown.OGTTEN.GFPTEN)<-c('logFC.OGTTEN.GFPTEN','adj.P.Val.OGTTEN.GFPTEN','Symbol')
head(updown.OGTTEN.GFPTEN)
rownames(updown.OGTTEN.GFPTEN)<-updown.OGTTEN.GFPTEN$Symbol
head(updown.OGTTEN.GFPTEN)
updown.OGTTEN.GFPTEN<-updown.OGTTEN.GFPTEN[,-3]
head(updown.OGTTEN.GFPTEN)
updown.OGTTEN.GFPTEN<-as.data.frame(updown.OGTTEN.GFPTEN)
head(updown.OGTTEN.GFPTEN)
#convert e to nmeric so can sort
#updown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN <- as.numeric(formatC(updown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN, format = "e", digits = 2))
#head(updown.OGTTEN.GFPTEN)
#sortupdown.OGTTEN.GFPTEN<-as.data.frame(updown.OGTTEN.GFPTEN[order(updown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN,decreasing = TRUE),]) 
#head(sortupdown.OGTTEN.GFPTEN) #ordered by decreasing p value
updown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN <- as.numeric(formatC(updown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN, format = "e", digits = 2))
head(updown.OGTTEN.GFPTEN)
type(updown.OGTTEN.GFPTEN)
updown.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN <- as.numeric(formatC(updown.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN, format = "e", digits = 2))
head(updown.OGTTEN.GFPTEN)
type(updown.OGTTEN.GFPTEN)
sortupdown.OGTTEN.GFPTEN<-as.data.frame(updown.OGTTEN.GFPTEN[order(updown.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN,decreasing = TRUE),]) 
head(sortupdown.OGTTEN.GFPTEN) #Ordered by decreasing LogFC
dim(sortupdown.OGTTEN.GFPTEN) 
#[1] 8477    2

#----------#
#for loop
#----------#

#========================================#
#Dataset: sortupdown.OGTTEN.GFPTEN 
#Dimensions (Rows x Columns): 8477 x 2 
#========================================#
#logFC.OGTTEN.GFPTEN adj.P.Val.OGTTEN.GFPTEN
#ATP23                 2.70                0.002280
#RNF5                  2.02                0.009050
#ZNF180                1.93                0.000395
#ZNF658                1.93                0.000259
#KAT6A                 1.92                0.000204
#KAT6B                 1.92                0.000204

#----------#
# long code
#----------#

#logFC.OGTTEN.GFPTEN adj.P.Val.OGTTEN.GFPTEN
#ATP23             2.700101            0.0022789629
#RNF5              2.020678            0.0090486275
#ZNF180            1.934382            0.0003946340
#ZNF658            1.932270            0.0002593913
#KAT6A             1.915693            0.0002041224
#KAT6B             1.915693            0.0002041224

#[1] 8477    2

#Bot match, so use for loop as it is more time efficient

#--------------------#
#2. OGTCTRL.GFPCTRL
#--------------------#

head(updown.OGTCTRL.GFPCTRL)
updown.OGTCTRL.GFPCTRL<-as.data.frame(updown.OGTCTRL.GFPCTRL)
colnames(updown.OGTCTRL.GFPCTRL)<-c('logFC.OGTCTRL.GFPCTRL','adj.P.Val.OGTCTRL.GFPCTRL','Symbol')
head(updown.OGTCTRL.GFPCTRL)
rownames(updown.OGTCTRL.GFPCTRL)<-updown.OGTCTRL.GFPCTRL$Symbol
head(updown.OGTCTRL.GFPCTRL)
updown.OGTCTRL.GFPCTRL<-updown.OGTCTRL.GFPCTRL[,-3]
head(updown.OGTCTRL.GFPCTRL)
updown.OGTCTRL.GFPCTRL<-as.data.frame(updown.OGTCTRL.GFPCTRL)
head(updown.OGTCTRL.GFPCTRL)
#convert e to nmeric so can sort
#updown.OGTCTRL.GFPCTRL$adj.P.Val.OGTCTRL.GFPCTRL <- as.numeric(formatC(updown.OGTCTRL.GFPCTRL$adj.P.Val.OGTCTRL.GFPCTRL, format = "e", digits = 2))
#head(updown.OGTCTRL.GFPCTRL)
#sortupdown.OGTCTRL.GFPCTRL<-as.data.frame(updown.OGTCTRL.GFPCTRL[order(updown.OGTCTRL.GFPCTRL$adj.P.Val.OGTCTRL.GFPCTRL,decreasing = TRUE),]) 
#head(sortupdown.OGTCTRL.GFPCTRL) #ordered by decreasing p value
updown.OGTCTRL.GFPCTRL$adj.P.Val.OGTCTRL.GFPCTRL <- as.numeric(formatC(updown.OGTCTRL.GFPCTRL$adj.P.Val.OGTCTRL.GFPCTRL, format = "e", digits = 2))
head(updown.OGTCTRL.GFPCTRL)
type(updown.OGTCTRL.GFPCTRL)
updown.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL <- as.numeric(formatC(updown.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL, format = "e", digits = 2))
head(updown.OGTCTRL.GFPCTRL)
type(updown.OGTCTRL.GFPCTRL)
sortupdown.OGTCTRL.GFPCTRL<-as.data.frame(updown.OGTCTRL.GFPCTRL[order(updown.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL,decreasing = TRUE),]) 
head(sortupdown.OGTCTRL.GFPCTRL) #Ordered by decreasing LogFC
dim(sortupdown.OGTCTRL.GFPCTRL) 

# Long code

#logFC.OGTCTRL.GFPCTRL adj.P.Val.OGTCTRL.GFPCTRL
#RNF5                 1.915524               0.033481287
#SHROOM2              1.498665               0.089624513
#ZNF658               1.393838               0.007173833
#ZNF180               1.357728               0.011540206
#ZBTB39               1.232913               0.005232161
#ZNF32                1.232913               0.005232161
#> dim(sortupdown.OGTCTRL.GFPCTRL) 
#[1] 8477    2

# For loop

#========================================#
#Dataset: sortupdown.OGTCTRL.GFPCTRL 
#Dimensions (Rows x Columns): 8477 x 2 
#========================================#
#logFC.OGTCTRL.GFPCTRL adj.P.Val.OGTCTRL.GFPCTRL
#RNF5                     1.92                   0.03350
#SHROOM2                  1.50                   0.08960
#ZNF658                   1.39                   0.00717
#ZNF180                   1.36                   0.01150
#ZBTB39                   1.23                   0.00523
#ZNF32                    1.23                   0.00523

#Both match, so use for loop as it is more time efficient

#-------------------#
#3. OGTTEN.OGTCTRL
#-------------------#

head(updown.OGTTEN.OGTCTRL)
updown.OGTTEN.OGTCTRL<-as.data.frame(updown.OGTTEN.OGTCTRL)
colnames(updown.OGTTEN.OGTCTRL)<-c('logFC.OGTTEN.OGTCTRL','adj.P.Val.OGTTEN.OGTCTRL','Symbol')
head(updown.OGTTEN.OGTCTRL)
rownames(updown.OGTTEN.OGTCTRL)<-updown.OGTTEN.OGTCTRL$Symbol
head(updown.OGTTEN.OGTCTRL)
updown.OGTTEN.OGTCTRL<-updown.OGTTEN.OGTCTRL[,-3]
head(updown.OGTTEN.OGTCTRL)
updown.OGTTEN.OGTCTRL<-as.data.frame(updown.OGTTEN.OGTCTRL)
head(updown.OGTTEN.OGTCTRL)
#convert e to nmeric so can sort
#updown.OGTTEN.OGTCTRL$adj.P.Val.OGTTEN.OGTCTRL<- as.numeric(formatC(updown.OGTTEN.OGTCTRL$adj.P.Val.OGTTEN.OGTCTRL, format = "e", digits = 2))
#head(updown.OGTTEN.OGTCTRL)
#sortupdown.OGTTEN.OGTCTRL<-as.data.frame(updown.OGTTEN.OGTCTRL[order(updown.OGTTEN.OGTCTRL$adj.P.Val.OGTTEN.OGTCTRL,decreasing = TRUE),]) 
#head(sortupdown.OGTTEN.OGTCTRL) #ordered by decreasing p value
updown.OGTTEN.OGTCTRL$adj.P.Val.OGTTEN.OGTCTRL <- as.numeric(formatC(updown.OGTTEN.OGTCTRL$adj.P.Val.OGTTEN.OGTCTRL, format = "e", digits = 2))
head(updown.OGTTEN.OGTCTRL)
type(updown.OGTTEN.OGTCTRL)
updown.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL <- as.numeric(formatC(updown.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL, format = "e", digits = 2))
head(updown.OGTTEN.OGTCTRL)
type(updown.OGTTEN.OGTCTRL)
sortupdown.OGTTEN.OGTCTRL<-as.data.frame(updown.OGTTEN.OGTCTRL[order(updown.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL,decreasing = TRUE),]) 
head(sortupdown.OGTTEN.OGTCTRL) #Ordered by decreasing LogFC

#------------------#
#4. ROGTTEN.GFPTEN
#------------------#

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
head(updown.ROGTTEN.GFPTEN)
head(dea2$tt)
head(updown.ROGTTEN.GFPTEN)
updown.ROGTTEN.GFPTEN<-as.data.frame(updown.ROGTTEN.GFPTEN)
head(updown.ROGTTEN.GFPTEN)
colnames(updown.ROGTTEN.GFPTEN)<-c('logFC.ROGTTEN.GFPTEN','adj.P.Val.ROGTTEN.GFPTEN','Symbol')
head(updown.ROGTTEN.GFPTEN)
rownames(updown.ROGTTEN.GFPTEN)<-updown.ROGTTEN.GFPTEN$Symbol
head(updown.ROGTTEN.GFPTEN)
updown.ROGTTEN.GFPTEN<-updown.ROGTTEN.GFPTEN[,-3]
head(updown.ROGTTEN.GFPTEN)
updown.ROGTTEN.GFPTEN<-as.data.frame(updown.ROGTTEN.GFPTEN)
head(updown.ROGTTEN.GFPTEN)
#convert e to nmeric so can sort
#updown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN<- as.numeric(formatC(updown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN, format = "e", digits = 2))
#head(updown.ROGTTEN.GFPTEN)
#sortupdown.ROGTTEN.GFPTEN<-as.data.frame(updown.ROGTTEN.GFPTEN[order(updown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN,decreasing = TRUE),]) 
#head(sortupdown.ROGTTEN.GFPTEN) #ordered by decreasing p value
updown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN <- as.numeric(formatC(updown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN, format = "e", digits = 2))
head(updown.ROGTTEN.GFPTEN)
type(updown.ROGTTEN.GFPTEN)
updown.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN <- as.numeric(formatC(updown.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN, format = "e", digits = 2))
head(updown.ROGTTEN.GFPTEN)
type(updown.ROGTTEN.GFPTEN)
sortupdown.ROGTTEN.GFPTEN<-as.data.frame(updown.ROGTTEN.GFPTEN[order(updown.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN,decreasing = TRUE),]) 
head(sortupdown.ROGTTEN.GFPTEN) #Ordered by decreasing LogFC

#----------------------#
#5. ROGTCTRL.GFPCTRL
#----------------------#

head(updown.ROGTCTRL.GFPCTRL)
updown.ROGTCTRL.GFPCTRL<-as.data.frame(updown.ROGTCTRL.GFPCTRL)
colnames(updown.ROGTCTRL.GFPCTRL)<-c('logFC.ROGTCTRL.GFPCTRL','adj.P.Val.ROGTCTRL.GFPCTRL','Symbol')
head(updown.ROGTCTRL.GFPCTRL)
rownames(updown.ROGTCTRL.GFPCTRL)<-updown.ROGTCTRL.GFPCTRL$Symbol
head(updown.ROGTCTRL.GFPCTRL)
updown.ROGTCTRL.GFPCTRL<-updown.ROGTCTRL.GFPCTRL[,-3]
head(updown.ROGTCTRL.GFPCTRL)
updown.ROGTCTRL.GFPCTRL<-as.data.frame(updown.ROGTCTRL.GFPCTRL)
head(updown.ROGTCTRL.GFPCTRL)
#convert e to nmeric so can sort
#updown.ROGTCTRL.GFPCTRL$adj.P.Val.ROGTCTRL.GFPCTRL<- as.numeric(formatC(updown.ROGTCTRL.GFPCTRL$adj.P.Val.ROGTCTRL.GFPCTRL, format = "e", digits = 2))
#head(updown.ROGTCTRL.GFPCTRL)
#sortupdown.ROGTCTRL.GFPCTRL<-as.data.frame(updown.ROGTCTRL.GFPCTRL[order(updown.ROGTCTRL.GFPCTRL$adj.P.Val.ROGTCTRL.GFPCTRL,decreasing = TRUE),]) 
#head(sortupdown.ROGTCTRL.GFPCTRL) #ordered by decreasing p value
updown.ROGTCTRL.GFPCTRL$adj.P.Val.ROGTCTRL.GFPCTRL <- as.numeric(formatC(updown.ROGTCTRL.GFPCTRL$adj.P.Val.ROGTCTRL.GFPCTRL, format = "e", digits = 2))
head(updown.ROGTCTRL.GFPCTRL)
type(updown.ROGTCTRL.GFPCTRL)
updown.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL <- as.numeric(formatC(updown.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL, format = "e", digits = 2))
head(updown.ROGTCTRL.GFPCTRL)
type(updown.ROGTCTRL.GFPCTRL)
sortupdown.ROGTCTRL.GFPCTRL<-as.data.frame(updown.ROGTCTRL.GFPCTRL[order(updown.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL,decreasing = TRUE),]) 
head(sortupdown.ROGTCTRL.GFPCTRL) #Ordered by decreasing LogFC

#----------------------#
#6. ROGTTEN.ROGTCTRL
#----------------------#

head(updown.ROGTTEN.ROGTCTRL)
updown.ROGTTEN.ROGTCTRL<-as.data.frame(updown.ROGTTEN.ROGTCTRL)
colnames(updown.ROGTTEN.ROGTCTRL)<-c('logFC.ROGTTEN.ROGTCTRL','adj.P.Val.ROGTTEN.ROGTCTRL','Symbol')
head(updown.ROGTTEN.ROGTCTRL)
rownames(updown.ROGTTEN.ROGTCTRL)<-updown.ROGTTEN.ROGTCTRL$Symbol
head(updown.ROGTTEN.ROGTCTRL)
updown.ROGTTEN.ROGTCTRL<-updown.ROGTTEN.ROGTCTRL[,-3]
head(updown.ROGTTEN.ROGTCTRL)
updown.ROGTTEN.ROGTCTRL<-as.data.frame(updown.ROGTTEN.ROGTCTRL)
head(updown.ROGTTEN.ROGTCTRL)
#convert e to nmeric so can sort
#updown.ROGTTEN.ROGTCTRL$adj.P.Val.ROGTTEN.ROGTCTRL<- as.numeric(formatC(updown.ROGTTEN.ROGTCTRL$adj.P.Val.ROGTTEN.ROGTCTRL, format = "e", digits = 2))
#head(updown.ROGTTEN.ROGTCTRL)
#sortupdown.ROGTTEN.ROGTCTRL<-as.data.frame(updown.ROGTTEN.ROGTCTRL[order(updown.ROGTTEN.ROGTCTRL$adj.P.Val.ROGTTEN.ROGTCTRL,decreasing = TRUE),]) 
#head(sortupdown.ROGTTEN.ROGTCTRL) #ordered by decreasing p value
updown.ROGTTEN.ROGTCTRL$adj.P.Val.ROGTTEN.ROGTCTRL <- as.numeric(formatC(updown.ROGTTEN.ROGTCTRL$adj.P.Val.ROGTTEN.ROGTCTRL, format = "e", digits = 2))
head(updown.ROGTTEN.ROGTCTRL)
type(updown.ROGTTEN.ROGTCTRL)
updown.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL <- as.numeric(formatC(updown.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL, format = "e", digits = 2))
head(updown.ROGTTEN.ROGTCTRL)
type(updown.ROGTTEN.ROGTCTRL)
sortupdown.ROGTTEN.ROGTCTRL<-as.data.frame(updown.ROGTTEN.ROGTCTRL[order(updown.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL,decreasing = TRUE),]) 
head(sortupdown.ROGTTEN.ROGTCTRL) #Ordered by decreasing LogFC

#----------------------#
#7. GFPTEN.GFPCTRL
#----------------------#

head(updown.GFPTEN.GFPCTRL)
updown.GFPTEN.GFPCTRL<-as.data.frame(updown.GFPTEN.GFPCTRL)
colnames(updown.GFPTEN.GFPCTRL)<-c('logFC.GFPTEN.GFPCTRL','adj.P.Val.GFPTEN.GFPCTRL','Symbol')
head(updown.GFPTEN.GFPCTRL)
rownames(updown.GFPTEN.GFPCTRL)<-updown.GFPTEN.GFPCTRL$Symbol
head(updown.GFPTEN.GFPCTRL)
updown.GFPTEN.GFPCTRL<-updown.GFPTEN.GFPCTRL[,-3]
head(updown.GFPTEN.GFPCTRL)
updown.GFPTEN.GFPCTRL<-as.data.frame(updown.GFPTEN.GFPCTRL)
head(updown.GFPTEN.GFPCTRL)
#convert e to nmeric so can sort
#updown.GFPTEN.GFPCTRL$adj.P.Val.GFPTEN.GFPCTRL<- as.numeric(formatC(updown.GFPTEN.GFPCTRL$adj.P.Val.GFPTEN.GFPCTRL, format = "e", digits = 2))
#head(updown.GFPTEN.GFPCTRL)
#sortupdown.GFPTEN.GFPCTRL<-as.data.frame(updown.GFPTEN.GFPCTRL[order(updown.GFPTEN.GFPCTRL$adj.P.Val.GFPTEN.GFPCTRL,decreasing = TRUE),]) 
#head(sortupdown.GFPTEN.GFPCTRL) #ordered by decreasing p value
updown.GFPTEN.GFPCTRL$adj.P.Val.GFPTEN.GFPCTRL <- as.numeric(formatC(updown.GFPTEN.GFPCTRL$adj.P.Val.GFPTEN.GFPCTRL, format = "e", digits = 2))
head(updown.GFPTEN.GFPCTRL)
type(updown.GFPTEN.GFPCTRL)
updown.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL <- as.numeric(formatC(updown.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL, format = "e", digits = 2))
head(updown.GFPTEN.GFPCTRL)
type(updown.GFPTEN.GFPCTRL)
sortupdown.GFPTEN.GFPCTRL<-as.data.frame(updown.GFPTEN.GFPCTRL[order(updown.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL,decreasing = TRUE),]) 
head(sortupdown.GFPTEN.GFPCTRL) #Ordered by decreasing LogFC

#===================================================================================================================================================#

#NOT NEEDED FOR GSEA BUT CAN DO FOR KEGG AND GO

#===========================================#
#Remove genes not significant p>0.05 
#===========================================#

#-------------#
#For loop
#-------------#

#----------------#
#OGTTEN.GFPTEN
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
head(sortupdown.OGTTEN.GFPTEN)
dim(sortupdown.OGTTEN.GFPTEN)
remsortupdown.OGTTEN.GFPTEN <- sortupdown.OGTTEN.GFPTEN[!(sortupdown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN>0.044), ]
head(remsortupdown.OGTTEN.GFPTEN)
remsortupdown.OGTTEN.GFPTEN<-as.data.frame(remsortupdown.OGTTEN.GFPTEN)
head(remsortupdown.OGTTEN.GFPTEN)
dim(remsortupdown.OGTTEN.GFPTEN)
#3265 proteins after < p=0.05 and FDR Filtration 

head(dea2$tt)

# Define your 7 exact treatment comparisons
datafile <- c("sortupdown.OGTTEN.GFPTEN", "sortupdown.OGTCTRL.GFPCTRL", "sortupdown.OGTTEN.OGTCTRL", 
               "sortupdown.ROGTTEN.GFPTEN", "sortupdown.ROGTCTRL.GFPCTRL", "sortupdown.ROGTTEN.ROGTCTRL", "sortupdown.GFPTEN.GFPCTRL")

for (b in datafile) {
  
  #b: This is a temporary "placeholder" variable.
  #The first time the loop runs, comp equals "OGTTEN.GFPTEN". 
  #The second time, comp automatically changes to "OGTCTRL.GFPCTRL", 
  #and so on, until it hits the end of the list.
  
  # 1. Construct the variable name string
  var_name <- paste0("rem", b)
  
  var_name <- b[!(b$adj.P.Val.b>0.044), ]
   
}


datafile <- c("sortupdown.OGTTEN.GFPTEN", "sortupdown.OGTCTRL.GFPCTRL", "sortupdown.OGTTEN.OGTCTRL", 
              "sortupdown.ROGTTEN.GFPTEN", "sortupdown.ROGTCTRL.GFPCTRL", "sortupdown.ROGTTEN.ROGTCTRL", "sortupdown.GFPTEN.GFPCTRL")

for (b in datafile) {
  
  # 1. Fetch the actual data frame from memory using get()
  current_df <- get(b)
  
  # 2. Extract the contrast suffix (e.g., "OGTTEN.GFPTEN") to find the right column name
  suffix <- sub("sortupdown.", "", b)
  p_val_col <- paste0("adj.P.Val.", suffix)
  
  # 3. Filter the data frame using your exact cutoff (keeping values <= 0.044)
  filtered_df <- current_df[!(current_df[[p_val_col]] > 0.044), ]
  filtered_df <- as.data.frame(filtered_df)
  
  # 4. Create the new variable name (e.g., "remsortupdown.OGTTEN.GFPTEN")
  new_var_name <- paste0("rem", b)
  
  # 5. Save it back to your R workspace
  assign(new_var_name, filtered_df, envir = .GlobalEnv)
}







#print head and dim of all using for loop

head(sortupdown.GFPTEN.GFPCTRL)
dim(sortupdown.GFPTEN.GFPCTRL)
contrasts <- c("OGTTEN.GFPTEN", "OGTCTRL.GFPCTRL", "OGTTEN.OGTCTRL", 
               "ROGTTEN.GFPTEN", "ROGTCTRL.GFPCTRL", "ROGTTEN.ROGTCTRL", "GFPTEN.GFPCTRL")

#for (a in contrasts) {
#print(head(sortupdown.a))
#  }

#for (a in contrasts) {
#  print(dim(sortupdown.a))
#}

for (a in contrasts) {
  # 1. Construct the variable name string
  var_name <- paste0("sortupdown.", a)
  
  # 2. Print a clean separator header
  cat("\n#========================================#\n")
  cat("Dataset:", var_name, "\n")
  
  # 3. Fetch the data frame using get()
  df_to_print <- get(var_name)
  
  # 4. Print the dimensions (Rows, Columns)
  cat("Dimensions (Rows x Columns):", paste(dim(df_to_print), collapse = " x "), "\n")
  cat("#========================================#\n")
  
  # 5. Print the top 6 rows
  print(head(df_to_print))
}

#msg <- "Analyzing data"
#print(msg)
#cat(msg) # cat = Clean print (cat stands for concatenate and print)

#----------------#
#OGTTEN.GFPTEN
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
head(sortupdown.OGTTEN.GFPTEN)
dim(sortupdown.OGTTEN.GFPTEN)
remsortupdown.OGTTEN.GFPTEN <- sortupdown.OGTTEN.GFPTEN[!(sortupdown.OGTTEN.GFPTEN$adj.P.Val.OGTTEN.GFPTEN>0.044), ]
head(remsortupdown.OGTTEN.GFPTEN)
remsortupdown.OGTTEN.GFPTEN<-as.data.frame(remsortupdown.OGTTEN.GFPTEN)
head(remsortupdown.OGTTEN.GFPTEN)
dim(remsortupdown.OGTTEN.GFPTEN)
#3265 proteins after < p=0.05 and FDR Filtration 

#----------------#
#OGTCTRL.GFPCTRL
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
colnames(dea$tt)
nrow(dea$tt)
head(sortupdown.OGTCTRL.GFPCTRL)
dim(sortupdown.OGTCTRL.GFPCTRL)
remsortupdown.OGTCTRL.GFPCTRL <- sortupdown.OGTCTRL.GFPCTRL[!(sortupdown.OGTCTRL.GFPCTRL$adj.P.Val.OGTCTRL.GFPCTRL>0.044), ]
head(remsortupdown.OGTCTRL.GFPCTRL)
remsortupdown.OGTCTRL.GFPCTRL<-as.data.frame(remsortupdown.OGTCTRL.GFPCTRL)
head(remsortupdown.OGTCTRL.GFPCTRL)
dim(remsortupdown.OGTCTRL.GFPCTRL)
#3265 proteins after < p=0.05 and FDR Filtration 

#----------------#
#OGTTEN.OGTCTRL
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
colnames(dea$tt)
nrow(dea$tt)
head(sortupdown.OGTTEN.OGTCTRL)
dim(sortupdown.OGTTEN.OGTCTRL)
remsortupdown.OGTTEN.OGTCTRL <- sortupdown.OGTTEN.OGTCTRL[!(sortupdown.OGTTEN.OGTCTRL$adj.P.Val.OGTTEN.OGTCTRL>0.044), ]
head(remsortupdown.OGTTEN.OGTCTRL)
remsortupdown.OGTTEN.OGTCTRL<-as.data.frame(remsortupdown.OGTTEN.OGTCTRL)
head(remsortupdown.OGTTEN.OGTCTRL)
dim(remsortupdown.OGTTEN.OGTCTRL)
#3265 proteins after < p=0.05 and FDR Filtration 

#--------------------#
#ROGTTEN.GFPTEN
#--------------------#
#OGTKD 606 10 MIN VS GFP 10 MIN
head(sortupdown.ROGTTEN.GFPTEN)
dim(sortupdown.ROGTTEN.GFPTEN)
remsortupdown.ROGTTEN.GFPTEN <- sortupdown.ROGTTEN.GFPTEN[!(sortupdown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN>0.045), ]
head(remsortupdown.ROGTTEN.GFPTEN)
remsortupdown.ROGTTEN.GFPTEN<-as.data.frame(remsortupdown.ROGTTEN.GFPTEN)
head(remsortupdown.ROGTTEN.GFPTEN)
dim(remsortupdown.ROGTTEN.GFPTEN)
#[1] 4636    2 < p=0.05 and FDR Filtration 

#----------------#
#ROGTCTRL.GFPCTRL
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
colnames(dea$tt)
nrow(dea$tt)
head(sortupdown.ROGTCTRL.GFPCTRL)
dim(sortupdown.ROGTCTRL.GFPCTRL)
remsortupdown.ROGTCTRL.GFPCTRL <- sortupdown.ROGTCTRL.GFPCTRL[!(sortupdown.ROGTCTRL.GFPCTRL$adj.P.Val.ROGTCTRL.GFPCTRL>0.044), ]
head(remsortupdown.ROGTCTRL.GFPCTRL)
remsortupdown.ROGTCTRL.GFPCTRL<-as.data.frame(remsortupdown.ROGTCTRL.GFPCTRL)
head(remsortupdown.ROGTCTRL.GFPCTRL)
dim(remsortupdown.ROGTCTRL.GFPCTRL)
#3265 proteins after < p=0.05 and FDR Filtration 

#----------------#
#ROGTTEN.ROGTCTRL
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
colnames(dea$tt)
nrow(dea$tt)
head(sortupdown.ROGTTEN.ROGTCTRL)
dim(sortupdown.ROGTTEN.ROGTCTRL)
remsortupdown.ROGTTEN.ROGTCTRL <- sortupdown.ROGTTEN.ROGTCTRL[!(sortupdown.ROGTTEN.ROGTCTRL$adj.P.Val.ROGTTEN.ROGTCTRL>0.044), ]
head(remsortupdown.ROGTTEN.ROGTCTRL)
remsortupdown.ROGTTEN.ROGTCTRL<-as.data.frame(remsortupdown.ROGTTEN.ROGTCTRL)
head(remsortupdown.ROGTTEN.ROGTCTRL)
dim(remsortupdown.ROGTTEN.ROGTCTRL)
#3265 proteins after < p=0.05 and FDR Filtration 

#----------------#
#GFPTEN.GFPCTRL
#----------------#
#OGTKD 605 10 MIN VS GFP 10 MIN
colnames(dea$tt)
nrow(dea$tt)
head(sortupdown.GFPTEN.GFPCTRL)
dim(sortupdown.GFPTEN.GFPCTRL)
remsortupdown.GFPTEN.GFPCTRL <- sortupdown.GFPTEN.GFPCTRL[!(sortupdown.GFPTEN.GFPCTRL$adj.P.Val.GFPTEN.GFPCTRL>0.044), ]
head(remsortupdown.GFPTEN.GFPCTRL)
remsortupdown.GFPTEN.GFPCTRL<-as.data.frame(remsortupdown.GFPTEN.GFPCTRL)
head(remsortupdown.GFPTEN.GFPCTRL)
dim(remsortupdown.GFPTEN.GFPCTRL)
#3265 proteins after < p=0.05 and FDR Filtration 

#==============================================#
# Filter log FC >1.5 (lfc = 0.58)
#==============================================#

#-------------------------#
#1. OGTTEN.GFPTEN
#-------------------------#

head(remsortupdown.OGTTEN.GFPTEN)
dim(remsortupdown.OGTTEN.GFPTEN)
logFCgreat.OGTTEN.GFPTEN<-remsortupdown.OGTTEN.GFPTEN
head(logFCgreat.OGTTEN.GFPTEN)
logFCgreat.OGTTEN.GFPTEN<-as.data.frame(logFCgreat.OGTTEN.GFPTEN)
#logFCgreat.OGTTEN.GFPTEN<-logFCgreat.OGTTEN.GFPTEN[,-2]
head(logFCgreat.OGTTEN.GFPTEN)
logFCgreat.OGTTEN.GFPTEN<-subset(logFCgreat.OGTTEN.GFPTEN,(logFCgreat.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN)>0.58)
head(logFCgreat.OGTTEN.GFPTEN)
dim(logFCgreat.OGTTEN.GFPTEN)
#351, 2

logFCless.OGTTEN.GFPTEN<-remsortupdown.OGTTEN.GFPTEN
head(logFCless.OGTTEN.GFPTEN)
logFCless.OGTTEN.GFPTEN<-as.data.frame(logFCless.OGTTEN.GFPTEN)
head(logFCless.OGTTEN.GFPTEN)
type(logFCless.OGTTEN.GFPTEN)
#logFCless.OGTTEN.GFPTEN<-as.numeric(logFCless.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN)
#head(logFCless.OGTTEN.GFPTEN)
logFCless.OGTTEN.GFPTEN <- transform(logFCless.OGTTEN.GFPTEN, logFC.OGTTEN.GFPTEN = as.numeric(logFC.OGTTEN.GFPTEN))
head(logFCless.OGTTEN.GFPTEN)
type(logFCless.OGTTEN.GFPTEN)
#logFCless.OGTTEN.GFPTEN<-subset(logFCless.OGTTEN.GFPTEN,(((logFCless.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN)<(-1.5))&((logFCless.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN)>(1.5))))
logFCless.OGTTEN.GFPTEN<-subset(logFCless.OGTTEN.GFPTEN,(((logFCless.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN)<(-0.58))))
head(logFCless.OGTTEN.GFPTEN)
dim(logFCless.OGTTEN.GFPTEN)
#335, 2


logFCgreat.less.OGTTEN.GFPTEN<-remsortupdown.OGTTEN.GFPTEN
head(logFCgreat.less.OGTTEN.GFPTEN)
dim(logFCgreat.less.OGTTEN.GFPTEN)
logFCgreat.less.OGTTEN.GFPTEN<-as.data.frame(logFCgreat.less.OGTTEN.GFPTEN)
head(logFCgreat.less.OGTTEN.GFPTEN)
type(logFCgreat.less.OGTTEN.GFPTEN)
#logFCgreat.less.OGTTEN.GFPTEN<-as.numeric(logFCgreat.less.OGTTEN.GFPTEN$logFCgreat.less.OGTTEN.GFPTEN)
#head(logFCgreat.less.OGTTEN.GFPTEN)
#logFCless.OGTTEN.GFPTEN <- transform(logFCless.OGTTEN.GFPTEN, logFC.OGTTEN.GFPTEN = as.numeric(logFC.OGTTEN.GFPTEN))
logFCgreat.less.OGTTEN.GFPTEN <- transform(logFCgreat.less.OGTTEN.GFPTEN, logFC.OGTTEN.GFPTEN = as.numeric(logFC.OGTTEN.GFPTEN))
head(logFCgreat.less.OGTTEN.GFPTEN)
type(logFCgreat.less.OGTTEN.GFPTEN)
#logFCgreat.less.OGTTEN.GFPTEN<-subset(logFCgreat.less.OGTTEN.GFPTEN,(((logFCgreat.less.OGTTEN.GFPTEN$logFCgreat.less.OGTTEN.GFPTEN)<(-1.5))&((logFCgreat.less.OGTTEN.GFPTEN$logFCgreat.less.OGTTEN.GFPTEN)>(1.5))))
#logFCgreat.less.OGTTEN.GFPTEN<-subset(logFCgreat.less.OGTTEN.GFPTEN,(((logFCgreat.less.OGTTEN.GFPTEN$logFCgreat.less.OGTTEN.GFPTEN)>(1.5))))
logFCgreat.less.OGTTEN.GFPTEN<-as.data.frame(logFCgreat.less.OGTTEN.GFPTEN)
head(logFCgreat.less.OGTTEN.GFPTEN)
type(logFCgreat.less.OGTTEN.GFPTEN)
logFCgreat.less.OGTTEN.GFPTEN<-subset(logFCgreat.less.OGTTEN.GFPTEN,logFCgreat.less.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN>0.58 | logFCgreat.less.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN<(-0.58))
#logFCgreat.OGTTEN.GFPTEN <- logFCgreat.OGTTEN.GFPTEN[logFCgreat.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN > 1.5 | logFCgreat.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN < -1.5, ]
head(logFCgreat.less.OGTTEN.GFPTEN)
dim(logFCgreat.less.OGTTEN.GFPTEN)
#686   2
#351+335
nrow(logFCgreat.OGTTEN.GFPTEN) #351
nrow(logFCless.OGTTEN.GFPTEN)  #335
nrow(logFCgreat.OGTTEN.GFPTEN)+nrow(logFCless.OGTTEN.GFPTEN)

#-------------------------#
#2. OGTCTRL.GFPCTRL
#-------------------------#
colnames(dea$tt)
dim(dea$tt)
head(remsortupdown.OGTCTRL.GFPCTRL)
dim(remsortupdown.OGTCTRL.GFPCTRL)
logFCgreat.OGTCTRL.GFPCTRL<-remsortupdown.OGTCTRL.GFPCTRL
head(logFCgreat.OGTCTRL.GFPCTRL)
logFCgreat.OGTCTRL.GFPCTRL<-as.data.frame(logFCgreat.OGTCTRL.GFPCTRL)
#logFCgreat.OGTCTRL.GFPCTRL<-logFCgreat.OGTCTRL.GFPCTRL[,-2]
head(logFCgreat.OGTCTRL.GFPCTRL)
logFCgreat.OGTCTRL.GFPCTRL<-subset(logFCgreat.OGTCTRL.GFPCTRL,(logFCgreat.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL)>0.58)
head(logFCgreat.OGTCTRL.GFPCTRL)
dim(logFCgreat.OGTCTRL.GFPCTRL)
#351, 2

logFCless.OGTCTRL.GFPCTRL<-remsortupdown.OGTCTRL.GFPCTRL
head(logFCless.OGTCTRL.GFPCTRL)
logFCless.OGTCTRL.GFPCTRL<-as.data.frame(logFCless.OGTCTRL.GFPCTRL)
head(logFCless.OGTCTRL.GFPCTRL)
type(logFCless.OGTCTRL.GFPCTRL)
#logFCless.OGTCTRL.GFPCTRL<-as.numeric(logFCless.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL)
#head(logFCless.OGTCTRL.GFPCTRL)
logFCless.OGTCTRL.GFPCTRL <- transform(logFCless.OGTCTRL.GFPCTRL, logFC.OGTCTRL.GFPCTRL = as.numeric(logFC.OGTCTRL.GFPCTRL))
head(logFCless.OGTCTRL.GFPCTRL)
type(logFCless.OGTCTRL.GFPCTRL)
#logFCless.OGTCTRL.GFPCTRL<-subset(logFCless.OGTCTRL.GFPCTRL,(((logFCless.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL)<(-1.5))&((logFCless.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL)>(1.5))))
logFCless.OGTCTRL.GFPCTRL<-subset(logFCless.OGTCTRL.GFPCTRL,(((logFCless.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL)<(-0.58))))
head(logFCless.OGTCTRL.GFPCTRL)
dim(logFCless.OGTCTRL.GFPCTRL)
#335, 2


logFCgreat.less.OGTCTRL.GFPCTRL<-remsortupdown.OGTCTRL.GFPCTRL
head(logFCgreat.less.OGTCTRL.GFPCTRL)
dim(logFCgreat.less.OGTCTRL.GFPCTRL)
logFCgreat.less.OGTCTRL.GFPCTRL<-as.data.frame(logFCgreat.less.OGTCTRL.GFPCTRL)
head(logFCgreat.less.OGTCTRL.GFPCTRL)
type(logFCgreat.less.OGTCTRL.GFPCTRL)
#logFCgreat.less.OGTCTRL.GFPCTRL<-as.numeric(logFCgreat.less.OGTCTRL.GFPCTRL$logFCgreat.less.OGTCTRL.GFPCTRL)
#head(logFCgreat.less.OGTCTRL.GFPCTRL)
#logFCless.OGTCTRL.GFPCTRL <- transform(logFCless.OGTCTRL.GFPCTRL, logFC.OGTCTRL.GFPCTRL = as.numeric(logFC.OGTCTRL.GFPCTRL))
logFCgreat.less.OGTCTRL.GFPCTRL <- transform(logFCgreat.less.OGTCTRL.GFPCTRL, logFC.OGTCTRL.GFPCTRL = as.numeric(logFC.OGTCTRL.GFPCTRL))
head(logFCgreat.less.OGTCTRL.GFPCTRL)
type(logFCgreat.less.OGTCTRL.GFPCTRL)
#logFCgreat.less.OGTCTRL.GFPCTRL<-subset(logFCgreat.less.OGTCTRL.GFPCTRL,(((logFCgreat.less.OGTCTRL.GFPCTRL$logFCgreat.less.OGTCTRL.GFPCTRL)<(-1.5))&((logFCgreat.less.OGTCTRL.GFPCTRL$logFCgreat.less.OGTCTRL.GFPCTRL)>(1.5))))
#logFCgreat.less.OGTCTRL.GFPCTRL<-subset(logFCgreat.less.OGTCTRL.GFPCTRL,(((logFCgreat.less.OGTCTRL.GFPCTRL$logFCgreat.less.OGTCTRL.GFPCTRL)>(1.5))))
logFCgreat.less.OGTCTRL.GFPCTRL<-as.data.frame(logFCgreat.less.OGTCTRL.GFPCTRL)
head(logFCgreat.less.OGTCTRL.GFPCTRL)
type(logFCgreat.less.OGTCTRL.GFPCTRL)
logFCgreat.less.OGTCTRL.GFPCTRL<-subset(logFCgreat.less.OGTCTRL.GFPCTRL,logFCgreat.less.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL>0.58 | logFCgreat.less.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL<(-0.58))
#logFCgreat.OGTCTRL.GFPCTRL <- logFCgreat.OGTCTRL.GFPCTRL[logFCgreat.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL > 1.5 | logFCgreat.OGTCTRL.GFPCTRL$logFC.OGTCTRL.GFPCTRL < -1.5, ]
head(logFCgreat.less.OGTCTRL.GFPCTRL)
dim(logFCgreat.less.OGTCTRL.GFPCTRL)
#686   2
#351+335
nrow(logFCgreat.OGTCTRL.GFPCTRL) #351
nrow(logFCless.OGTCTRL.GFPCTRL)  #335
nrow(logFCgreat.OGTCTRL.GFPCTRL)+nrow(logFCless.OGTCTRL.GFPCTRL)


#-------------------------#
#3. OGTTEN.OGTCTRL
#-------------------------#

head(remsortupdown.OGTTEN.OGTCTRL)
dim(remsortupdown.OGTTEN.OGTCTRL)
logFCgreat.OGTTEN.OGTCTRL<-remsortupdown.OGTTEN.OGTCTRL
head(logFCgreat.OGTTEN.OGTCTRL)
logFCgreat.OGTTEN.OGTCTRL<-as.data.frame(logFCgreat.OGTTEN.OGTCTRL)
#logFCgreat.OGTTEN.OGTCTRL<-logFCgreat.OGTTEN.OGTCTRL[,-2]
head(logFCgreat.OGTTEN.OGTCTRL)
logFCgreat.OGTTEN.OGTCTRL<-subset(logFCgreat.OGTTEN.OGTCTRL,(logFCgreat.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL)>0.58)
head(logFCgreat.OGTTEN.OGTCTRL)
dim(logFCgreat.OGTTEN.OGTCTRL)
#351, 2

logFCless.OGTTEN.OGTCTRL<-remsortupdown.OGTTEN.OGTCTRL
head(logFCless.OGTTEN.OGTCTRL)
logFCless.OGTTEN.OGTCTRL<-as.data.frame(logFCless.OGTTEN.OGTCTRL)
head(logFCless.OGTTEN.OGTCTRL)
type(logFCless.OGTTEN.OGTCTRL)
#logFCless.OGTTEN.OGTCTRL<-as.numeric(logFCless.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL)
#head(logFCless.OGTTEN.OGTCTRL)
logFCless.OGTTEN.OGTCTRL <- transform(logFCless.OGTTEN.OGTCTRL, logFC.OGTTEN.OGTCTRL = as.numeric(logFC.OGTTEN.OGTCTRL))
head(logFCless.OGTTEN.OGTCTRL)
type(logFCless.OGTTEN.OGTCTRL)
#logFCless.OGTTEN.OGTCTRL<-subset(logFCless.OGTTEN.OGTCTRL,(((logFCless.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL)<(-1.5))&((logFCless.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL)>(1.5))))
logFCless.OGTTEN.OGTCTRL<-subset(logFCless.OGTTEN.OGTCTRL,(((logFCless.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL)<(-0.58))))
head(logFCless.OGTTEN.OGTCTRL)
dim(logFCless.OGTTEN.OGTCTRL)
#335, 2


logFCgreat.less.OGTTEN.OGTCTRL<-remsortupdown.OGTTEN.OGTCTRL
head(logFCgreat.less.OGTTEN.OGTCTRL)
dim(logFCgreat.less.OGTTEN.OGTCTRL)
logFCgreat.less.OGTTEN.OGTCTRL<-as.data.frame(logFCgreat.less.OGTTEN.OGTCTRL)
head(logFCgreat.less.OGTTEN.OGTCTRL)
type(logFCgreat.less.OGTTEN.OGTCTRL)
#logFCgreat.less.OGTTEN.OGTCTRL<-as.numeric(logFCgreat.less.OGTTEN.OGTCTRL$logFCgreat.less.OGTTEN.OGTCTRL)
#head(logFCgreat.less.OGTTEN.OGTCTRL)
#logFCless.OGTTEN.OGTCTRL <- transform(logFCless.OGTTEN.OGTCTRL, logFC.OGTTEN.OGTCTRL = as.numeric(logFC.OGTTEN.OGTCTRL))
logFCgreat.less.OGTTEN.OGTCTRL <- transform(logFCgreat.less.OGTTEN.OGTCTRL, logFC.OGTTEN.OGTCTRL = as.numeric(logFC.OGTTEN.OGTCTRL))
head(logFCgreat.less.OGTTEN.OGTCTRL)
type(logFCgreat.less.OGTTEN.OGTCTRL)
#logFCgreat.less.OGTTEN.OGTCTRL<-subset(logFCgreat.less.OGTTEN.OGTCTRL,(((logFCgreat.less.OGTTEN.OGTCTRL$logFCgreat.less.OGTTEN.OGTCTRL)<(-1.5))&((logFCgreat.less.OGTTEN.OGTCTRL$logFCgreat.less.OGTTEN.OGTCTRL)>(1.5))))
#logFCgreat.less.OGTTEN.OGTCTRL<-subset(logFCgreat.less.OGTTEN.OGTCTRL,(((logFCgreat.less.OGTTEN.OGTCTRL$logFCgreat.less.OGTTEN.OGTCTRL)>(1.5))))
logFCgreat.less.OGTTEN.OGTCTRL<-as.data.frame(logFCgreat.less.OGTTEN.OGTCTRL)
head(logFCgreat.less.OGTTEN.OGTCTRL)
type(logFCgreat.less.OGTTEN.OGTCTRL)
logFCgreat.less.OGTTEN.OGTCTRL<-subset(logFCgreat.less.OGTTEN.OGTCTRL,logFCgreat.less.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL>0.58 | logFCgreat.less.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL<(-0.58))
#logFCgreat.OGTTEN.OGTCTRL <- logFCgreat.OGTTEN.OGTCTRL[logFCgreat.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL > 1.5 | logFCgreat.OGTTEN.OGTCTRL$logFC.OGTTEN.OGTCTRL < -1.5, ]
head(logFCgreat.less.OGTTEN.OGTCTRL)
dim(logFCgreat.less.OGTTEN.OGTCTRL)
#686   2
#351+335
nrow(logFCgreat.OGTTEN.OGTCTRL) #351
nrow(logFCless.OGTTEN.OGTCTRL)  #335
nrow(logFCgreat.OGTTEN.OGTCTRL)+nrow(logFCless.OGTTEN.OGTCTRL)


#-------------------------#
#4. ROGTTEN.GFPTEN
#-------------------------#

head(remsortupdown.ROGTTEN.GFPTEN)
dim(remsortupdown.ROGTTEN.GFPTEN)
logFCgreat.ROGTTEN.GFPTEN<-remsortupdown.ROGTTEN.GFPTEN
head(logFCgreat.ROGTTEN.GFPTEN)
logFCgreat.ROGTTEN.GFPTEN<-as.data.frame(logFCgreat.ROGTTEN.GFPTEN)
#logFCgreat.ROGTTEN.GFPTEN<-logFCgreat.ROGTTEN.GFPTEN[,-2]
head(logFCgreat.ROGTTEN.GFPTEN)
logFCgreat.ROGTTEN.GFPTEN<-subset(logFCgreat.ROGTTEN.GFPTEN,(logFCgreat.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN)>0.58)
head(logFCgreat.ROGTTEN.GFPTEN)
dim(logFCgreat.ROGTTEN.GFPTEN)
#351, 2

logFCless.ROGTTEN.GFPTEN<-remsortupdown.ROGTTEN.GFPTEN
head(logFCless.ROGTTEN.GFPTEN)
logFCless.ROGTTEN.GFPTEN<-as.data.frame(logFCless.ROGTTEN.GFPTEN)
head(logFCless.ROGTTEN.GFPTEN)
type(logFCless.ROGTTEN.GFPTEN)
#logFCless.ROGTTEN.GFPTEN<-as.numeric(logFCless.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN)
#head(logFCless.ROGTTEN.GFPTEN)
logFCless.ROGTTEN.GFPTEN <- transform(logFCless.ROGTTEN.GFPTEN, logFC.ROGTTEN.GFPTEN = as.numeric(logFC.ROGTTEN.GFPTEN))
head(logFCless.ROGTTEN.GFPTEN)
type(logFCless.ROGTTEN.GFPTEN)
#logFCless.ROGTTEN.GFPTEN<-subset(logFCless.ROGTTEN.GFPTEN,(((logFCless.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN)<(-1.5))&((logFCless.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN)>(1.5))))
logFCless.ROGTTEN.GFPTEN<-subset(logFCless.ROGTTEN.GFPTEN,(((logFCless.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN)<(-0.58))))
head(logFCless.ROGTTEN.GFPTEN)
dim(logFCless.ROGTTEN.GFPTEN)
#335, 2


logFCgreat.less.ROGTTEN.GFPTEN<-remsortupdown.ROGTTEN.GFPTEN
head(logFCgreat.less.ROGTTEN.GFPTEN)
dim(logFCgreat.less.ROGTTEN.GFPTEN)
logFCgreat.less.ROGTTEN.GFPTEN<-as.data.frame(logFCgreat.less.ROGTTEN.GFPTEN)
head(logFCgreat.less.ROGTTEN.GFPTEN)
type(logFCgreat.less.ROGTTEN.GFPTEN)
#logFCgreat.less.ROGTTEN.GFPTEN<-as.numeric(logFCgreat.less.ROGTTEN.GFPTEN$logFCgreat.less.ROGTTEN.GFPTEN)
#head(logFCgreat.less.ROGTTEN.GFPTEN)
#logFCless.ROGTTEN.GFPTEN <- transform(logFCless.ROGTTEN.GFPTEN, logFC.ROGTTEN.GFPTEN = as.numeric(logFC.ROGTTEN.GFPTEN))
logFCgreat.less.ROGTTEN.GFPTEN <- transform(logFCgreat.less.ROGTTEN.GFPTEN, logFC.ROGTTEN.GFPTEN = as.numeric(logFC.ROGTTEN.GFPTEN))
head(logFCgreat.less.ROGTTEN.GFPTEN)
type(logFCgreat.less.ROGTTEN.GFPTEN)
#logFCgreat.less.ROGTTEN.GFPTEN<-subset(logFCgreat.less.ROGTTEN.GFPTEN,(((logFCgreat.less.ROGTTEN.GFPTEN$logFCgreat.less.ROGTTEN.GFPTEN)<(-1.5))&((logFCgreat.less.ROGTTEN.GFPTEN$logFCgreat.less.ROGTTEN.GFPTEN)>(1.5))))
#logFCgreat.less.ROGTTEN.GFPTEN<-subset(logFCgreat.less.ROGTTEN.GFPTEN,(((logFCgreat.less.ROGTTEN.GFPTEN$logFCgreat.less.ROGTTEN.GFPTEN)>(1.5))))
logFCgreat.less.ROGTTEN.GFPTEN<-as.data.frame(logFCgreat.less.ROGTTEN.GFPTEN)
head(logFCgreat.less.ROGTTEN.GFPTEN)
type(logFCgreat.less.ROGTTEN.GFPTEN)
logFCgreat.less.ROGTTEN.GFPTEN<-subset(logFCgreat.less.ROGTTEN.GFPTEN,logFCgreat.less.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN>0.58 | logFCgreat.less.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN<(-0.58))
#logFCgreat.ROGTTEN.GFPTEN <- logFCgreat.ROGTTEN.GFPTEN[logFCgreat.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN > 1.5 | logFCgreat.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN < -1.5, ]
head(logFCgreat.less.ROGTTEN.GFPTEN)
dim(logFCgreat.less.ROGTTEN.GFPTEN)
#686   2
#351+335
nrow(logFCgreat.ROGTTEN.GFPTEN) #351
nrow(logFCless.ROGTTEN.GFPTEN)  #335
nrow(logFCgreat.ROGTTEN.GFPTEN)+nrow(logFCless.ROGTTEN.GFPTEN)


#-------------------------#
#5. ROGTCTRL.GFPCTRL
#-------------------------#

head(remsortupdown.ROGTCTRL.GFPCTRL)
dim(remsortupdown.ROGTCTRL.GFPCTRL)
logFCgreat.ROGTCTRL.GFPCTRL<-remsortupdown.ROGTCTRL.GFPCTRL
head(logFCgreat.ROGTCTRL.GFPCTRL)
logFCgreat.ROGTCTRL.GFPCTRL<-as.data.frame(logFCgreat.ROGTCTRL.GFPCTRL)
#logFCgreat.ROGTCTRL.GFPCTRL<-logFCgreat.ROGTCTRL.GFPCTRL[,-2]
head(logFCgreat.ROGTCTRL.GFPCTRL)
logFCgreat.ROGTCTRL.GFPCTRL<-subset(logFCgreat.ROGTCTRL.GFPCTRL,(logFCgreat.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL)>0.58)
head(logFCgreat.ROGTCTRL.GFPCTRL)
dim(logFCgreat.ROGTCTRL.GFPCTRL)
#351, 2

logFCless.ROGTCTRL.GFPCTRL<-remsortupdown.ROGTCTRL.GFPCTRL
head(logFCless.ROGTCTRL.GFPCTRL)
logFCless.ROGTCTRL.GFPCTRL<-as.data.frame(logFCless.ROGTCTRL.GFPCTRL)
head(logFCless.ROGTCTRL.GFPCTRL)
type(logFCless.ROGTCTRL.GFPCTRL)
#logFCless.ROGTCTRL.GFPCTRL<-as.numeric(logFCless.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL)
#head(logFCless.ROGTCTRL.GFPCTRL)
logFCless.ROGTCTRL.GFPCTRL <- transform(logFCless.ROGTCTRL.GFPCTRL, logFC.ROGTCTRL.GFPCTRL = as.numeric(logFC.ROGTCTRL.GFPCTRL))
head(logFCless.ROGTCTRL.GFPCTRL)
type(logFCless.ROGTCTRL.GFPCTRL)
#logFCless.ROGTCTRL.GFPCTRL<-subset(logFCless.ROGTCTRL.GFPCTRL,(((logFCless.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL)<(-1.5))&((logFCless.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL)>(1.5))))
logFCless.ROGTCTRL.GFPCTRL<-subset(logFCless.ROGTCTRL.GFPCTRL,(((logFCless.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL)<(-0.58))))
head(logFCless.ROGTCTRL.GFPCTRL)
dim(logFCless.ROGTCTRL.GFPCTRL)
#335, 2


logFCgreat.less.ROGTCTRL.GFPCTRL<-remsortupdown.ROGTCTRL.GFPCTRL
head(logFCgreat.less.ROGTCTRL.GFPCTRL)
dim(logFCgreat.less.ROGTCTRL.GFPCTRL)
logFCgreat.less.ROGTCTRL.GFPCTRL<-as.data.frame(logFCgreat.less.ROGTCTRL.GFPCTRL)
head(logFCgreat.less.ROGTCTRL.GFPCTRL)
type(logFCgreat.less.ROGTCTRL.GFPCTRL)
#logFCgreat.less.ROGTCTRL.GFPCTRL<-as.numeric(logFCgreat.less.ROGTCTRL.GFPCTRL$logFCgreat.less.ROGTCTRL.GFPCTRL)
#head(logFCgreat.less.ROGTCTRL.GFPCTRL)
#logFCless.ROGTCTRL.GFPCTRL <- transform(logFCless.ROGTCTRL.GFPCTRL, logFC.ROGTCTRL.GFPCTRL = as.numeric(logFC.ROGTCTRL.GFPCTRL))
logFCgreat.less.ROGTCTRL.GFPCTRL <- transform(logFCgreat.less.ROGTCTRL.GFPCTRL, logFC.ROGTCTRL.GFPCTRL = as.numeric(logFC.ROGTCTRL.GFPCTRL))
head(logFCgreat.less.ROGTCTRL.GFPCTRL)
type(logFCgreat.less.ROGTCTRL.GFPCTRL)
#logFCgreat.less.ROGTCTRL.GFPCTRL<-subset(logFCgreat.less.ROGTCTRL.GFPCTRL,(((logFCgreat.less.ROGTCTRL.GFPCTRL$logFCgreat.less.ROGTCTRL.GFPCTRL)<(-1.5))&((logFCgreat.less.ROGTCTRL.GFPCTRL$logFCgreat.less.ROGTCTRL.GFPCTRL)>(1.5))))
#logFCgreat.less.ROGTCTRL.GFPCTRL<-subset(logFCgreat.less.ROGTCTRL.GFPCTRL,(((logFCgreat.less.ROGTCTRL.GFPCTRL$logFCgreat.less.ROGTCTRL.GFPCTRL)>(1.5))))
logFCgreat.less.ROGTCTRL.GFPCTRL<-as.data.frame(logFCgreat.less.ROGTCTRL.GFPCTRL)
head(logFCgreat.less.ROGTCTRL.GFPCTRL)
type(logFCgreat.less.ROGTCTRL.GFPCTRL)
logFCgreat.less.ROGTCTRL.GFPCTRL<-subset(logFCgreat.less.ROGTCTRL.GFPCTRL,logFCgreat.less.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL>0.58 | logFCgreat.less.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL<(-0.58))
#logFCgreat.ROGTCTRL.GFPCTRL <- logFCgreat.ROGTCTRL.GFPCTRL[logFCgreat.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL > 1.5 | logFCgreat.ROGTCTRL.GFPCTRL$logFC.ROGTCTRL.GFPCTRL < -1.5, ]
head(logFCgreat.less.ROGTCTRL.GFPCTRL)
dim(logFCgreat.less.ROGTCTRL.GFPCTRL)
#686   2
#351+335
nrow(logFCgreat.ROGTCTRL.GFPCTRL) #351
nrow(logFCless.ROGTCTRL.GFPCTRL)  #335
nrow(logFCgreat.ROGTCTRL.GFPCTRL)+nrow(logFCless.ROGTCTRL.GFPCTRL)


#-------------------------#
#6. ROGTTEN.ROGTCTRL
#-------------------------#

head(remsortupdown.ROGTTEN.ROGTCTRL)
dim(remsortupdown.ROGTTEN.ROGTCTRL)
logFCgreat.ROGTTEN.ROGTCTRL<-remsortupdown.ROGTTEN.ROGTCTRL
head(logFCgreat.ROGTTEN.ROGTCTRL)
logFCgreat.ROGTTEN.ROGTCTRL<-as.data.frame(logFCgreat.ROGTTEN.ROGTCTRL)
#logFCgreat.ROGTTEN.ROGTCTRL<-logFCgreat.ROGTTEN.ROGTCTRL[,-2]
head(logFCgreat.ROGTTEN.ROGTCTRL)
logFCgreat.ROGTTEN.ROGTCTRL<-subset(logFCgreat.ROGTTEN.ROGTCTRL,(logFCgreat.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL)>0.58)
head(logFCgreat.ROGTTEN.ROGTCTRL)
dim(logFCgreat.ROGTTEN.ROGTCTRL)
#351, 2

logFCless.ROGTTEN.ROGTCTRL<-remsortupdown.ROGTTEN.ROGTCTRL
head(logFCless.ROGTTEN.ROGTCTRL)
logFCless.ROGTTEN.ROGTCTRL<-as.data.frame(logFCless.ROGTTEN.ROGTCTRL)
head(logFCless.ROGTTEN.ROGTCTRL)
type(logFCless.ROGTTEN.ROGTCTRL)
#logFCless.ROGTTEN.ROGTCTRL<-as.numeric(logFCless.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL)
#head(logFCless.ROGTTEN.ROGTCTRL)
logFCless.ROGTTEN.ROGTCTRL <- transform(logFCless.ROGTTEN.ROGTCTRL, logFC.ROGTTEN.ROGTCTRL = as.numeric(logFC.ROGTTEN.ROGTCTRL))
head(logFCless.ROGTTEN.ROGTCTRL)
type(logFCless.ROGTTEN.ROGTCTRL)
#logFCless.ROGTTEN.ROGTCTRL<-subset(logFCless.ROGTTEN.ROGTCTRL,(((logFCless.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL)<(-1.5))&((logFCless.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL)>(1.5))))
logFCless.ROGTTEN.ROGTCTRL<-subset(logFCless.ROGTTEN.ROGTCTRL,(((logFCless.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL)<(-0.58))))
head(logFCless.ROGTTEN.ROGTCTRL)
dim(logFCless.ROGTTEN.ROGTCTRL)
#335, 2


logFCgreat.less.ROGTTEN.ROGTCTRL<-remsortupdown.ROGTTEN.ROGTCTRL
head(logFCgreat.less.ROGTTEN.ROGTCTRL)
dim(logFCgreat.less.ROGTTEN.ROGTCTRL)
logFCgreat.less.ROGTTEN.ROGTCTRL<-as.data.frame(logFCgreat.less.ROGTTEN.ROGTCTRL)
head(logFCgreat.less.ROGTTEN.ROGTCTRL)
type(logFCgreat.less.ROGTTEN.ROGTCTRL)
#logFCgreat.less.ROGTTEN.ROGTCTRL<-as.numeric(logFCgreat.less.ROGTTEN.ROGTCTRL$logFCgreat.less.ROGTTEN.ROGTCTRL)
#head(logFCgreat.less.ROGTTEN.ROGTCTRL)
#logFCless.ROGTTEN.ROGTCTRL <- transform(logFCless.ROGTTEN.ROGTCTRL, logFC.ROGTTEN.ROGTCTRL = as.numeric(logFC.ROGTTEN.ROGTCTRL))
logFCgreat.less.ROGTTEN.ROGTCTRL <- transform(logFCgreat.less.ROGTTEN.ROGTCTRL, logFC.ROGTTEN.ROGTCTRL = as.numeric(logFC.ROGTTEN.ROGTCTRL))
head(logFCgreat.less.ROGTTEN.ROGTCTRL)
type(logFCgreat.less.ROGTTEN.ROGTCTRL)
#logFCgreat.less.ROGTTEN.ROGTCTRL<-subset(logFCgreat.less.ROGTTEN.ROGTCTRL,(((logFCgreat.less.ROGTTEN.ROGTCTRL$logFCgreat.less.ROGTTEN.ROGTCTRL)<(-1.5))&((logFCgreat.less.ROGTTEN.ROGTCTRL$logFCgreat.less.ROGTTEN.ROGTCTRL)>(1.5))))
#logFCgreat.less.ROGTTEN.ROGTCTRL<-subset(logFCgreat.less.ROGTTEN.ROGTCTRL,(((logFCgreat.less.ROGTTEN.ROGTCTRL$logFCgreat.less.ROGTTEN.ROGTCTRL)>(1.5))))
logFCgreat.less.ROGTTEN.ROGTCTRL<-as.data.frame(logFCgreat.less.ROGTTEN.ROGTCTRL)
head(logFCgreat.less.ROGTTEN.ROGTCTRL)
type(logFCgreat.less.ROGTTEN.ROGTCTRL)
logFCgreat.less.ROGTTEN.ROGTCTRL<-subset(logFCgreat.less.ROGTTEN.ROGTCTRL,logFCgreat.less.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL>0.58 | logFCgreat.less.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL<(-0.58))
#logFCgreat.ROGTTEN.ROGTCTRL <- logFCgreat.ROGTTEN.ROGTCTRL[logFCgreat.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL > 1.5 | logFCgreat.ROGTTEN.ROGTCTRL$logFC.ROGTTEN.ROGTCTRL < -1.5, ]
head(logFCgreat.less.ROGTTEN.ROGTCTRL)
dim(logFCgreat.less.ROGTTEN.ROGTCTRL)
#686   2
#351+335
nrow(logFCgreat.ROGTTEN.ROGTCTRL) #351
nrow(logFCless.ROGTTEN.ROGTCTRL)  #335
nrow(logFCgreat.ROGTTEN.ROGTCTRL)+nrow(logFCless.ROGTTEN.ROGTCTRL)


#-------------------------#
#7. GFPTEN.GFPCTRL
#-------------------------#

head(remsortupdown.GFPTEN.GFPCTRL)
dim(remsortupdown.GFPTEN.GFPCTRL)
logFCgreat.GFPTEN.GFPCTRL<-remsortupdown.GFPTEN.GFPCTRL
head(logFCgreat.GFPTEN.GFPCTRL)
logFCgreat.GFPTEN.GFPCTRL<-as.data.frame(logFCgreat.GFPTEN.GFPCTRL)
#logFCgreat.GFPTEN.GFPCTRL<-logFCgreat.GFPTEN.GFPCTRL[,-2]
head(logFCgreat.GFPTEN.GFPCTRL)
logFCgreat.GFPTEN.GFPCTRL<-subset(logFCgreat.GFPTEN.GFPCTRL,(logFCgreat.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL)>0.58)
head(logFCgreat.GFPTEN.GFPCTRL)
dim(logFCgreat.GFPTEN.GFPCTRL)
#351, 2

logFCless.GFPTEN.GFPCTRL<-remsortupdown.GFPTEN.GFPCTRL
head(logFCless.GFPTEN.GFPCTRL)
logFCless.GFPTEN.GFPCTRL<-as.data.frame(logFCless.GFPTEN.GFPCTRL)
head(logFCless.GFPTEN.GFPCTRL)
type(logFCless.GFPTEN.GFPCTRL)
#logFCless.GFPTEN.GFPCTRL<-as.numeric(logFCless.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL)
#head(logFCless.GFPTEN.GFPCTRL)
logFCless.GFPTEN.GFPCTRL <- transform(logFCless.GFPTEN.GFPCTRL, logFC.GFPTEN.GFPCTRL = as.numeric(logFC.GFPTEN.GFPCTRL))
head(logFCless.GFPTEN.GFPCTRL)
type(logFCless.GFPTEN.GFPCTRL)
#logFCless.GFPTEN.GFPCTRL<-subset(logFCless.GFPTEN.GFPCTRL,(((logFCless.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL)<(-1.5))&((logFCless.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL)>(1.5))))
logFCless.GFPTEN.GFPCTRL<-subset(logFCless.GFPTEN.GFPCTRL,(((logFCless.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL)<(-0.58))))
head(logFCless.GFPTEN.GFPCTRL)
dim(logFCless.GFPTEN.GFPCTRL)
#335, 2


logFCgreat.less.GFPTEN.GFPCTRL<-remsortupdown.GFPTEN.GFPCTRL
head(logFCgreat.less.GFPTEN.GFPCTRL)
dim(logFCgreat.less.GFPTEN.GFPCTRL)
logFCgreat.less.GFPTEN.GFPCTRL<-as.data.frame(logFCgreat.less.GFPTEN.GFPCTRL)
head(logFCgreat.less.GFPTEN.GFPCTRL)
type(logFCgreat.less.GFPTEN.GFPCTRL)
#logFCgreat.less.GFPTEN.GFPCTRL<-as.numeric(logFCgreat.less.GFPTEN.GFPCTRL$logFCgreat.less.GFPTEN.GFPCTRL)
#head(logFCgreat.less.GFPTEN.GFPCTRL)
#logFCless.GFPTEN.GFPCTRL <- transform(logFCless.GFPTEN.GFPCTRL, logFC.GFPTEN.GFPCTRL = as.numeric(logFC.GFPTEN.GFPCTRL))
logFCgreat.less.GFPTEN.GFPCTRL <- transform(logFCgreat.less.GFPTEN.GFPCTRL, logFC.GFPTEN.GFPCTRL = as.numeric(logFC.GFPTEN.GFPCTRL))
head(logFCgreat.less.GFPTEN.GFPCTRL)
type(logFCgreat.less.GFPTEN.GFPCTRL)
#logFCgreat.less.GFPTEN.GFPCTRL<-subset(logFCgreat.less.GFPTEN.GFPCTRL,(((logFCgreat.less.GFPTEN.GFPCTRL$logFCgreat.less.GFPTEN.GFPCTRL)<(-1.5))&((logFCgreat.less.GFPTEN.GFPCTRL$logFCgreat.less.GFPTEN.GFPCTRL)>(1.5))))
#logFCgreat.less.GFPTEN.GFPCTRL<-subset(logFCgreat.less.GFPTEN.GFPCTRL,(((logFCgreat.less.GFPTEN.GFPCTRL$logFCgreat.less.GFPTEN.GFPCTRL)>(1.5))))
logFCgreat.less.GFPTEN.GFPCTRL<-as.data.frame(logFCgreat.less.GFPTEN.GFPCTRL)
head(logFCgreat.less.GFPTEN.GFPCTRL)
type(logFCgreat.less.GFPTEN.GFPCTRL)
logFCgreat.less.GFPTEN.GFPCTRL<-subset(logFCgreat.less.GFPTEN.GFPCTRL,logFCgreat.less.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL>0.58 | logFCgreat.less.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL<(-0.58))
#logFCgreat.GFPTEN.GFPCTRL <- logFCgreat.GFPTEN.GFPCTRL[logFCgreat.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL > 1.5 | logFCgreat.GFPTEN.GFPCTRL$logFC.GFPTEN.GFPCTRL < -1.5, ]
head(logFCgreat.less.GFPTEN.GFPCTRL)
dim(logFCgreat.less.GFPTEN.GFPCTRL)
#686   2
#351+335
nrow(logFCgreat.GFPTEN.GFPCTRL) #351
nrow(logFCless.GFPTEN.GFPCTRL)  #335
nrow(logFCgreat.GFPTEN.GFPCTRL)+nrow(logFCless.GFPTEN.GFPCTRL)

#==================================================================================================================================#

#============================================================================#
#Prepare protlist for waterfall plot (Plot of ranked Fold Changes) 
   # - rank-ordered protein (gene) list
#============================================================================#

# Prepare the gene list
#gene.list <- d.tmp$logFC             # rank-ordered gene list
#head(dea$tt)
#names(gene.list) <- d.tmp$hgnc_symbol
#head(gene.list)

#-------------------------#
#ROGTTEN.GFPTEN
#-------------------------#
#GSEA 606 PROT LIST
head(sortupdown.ROGTTEN.GFPTEN)
dim(sortupdown.ROGTTEN.GFPTEN)
#8477   2
prot.list.ROGTTEN.GFPTEN <- sortupdown.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN    # rank-ordered gene list
head(prot.list.ROGTTEN.GFPTEN)
prot.list.ROGTTEN.GFPTEN<-as.numeric(prot.list.ROGTTEN.GFPTEN)
head(prot.list.ROGTTEN.GFPTEN)
head(sortupdown.ROGTTEN.GFPTEN)
names(prot.list.ROGTTEN.GFPTEN) <- rownames(sortupdown.ROGTTEN.GFPTEN)
head(prot.list.ROGTTEN.GFPTEN)

#-------------------------#
#OGTTEN.GFPTEN
#-------------------------#
#GSEA 605 PROT LIST
head(sortupdown.OGTTEN.GFPTEN)
dim(sortupdown.OGTTEN.GFPTEN)
#[1] 8477    2
prot.list.OGTTEN.GFPTEN <- sortupdown.OGTTEN.GFPTEN$logFC.OGTTEN.GFPTEN             # rank-ordered gene list
head(prot.list.OGTTEN.GFPTEN)
prot.list.OGTTEN.GFPTEN<-as.numeric(prot.list.OGTTEN.GFPTEN)
head(prot.list.OGTTEN.GFPTEN)
head(sortupdown.OGTTEN.GFPTEN)
names(prot.list.OGTTEN.GFPTEN) <- rownames(sortupdown.OGTTEN.GFPTEN)
head(prot.list.OGTTEN.GFPTEN)

##############################################################################################################################
#head(remsortupdownOGTTENGFPTEN) #This is sorted based on p value and filtered for FC and p value so dont use it.
#prot.list <- remsortupdownOGTTENGFPTEN$logFC.OGTTEN.GFPTEN             # rank-ordered gene list
#head(prot.list)
#prot.list<-as.numeric(prot.list)
#head(prot.list)
#names(prot.list) <- rownames(remsortupdownOGTTENGFPTEN)
#head(prot.list)

#head(logFCgreat)
#dim(logFCgreat)
#prot.list2 <- logFCgreat$logFC.OGTTEN.GFPTEN             # rank-ordered gene list
#head(prot.list2)
#prot.list2<-as.numeric(prot.list2)
#head(prot.list2)
#names(prot.list2) <- rownames(logFCgreat)
#head(prot.list2)
##############################################################################################################################

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
head(remsortupdown.ROGTTEN.GFPTEN)
d.go2 <- remsortupdown.ROGTTEN.GFPTEN
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

###KEGG analysis
d.kegg <- d.tmp
#d.kegg.DE <- subset(d.kegg,FDR<0.05) 
d.kegg.DE <- subset(d.kegg,PValue<0.05)
head(d.kegg.DE)
all(d.kegg.DE$hgnc_symbol==names(d.entrez.id)) 

kegg.test <- kegga(d.entrez.id,species="Hs")
kegg.results <- topKEGG(kegg.test, sort = "DE", number = Inf)
head(kegg.results)
sum(kegg.results$P.DE<10^(-5))
sum(kegg.results$P.DE<0.05)

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

### GSEA analysis (Template)
# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
#all.gene.sets <- gmtPathways("C:\\Users\\Sophiya John\\Desktop\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
#class(all.gene.sets)
#length(all.gene.sets)
#all.gene.sets[1:2]
#Show first a few pathways, and within those, show only the first few genes. 
#library(tidyverse)
#all.gene.sets %>% head() %>% lapply(head)

### Now run fgsea 
#fgseaRes <- fgsea(pathways = all.gene.sets,stats = gene.list,minSize=15,maxSize=500,eps=0)
#head(fgseaRes)
#head(fgseaRes[order(pval), ])
#sum(fgseaRes[, padj < 0.05])
# Make a few Enrichment Plots
#Plot1
#plotEnrichment(all.gene.sets[["GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS"]],prot.list) + labs(title="GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS")
#Plot2
#plotEnrichment(all.gene.sets[["CHEN_METABOLIC_SYNDROM_NETWORK"]],gene.list) + labs(title="CHEN_METABOLIC_SYNDROM_NETWORK")
#tail(fgseaRes[order(pval), ])
#Plot3
#plotEnrichment(all.gene.sets[["ZNF610_TARGET_GENES"]],gene.list) + labs(title="ZNF610_TARGET_GENES")
#Plot4
#plotEnrichment(all.gene.sets[["ZSCAN5DP_TARGET_GENES"]],gene.list) + labs(title="ZSCAN5DP_TARGET_GENES")
#fgseaRes[order(pval), ][15:20,]
#Plot5
#plotEnrichment(all.gene.sets[["GOBP_RESPONSE_TO_BIOTIC_STIMULUS"]],gene.list) + labs(title="GOBP_RESPONSE_TO_BIOTIC_STIMULUS")

# Make a table plot for a bunch of selected pathways:
#topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
#topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
#topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#plotGseaTable(all.gene.sets[topPathways], gene.list, fgseaRes,gseaParam = 0.5)

###################################################################################################################

######################################################################################################################

### J_GSEA analysis

# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
library(fgsea)

# Install and load the library
# install.packages("msigdbr")
library(msigdbr)

# Pull the latest Human (Hs) gene sets (e.g., Hallmarks)
#h_df <- msigdbr(species = "Homo sapiens", category = "H") #Hallmark list

# Convert to a list format for fgsea
#all_gene_sets <- split(x = h_df$gene_symbol, f = h_df$gs_name)
#head(h_df)
#colnames(h_df)
# Check the version being used
#unique(h_df$db_version)
#class(all_gene_sets)
#length(all_gene_sets)

# Pull all MSigDB collections (caution: this will be a huge list)
all_msigdb_df <- msigdbr(species = "Homo sapiens")
all_gene_sets2 <- split(x = all_msigdb_df$gene_symbol, f = all_msigdb_df$gs_name)

length(all_gene_sets2) # Should be ~33,000+
#[1] 35134
#names(all_gene_sets2) <- gsub("^[^_]*_", "", names(all_gene_sets2))
#head(all_gene_sets2)

#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v2023.2.Hs.symbols.gmt") #got lot of warnings
#all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_DESKTOP 2025\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
#class(all.gene.sets)
#length(all.gene.sets)
#all.gene.sets[1:2]
#all_gene_sets2[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
library(tidyverse)

#The following command will give a quick "peek" at the first 6 pathways 
#and the first 6 genes within each.

#all.gene.sets %>% head() %>% lapply(head) #2023 list
all_gene_sets2 %>% head() %>% lapply(head) #updated 2025 list
#all_gene_sets %>% head() %>% lapply(head) #Hallmark list

#all.gene.sets1 <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v2023.2.Hs.symbols.gmt")
#class(all.gene.sets1)
#length(all.gene.sets1)
#all.gene.sets1[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
#library(tidyverse)
#all.gene.sets1 %>% head() %>% lapply(head)

head(prot.list) ##GSEA 606 PROT LIST
head(prot.list2)#GSEA 605 PROT LIST

#head(d.entrez.id1)
#d.entrez.id2<-
### Now run fgsea 
#names(all_gene_sets2) <- gsub("^[^_]*_", "", names(all_gene_sets2))
#head(all_gene_sets2)

#fgseaRes1 <- fgsea(pathways = all.gene.sets, stats = prot.list, minSize=15, maxSize=500,eps=0)
fgseaRes2 <- fgsea(pathways = all_gene_sets2, stats = prot.list, minSize=15, maxSize=500,eps=0)
#fgseaRes3 <- fgsea(pathways = all_gene_sets, stats = prot.list, minSize=15, maxSize=500,eps=0)

#head(fgseaRes1)#2023 list
head(fgseaRes2) # updated 2025 list
#head(fgseaRes3) # Hallmark list
#head(fgseaRes1[order(pval), ])
head(fgseaRes2[order(pval), ])
#head(fgseaRes3[order(pval), ])

#sum(fgseaRes1[, padj < 0.05]) #261 #2023 list
sum(fgseaRes2[, padj < 0.05]) #289 # updated 2025 list - so better to use
#sum(fgseaRes3[, padj < 0.05]) #1   # Hallmark list

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

# This removes the prefix (e.g., "KEGG_") from the names in the list

p<-plotGseaTable(all_gene_sets2[topPathways], prot.list, fgseaRes2,gseaParam = 0.5, pathwayLabelStyle=list(size=7, fontface = "bold"), 
                 valueStyle = list(size=7),
                 axisLabelStyle = list(size=4))

# 2. Add the title
grid.arrange(p, top = textGrob("GSEA - OGT KD 606 10 min vs GFP 10 min", 
                               x = 0.62, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))

# 1. Clean the names (remove prefixes)
fgseaRes2[, pathway := gsub("^[^_]*_", "", pathway)]

# 2. Sort by significance (padj) so the best version of a pathway is on top
fgseaRes2 <- fgseaRes2[order(padj)]

# 3. Remove duplicates based on the pathway name
# This keeps ONLY the most significant version of "OXIDATIVE_PHOSPHORYLATION"
fgseaRes_unique <- fgseaRes2[!duplicated(pathway)]

# 4. Now select your Top 10 from the UNIQUE list
topPathwaysUp <- fgseaRes_unique[ES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes_unique[ES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# 5. Make sure the pathway list names match the cleaned names
names(all_gene_sets2) <- gsub("^[^_]*_", "", names(all_gene_sets2))

# 6. Re-generate the plot
p <- plotGseaTable(all_gene_sets2[topPathways], 
                   prot.list, 
                   fgseaRes_unique, 
                   gseaParam = 0.5, 
                   pathwayLabelStyle = list(size=7, fontface = "bold"), 
                   valueStyle = list(size=7),
                   axisLabelStyle = list(size=4))

# 7. Draw with the title
grid.arrange(p, top = textGrob("GSEA - OGT KD 606 10 min vs GFP 10 min", 
                               x = 0.65, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))


#ggsave(p, width=10, height=8, file="table.png")

#==================================================================#
### J_GSEA analysis - Final code - OGT KD 606 10MIN VS GFP 10 MIN
#==================================================================#

# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
  
library(fgsea)
library(msigdbr)

# Pull all MSigDB collections (caution: this will be a huge list)
all_msigdb_df <- msigdbr(species = "Homo sapiens")
all_gene_sets2 <- split(x = all_msigdb_df$gene_symbol, f = all_msigdb_df$gs_name)

length(all_gene_sets2) # Should be ~33,000+
#[1] 35134

library(tidyverse)

#The following command will give a quick "peek" at the first 6 pathways 
#and the first 6 genes within each.

#all.gene.sets %>% head() %>% lapply(head) #2023 list
all_gene_sets2 %>% head() %>% lapply(head) #updated 2025 list

head(prot.list) ##GSEA 606 PROT LIST
head(prot.list2)#GSEA 605 PROT LIST

fgseaRes2 <- fgsea(pathways = all_gene_sets2, stats = prot.list, minSize=15, maxSize=500,eps=0)
head(fgseaRes2) # updated 2025 list
head(fgseaRes2[order(pval), ])
sum(fgseaRes2[, padj < 0.05]) #289 # updated 2025 list - so better to use

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
grid.arrange(p, top = textGrob("GSEA - OGT KD 606 10 min vs GFP 10 min", 
                               x = 0.62, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))

# 1. Clean the names (remove prefixes)
fgseaRes2[, pathway := gsub("^[^_]*_", "", pathway)]

# 2. Sort by significance (padj) so the best version of a pathway is on top
fgseaRes2 <- fgseaRes2[order(padj)]

# 3. Remove duplicates based on the pathway name
# This keeps ONLY the most significant version of "OXIDATIVE_PHOSPHORYLATION"
fgseaRes_unique <- fgseaRes2[!duplicated(pathway)]

# 4. Now select your Top 10 from the UNIQUE list
topPathwaysUp <- fgseaRes_unique[ES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes_unique[ES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# 5. Make sure the pathway list names match the cleaned names
names(all_gene_sets2) <- gsub("^[^_]*_", "", names(all_gene_sets2))

# 6. Re-generate the plot
p <- plotGseaTable(all_gene_sets2[topPathways], 
                   prot.list, 
                   fgseaRes_unique, 
                   gseaParam = 0.5, 
                   pathwayLabelStyle = list(size=7, fontface = "bold"), 
                   valueStyle = list(size=7),
                   axisLabelStyle = list(size=4))

# 7. Draw with the title
grid.arrange(p, top = textGrob("GSEA - OGT KD 606 10 min vs GFP 10 min", 
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
plotEnrichment(all_gene_sets2[["OXIDATIVE_PHOSPHORYLATION"]],prot.list) + labs(title="OXIDATIVE_PHOSPHORYLATION")
plotEnrichment(all_gene_sets2[["OXIDATIVE_PHOSPHORYLATION"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n OXIDATIVE_PHOSPHORYLATION",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot2
plotEnrichment(all_gene_sets2[["ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]],prot.list) + labs(title="ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA")
plotEnrichment(all_gene_sets2[["ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot3
plotEnrichment(all_gene_sets2[["ORGANELLE_INNER_MEMBRANE"]],prot.list) + labs(title="ORGANELLE_INNER_MEMBRANE")
plotEnrichment(all_gene_sets2[["ORGANELLE_INNER_MEMBRANE"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n ORGANELLE_INNER_MEMBRANE",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot4
plotEnrichment(all_gene_sets2[["AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT"]],prot.list) + labs(title="AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT")
plotEnrichment(all_gene_sets2[["AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Top 4 down

#Plot1
plotEnrichment(all.gene.sets[["EUKARYOTIC_TRANSLATION_INITIATION"]],prot.list) + labs(title="EUKARYOTIC_TRANSLATION_INITIATION")
plotEnrichment(all_gene_sets2[["EUKARYOTIC_TRANSLATION_INITIATION"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n EUKARYOTIC_TRANSLATION_INITIATION",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )
#gseaplot(all.gene.sets[["REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"]],prot.list) + labs(title="REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")

#Plot2
plotEnrichment(all.gene.sets[["NONSENSE_MEDIATED_DECAY_NMD"]],prot.list) + labs(title="NONSENSE_MEDIATED_DECAY_NMD")
plotEnrichment(all_gene_sets2[["NONSENSE_MEDIATED_DECAY_NMD"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n NONSENSE_MEDIATED_DECAY_NMD",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot3
plotEnrichment(all.gene.sets[["EUKARYOTIC_TRANSLATION_ELONGATION"]],prot.list) + labs(title="EUKARYOTIC_TRANSLATION_ELONGATION")
plotEnrichment(all_gene_sets2[["EUKARYOTIC_TRANSLATION_ELONGATION"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n EUKARYOTIC_TRANSLATION_ELONGATION",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot4
plotEnrichment(all.gene.sets[["SELENOAMINO_ACID_METABOLISM"]],prot.list) + labs(title="SELENOAMINO_ACID_METABOLISM")
plotEnrichment(all_gene_sets2[["SELENOAMINO_ACID_METABOLISM"]], prot.list) + 
  labs(
    title = "GSEA_OGT KD 606 10 MIN vs. GFP 10 MIN \n SELENOAMINO_ACID_METABOLISM",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#----------------------------------------------------------------------------------------------------

#Pretty gsea visualizations    #07152024 start here

library(org.Hs.eg.db)
#install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
head(prot.list)

head(prot.list)
type(prot.list)


#GSEA 606 PROT LIST


library(enrichplot)

head(fgseaRes2)
#listpathway<-fgseaRes[,c("pathway","padj","ES")]
#head(listpathway)
#topPathwaysUp1 <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
#listpathway <- listpathway[order(listpathway$padj),] #to sort based on padj
#head(listpathway)
#listpathway<-as.matrix(listpathway)

topPathwaysUp1 <- fgseaRes_unique[ES > 0][head(order(padj), n=10)]
topPathwaysUp1<-topPathwaysUp1[,-2]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-(4:5)]
topPathwaysUp1
topPathwaysUp1<-as.data.frame(topPathwaysUp1)

topPathwaysDown1 <- fgseaRes_unique[ES < 0][head(order(padj), n=10)]
topPathwaysDown1<-topPathwaysDown1[,-2]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-(4:5)]
topPathwaysDown1
topPathwaysDown1<-as.data.frame(topPathwaysDown1)


topPathways1 <- rbind(topPathwaysUp1, rev(topPathwaysDown1)) #merge rows
topPathways1
topPathways1<-as.data.frame(topPathways1)
topPathways1

#https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package
#ggplot(topPathways1, aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=padj<0.05&NES>0)) +
#  coord_flip() +
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="GSEA TOP 10 UP AND DOWN REGULATED PATHWAYS \n OGT KD 606 10 MIN. vs GFP 10 MIN.") + 
#  theme_minimal()+
#  theme(text = element_text(face = "bold"))

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
    title = "GSEA: Top 10 Up and Down Regulated Pathways \n OGT KD 606 10 MIN. vs GFP 10 MIN."
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
    title = "GSEA: Top 10 Up and Down Regulated Pathways \n OGT KD 606 10 MIN. vs GFP 10 MIN."
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


#==================================================================#
### J_GSEA analysis - Final code - OGT KD 605 10MIN VS GFP 10 MIN
#==================================================================#


# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb

library(fgsea)
library(msigdbr)

# Pull all MSigDB collections (caution: this will be a huge list)
all_msigdb_df <- msigdbr(species = "Homo sapiens")
all_gene_sets2 <- split(x = all_msigdb_df$gene_symbol, f = all_msigdb_df$gs_name)

length(all_gene_sets2) # Should be ~33,000+
#[1] 35134

library(tidyverse)

#The following command will give a quick "peek" at the first 6 pathways 
#and the first 6 genes within each.

#all.gene.sets %>% head() %>% lapply(head) #2023 list
all_gene_sets2 %>% head() %>% lapply(head) #updated 2025 list

head(prot.list) ##GSEA 606 PROT LIST
head(prot.list2)#GSEA 605 PROT LIST

fgseaRes3 <- fgsea(pathways = all_gene_sets2, stats = prot.list2, minSize=15, maxSize=500,eps=0)
head(fgseaRes3) # updated 2025 list
head(fgseaRes3[order(pval), ])
sum(fgseaRes3[, padj < 0.05]) #278 # updated 2025 list - so better to use
head(fgseaRes3)
fgseaRes5<-fgseaRes3[order(padj), ]
head(fgseaRes5)

# Make a table plot for a bunch of selected pathways:
#if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)
library(grid)

topPathwaysUp2 <- fgseaRes3[ES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
#topPathwaysUp <- fgseaRes2[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown2 <- fgseaRes3[ES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways2 <- c(topPathwaysUp2, rev(topPathwaysDown2))
head(topPathways2)
#with prefixes 
p2<-plotGseaTable(all_gene_sets2[topPathways2], prot.list2, fgseaRes3,gseaParam = 0.5, pathwayLabelStyle=list(size=7, fontface = "bold"), 
                 valueStyle = list(size=7),
                 axisLabelStyle = list(size=4))

# 2. Add the title
grid.arrange(p2, top = textGrob("GSEA - OGT KD 605 10 min vs GFP 10 min", 
                               x = 0.62, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))

# 1. Clean the names (remove prefixes)
fgseaRes3[, pathway := gsub("^[^_]*_", "", pathway)]

# 2. Sort by significance (padj) so the best version of a pathway is on top
fgseaRes3 <- fgseaRes3[order(padj)]

# 3. Remove duplicates based on the pathway name
# This keeps ONLY the most significant version of "OXIDATIVE_PHOSPHORYLATION"
fgseaRes_unique2 <- fgseaRes3[!duplicated(pathway)]

# 4. Now select your Top 10 from the UNIQUE list
topPathwaysUp2 <- fgseaRes_unique2[ES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysDown2 <- fgseaRes_unique2[ES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways2 <- c(topPathwaysUp2, rev(topPathwaysDown2))

# 5. Make sure the pathway list names match the cleaned names
names(all_gene_sets2) <- gsub("^[^_]*_", "", names(all_gene_sets2))

# 6. Re-generate the plot
p2 <- plotGseaTable(all_gene_sets2[topPathways2], 
                   prot.list2, 
                   fgseaRes_unique2, 
                   gseaParam = 0.5, 
                   pathwayLabelStyle = list(size=7, fontface = "bold"), 
                   valueStyle = list(size=7),
                   axisLabelStyle = list(size=4))

# 7. Draw with the title
grid.arrange(p2, top = textGrob("GSEA - OGT KD 605 10 min vs GFP 10 min", 
                               x = 0.65, 
                               just = "center",
                               gp = gpar(fontsize = 12, fontface = "bold")))

head(topPathways2)

library(ggplot2)
topPathways2
# Make a few Enrichment Plots

#Top 6 up
head(prot.list2)
#Plot6
plotEnrichment(all_gene_sets2[["OXIDATIVE_PHOSPHORYLATION"]],prot.list2) + labs(title="OXIDATIVE_PHOSPHORYLATION")
plotEnrichment(all_gene_sets2[["OXIDATIVE_PHOSPHORYLATION"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n OXIDATIVE_PHOSPHORYLATION (PATHWAY 6)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot5
plotEnrichment(all_gene_sets2[["ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]],prot.list2) + labs(title="ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA")
plotEnrichment(all_gene_sets2[["ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA (PATHWAY 5)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot1
plotEnrichment(all_gene_sets2[["RNA_POLYMERASE_II_TRANSCRIPTION_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING"]],prot.list2) + labs(title="RNA_POLYMERASE_II_TRANSCRIPTION_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING")
plotEnrichment(all_gene_sets2[["RNA_POLYMERASE_II_TRANSCRIPTION_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n RNA_POLYMERASE_II_TRANSCRIPTION_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING (PATHWAY 1)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot2
plotEnrichment(all_gene_sets2[["DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY"]],prot.list2) + labs(title="DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY")
plotEnrichment(all_gene_sets2[["DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY (PATHWAY 2)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot3
plotEnrichment(all_gene_sets2[["CIS_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING"]],prot.list2) + labs(title="CIS_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING")
plotEnrichment(all_gene_sets2[["CIS_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n CIS_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING (PATHWAY 3)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot4
plotEnrichment(all_gene_sets2[["PARKINSONS_DISEASE"]],prot.list2) + labs(title="PARKINSONS_DISEASE")
plotEnrichment(all_gene_sets2[["PARKINSONS_DISEASE"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n PARKINSONS_DISEASE (PATHWAY 4)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Top 4 down

#Plot20
plotEnrichment(all_gene_sets2[["MATRISOME"]],prot.list2) + labs(title="MATRISOME")
plotEnrichment(all_gene_sets2[["MATRISOME"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n MATRISOME (PATHWAY 20 - TOP DOWN REGULATED)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot19
plotEnrichment(all_gene_sets2[["SELENOAMINO_ACID_METABOLISM"]],prot.list2) + labs(title="SELENOAMINO_ACID_METABOLISM")
plotEnrichment(all_gene_sets2[["SELENOAMINO_ACID_METABOLISM"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n SELENOAMINO_ACID_METABOLISM (PATHWAY 19 - TOP DOWN REGULATED)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot17
plotEnrichment(all_gene_sets2[["APOPTOSIS_BY_EPOXOMICIN_UP"]],prot.list2) + labs(title="APOPTOSIS_BY_EPOXOMICIN_UP")
plotEnrichment(all_gene_sets2[["APOPTOSIS_BY_EPOXOMICIN_UP"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n APOPTOSIS_BY_EPOXOMICIN_UP (PATHWAY 17 - TOP DOWN REGULATED)",
    x = "Rank", 
    y = "Enrichment Score"
  ) +
  theme(
    # hjust=0.5 centers, face="bold" bolds, lineheight adds the space
    plot.title = element_text(hjust = 0.5, face = "bold", lineheight = 1.5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

#Plot13
plotEnrichment(all_gene_sets2[["EUKARYOTIC_TRANSLATION_ELONGATION"]],prot.list2) + labs(title="EUKARYOTIC_TRANSLATION_ELONGATION")
plotEnrichment(all_gene_sets2[["EUKARYOTIC_TRANSLATION_ELONGATION"]], prot.list2) + 
  labs(
    title = "GSEA_OGT KD 605 10 MIN vs. GFP 10 MIN \n EUKARYOTIC_TRANSLATION_ELONGATION (PATHWAY 13 - TOP DOWN REGULATED)",
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

library(org.Hs.eg.db)
#install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
head(prot.list2)

head(prot.list2)
type(prot.list2)


#GSEA 605 PROT LIST


library(enrichplot)

head(fgseaRes3)
head(fgseaRes_unique2)
#listpathway<-fgseaRes[,c("pathway","padj","ES")]
#head(listpathway)
#topPathwaysUp1 <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
#listpathway <- listpathway[order(listpathway$padj),] #to sort based on padj
#head(listpathway)
#listpathway<-as.matrix(listpathway)

topPathwaysUp1 <- fgseaRes_unique2[ES > 0][head(order(padj), n=10)]
topPathwaysUp1<-topPathwaysUp1[,-2]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-3]
topPathwaysUp1
topPathwaysUp1<-topPathwaysUp1[,-(4:5)]
topPathwaysUp1
topPathwaysUp1<-as.data.frame(topPathwaysUp1)

topPathwaysDown1 <- fgseaRes_unique2[ES < 0][head(order(padj), n=10)]
topPathwaysDown1<-topPathwaysDown1[,-2]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-3]
topPathwaysDown1
topPathwaysDown1<-topPathwaysDown1[,-(4:5)]
topPathwaysDown1
topPathwaysDown1<-as.data.frame(topPathwaysDown1)


topPathways1 <- rbind(topPathwaysUp1, rev(topPathwaysDown1)) #merge rows
topPathways1
topPathways1<-as.data.frame(topPathways1)
topPathways1

#https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package
#ggplot(topPathways1, aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=padj<0.05&NES>0)) +
#  coord_flip() +
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="GSEA TOP 10 UP AND DOWN REGULATED PATHWAYS \n OGT KD 606 10 MIN. vs GFP 10 MIN.") + 
#  theme_minimal()+
#  theme(text = element_text(face = "bold"))

#Pathways with underscore

# 1. Clean the data frame FIRST to fix the "2" and underscores
topPathways1 <- topPathways1 %>%
  mutate(
    # Fix the specific "2" pathway to "Module 2"
    pathway = ifelse(pathway == "2" | pathway == "module_2", "Module 2", pathway),
    # Remove all underscores from all pathway names
    pathway = gsub("_", " ", pathway)
  )

topPathways1
ggplot(topPathways1, aes(reorder(pathway, NES), NES)) +
  # Removed color="black" to take away the borders around the bars
  geom_col(aes(fill = NES > 0)) + 
  scale_fill_manual(
    values = c("TRUE" = "#C21E56", "FALSE" = "#008080"), 
    labels = c("TRUE" = "Upregulated", "FALSE" = "Downregulated"),
    name = "Direction"
  ) +
  coord_flip() +
  labs(
    x = "Pathway", 
    y = "Normalized Enrichment Score (NES)",
    title = "GSEA: Top 10 Up and Down Regulated Pathways \n OGT KD 605 10 MIN. vs GFP 10 MIN."
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
##########################################################################################

#========================#
#J_HEATMAP
#========================#

#--------------------------------#
#OGT KD 606 10 MIN VS GFP 10 MIN
#--------------------------------#

head(prot.list)
head(sortupdown.ROGTTEN.GFPTEN)

# 1. Filter for adjusted p-value < 0.05
sig_only <- sortupdown.ROGTTEN.GFPTEN[sortupdown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN < 0.045, ]
head(sig_only)
tail(sig_only)
# 2. Top 20 Up of sorted list
top_20_up_all <- sig_only[order(as.numeric(sig_only$logFC.ROGTTEN.GFPTEN), decreasing = TRUE), ]
top_20_up <- head(top_20_up_all, 20)
top_20_up

# 3. Top 20 Down (Bottom of the list)
top_20_down_all <- sig_only[order(as.numeric(sig_only$logFC.ROGTTEN.GFPTEN), decreasing = FALSE), ]
top_20_down <- head(top_20_down_all, 20)
top_20_down
# 4. Combining them
heatmap_df <- rbind(top_20_up, top_20_down)

plot_mat <- cbind(
  LogFC = heatmap_df$logFC.ROGTTEN.GFPTEN,
  Significance = heatmap_df$adj.P.Val.ROGTTEN.GFPTEN
)
head(plot_mat)
rownames(plot_mat) <- rownames(heatmap_df)
plot_mat

# 1. Create an empty matrix with the correct dimensions
display_mat <- matrix("", nrow = nrow(plot_mat), ncol = ncol(plot_mat))
display_mat
# 2. Fill the first column (LogFC)
display_mat[, 1] <- as.character(round(plot_mat[, 1], 2))
#plot_mat[, 1]: This looks at the first column of your numeric data matrix (plot_mat).
#round(..., 2): It rounds those numbers to 2 decimal places (e.g., 0.12345 becomes 0.12).
#as.character(...): It converts those numbers into text (strings). In R, you can't display actual numbers inside a pheatmap cell unless they are treated as text labels.
#display_mat[, 1] <-: It saves these rounded text labels into the first column of a new matrix called display_mat.

# 3. Fill the second column (Significance) 
display_mat[, 2] <- formatC(plot_mat[, 2], format = "e", digits = 2)

display_mat

# 1. Create a matrix for colors where both columns ARE the LogFC
color_mat <- cbind(plot_mat[, 1], plot_mat[, 1]) 
color_mat

# 2. Creating the Heatmap

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
         main = ("OGT KD 606 10 MIN VD GFP 10 MIN"))

pheatmap(color_mat, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE, 
         # This function reorders the branches by LogFC values
         clustering_callback = function(hc, mat){
           # Calculate weights based on the first column (LogFC)
           wts <- mat[,1] 
           # Reorder the dendrogram: decreasing = TRUE puts high values at the top
           dend <- stats::reorder(stats::as.dendrogram(hc), wts = wts, agglo.FUN = mean)
           # Reversing may be needed depending on default orientation
           stats::as.hclust(rev(dend)) 
         },
         display_numbers = display_mat, 
         number_color = "black", 
         color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(256), 
         breaks = seq(-max(abs(plot_mat[,1])), max(abs(plot_mat[,1])), length.out = 257), 
         labels_col = c(expression(bold("LogFC")), expression(bold("Adj. P-Value"))), 
         angle_col = 0, 
         main = "Total-Proteomics Heatmap - OGT KD 606 10 min. vs. GFP 10 min.")

#--------------------------------#
#OGT KD 605 10 MIN VS GFP 10 MIN
#--------------------------------#

head(prot.list2)
head(sortupdownOGTTENGFPTEN)

# 1. Filter for adjusted p-value < 0.05
sig_only <- sortupdownOGTTENGFPTEN[sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.OGTCTRL< 0.045, ]
head(sig_only)
# 2. Get top 20 Up of sorted list
top_20_up_all <- sig_only[order(as.numeric(sig_only$logFC.OGTTEN.GFPTEN), decreasing = TRUE), ]
top_20_up <- head(top_20_up_all, 20)
top_20_up
# 3. Get top 20 Down (Bottom of the list)
top_20_down_all <- sig_only[order(as.numeric(sig_only$logFC.OGTTEN.GFPTEN), decreasing = FALSE), ]
top_20_down <- head(top_20_down_all, 20)
top_20_down
# 4. Combine them
heatmap_df <- rbind(top_20_up, top_20_down)
heatmap_df
plot_mat <- matrix(c(as.numeric(heatmap_df$logFC.OGTTEN.GFPTEN), 
                     as.numeric(heatmap_df$adj.P.Val.OGTTEN.GFPTEN)), 
                   ncol = 2)
head(plot_mat)
colnames(plot_mat) <- c("LogFC", "Significance")
rownames(plot_mat) <- rownames(heatmap_df)
plot_mat

# 1. Create an empty matrix with the correct dimensions
display_mat <- matrix("", nrow = nrow(plot_mat), ncol = ncol(plot_mat))
display_mat
# 2. Fill the first column (LogFC)
display_mat[, 1] <- as.character(round(plot_mat[, 1], 2))

#plot_mat[, 1]: This looks at the first column of your numeric data matrix (plot_mat).
#round(..., 2): It rounds those numbers to 2 decimal places (e.g., 0.12345 becomes 0.12).
#as.character(...): It converts those numbers into text (strings). In R, you can't display actual numbers inside a pheatmap cell unless they are treated as text labels.
#display_mat[, 1] <-: It saves these rounded text labels into the first column of a new matrix called display_mat.

# 3. Fill the second column (Significance) 
display_mat[, 2] <- formatC(plot_mat[, 2], format = "e", digits = 2)

display_mat

# 1. Create a matrix for colors where both columns ARE the LogFC
color_mat <- cbind(plot_mat[, 1], plot_mat[, 1]) 
color_mat

# 2. Creating the Heatmap

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
         main = ("OGT KD 605 10 MIN VD GFP 10 MIN"))

pheatmap(color_mat, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE, 
         # This function reorders the branches by LogFC values
         clustering_callback = function(hc, mat){
           # Calculate weights based on the first column (LogFC)
           wts <- mat[,1] 
           # Reorder the dendrogram: decreasing = TRUE puts high values at the top
           dend <- stats::reorder(stats::as.dendrogram(hc), wts = wts, agglo.FUN = mean)
           # Reversing may be needed depending on default orientation
           stats::as.hclust(rev(dend)) 
         },
         display_numbers = display_mat, 
         number_color = "black", 
         color = colorRampPalette(c("#0072B2", "white", "#E69F00"))(256), 
         breaks = seq(-max(abs(plot_mat[,1])), max(abs(plot_mat[,1])), length.out = 257), 
         labels_col = c(expression(bold("LogFC")), expression(bold("Adj. P-Value"))), 
         angle_col = 0, 
         main = "Total-proteomics Heatmap - OGT KD 605 10 min. vs. GFP 10 min.")



##########################################################################################
#J_PPI NETWROK ANALYSIS
#=========================#

head(sortupdownOGTTENGFPTEN)

# 1. Filter for adjusted p-value < 0.05
sig_only <- sortupdownOGTTENGFPTEN[sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN < 0.044, ]
head(sig_only)

# 1. Force the logFC column to be numeric (this fixes the error)
sortupdownOGTTENGFPTEN$logFC.numeric <- as.numeric(as.character(sortupdownOGTTENGFPTEN$logFC.OGTTEN.GFPTEN))

# 2. Apply both filters (Adjusted P < 0.044 AND absolute LogFC > 0.58)
sig_only <- sortupdownOGTTENGFPTEN[
  !is.na(sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN) & # Remove NAs
    sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN < 0.044 & 
    abs(sortupdownOGTTENGFPTEN$logFC.numeric) > 0.58, 
]

# 3. Check how many proteins are left
nrow(sig_only)
head(sig_only)
# Check the count of proteins remaining
nrow(sig_only)
head(sig_only)
# 2. Get the list of significant gene names
genes_for_ppi <- rownames(sig_only)
head(genes_for_ppi)

library(STRINGdb)
# Initialize for Human (9606) with a medium confidence score (400)
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=700)
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
string_db$plot_network(hits)

#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("STRINGdb", update = TRUE, ask = FALSE)

# 1. Get the graph object locally
graph <- string_db$get_graph()
library(igraph)
# 2. Extract the subgraph for your 'hits'
# Ensure 'hits' is a character vector of STRING IDs
subgraph <- igraph::induced_subgraph(graph, vids = V(graph)[name %in% hits])

# 3. Plot it using standard R graphics
plot(subgraph, vertex.label = V(subgraph)$name, vertex.size = 5)

# 1. Map the Gene Symbols to the graph nodes
# This creates a named vector: names = STRING IDs, values = Symbols
symbol_map <- setNames(prot_mapped$Symbol, prot_mapped$STRING_id)

# 2. Assign these symbols as a new attribute to your subgraph
V(subgraph)$label <- symbol_map[V(subgraph)$name]

# 3. Clean up the plot settings for a more professional look

  
library(igraph)

# 1. Use the 'Fruchterman-Reingold' layout with more iterations for better spacing
# niter = 1000 gives the algorithm more time to push nodes apart
l <- layout_with_fr(subgraph, niter = 1000)

# 2. Plot with ultra-fine settings for a professional look
plot(subgraph, 
     layout = l,
     vertex.label = NA,             # No text
     vertex.size = 1.0,             # Very small dots
     vertex.color = "#1f77b4",      # Professional blue
     vertex.frame.color = NA,       # No white borders
     edge.width = 0.1,              # Very thin lines
     edge.color = adjustcolor("gray", alpha.f = 0.2), # Transparent lines
     main = "Systems-Level Interaction Network: OGT Knockdown 605"
)

# 1. Create a color mapping based on LogFC
# Everything > 0.58 is Red, everything < -0.58 is Blue
node_colors <- ifelse(sig_only$logFC.numeric > 0, "#d7191c", "#2c7bb6")

# 2. Assign these colors to the graph nodes
# Match the order of your 'hits' with the 'sig_only' dataframe
V(subgraph)$color <- node_colors[match(V(subgraph)$name, prot_mapped$STRING_id)]

# 3. Plot with the new color scheme
library(igraph)
l <- layout_with_fr(subgraph, niter = 1000)

plot(subgraph, 
     layout = l,
     vertex.label = NA,             # Keep it clean
     vertex.size = 2,               # Slightly larger to see the colors
     vertex.color = V(subgraph)$color, 
     vertex.frame.color = NA,       # Cleaner look
     edge.width = 0.15, 
     edge.color = "gray85",         # Light gray so the colors pop
     main = "OGT Knockdown 605: Directional Protein Interactions",
     sub = "Red = Up-regulated (LogFC > 0.58), Blue = Down-regulated (LogFC < -0.58)")

library(igraph)

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

# 4. Plot the cleaned, labeled central hub
l <- layout_with_fr(giant_component, niter = 1000)

plot(giant_component, 
     layout = l,
     vertex.label = V(giant_component)$label, 
     vertex.label.cex = 0.6,              # Small, readable text
     vertex.label.color = "black",        # High contrast
     vertex.label.dist = 0.6,             # Push labels slightly away from dots
     vertex.size = 3.5,                   # Distinct nodes
     vertex.color = V(giant_component)$color, # Red (Up) / Blue (Down)
     vertex.frame.color = "gray30",       # Clean node borders
     edge.width = 0.3, 
     edge.color = "gray85",               # Light edges to keep the focus on nodes
     main = "Central Interaction Hub: OGT Knockdown",
     sub = "Largest Connected Component (Outliers Removed)")

# 2. Map Gene Symbols to the remaining nodes
# We use the 'Symbol' column we created earlier in 'prot_mapped'
V(giant_component)$label <- prot_mapped$Symbol[match(V(giant_component)$name, prot_mapped$STRING_id)]

# 3. Plot the final, labeled network
# Fruchterman-Reingold layout helps spread the labels for readability
l <- layout_with_fr(giant_component, niter = 1000)

plot(giant_component, 
     layout = l,
     vertex.label = V(giant_component)$label, 
     vertex.label.cex = 0.6,             # Shrink text to prevent overcrowding
     vertex.label.color = "black",       # High contrast for symbols
     vertex.label.font = 2,              # Bold symbols for better visibility
     vertex.size = 3,                    # Slightly larger dots for the main hub
     vertex.color = V(giant_component)$color, # Keeps your Red/Blue directionality
     vertex.frame.color = "gray40",      # Soft border for nodes
     edge.width = 0.2, 
     edge.color = "gray80", 
     main = "Central Interaction Hub: OGT Knockdown",
     sub = "Outliers removed; Red = Up-regulated, Blue = Down-regulated")

library(igraph)

# Spacing out the nodes for a cleaner look
#l <- layout_with_kk(giant_component)

# Plotting with optimized aesthetics
#plot(giant_component, 
#    layout = l,
#    vertex.label = V(giant_component)$label, 
#    vertex.label.cex = 0.5,             # Smaller text for high density
#    vertex.label.font = 2,              # Bold symbols
#    vertex.label.color = "black", 
#    vertex.label.dist = 0.4,            # Pushes text away from the dot
#    vertex.size = 4,                    # Larger dots for visibility
#    vertex.color = V(giant_component)$color, 
#    vertex.frame.color = "white",       # White border makes dots "pop"
#    vertex.frame.width = 0.5,
#    edge.width = 0.2, 
#    edge.color = adjustcolor("gray", alpha.f = 0.3), # Subtle edges
#    main = "Central Hub: OGT Knockdown (Zoomed)",
#    margin = -0.1)                      # Negative margin "zooms in" the plot

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("RCy3")
library(RCy3)

cytoscapePing() # Confirms R can "see" the running Cytoscape app


createNetworkFromIgraph(giant_component, 
                        title = "OGT_Knockdown_Giant_Hub", 
                        collection = "STRING_Analysis")
# 1. Set the background color
setCanvasBackgroundColor("#FFFFFF")
library(RCy3)
setCanvasBackgroundColor('#FFFFFF')

setVisualPropertyDefault(list(visualProperty = "NETWORK_BACKGROUND_PAINT", value = "#FFFFFF"))

# 2. Map Node Colors (based on the 'color' attribute you created in igraph)
# Note: If 'color' didn't carry over perfectly, you can map by logFC.numeric
setNodeColorMapping('color', colors = unique(V(giant_component)$color), mapping.type = 'passthrough')

# 3. Adjust Node and Label sizes
setNodeSizeDefault(35)
setNodeLabelMapping('label')
setNodeFontSizeDefault(12)

# 4. Use a better layout (Cytoscape's Force-Directed is usually best)
layoutNetwork('force-directed')

#==============================#
#OGTKD 606 - Netwrok Analysis
#==============================#

head(sortupdown.ROGTTEN.GFPTEN)

# 1. Filter for adjusted p-value < 0.05
sig_only <- sortupdown.ROGTTEN.GFPTEN[sortupdown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN < 0.044, ]
head(sig_only)

# 1. Force the logFC column to be numeric (this fixes the error)
sortupdown.ROGTTEN.GFPTEN$logFC.numeric <- as.numeric(as.character(sortupdown.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN))
head(sortupdown.ROGTTEN.GFPTEN)

# 2. Apply both filters (Adjusted P < 0.044 AND absolute LogFC > 0.58)
sig_only <- sortupdown.ROGTTEN.GFPTEN[
  !is.na(sortupdown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN) & # Remove NAs
    sortupdown.ROGTTEN.GFPTEN$adj.P.Val.ROGTTEN.GFPTEN < 0.044 & 
    abs(sortupdown.ROGTTEN.GFPTEN$logFC.numeric) > 0.58, 
]

# 3. Check how many proteins are left
nrow(sig_only)
head(sig_only)

# Get the list of significant gene names
genes_for_ppi <- rownames(sig_only)
head(genes_for_ppi)

library(STRINGdb)
# Initialize for Human (9606) with a medium confidence score (400)
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=700)
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
                        title = "OGT_Knockdown 606_Giant_Hub v3", 
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

# 2. Make the labels bold to make them stand out
setNodeFontFaceDefault("Arial-Bold")

# 1. Change the global edge color to a light-medium gray
setEdgeColorDefault('#000000')

setVisualPropertyDefault(list(visualProperty="NODE_LABEL_FONT_SIZE", value=20))

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


#==============================#
#OGTKD 605 - Netwrok Analysis
#==============================#
head(sortupdown.ROGTTEN.GFPTEN)
head(sortupdownOGTTENGFPTEN)

# 1. Filter for adjusted p-value < 0.05
sig_only <- sortupdownOGTTENGFPTEN[sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN < 0.044, ]
head(sig_only)
tail(sig_only)
# 1. Force the logFC column to be numeric (this fixes the error)
sortupdownOGTTENGFPTEN$logFC.numeric <- as.numeric(as.character(sortupdownOGTTENGFPTEN$logFC.OGTTEN.GFPTEN))
head(sortupdownOGTTENGFPTEN)

# 2. Apply both filters (Adjusted P < 0.044 AND absolute LogFC > 0.58)
sig_only <- sortupdownOGTTENGFPTEN[
  !is.na(sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN) & # Remove NAs
    sortupdownOGTTENGFPTEN$adj.P.Val.OGTTEN.GFPTEN < 0.044 & 
    abs(sortupdownOGTTENGFPTEN$logFC.numeric) > 0.58, 
]

head(sig_only)
tail(sig_only)

# 3. Check how many proteins are left
nrow(sig_only)
head(sig_only)

# Get the list of significant gene names
genes_for_ppi <- rownames(sig_only)
head(genes_for_ppi)

library(STRINGdb)
# Initialize for Human (9606) with a medium confidence score (400)
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=700)
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

cytoscapePing() # Confirms R can "see" the running Cytoscape app # Should return "You are connected to Cytoscape!"
#layoutNetwork('force-directed')
#fitContent()

createNetworkFromIgraph(giant_component, 
                        title = "OGT_Knockdown 605_Giant_Hub v1", 
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

# 2. Make the labels bold to make them stand out
setNodeFontFaceDefault("Arial-Bold")

# 1. Change the global edge color to a light-medium gray
setEdgeColorDefault('#000000')

setVisualPropertyDefault(list(visualProperty="NODE_LABEL_FONT_SIZE", value=20))

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
                   sizes = c(100, 20, 100), 
                   mapping.type = 'c')

setNodeShapeDefault("ELLIPSE")

# Fit everything back into the window view
fitContent()







#===================
# 2. Map Node Colors (based on the 'color' attribute you created in igraph)
# Note: If 'color' didn't carry over perfectly, you can map by logFC.numeric
setNodeColorMapping('color', colors = unique(V(giant_component)$color), mapping.type = 'passthrough')

# 3. Adjust Node and Label sizes
setNodeSizeDefault(35)
setNodeLabelMapping('label')
setNodeFontSizeDefault(12)

# 4. Use a better layout (Cytoscape's Force-Directed is usually best)
layoutNetwork('force-directed')


#GSEA 606 PROT LIST
head(sortupdown.ROGTTEN.GFPTEN)
prot.list <- sortupdown.ROGTTEN.GFPTEN$logFC.ROGTTEN.GFPTEN             # rank-ordered gene list
head(prot.list)
prot.list<-as.numeric(prot.list)
head(prot.list)
head(sortupdown.ROGTTEN.GFPTEN)
names(prot.list) <- rownames(sortupdown.ROGTTEN.GFPTEN)
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















#J_AMEND

#============#
# PPI Network ----
#============#
library(igraph)
library(data.table)

# Download and pre-process protein-protein interaction (PPI) network from STRING
# Website: https://string-db.org/ 
#   - Go to Downloads, then type in an organism (e.g., Homo sapiens, Mus musculus)
#   - Choose either 'protein.links' file (includes both physical and functional PPIs) or 'protein.physical.links' (includes only physical PPIs)

#######################################################################
#sample code
x <- data.frame(
  from = c("A", "A", "B", "C"),
  to = c("B", "C", "C", "D")
)
x

library(igraph)

# Convert the data frame to a matrix (edges are in columns 1 and 2)
edge_matrix <- as.matrix(x[, 1:2])
edge_matrix
# Create an undirected graph
g <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)

# Print the graph
print(g)
#########################################################################

# Helper functions
remove_duplicate_edges = function(x){
  
  ## Create a graph from the edge list (x) without directionality
 
   g = igraph::graph_from_edgelist(as.matrix(x[,1:2]), directed = FALSE)
  
  #edgelist: A matrix or data frame where each row represents an edge. The first 
  #column is the source vertex, and the second column is the target vertex. 
  #If the graph is directed, edges will point from the source to the target.
  #A graph is a slightly abstract represention of some objects that are related 
  #to each other in some way, and of these relationships. The objects are drawn 
  #as dots called vertices (nodes); a line (or curve) connects any two vertices 
  #representing objects that are related or adjacent; such a line is called an edge.
  
  #Vertex=dots = nodes; lines = edges
  
   # Create an edge list
   edgelist <- matrix(c(1, 2,
                        2, 3,
                        3, 4,
                        4, 1), 
                      ncol = 2, byrow = TRUE)
   edgelist
   #row = edges, col=nodes/vertices
   #edgelist: A matrix or data frame where each row represents an edge. 
   #The first column is the source vertex, and the second column is the target vertex.
   #If the graph is directed, edges will point from the source to the target.
   
   # Create an undirected graph from the edge list
   g1 <- graph_from_edgelist(edgelist, directed = FALSE)
   
   # Plot the graph
   plot(g1)
   
  # Assign weights to the edges
  igraph::E(g)$weight = as.numeric(x[,3])
  ## Simplify the graph: merge duplicate edges and use the minimum weight
  g = igraph::simplify(graph = g, edge.attr.comb = list(weight = "min"))
  # Extract the simplified edge list
  el = igraph::as_edgelist(graph = g)
  # Add weights to the edge list
  el = cbind(el, igraph::E(g)$weight)
  return(el)
}
largest_connected_component = function(g){
  # Get the largest connected component - check if the graph is connected. If the graph is not connected, issue a warning.
  if(!igraph::is_connected(g)){
    message("Taking largest connected component.")
    #Computes the connected components of the graph.
    comps = igraph::components(g)
    largest_comp_id = which.max(comps$csize)
    # Get the largest connected component
    g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  g
}

# Load in full PPIN as an edge list
edge.threshold = 0.7
#path.2.network = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\J_AMEND\\9606.protein.links.v12.0.txt"
#path.2.network = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\J_AMEND\\10090.protein.links.v12.0.txt"
path.2.network = "C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\J_AMEND\\New folder\\9606.protein.links.v12.0.txt"

ppi = data.table::fread(file = path.2.network, header = TRUE) %>%
  dplyr::mutate(combined_score = combined_score / 1000) %>% # Divide edge scores by 1000 so that they are between 0 and 1
  dplyr::filter(combined_score >= edge.threshold) %>% # Threshold edge scores 
  dplyr::mutate(protein1 = extract_string(protein1, "\\.", 2), # Remove the taxonomy ID '10090'
                protein2 = extract_string(protein2, "\\.", 2)) %>% 
  as.matrix() %>%
  remove_duplicate_edges()

head(ppi)

# Making igraph object, taking largest connected component
g = graph_from_edgelist(ppi[,1:2], directed = FALSE)
E(g)$weight = as.numeric(ppi[,3])
g = largest_connected_component(g)

head(g)

# Mapping ENSEMBL peptide IDs from graph to gene symbol 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
# Pick which species you need
#species = c('mmusculus', 'hsapiens')[1]
species = c('hsapiens')[1]
mart_data <- useDataset(paste0(species, "_gene_ensembl"), mart = mart)
return.type = ifelse(species == 'hsapiens', 'hgnc_symbol', 'mgi_symbol') # If you want to map ensembl peptide IDs to Gene Symbols
mapping <- getBM(attributes = c("ensembl_peptide_id", return.type),
                 filters = "ensembl_peptide_id",
                 values = unique(V(g)$name),
                 mart = mart_data)
mapping = mapping[apply(mapping,1,function(x) all(x != "")), ] 

head(mapping)

# Check relationship of mapping (e.g., many-to-one, many-to-many, etc.)
table(table(mapping[,1]))
table(table(mapping[,2]))

# Check quality of mapping
#mean(de$tt$Symbol %in% mapping[,2]) # 91 % 
mean(dea$tt$Symbol %in% mapping[,2]) # 89.67 % 
mean(V(g)$name %in% mapping[,1]) # 97 %

# Match vertex names to the Ensembl peptide IDs in the mapping
id = match(V(g)$name, mapping$ensembl_peptide_id)

# Assign 'Symbol' attribute to vertices with valid matches
# which(!is.na(id)) gives the indices of vertices with valid IDs
# mapping$mgi_symbol[id[!is.na(id)]] gives the corresponding symbols from the mapping data frame
vertex_attr(g, "Symbol", which(!is.na(id))) = mapping$hgnc_symbol[id[!is.na(id)]]

head(id)




#====================#
# Map Data to Network ----
#====================#
# Map omics data to nodes in the network
#col.ids.to.network = c(2,3,4,6,8,10)
col.ids.to.network = c(2,3,4,5,6,7,8,10,12,14,16,18,20,22)
for(i in seq_along(col.ids.to.network)){
  id = match(V(g)$Symbol, dea$tt$Symbol)
  igraph::vertex_attr(g, colnames(dea$tt)[col.ids.to.network[i]], which(!is.na(id))) = dea$tt[id[!is.na(id)],col.ids.to.network[i]]
}

head(g)
head(de$tt)
head(dea$tt)
#=================#
# Network Analysis ----
#=================#
#install.packages("devtools")
#devtools::install_github("samboyd0/AMEND")
library(AMEND)
#devtools::install_github("samboyd0/AMEND", build_vignettes = TRUE)

# Perform AMEND to get subnetworks of nodes with large experimental values
# ?run_AMEND
# FUN:
#   shift_scale: For values centered about 0. Takes absolute value then down-weights values not in direction of interest (DOI) by w. 
#       Parameters: DOI, w
#   exp: For log fold changes. exp(DOI * x) for DOI=1,-1, and exp(abs(DOI) * x) for DOI = 0
#       Parameters: DOI
#   p_value: -log10(x + 1e-6)

# Vertex attributes of PPIN
igraph::vertex_attr_names(g)
# Pick 'data.name' to be one of these vertex attributes

#== Log fold change example ==#
data.name = "ROGTTEN.GFPTEN"
head(data.name)
# Searching for both up- and down-regulated nodes
subnet = run_AMEND(graph = g, n = 50, FUN = "exp", FUN.params = list(DOI = 0), data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)
# Searching for only up-regulated nodes, in 2 different ways
subnet = run_AMEND(graph = g, n = 50, FUN = "exp", FUN.params = list(DOI = 1), data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)
subnet = run_AMEND(graph = g, n = 50, FUN = "shift_scale", FUN.params = list(DOI = 1, w = 0.5), data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)
head(subnet)
## Inspect output from run_AMEND
# igraph object
subnet$module
# size
vcount(subnet$module) #57
# mean data value
mean(vertex_attr(subnet$module, data.name))
# Ensembl peptide IDs of nodes in module
V(subnet$module)$name
# Symbols of nodes in module
V(subnet$module)$Symbol

head(dea$tt)
#== P-value example ==#
#data.name = "P.Value.ROGTTEN.GFPTEN"
#subnet = run_AMEND(graph = g, n = 50, FUN = "p_value", data = data.name, 
      #             normalize = "modified_degree", verbose = TRUE)

## Inspect output from run_AMEND
# igraph object
subnet$module
# size
vcount(subnet$module)
# mean data value
mean(vertex_attr(subnet$module, data.name))
# Ensembl peptide IDs of nodes in module
V(subnet$module)$name
# Symbols of nodes in module
V(subnet$module)$Symbol

#=============================#
# Over-representation Analysis ----
#=============================#
library(fgsea)

#https://rest.kegg.jp/list/pathway/hsa

hsa1 = read.delim("C:\\Users\\sophi\\OneDrive\\Desktop\\J_Dr. Slawson projects_ 2023\\J_ERK MS\\J_TOTAL PROTEOME_R\\J_PROTEOMICS R CODE\\J_AMEND\\hsakeggMetabolicpathways.txt")
head(hsa1)
saveRDS(hsa, "hsa1.rds")
#hsa <- readRDS("hsa1.rds")
#head(hsa)

#path.to.pathways = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/Subnetwork Data/m_pathways.RDS"
#pathway.list = readRDS(path.to.pathways)

pathway.list <- readRDS("hsa1.rds")
head(pathway.list)
#genekegg_6 <- genekegg_6[!duplicated(names(genekegg_6))]
pathway.list <- pathway.list[!duplicated(names(pathway.list))]
head(pathway.list)
genes<-V(subnet$module)$Symbol
genes
genes <- genes[!duplicated(names(genes))]
head(genes)

ora.paths = fora(pathways = pathway.list, # List of pathways
                 genes = V(subnet$module)$Symbol, # Gene symbols from AMEND module
                 universe = V(g)$Symbol) %>% # Gene symbols from full graph
  dplyr::mutate(overlapGenes = unlist(lapply(overlapGenes, function(x) paste(x, collapse = ", ")))) %>%
  dplyr::filter(padj <= 0.05) # Filter out pathways with adjusted p-value greater than 0.05

head(ora.paths)

#=====================================#
# Visualize AMEND Results in Cytoscape ----
#=====================================#
library(RCy3)
library(RColorBrewer)

#' @title Construct a graph of nested pathways
#' 
#' @description
#' Given a list of pathways resulting from functional enrichment analysis, construct a directed, acyclic graph (DAG) where nodes represent pathways and directed edges go from smaller pathways to larger pathways in which the smaller pathway is partially nested.
#' 
#' @details
#' For two sets A and B, the Nested Index is given by f(A,B) = |intersect(A,B)| / min(|A|,|B|)
#' 
#' In the context of AMEND analysis, which returns a module, we can identify 'representative pathways' that are associated with clusters in the module. First, we run a topological clustering algorithm on the AMEND module to get clusters. Then, for each cluster, we determine a representative pathway by choosing the pathway containing the most nodes in that cluster.
#' 
#' By using 'rep.path.ids', we assign colors to nodes based on whether they are a representative pathway or not.
#' 
#' This may take a very long time to run for inputs of greater than ~30 pathways. In these scenarios, use first.n, or just modify pathway data.frame to contain fewer pathways prior to input. Also, running with parallel=TRUE will help.
#' 
#' @param enriched.pathways data.frame of pathway names along with adjusted and un-adjusted p-values
#' @param rep.path.ids Row IDs of 'enriched.pathways' corresponding to the representative pathways
#' @param all.pathways List of all pathways used in Pathway Analysis
#' @param threshold numeric value between 0 and 1. Threshold for considering one pathway nested within another
#' @param first.n Integer value or NULL. Set to non-NULL if you only want to consider the top n pathways (using the ordering given in 'enriched.pathways')
#' @param parallel logical. Set to TRUE if you want to run in parallel
#' @param nclust Integer value. Number of clusters to use for parallel computing
#' @param padj Logical. Whether to use adjusted p-value for 'pval' vertex attribute of return graph
#' 
#' @return igraph object
#' 
nested_pathways = function(enriched.pathways, rep.path.ids, all.pathways, threshold = 0.8, first.n = NULL, parallel = FALSE, nclust, padj = TRUE){
  require(doParallel)
  require(foreach)
  
  if(nrow(enriched.pathways) == 0){
    stop("nrow of \'enriched.pathways\' is zero")
  }else if(nrow(enriched.pathways) == 1){
    ora_path = enriched.pathways$pathway
    g = igraph::make_empty_graph(n = length(ora_path))
    vlabs = str_replace(ora_path, "Mus musculus: ", "") 
    V(g)$name = vlabs
    temp = numeric(vcount(g))
    for(i in seq_along(temp)){
      temp[i] = length(all.pathways[[ora_path[i]]])
    }
    V(g)$size = temp
    V(g)$pval = ifelse(padj, enriched.pathways$padj, enriched.pathways$pval)
    col.tmp = rep(0, igraph::vcount(g))
    col.tmp[rep.path.ids] = seq_along(rep.path.ids)
    V(g)$color = col.tmp
    return(g)
  }
  
  if(!is.null(first.n)) enriched.pathways = enriched.pathways[1:min(first.n, nrow(enriched.pathways)),]
  
  nested_index = function(a, b){
    intersection = length(intersect(a, b))
    denom = min(length(a), length(b))
    return(intersection / denom)
  }
  
  ora_path = enriched.pathways$pathway
  p.id1 = c()
  p.id2 = c()
  sim.mat = matrix(0, nrow = length(ora_path), ncol = length(ora_path), dimnames = list(ora_path, ora_path))
  for(i in 1:(length(ora_path)-1)){
    for(j in (i+1):length(ora_path)){
      sim.mat[i,j] = nested_index(all.pathways[[ora_path[i]]], all.pathways[[ora_path[j]]])
      sim.mat[j,i] = sim.mat[i,j]
      if(sim.mat[i,j] > threshold){ 
        p.id1 = c(p.id1, ifelse(min(length(all.pathways[[ora_path[i]]]), length(all.pathways[[ora_path[j]]])) == length(all.pathways[[ora_path[i]]]), i, j))
        p.id2 = c(p.id2, ifelse(min(length(all.pathways[[ora_path[i]]]), length(all.pathways[[ora_path[j]]])) == length(all.pathways[[ora_path[i]]]), j, i))
      }
    }
  }
  
  # Finding nested pathways 
  # nest.mat is an edgelist, where a vertex in col 1 is nested within a vertex of col 2 (vertex = pathway)
  if(is.null(p.id1) || is.null(p.id2)){
    g = igraph::make_empty_graph(n = length(ora_path))
    vlabs = str_replace(ora_path, "Mus musculus: ", "") 
    V(g)$name = vlabs
    temp = numeric(vcount(g))
    for(i in seq_along(temp)){
      temp[i] = length(all.pathways[[ora_path[i]]])
    }
    V(g)$size = temp
    V(g)$pval = ifelse(padj, enriched.pathways$padj, enriched.pathways$pval)
    col.tmp = rep(0, igraph::vcount(g))
    col.tmp[rep.path.ids] = seq_along(rep.path.ids)
    V(g)$color = col.tmp
    return(g)
  }else{
    nest.mat = matrix(c(p.id1, p.id2), ncol = 2)
    g = igraph::graph_from_edgelist(nest.mat, directed = TRUE)
    for(i in seq_along(setdiff(1:length(ora_path), V(g)))) g = igraph::add_vertices(g, 1)
  }
  
  # Getting longest simple, directed paths between pair of vertices in edgelist
  # Rationale: If a pathway is nested w/in a larger pathway, it may be nested in other pathways that are themselves nested in that larger pathway.
  #   I want to see the hierarchy of nested-ness
  if(parallel){ # parallel FOR loop
    cl = makeForkCluster(nclust, outfile = "")
    registerDoParallel(cl)
    nest.id0 = foreach(i = 1:nrow(nest.mat), .packages = "igraph") %dopar% {
      message(paste0(i, ":", nrow(nest.mat)))
      asp = all_simple_paths(graph = g, from = nest.mat[i,1], to = nest.mat[i,2], mode = "out")
      asp[[which.max(do.call("c", lapply(asp, length)))]]
    }
    stopCluster(cl)
  }else{ # normal FOR loop
    nest.id0 = vector("list", nrow(nest.mat))
    for(i in 1:nrow(nest.mat)){ # This for loop may take a long time
      message(paste0(i, ":", nrow(nest.mat)))
      asp = all_simple_paths(graph = g, from = nest.mat[i,1], to = nest.mat[i,2], mode = "out")
      nest.id0[[i]] = asp[[which.max(do.call("c", lapply(asp, length)))]]
      # if(i %% 50 == 0) message(paste0(i, ":", nrow(nest.mat)))
    }
  }
  
  # This gets rid of redundant paths
  nest.id = list()
  for(i in 1:length(nest.id0)){
    cond = T
    for(j in (1:length(nest.id0))[(1:length(nest.id0)) != i]){
      if(all(nest.id0[[i]] %in% nest.id0[[j]])) cond = F
    }
    if(cond) nest.id[[length(nest.id) + 1]] = nest.id0[[i]]
  }
  
  # Converting vertex ids to pathway names
  nested.paths = vector(mode = "list", length = length(nest.id))
  for(i in 1:length(nest.id)){
    nested.paths[[i]] = ora_path[nest.id[[i]]]
  }
  # Going from most nested to least nested (smallest to largest)
  
  # Getting graph showing nested hierarchy structure
  el = matrix(nrow = 1000, ncol = 2)
  r.id = 1
  for(i in seq_along(nest.id)){
    x = as.numeric(nest.id[[i]])
    for(j in 1:(length(x) - 1)){
      new.row = x[j:(j+1)]
      if(any(apply(el, 1, function(y) all(new.row %in% y)))) next
      el[r.id,] = new.row 
      r.id = r.id + 1
    }
  }
  el = na.omit(el)
  g = igraph::graph_from_edgelist(el = el, directed = TRUE)
  for(i in seq_along(setdiff(1:length(ora_path), V(g)))) g = igraph::add_vertices(g, 1)
  vlabs = str_replace(ora_path, "Mus musculus: ", "") 
  V(g)$name = vlabs
  temp = numeric(vcount(g))
  for(i in seq_along(temp)){
    temp[i] = length(all.pathways[[ora_path[i]]])
  }
  V(g)$size = temp
  V(g)$pval = enriched.pathways$padj
  col.tmp = rep(0, igraph::vcount(g))
  col.tmp[rep.path.ids] = seq_along(rep.path.ids)
  V(g)$color = col.tmp
  return(g)
}

#=== Upload AMEND module to Cytoscape ===#
createNetworkFromIgraph(subnet$module, title = "AMEND Module", collection = "")
# net.id = getNetworkSuid()

### Set the visual style
# Looking at color schemes
if(0){ 
  display.brewer.all()
  display.brewer.pal(n = 9, name = "RdBu")
  brewer.pal(n = 11, name = "RdBu")
}
color.theme = "RdBu"
# The above may need to change depending on if you want divergent or gradient color themes, which will depend on the data type used in AMEND

style.name = "AMEND_module"
createVisualStyle(style.name)

# Set node shape and size
lockNodeDimensions(new.state = TRUE, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SIZE", value = 40), style.name = style.name)

# Set node color as a function of Experimental data used in AMEND. This will need to be modified for different data types
color.scheme = brewer.pal(n = 3, name = color.theme)
v.attr.name = "TMG.saline" # The name of the vertex attribute used as seed values in AMEND
v.range = c(min(vertex_attr(g, v.attr.name)), max(vertex_attr(g, v.attr.name)))
setNodeColorMapping(table.column = v.attr.name, table.column.values = c(v.range[1], 0, v.range[2]), colors = color.scheme, style.name = style.name)
# Node labels, position, and font size
setNodeLabelMapping(table.column = "Symbol", style.name = style.name)
setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 28), style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = "AlBayan-Bold"), style.name = style.name)
# Set edge style and width
setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
setEdgeLineWidthDefault(new.width = 3, style.name = style.name)

setVisualStyle(style.name)

#=====================================#
# Visualizing ORA Results in Cytoscape ----
#=====================================#
# Get that pathways that each node in module is associated with
node.pathway = vector("list", vcount(subnet$module)); names(node.pathway) = V(subnet$module)$Symbol
for(i in seq_along(node.pathway)){
  node.pathway[[i]] = ora.paths$pathway[unlist(lapply(strsplit(ora.paths$overlapGenes, ", "), function(x) V(subnet$module)$Symbol[i] %in% x))]
}

## Getting cellular functions associated with each cluster in module
set.seed(2088) # Louvain is stochastic, so we set a seed value for reproducibility 
# Running Louvain clustering algorithm on AMEND module to get clusters
cl = igraph::cluster_louvain(graph = subnet$module, resolution = 1)
cl.id = unique(cl$membership)
module.cluster = character(length(cl.id))
# For each cluster, identify pathway that contains the most nodes of that cluster
for(j in seq_along(cl.id)){
  node.names = cl$names[which(cl$membership == cl.id[j])]
  tmp = unlist(node.pathway[which(V(subnet$module)$name %in% node.names)])
  if(length(tmp) == 0){ # If the selected nodes aren't part of any pathway...
    module.cluster[j] = "unknown"
  }else{
    p.count = table(tmp)
    tmp.id = which(p.count == max(p.count))
    # In case of ties for most occurring pathway, take largest pathway
    if(length(tmp.id) > 1) tmp.id = tmp.id[which.max(unlist(lapply(pathway.list[names(p.count)[tmp.id]], length)))]
    module.cluster[j] = names(p.count)[tmp.id]
  }
  names(module.cluster)[j] = paste(V(subnet$module)$Symbol[V(subnet$module)$name %in% node.names], collapse = ", ")
}
if(any(table(module.cluster) > 1)){ # Merge clusters that have same Representative Pathway
  tmp = table(module.cluster)
  tmp.id = which(tmp > 1)
  for(j in seq_along(tmp.id)){
    tmp.id2 = which(module.cluster %in% names(tmp)[tmp.id[j]])
    keep.id = tmp.id2[1]
    rm.id = tmp.id2[-1]
    names(module.cluster)[keep.id] = paste(c(names(module.cluster)[keep.id], names(module.cluster)[rm.id]), collapse = ", ")
    module.cluster = module.cluster[-rm.id]
  }
}
tmp = data.frame(.function = module.cluster, genes = names(module.cluster))
tmp = tmp[tmp$.function != "unknown",]
module.cluster = tmp

## Visualize pathways in Cytoscape
top.k.pathways = 15
threshold = 0.8
tkp = min(top.k.pathways, nrow(ora.paths))
p.tmp = ora.paths[1:tkp,]
rp.tmp = module.cluster$.function[!module.cluster$.function %in% p.tmp$pathway]
if(length(rp.tmp) == 0){
  rp.ids = which(p.tmp$pathway %in% module.cluster$.function)
}else{
  id.tmp = which(!p.tmp$pathway %in% module.cluster$.function)
  p.tmp = rbind2(p.tmp[-tail(id.tmp, length(rp.tmp)),], ora.paths[ora.paths$pathway %in% rp.tmp,])
  rp.ids = which(p.tmp$pathway %in% module.cluster$.function)
}
g.path = nested_pathways(enriched.pathways = p.tmp, rep.path.ids = rp.ids, all.pathways = pathway.list, threshold = threshold, first.n = NULL, parallel = TRUE, nclust = 4)
createNetworkFromIgraph(g.path, title = "", collection = "Nested Pathways Hierarchy")
layoutNetwork("hierarchical")

## Set Visual Style
style.name = "ORA_viz"
createVisualStyle(style.name)

# Set node shape and size
lockNodeDimensions(new.state = TRUE, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
setNodeSizeMapping(table.column = "size", table.column.values = c(10, 500, 1000), 
                   sizes = c(15, 55, 100), style.name = style.name)

# Set node color
gray = brewer.pal(n = 9, name = "Set1")[9]
col.tmp = rep(gray, igraph::vcount(g.path))
col.tmp[V(g.path)$color != 0] = brewer.pal(n = max(sum(V(g.path)$color != 0),3), name = "Set1")[1:sum(V(g.path)$color != 0)]

setNodeColorMapping(table.column = "color", table.column.values = V(g.path)$color, colors = col.tmp, mapping.type = "d", style.name = style.name)

# Node labels, position, and font size
setNodeLabelMapping(table.column = "name", style.name = style.name)
setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 20), style.name = style.name)
font_style = c("Rockwell-Bold", "ArialNarrow-Bold", "Arial-BoldMT", "Arial-Black")[4]
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = font_style), style.name = style.name)
# Set edge style and width
setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
setEdgeLineWidthDefault(new.width = 4, style.name = style.name)
setEdgeTargetArrowShapeDefault(new.shape = "DELTA", style.name = style.name)
setVisualStyle(style.name)

### Select Nodes in AMEND module corresponding to Representative pathways
module.cluster$.function # Representative pathways
n.cl = 1 # Pick one
module.cluster$.function[n.cl]
# Select nodes that are in the cluster corresponding to this representative pathway
selectNodes(nodes = unlist(strsplit(module.cluster$genes[n.cl], ", ")), by.col = "Symbol",
            preserve.current.selection = FALSE)
# To link the ORA graph with the AMEND module, you'll need to export AMEND module into Powerpoint 
# and insert transparent shapes over module clusters whose color corresponds to color of representative pathway in ORA graph


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#J_AMEND SAM CODE

#============#
# PPI Network ----
#============#
library(igraph)
library(data.table)

# Download and pre-process protein-protein interaction (PPI) network from STRING
# Website: https://string-db.org/ 
#   - Go to Downloads, then type in an organism (e.g., Homo sapiens, Mus musculus)
#   - Choose either 'protein.links' file (includes both physical and functional PPIs) or 'protein.physical.links' (includes only physical PPIs)

# Helper functions
remove_duplicate_edges = function(x){
  g = igraph::graph_from_edgelist(as.matrix(x[,1:2]), directed = FALSE)
  igraph::E(g)$weight = as.numeric(x[,3])
  g = igraph::simplify(graph = g, edge.attr.comb = list(weight = "min"))
  el = igraph::as_edgelist(graph = g)
  el = cbind(el, igraph::E(g)$weight)
  return(el)
}
largest_connected_component = function(g){
  if(!igraph::is_connected(g)){
    message("Taking largest connected component.")
    comps = igraph::components(g)
    largest_comp_id = which.max(comps$csize)
    g = igraph::induced_subgraph(g, which(comps$membership == largest_comp_id))
  }
  g
}

# Load in full PPIN as an edge list
edge.threshold = 0.7
path.2.network = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/R_package/10090.protein.links.v11.0.txt.gz"
ppi = data.table::fread(file = path.2.network, header = TRUE) %>%
  dplyr::mutate(combined_score = combined_score / 1000) %>% # Divide edge scores by 1000 so that they are between 0 and 1
  dplyr::filter(combined_score >= edge.threshold) %>% # Threshold edge scores 
  dplyr::mutate(protein1 = extract_string(protein1, "\\.", 2), # Remove the taxonomy ID '10090'
                protein2 = extract_string(protein2, "\\.", 2)) %>% 
  as.matrix() %>%
  remove_duplicate_edges()

head(ppi)

# Making igraph object, taking largest connected component
g = graph_from_edgelist(ppi[,1:2], directed = FALSE)
E(g)$weight = as.numeric(ppi[,3])
g = largest_connected_component(g)

# Mapping ENSEMBL peptide IDs from graph to gene symbol 
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
# Pick which species you need
species = c('mmusculus', 'hsapiens')[1]
mart_data <- useDataset(paste0(species, "_gene_ensembl"), mart = mart)
return.type = ifelse(species == 'hsapiens', 'hgnc_symbol', 'mgi_symbol') # If you want to map ensembl peptide IDs to Gene Symbols
mapping <- getBM(attributes = c("ensembl_peptide_id", return.type),
                 filters = "ensembl_peptide_id",
                 values = unique(V(g)$name),
                 mart = mart_data)
mapping = mapping[apply(mapping,1,function(x) all(x != "")), ] 

head(mapping)

# Check relationship of mapping (e.g., many-to-one, many-to-many, etc.)
table(table(mapping[,1]))
table(table(mapping[,2]))

# Check quality of mapping
mean(de$tt$Symbol %in% mapping[,2]) # 91 % 
mean(V(g)$name %in% mapping[,1]) # 97 %

id = match(V(g)$name, mapping$ensembl_peptide_id)
vertex_attr(g, "Symbol", which(!is.na(id))) = mapping$mgi_symbol[id[!is.na(id)]]

#====================#
# Map Data to Network ----
#====================#
# Map omics data to nodes in the network
col.ids.to.network = c(2,3,4,6,8,10)
for(i in seq_along(col.ids.to.network)){
  id = match(V(g)$Symbol, de$tt$Symbol)
  igraph::vertex_attr(g, colnames(de$tt)[col.ids.to.network[i]], which(!is.na(id))) = de$tt[id[!is.na(id)],col.ids.to.network[i]]
}

#=================#
# Network Analysis ----
#=================#
devtools::install_github("samboyd0/AMEND")
library(AMEND)

# Perform AMEND to get subnetworks of nodes with large experimental values
# ?run_AMEND
# FUN:
#   shift_scale: For values centered about 0. Takes absolute value then down-weights values not in direction of interest (DOI) by w. 
#       Parameters: DOI, w
#   exp: For log fold changes. exp(DOI * x) for DOI=1,-1, and exp(abs(DOI) * x) for DOI = 0
#       Parameters: DOI
#   p_value: -log10(x + 1e-6)

# Vertex attributes of PPIN
igraph::vertex_attr_names(g)
# Pick 'data.name' to be one of these vertex attributes

#== Log fold change example ==#
data.name = "TMG.saline"
# Searching for both up- and down-regulated nodes
subnet = run_AMEND(graph = g, n = 50, FUN = "exp", FUN.params = list(DOI = 0), data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)
# Searching for only up-regulated nodes, in 2 different ways
subnet = run_AMEND(graph = g, n = 50, FUN = "exp", FUN.params = list(DOI = 1), data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)
subnet = run_AMEND(graph = g, n = 50, FUN = "shift_scale", FUN.params = list(DOI = 1, w = 0.5), data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)

## Inspect output from run_AMEND
# igraph object
subnet$module
# size
vcount(subnet$module)
# mean data value
mean(vertex_attr(subnet$module, data.name))
# Ensembl peptide IDs of nodes in module
V(subnet$module)$name
# Symbols of nodes in module
V(subnet$module)$Symbol

#== P-value example ==#
data.name = "P.Value.TMG.saline"
subnet = run_AMEND(graph = g, n = 50, FUN = "p_value", data = data.name, 
                   normalize = "modified_degree", verbose = TRUE)

## Inspect output from run_AMEND
# igraph object
subnet$module
# size
vcount(subnet$module)
# mean data value
mean(vertex_attr(subnet$module, data.name))
# Ensembl peptide IDs of nodes in module
V(subnet$module)$name
# Symbols of nodes in module
V(subnet$module)$Symbol

#=============================#
# Over-representation Analysis ----
#=============================#
library(fgsea)

path.to.pathways = "/Users/samboyd/Documents/GRA/Network Analysis/AMEND/Subnetwork Data/m_pathways.RDS"
pathway.list = readRDS(path.to.pathways)

ora.paths = fora(pathways = pathway.list, # List of pathways
                 genes = V(subnet$module)$Symbol, # Gene symbols from AMEND module
                 universe = V(g)$Symbol) %>% # Gene symbols from full graph
  dplyr::mutate(overlapGenes = unlist(lapply(overlapGenes, function(x) paste(x, collapse = ", ")))) %>%
  dplyr::filter(padj <= 0.05) # Filter out pathways with adjusted p-value greater than 0.05

head(ora.paths)

#=====================================#
# Visualize AMEND Results in Cytoscape ----
#=====================================#
library(RCy3)
library(RColorBrewer)

#' @title Construct a graph of nested pathways
#' 
#' @description
#' Given a list of pathways resulting from functional enrichment analysis, construct a directed, acyclic graph (DAG) where nodes represent pathways and directed edges go from smaller pathways to larger pathways in which the smaller pathway is partially nested.
#' 
#' @details
#' For two sets A and B, the Nested Index is given by f(A,B) = |intersect(A,B)| / min(|A|,|B|)
#' 
#' In the context of AMEND analysis, which returns a module, we can identify 'representative pathways' that are associated with clusters in the module. First, we run a topological clustering algorithm on the AMEND module to get clusters. Then, for each cluster, we determine a representative pathway by choosing the pathway containing the most nodes in that cluster.
#' 
#' By using 'rep.path.ids', we assign colors to nodes based on whether they are a representative pathway or not.
#' 
#' This may take a very long time to run for inputs of greater than ~30 pathways. In these scenarios, use first.n, or just modify pathway data.frame to contain fewer pathways prior to input. Also, running with parallel=TRUE will help.
#' 
#' @param enriched.pathways data.frame of pathway names along with adjusted and un-adjusted p-values
#' @param rep.path.ids Row IDs of 'enriched.pathways' corresponding to the representative pathways
#' @param all.pathways List of all pathways used in Pathway Analysis
#' @param threshold numeric value between 0 and 1. Threshold for considering one pathway nested within another
#' @param first.n Integer value or NULL. Set to non-NULL if you only want to consider the top n pathways (using the ordering given in 'enriched.pathways')
#' @param parallel logical. Set to TRUE if you want to run in parallel
#' @param nclust Integer value. Number of clusters to use for parallel computing
#' @param padj Logical. Whether to use adjusted p-value for 'pval' vertex attribute of return graph
#' 
#' @return igraph object
#' 
nested_pathways = function(enriched.pathways, rep.path.ids, all.pathways, threshold = 0.8, first.n = NULL, parallel = FALSE, nclust, padj = TRUE){
  require(doParallel)
  require(foreach)
  
  if(nrow(enriched.pathways) == 0){
    stop("nrow of \'enriched.pathways\' is zero")
  }else if(nrow(enriched.pathways) == 1){
    ora_path = enriched.pathways$pathway
    g = igraph::make_empty_graph(n = length(ora_path))
    vlabs = str_replace(ora_path, "Mus musculus: ", "") 
    V(g)$name = vlabs
    temp = numeric(vcount(g))
    for(i in seq_along(temp)){
      temp[i] = length(all.pathways[[ora_path[i]]])
    }
    V(g)$size = temp
    V(g)$pval = ifelse(padj, enriched.pathways$padj, enriched.pathways$pval)
    col.tmp = rep(0, igraph::vcount(g))
    col.tmp[rep.path.ids] = seq_along(rep.path.ids)
    V(g)$color = col.tmp
    return(g)
  }
  
  if(!is.null(first.n)) enriched.pathways = enriched.pathways[1:min(first.n, nrow(enriched.pathways)),]
  
  nested_index = function(a, b){
    intersection = length(intersect(a, b))
    denom = min(length(a), length(b))
    return(intersection / denom)
  }
  
  ora_path = enriched.pathways$pathway
  p.id1 = c()
  p.id2 = c()
  sim.mat = matrix(0, nrow = length(ora_path), ncol = length(ora_path), dimnames = list(ora_path, ora_path))
  for(i in 1:(length(ora_path)-1)){
    for(j in (i+1):length(ora_path)){
      sim.mat[i,j] = nested_index(all.pathways[[ora_path[i]]], all.pathways[[ora_path[j]]])
      sim.mat[j,i] = sim.mat[i,j]
      if(sim.mat[i,j] > threshold){ 
        p.id1 = c(p.id1, ifelse(min(length(all.pathways[[ora_path[i]]]), length(all.pathways[[ora_path[j]]])) == length(all.pathways[[ora_path[i]]]), i, j))
        p.id2 = c(p.id2, ifelse(min(length(all.pathways[[ora_path[i]]]), length(all.pathways[[ora_path[j]]])) == length(all.pathways[[ora_path[i]]]), j, i))
      }
    }
  }
  
  # Finding nested pathways 
  # nest.mat is an edgelist, where a vertex in col 1 is nested within a vertex of col 2 (vertex = pathway)
  if(is.null(p.id1) || is.null(p.id2)){
    g = igraph::make_empty_graph(n = length(ora_path))
    vlabs = str_replace(ora_path, "Mus musculus: ", "") 
    V(g)$name = vlabs
    temp = numeric(vcount(g))
    for(i in seq_along(temp)){
      temp[i] = length(all.pathways[[ora_path[i]]])
    }
    V(g)$size = temp
    V(g)$pval = ifelse(padj, enriched.pathways$padj, enriched.pathways$pval)
    col.tmp = rep(0, igraph::vcount(g))
    col.tmp[rep.path.ids] = seq_along(rep.path.ids)
    V(g)$color = col.tmp
    return(g)
  }else{
    nest.mat = matrix(c(p.id1, p.id2), ncol = 2)
    g = igraph::graph_from_edgelist(nest.mat, directed = TRUE)
    for(i in seq_along(setdiff(1:length(ora_path), V(g)))) g = igraph::add_vertices(g, 1)
  }
  
  # Getting longest simple, directed paths between pair of vertices in edgelist
  # Rationale: If a pathway is nested w/in a larger pathway, it may be nested in other pathways that are themselves nested in that larger pathway.
  #   I want to see the hierarchy of nested-ness
  if(parallel){ # parallel FOR loop
    cl = makeForkCluster(nclust, outfile = "")
    registerDoParallel(cl)
    nest.id0 = foreach(i = 1:nrow(nest.mat), .packages = "igraph") %dopar% {
      message(paste0(i, ":", nrow(nest.mat)))
      asp = all_simple_paths(graph = g, from = nest.mat[i,1], to = nest.mat[i,2], mode = "out")
      asp[[which.max(do.call("c", lapply(asp, length)))]]
    }
    stopCluster(cl)
  }else{ # normal FOR loop
    nest.id0 = vector("list", nrow(nest.mat))
    for(i in 1:nrow(nest.mat)){ # This for loop may take a long time
      message(paste0(i, ":", nrow(nest.mat)))
      asp = all_simple_paths(graph = g, from = nest.mat[i,1], to = nest.mat[i,2], mode = "out")
      nest.id0[[i]] = asp[[which.max(do.call("c", lapply(asp, length)))]]
      # if(i %% 50 == 0) message(paste0(i, ":", nrow(nest.mat)))
    }
  }
  
  # This gets rid of redundant paths
  nest.id = list()
  for(i in 1:length(nest.id0)){
    cond = T
    for(j in (1:length(nest.id0))[(1:length(nest.id0)) != i]){
      if(all(nest.id0[[i]] %in% nest.id0[[j]])) cond = F
    }
    if(cond) nest.id[[length(nest.id) + 1]] = nest.id0[[i]]
  }
  
  # Converting vertex ids to pathway names
  nested.paths = vector(mode = "list", length = length(nest.id))
  for(i in 1:length(nest.id)){
    nested.paths[[i]] = ora_path[nest.id[[i]]]
  }
  # Going from most nested to least nested (smallest to largest)
  
  # Getting graph showing nested hierarchy structure
  el = matrix(nrow = 1000, ncol = 2)
  r.id = 1
  for(i in seq_along(nest.id)){
    x = as.numeric(nest.id[[i]])
    for(j in 1:(length(x) - 1)){
      new.row = x[j:(j+1)]
      if(any(apply(el, 1, function(y) all(new.row %in% y)))) next
      el[r.id,] = new.row 
      r.id = r.id + 1
    }
  }
  el = na.omit(el)
  g = igraph::graph_from_edgelist(el = el, directed = TRUE)
  for(i in seq_along(setdiff(1:length(ora_path), V(g)))) g = igraph::add_vertices(g, 1)
  vlabs = str_replace(ora_path, "Mus musculus: ", "") 
  V(g)$name = vlabs
  temp = numeric(vcount(g))
  for(i in seq_along(temp)){
    temp[i] = length(all.pathways[[ora_path[i]]])
  }
  V(g)$size = temp
  V(g)$pval = enriched.pathways$padj
  col.tmp = rep(0, igraph::vcount(g))
  col.tmp[rep.path.ids] = seq_along(rep.path.ids)
  V(g)$color = col.tmp
  return(g)
}

#=== Upload AMEND module to Cytoscape ===#
createNetworkFromIgraph(subnet$module, title = "AMEND Module", collection = "")
# net.id = getNetworkSuid()

### Set the visual style
# Looking at color schemes
if(0){ 
  display.brewer.all()
  display.brewer.pal(n = 9, name = "RdBu")
  brewer.pal(n = 11, name = "RdBu")
}
color.theme = "RdBu"
# The above may need to change depending on if you want divergent or gradient color themes, which will depend on the data type used in AMEND

style.name = "AMEND_module"
createVisualStyle(style.name)

# Set node shape and size
lockNodeDimensions(new.state = TRUE, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SIZE", value = 40), style.name = style.name)

# Set node color as a function of Experimental data used in AMEND. This will need to be modified for different data types
color.scheme = brewer.pal(n = 3, name = color.theme)
v.attr.name = "TMG.saline" # The name of the vertex attribute used as seed values in AMEND
v.range = c(min(vertex_attr(g, v.attr.name)), max(vertex_attr(g, v.attr.name)))
setNodeColorMapping(table.column = v.attr.name, table.column.values = c(v.range[1], 0, v.range[2]), colors = color.scheme, style.name = style.name)
# Node labels, position, and font size
setNodeLabelMapping(table.column = "Symbol", style.name = style.name)
setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 28), style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = "AlBayan-Bold"), style.name = style.name)
# Set edge style and width
setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
setEdgeLineWidthDefault(new.width = 3, style.name = style.name)

setVisualStyle(style.name)

#=====================================#
# Visualizing ORA Results in Cytoscape ----
#=====================================#
# Get that pathways that each node in module is associated with
node.pathway = vector("list", vcount(subnet$module)); names(node.pathway) = V(subnet$module)$Symbol
for(i in seq_along(node.pathway)){
  node.pathway[[i]] = ora.paths$pathway[unlist(lapply(strsplit(ora.paths$overlapGenes, ", "), function(x) V(subnet$module)$Symbol[i] %in% x))]
}

## Getting cellular functions associated with each cluster in module
set.seed(2088) # Louvain is stochastic, so we set a seed value for reproducibility 
# Running Louvain clustering algorithm on AMEND module to get clusters
cl = igraph::cluster_louvain(graph = subnet$module, resolution = 1)
cl.id = unique(cl$membership)
module.cluster = character(length(cl.id))
# For each cluster, identify pathway that contains the most nodes of that cluster
for(j in seq_along(cl.id)){
  node.names = cl$names[which(cl$membership == cl.id[j])]
  tmp = unlist(node.pathway[which(V(subnet$module)$name %in% node.names)])
  if(length(tmp) == 0){ # If the selected nodes aren't part of any pathway...
    module.cluster[j] = "unknown"
  }else{
    p.count = table(tmp)
    tmp.id = which(p.count == max(p.count))
    # In case of ties for most occurring pathway, take largest pathway
    if(length(tmp.id) > 1) tmp.id = tmp.id[which.max(unlist(lapply(pathway.list[names(p.count)[tmp.id]], length)))]
    module.cluster[j] = names(p.count)[tmp.id]
  }
  names(module.cluster)[j] = paste(V(subnet$module)$Symbol[V(subnet$module)$name %in% node.names], collapse = ", ")
}
if(any(table(module.cluster) > 1)){ # Merge clusters that have same Representative Pathway
  tmp = table(module.cluster)
  tmp.id = which(tmp > 1)
  for(j in seq_along(tmp.id)){
    tmp.id2 = which(module.cluster %in% names(tmp)[tmp.id[j]])
    keep.id = tmp.id2[1]
    rm.id = tmp.id2[-1]
    names(module.cluster)[keep.id] = paste(c(names(module.cluster)[keep.id], names(module.cluster)[rm.id]), collapse = ", ")
    module.cluster = module.cluster[-rm.id]
  }
}
tmp = data.frame(.function = module.cluster, genes = names(module.cluster))
tmp = tmp[tmp$.function != "unknown",]
module.cluster = tmp

## Visualize pathways in Cytoscape
top.k.pathways = 15
threshold = 0.8
tkp = min(top.k.pathways, nrow(ora.paths))
p.tmp = ora.paths[1:tkp,]
rp.tmp = module.cluster$.function[!module.cluster$.function %in% p.tmp$pathway]
if(length(rp.tmp) == 0){
  rp.ids = which(p.tmp$pathway %in% module.cluster$.function)
}else{
  id.tmp = which(!p.tmp$pathway %in% module.cluster$.function)
  p.tmp = rbind2(p.tmp[-tail(id.tmp, length(rp.tmp)),], ora.paths[ora.paths$pathway %in% rp.tmp,])
  rp.ids = which(p.tmp$pathway %in% module.cluster$.function)
}
g.path = nested_pathways(enriched.pathways = p.tmp, rep.path.ids = rp.ids, all.pathways = pathway.list, threshold = threshold, first.n = NULL, parallel = TRUE, nclust = 4)
createNetworkFromIgraph(g.path, title = "", collection = "Nested Pathways Hierarchy")
layoutNetwork("hierarchical")

## Set Visual Style
style.name = "ORA_viz"
createVisualStyle(style.name)

# Set node shape and size
lockNodeDimensions(new.state = TRUE, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_SHAPE", value = "ellipse"), style.name = style.name)
setNodeSizeMapping(table.column = "size", table.column.values = c(10, 500, 1000), 
                   sizes = c(15, 55, 100), style.name = style.name)

# Set node color
gray = brewer.pal(n = 9, name = "Set1")[9]
col.tmp = rep(gray, igraph::vcount(g.path))
col.tmp[V(g.path)$color != 0] = brewer.pal(n = max(sum(V(g.path)$color != 0),3), name = "Set1")[1:sum(V(g.path)$color != 0)]

setNodeColorMapping(table.column = "color", table.column.values = V(g.path)$color, colors = col.tmp, mapping.type = "d", style.name = style.name)

# Node labels, position, and font size
setNodeLabelMapping(table.column = "name", style.name = style.name)
setNodeLabelPositionDefault(new.nodeAnchor = "N", new.graphicAnchor = "S", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = style.name)
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 20), style.name = style.name)
font_style = c("Rockwell-Bold", "ArialNarrow-Bold", "Arial-BoldMT", "Arial-Black")[4]
setVisualPropertyDefault(style.string = list(visualProperty = "NODE_LABEL_FONT_FACE", value = font_style), style.name = style.name)
# Set edge style and width
setEdgeLineStyleDefault(new.line.style = "DOT", style.name = style.name)
setEdgeLineWidthDefault(new.width = 4, style.name = style.name)
setEdgeTargetArrowShapeDefault(new.shape = "DELTA", style.name = style.name)
setVisualStyle(style.name)

### Select Nodes in AMEND module corresponding to Representative pathways
module.cluster$.function # Representative pathways
n.cl = 1 # Pick one
module.cluster$.function[n.cl]
# Select nodes that are in the cluster corresponding to this representative pathway
selectNodes(nodes = unlist(strsplit(module.cluster$genes[n.cl], ", ")), by.col = "Symbol",
            preserve.current.selection = FALSE)
# To link the ORA graph with the AMEND module, you'll need to export AMEND module into Powerpoint 
# and insert transparent shapes over module clusters whose color corresponds to color of representative pathway in ORA graph







mutate(listpathway, qscore = -log(padj, base=10)) %>% 
  barplot(x="qscore")

barplot(listpathway, showCategory=20) 

barplot(
  height,
  x = "Count",
  color = "p.adjust",
  showCategory = 8,
  font.size = 12,
  title = "",
  label_format = 30,
  ...
)

require(clusterProfiler)

require(DOSE)

















if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")

library(DOSE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("enrichplot")

library(enrichplot)

head(fgseaRes)
topPathwaysUp1 <- fgseaRes[ES > 0][head(order(padj), n=10), pathway:size]
topPathwaysDown1 <- fgseaRes[ES < 0][head(order(padj), n=10), pathway:size]
topPathways1 <- rbind(topPathwaysUp1, rev(topPathwaysDown1))

topPathways1<-as.data.frame(topPathways1)
#dotplot(topPathwaysUp1, showCategory=30) + ggtitle("dotplot for GSEA")

library(ggplot2)

ggplot(topPathwaysUp1, aes(x = padj, y = pathway, color = padj)) + 
  geom_bar(stat = 'identity') + 
  xlab("padj") + ylab("pathway") + ggtitle("topPathwaysUp1") + 
  theme_bw()



BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
# For visualisation
install.packages('pheatmap')
#install.packages("DOSE")
#install.packages("enrichplot")
install.packages("ggupset")

enrichres <- new("topPathways1",
                 readable = FALSE,
                 result = res_df,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = df$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)

# Barplot
barplot(topPathways1, showCategory = 20) 
mutate(topPathways1, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")





mutate(topPathways, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

#edo2 <- gseDO(prot.list, exponent = 1,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              verbose = TRUE,
              seed = FALSE,
              by = "fgsea",
)
#dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")

edo2 <- gseDO(fgseaRes)
head(fgseaRes)
dotplot(fgseaRes, showCategory=30) + ggtitle("dotplot for GSEA")

ggplot(topPathways1, aes(x=NES, y=Pathway, size=FDR,color=-1*log10(FDR))) + geom_point()+scale_size_area(max_size = 5)+scale_colour_gradient(low="green",high="red")
library(ggplot2)


topPathways
topPathways1<-as.data.frame(topPathways)

fgseaRes[order(pval), ]
#Plot5
plotEnrichment(all.gene.sets[["THUM_SYSTOLIC_HEART_FAILURE_UP"]],prot.list) + labs(title="THUM_SYSTOLIC_HEART_FAILURE_UP")








































#Plot4
plotEnrichment(all.gene.sets[["WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS"]],prot.list) + labs(title="WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS")


#Plot3
plotEnrichment(all.gene.sets[["GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY"]],prot.list) + labs(title="GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY")
#Plot4
plotEnrichment(all.gene.sets[["GOMF_CIS_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING"]],prot.list) + labs(title="GOMF_CIS_REGULATORY_REGION_SEQUENCE_SPECIFIC_DNA_BINDING")
#Plot5
plotEnrichment(all.gene.sets[["WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA"]],prot.list) + labs(title="WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA")
#Plot6
plotEnrichment(all.gene.sets[["REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION"]],prot.list) + labs(title="REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION")

#Plot7
plotEnrichment(all.gene.sets[["KEGG_PARKINSONS_DISEASE"]],prot.list) + labs(title="KEGG_PARKINSONS_DISEASE")
fgseaRes[order(pval), ][15:20,]
#Plot34
plotEnrichment(all.gene.sets[["KEGG_ALZHEIMERS_DISEASE"]],prot.list) + labs(title="KEGG_ALZHEIMERS_DISEASE")

#Plot8
plotEnrichment(all.gene.sets[["GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT"]],prot.list) + labs(title="GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT")
#Plot9
plotEnrichment(all.gene.sets[["GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN"]],prot.list) + labs(title="GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN")
#Plot10
plotEnrichment(all.gene.sets[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],prot.list) + labs(title="KEGG_OXIDATIVE_PHOSPHORYLATION")
#Plot11
plotEnrichment(all.gene.sets[["GOBP_MEMORY"]],prot.list) + labs(title="GOBP_MEMORY")
#Plot12
plotEnrichment(all.gene.sets[["GOBP_ACTIN_FILAMENT_BASED_PROCESS"]],prot.list) + labs(title="GOBP_ACTIN_FILAMENT_BASED_PROCESS")
#Plot13
plotEnrichment(all.gene.sets[["MCBRYAN_PUBERTAL_BREAST_4_5WK_UP"]],prot.list) + labs(title="MCBRYAN_PUBERTAL_BREAST_4_5WK_UP")
#Plot14
plotEnrichment(all.gene.sets[["GOBP_COGNITION"]],prot.list) + labs(title="GOBP_COGNITION")
#Plot15
plotEnrichment(all.gene.sets[["GOCC_9PLUS0_NON_MOTILE_CILIUM"]],prot.list) + labs(title="GOCC_9PLUS0_NON_MOTILE_CILIUM")
#Plot16
plotEnrichment(all.gene.sets[["GOBP_MICROTUBULE_ORGANIZING_CENTER_ORGANIZATION"]],prot.list) + labs(title="GOBP_MICROTUBULE_ORGANIZING_CENTER_ORGANIZATION")
#Plot17
plotEnrichment(all.gene.sets[["GOBP_SYNAPTIC_SIGNALING"]],prot.list) + labs(title="GOBP_SYNAPTIC_SIGNALING")
#Plot18
plotEnrichment(all.gene.sets[["GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_UP"]],prot.list) + labs(title="GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_UP")
#Plot19
plotEnrichment(all.gene.sets[["GOBP_CYTOSKELETON_ORGANIZATION"]],prot.list) + labs(title="GOBP_CYTOSKELETON_ORGANIZATION")
#Plot20
plotEnrichment(all.gene.sets[["GOBP_CELL_PROJECTION_ORGANIZATION"]],prot.list) + labs(title="GOBP_CELL_PROJECTION_ORGANIZATION")


### GSEA analysis logFCgreat
# Load All gene sets file downloaded from Broad Institute 
# The following website contains the gene set collection or the complete Molecular Signatures Database (MSigDB)  
# http://software.broadinstitute.org/gsea/downloads.jsp#msigdb
#  
library(fgsea)
all.gene.sets <- gmtPathways("C:\\Users\\sophi\\OneDrive\\Desktop\\J_POSTDOC AND JOBS\\J_APPLIED DATA SCIENCE\\J_STATISTICS\\J_STATISTICAL GENOMICS BIOS 799\\msigdb.v7.4.symbols.gmt")
class(all.gene.sets)
length(all.gene.sets)
all.gene.sets[1:2]
# Show first a few pathways, and within those, show only the first few genes. 
library(tidyverse)
all.gene.sets %>% head() %>% lapply(head)

head(prot.list2)
head(d.entrez.id2)  
### Now run fgsea 
library(fgsea)
fgseaRes2 <- fgsea(pathways = all.gene.sets, stats = prot.list2, minSize=15, maxSize=500,eps=0)
head(fgseaRes2)
head(fgseaRes2[order(pval), ])
sum(fgseaRes2[, padj < 0.05])
fgseaRes3<-fgseaRes2[order(pval), ]
head(fgseaRes3)
library(ggplot2)
# Make a few Enrichment Plots
#Plot1
#plotEnrichment(all.gene.sets[["GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS"]],gene.list) + labs(title="GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS")

#Plot1
plotEnrichment(all.gene.sets[["GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY"]],prot.list2) + labs(title="GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY")
#Plot2
plotEnrichment(all.gene.sets[["GOMF_SEQUENCE_SPECIFIC_DNA_BINDING"]],prot.list2) + labs(title="GOMF_SEQUENCE_SPECIFIC_DNA_BINDING")
#tail(fgseaRes[order(pval), ])

# Make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes2[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes2[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(all.gene.sets[topPathways], prot.list, fgseaRes,gseaParam = 0.5)

head(topPathways)


########################################################################################################
#Heatmap
dd
head(dd)
heatmap(dd, main = "SALINE Vs TMG", xlab = "Samples", ylab = "DEG")
#https://www.biostars.org/p/374551/

head(d.go.DE1)
head(prot.list)
dd1<-d.go.DE1
head(dd1)
dd1<-as.data.frame(dd1)
head(dd1)
dd1<-dd1[,c(-2)]
dd1<-as.data.frame(dd1)
head(dd1)

head(mapping1) 
dim(mapping1)#[1] 8334   19
head(dat1)
dim(dat1) #8334   18
type(dat1)
dd1<-dat1
head(dd1)
type(dd1)


#dd1<-cbind(dd1,mapping1$Symbol)
colnames(dd1)
#colnames(dd1)<-c('GFPCTRL_1','GFPCTRL_2','GFPCTRL_3','GFPTEN_1','GFPTEN_2','GFPTEN_3','OGTCTRL_1','OGTCTRL_2','OGTCTRL_3','OGTTEN_1','OGTTEN_2','OGTTEN_3','ROGTCTRL_1','ROGTCTRL_2','ROGTCTRL_3','ROGTTEN_1','ROGTTEN_2','ROGTTEN_3','Symbol')
#colnames(dd1) <- c('GFP0_1','GFP0_2','GFP0_3','GFP10_1','GFP10_2','GFP10_3','T605_0_1','T605_0_2','T605_0_3','T605_10_1','T605_10_2','T605_10_3','T606_0_1','T606_0_2','T606_0_3','T606_10_1','T606_10_2','T606_10_3','Symbol')

head(dd)
heatmap(dd, main = "SALINE Vs TMG", xlab = "Samples", ylab = "DEG")

head(mapping1) 
dim(mapping1)#[1] 8334   19
head(dat1)
dim(dat1) #8334   18
type(dat1)
dd1<-dat1
head(dd1)
type(dd1)
mat <- as.matrix(dd1)
type(dd1)
type(mat)
mat
type(mat)
rownames(mat) <- mapping1$Symbol
mat[1:3,1:3]
class(mat)
typeof(mat)

heatmap(mat, main = "Heatmap", xlab = "Samples", ylab = "Proteins")

#https://www.biostars.org/p/374551/

#make heatmap of log FC
dea$tt
dd1<-dea$tt
dd2<-dea$tt[,-1]
head(dd2)
dd3<-dd2[,1:7]
head(dd3)
type(dd3)
mata <- as.matrix(dd3)
type(dd3)
type(mata)
head(dd1)
rownames(mata) <- dd1$Symbol
mata[1:3,1:3]
class(mata)
typeof(mata)
heatmap(mata, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
heatmap(mata, main = "Heatmap", ylab = "Proteins")

#make heatmap of log FC
dea$tt
dd1<-dea$tt
dd2<-dea$tt[,-1]
head(dd2)
dd3<-dd2[,1:7]
head(dd3)
type(dd3)
dd4<-dd3
type(dd4)
dd4<-dd4[,1:3]
head(dd4)
type(dd4)

matb <- as.matrix(dd4)
type(dd4)
type(matb)
head(dd4)
rownames(matb) <- dd1$Symbol
matb[1:3,1:3]
class(matb)
typeof(matb)
heatmap(matb, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
heatmap(matb, main = "Heatmap", ylab = "Proteins", cexCol = 0.8)

head(dd3)
dd5<-dd3
head(dd5)
type(dd5)
dd5<-dd5[,-(1:3)]
head(dd5)
type(dd5)
matc <- as.matrix(dd5)

type(matc)
head(matc)
rownames(matc) <- dd1$Symbol
head(matc)
class(matc)
typeof(matc)
#heatmap(matc, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
heatmap(matc, main = "Heatmap", ylab = "Proteins", cexCol = 0.8)

#Remove all ROGT
head(dd3)
dd6<-dd3
head(dd6)
type(dd6)
dd6<-dd6[,-(3)]
head(dd6)
dd6<-dd6[,-(5:6)]
type(dd6)
head(dd6)
matd <- as.matrix(dd6)

type(matd)
head(matd)
rownames(matd) <- dd1$Symbol
head(matd)
class(matd)
typeof(matd)
#heatmap(matc, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
head(matd)
heatmap(matd, main = "Heatmap", ylab = "Proteins", cexCol = 0.8)

#GFPTEN.GFPCTRL OGTTEN.OGTCTRL 
head(dd3)
dd7<-dd3
head(dd7)
type(dd7)
dd7<-dd7[,-(3)]
head(dd7)
dd7<-dd7[,-(3:6)]
type(dd7)
head(dd7)
mate <- as.matrix(dd7)

type(mate)
head(mate)
rownames(mate) <- dd1$Symbol
head(mate)
class(mate)
typeof(mate)
#heatmap(matc, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
head(mate)
heatmap(mate, main = "Heatmap", ylab = "Proteins", cexCol = 0.8)

#OGTCTRL.GFPCTRL OGTTEN.GFPTEN; ROGTCTRL.GFPCTRL ROGTTEN.GFPTEN
head(dd3)
dd8<-dd3
head(dd8)
type(dd8)
dd8<-dd8[,-(1:3)]
head(dd8)
#dd8<-dd8[,-(3:6)]
type(dd8)
head(dd8)
matf <- as.matrix(dd8)

type(matf)
head(matf)
rownames(matf) <- dd1$Symbol
head(matf)
class(matf)
typeof(matf)
#heatmap(matc, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
head(matf)
heatmap(matf, main = "Heatmap", ylab = "Proteins", cexCol = 0.8)

#OGTCTRL.GFPCTRL OGTTEN.GFPTEN; 
head(dd3)
dd9<-dd3
head(dd9)
type(dd9)
dd9<-dd9[,-(1:3)]
head(dd9)
dd9<-dd9[,-(3:4)]
type(dd9)
head(dd9)
matg <- as.matrix(dd9)

type(matg)
head(matg)
rownames(matg) <- dd1$Symbol
head(matg)
class(matg)
typeof(matg)
#heatmap(matc, main = "Heatmap", xlab = "Samples", ylab = "Proteins")
head(matg)
heatmap(matg, main = "Heatmap", ylab = "Proteins", cexCol = 0.8)

