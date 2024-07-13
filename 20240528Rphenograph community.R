
# Rphenograph clustering of normalised, concatenated and scaled data 

## Scripts required to run previous to this:
# -Normalisation, concatenation and scaling of single cell data

## Megan Cole 

#Load packages
library(ggplot2)
library(plyr)
#library(cytofkit)
library(Rphenograph)
print("Phenograph package loaded")

##############################
#### SET GLOBAL VARIABLES ####
##############################

# Name of data file to read in 
path = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile"

filename = "frequencies_25px_20240626.csv"

experiment_path = "/Results/Figures_output/"

output_path = "/Results/Figures_output/"
Sys.setenv(OMP_NUM_THREADS="8") # inside RStudio
Sys.getenv("OMP_NUM_THREADS")
# Load data 
# celldata = read.csv(paste(path, experiment_path, filename, sep = ""))
celldata = read.csv(paste(path, experiment_path, filename, sep = ""), check.names=FALSE)
# celldata = tail(celldata, n = 2000) 
# celldata <- celldata_o[1:2000,]
####################################################################

# Obtain information about the loaded dataset
dim(celldata)
names(celldata)
unique(celldata$filename)



# Indicate expected length of selectcolumns and select columns for clustering
#Alt MAC	B cell	Cancer	Cl MAC	Cl Mo	DCs cell	Endothelial cell	Int Mo	Mast cell	NK cell	Neutrophils	Non-Cl Mo	T other	Tc	Th	Treg
#check the columns carefully 
relevant_columns = c()
# 13 because note below`
columns_of_interest = 13
selectcolumns = c(grep("B cells", names(celldata), ignore.case = TRUE),
                  grep("Dendritic cells", names(celldata), ignore.case = TRUE), 
                  grep("Endothelium", names(celldata), ignore.case = TRUE),
                  grep("Epithelium", names(celldata), ignore.case = TRUE),
                  grep("Fibroblasts", names(celldata), ignore.case = TRUE),
                  grep("Macrophages type 1", names(celldata), ignore.case = TRUE),
                  grep("Macrophages type 2", names(celldata), ignore.case = TRUE),
                  grep("Neutrophils", names(celldata), ignore.case = TRUE),
                  grep("T cells CD4", names(celldata), ignore.case = TRUE),
                  grep("T cells CD8", names(celldata), ignore.case = TRUE),
                  grep("T reg cells", names(celldata), ignore.case = TRUE), 
                  grep("Tumour", names(celldata), ignore.case = TRUE))
# NOTA BENE! WE take both DC classes together with this but this could cause problems down the line. Make sure you select all the right columns here.
                  # grep("Dendritic cells CD103", names(celldata), ignore.case = TRUE))

selectcolumns
length(selectcolumns)
# Check if correct columns have been selected
if(length(selectcolumns) < columns_of_interest){
  print("WARNING: there is a problem with the data selection - missing column name(s)")
} else if(length(selectcolumns) > columns_of_interest){
  # Remove repeat selection of columns 
  selectcolumns = unique(selectcolumns)
  length(unique(selectcolumns))
  
  if(length(selectcolumns) > columns_of_interest){
    # If the length is still above 17, print a warning message that data selection is wrong
    stop(paste("There is a problem wtith data selection - after an attempt to fix, there are still too many columnns selected, expected:", 
               columns_of_interest, "found:", length(selectcolumns)))
  }
  else if(length(selectcolumns) < columns_of_interest){
    stop(paste("There is a problmen with the data selection - missing column name(s), expected:", columns_of_interest,
               "found:", length(selectcolumns)))
  } else{
    print(paste("The list of selected columns has been corrected - there are now", columns_of_interest, "selected columns"))
  }
} else {
  print(paste("The correct number of columns were found - ", columns_of_interest))
}

# Retrieve names of select columns 
names(celldata[,selectcolumns])

# Run Rphenograph on data of select columns
Rphenograph_out <- Rphenograph(celldata[,selectcolumns], k=250)

pheno_temp <- length(as.numeric(membership(Rphenograph_out[[2]])))
pheno_temp
dim(celldata)

celldata$cluster <- as.numeric(membership(Rphenograph_out[[2]]))

# Save data with clustering information added
write.csv(celldata, paste(path, output_path, "Rphenograph_Megan_26June_output_", max(celldata$cluster), "clusters_k250_", columns_of_interest, "ct_fractions.csv", sep = ""))
print("file saved")


############################################################################################################
