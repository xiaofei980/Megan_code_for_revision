print(Sys.time())
print(paste("C stack size:", Cstack_info()["size"]))
print("Loading packages...")
library(stringr)
library(dplyr)
library(Rtsne)
print("Packages loaded")





# Path for loading data and saving results 
# Project folder path:
base = "/mnt/Data1/imcyto/runs/Xiaofei/Megan_IMC"

# Name of input and output directories:
input_path = file.path(base, "Projectfile/Data/tSNE_input")
output_path = file.path(base, "Projectfile/Results/tSNE_output/")

# Name of the file wi th normalised and scaled data:
input_file = "240523_normalised_concatenated_scaled_data_001threshold.csv"


print(Sys.time())
print("Loading celldata...")
# Load celldata excluding tSNE coordinates
celldata = read.csv(file.path(input_path, "240523_normalised_concatenated_scaled_data_001threshold.csv"),
                    stringsAsFactors = F)


print("Celldata loaded")

# Select columns to be used for t-SNE generation:
# All markers except BG, DNA and domain
tSNEcolumns = grep("MI_", names(celldata), value = TRUE)
exclCOLs = c("MI_Argon80", "MI_Xenon131", "MI_Xenon132",  "MI_DNA1",
             "MI_DNA2")
tSNEcolumns = tSNEcolumns[!tSNEcolumns %in% exclCOLs]

print(Sys.time())
print("Running t-SNE")
print(Sys.time())
tSNE_results = Rtsne(celldata[,tSNEcolumns], num_threads = 4,
                     perplexity = 100, verbose = TRUE,
                     check_duplicates = FALSE, max_iter = 1500)
print(Sys.time())
print("t-SNE completed")

tSNE1 = tSNE_results$Y[,1]
tSNE2 = tSNE_results$Y[,2]

tSNE = cbind(tSNE1, tSNE2)

#date = str_sub(gsub("-", "", Sys.Date()), -6, -1)
#output_file = paste0(output_path, date, "_tSNE.csv")
#write.csv(tSNE, file = output_file, row.names = F)
#print(Sys.time())
#print("finished!")
date = str_sub(gsub("-", "", Sys.Date()), -6, -1)
output_file = paste0(output_path, date, "_tSNE.csv")
write.csv(tSNE, file = output_file, row.names = F)
print(Sys.time())
print("finished!")
# Save data with clustering information added 

