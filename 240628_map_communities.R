

#### Scripts for Treg paper - Figure 3 & Supplementary figure 3 ####

## Megan Cole (megan.cole@crick.ac.uk)

rm(list=ls())

# Load relevant libraries 
library(ggplot2) 
library(dplyr) 
library(tiff) 
library(reshape) 
library(tibble)
library(tidyverse)
library(conflicted)


date = format(Sys.Date(), "%Y%m%d")

## Global variables ##

input_path = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Data/Figures_input/"
output_path = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Results/Figures_output/"
out_fig_test = paste0(output_path,'Organized_Figures/map_community/')
dir.create(out_fig_test)
# Rphenograph20_240626 <- Rphenograph20_240626
# read.csv(paste0(input_path,"Rphenograph20_240626.csv"))

#Rphenograph20_240626= subset(Rphenograph20_240626,select = -c(annotation.y,ROI_ID.y))
# Rphenograph20_240626%>%
#   rename(ROI_ID=ROI_ID.x,
#          annotation=annotation.x)
names(Rphenograph20_240626)[names(Rphenograph20_240626)=='ROI_ID.x']<-'ROI_ID'
#names(Rphenograph20_240626)[names(Rphenograph20_240626)=='annotation.x']<-'annotation'
colours =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")
treat_col = c("MRTX+PD1" = "#F8766D", "MRTX+PD1+CTLA4" = "#00BFC4")


###################
# Spatial location of communities on each ROI
images1 = c("ROI1_BRAC8816_5C_MRTX+PD1_ROI_001_7",         
            "ROI1_BRAC9050_4B_MRTX+PD1_ROI_001_1" ,         
            "ROI1_BRAC9091_3B_MRTX+PD1+CTLA-4_ROI_001_3",   
            "ROI1_BRAC9425_1A_MRTX+PD1+CTLA-4_ROI_001_1_x1",
            "ROI1_BRAC9437_1A_MRTX+PD1+CTLA-4_ROI_001_1",  
            "ROI2_BRAC8965_6F_MRTX+PD1+CTLA-4_ROI_002_9_x1",
            "ROI2_BRAC9252_2A_MRTX+PD1_ROI_002_2" )

out_dir = paste0(out_fig_test, "20community/")
dir.create(out_dir)


# Table of comparable communities and labels - created in powerpoint 

Rphenograph20_240626$top5 = "Other"
Rphenograph20_240626[which(Rphenograph20_240626$community == 7), "top5"] <- "T/DC1"
Rphenograph20_240626[which(Rphenograph20_240626$community == 4), "top5"] <- "T/M1"
Rphenograph20_240626[which(Rphenograph20_240626$community == 5), "top5"] <- "T/reg"
Rphenograph20_240626[which(Rphenograph20_240626$community == 2), "top5"] <- "T/M2"
Rphenograph20_240626[which(Rphenograph20_240626$community == 9), "top5"] <- "T/NA"

unique(Rphenograph20_240626$top5)

#Rphenograph20_240626$top6 = factor(Rphenograph20_240626$top6)

# Define variables
date <- "2024_06_26"
top5 <- c("T/DC1", "T/M1", "T/reg", "T/M2", "T/NA")
top5_col <- c("T/DC1" = "green", "T/M1" = "red", "T/reg" = "#FFEA42", "T/M2" = "#4E79A7FF", "T/NA" = "#CC6666", "Other" = "#E1E1E1")

# Write CSV
write.csv(Rphenograph20_240626, paste0(output_path, date, "_neighbor_top5.csv"), row.names = FALSE)

# Ensure the levels of the factor are set correctly
Rphenograph20_240626$top5 <- factor(Rphenograph20_240626$top5, levels = top5)
unique(Rphenograph20_240626$top5)
 ## here others -->"NA"


# Prepare labels for plotting
nbclusters <- top5
labels <- as.data.frame(top5_col[1:5], stringsAsFactors = FALSE)
labels <- rownames_to_column(labels, "cluster")
colnames(labels)[2] <- 'colour'
add <- data.frame(cluster = c("101", "102"), colour = c("black", "white"), stringsAsFactors = FALSE)
labels <- rbind(add, labels)

# Define the cell_map function
cell_map <- function(ROI, files, data, clusters, path, labels) {
  
  # Filter the data for the specific ROI
  d <- data[data[[files]] == ROI, ]
  
  # Construct file paths for TIFF images
  base_path <- "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Data/Figures_input/"
  TIFFnm <- file.path(base_path, "all_cells_mask", paste0(ROI, "_all_cells_mask.tiff"))
  TIFFol <- file.path(base_path, "Cells_outline", paste0(ROI, "_Cells_outline.tiff"))
  
  # Print file paths for debugging
  print(paste("TIFFnm:", TIFFnm))
  print(paste("TIFFol:", TIFFol))
  
  # Read TIFF images
  if (!file.exists(TIFFnm)) {
    stop(paste("File not found:", TIFFnm))
  }
  if (!file.exists(TIFFol)) {
    stop(paste("File not found:", TIFFol))
  }
  
  TIFF <- readTIFF(TIFFnm)
  TIFF2 <- readTIFF(TIFFol)
  
  # Convert TIFF images to long format data frames
  to_long_df <- function(TIFF) {
    names(TIFF) <- seq_len(length(TIFF))
    TIFF <- melt(TIFF, id.vars = c(row(TIFF), names(TIFF)))
    names(TIFF)[1] <- "y"
    names(TIFF)[2] <- "x"
    TIFF
  }
  
  Outline <- TIFF2
  ROI1 <- TIFF
  ROI1 <- to_long_df(ROI1)
  Outline <- to_long_df(Outline)
  ROI1$unique_px_ID <- seq_len(nrow(ROI1))
  ROI2 <- ROI1
  n <- 2
  
  # Choose the plotting of neighbor clusters
  for (cl in nbclusters){
    print(cl)
    cluster_xy = d[which(d$top5 == cl), c("ObjectNumber","Location_Center_X","Location_Center_Y", "top5")]
    names(cluster_xy)[names(cluster_xy) == 'top5'] <- 'cluster' 
    #If cluster cl is not found in the image/ROI, remove it from the labels list
    if (dim(cluster_xy)[1] == 0){
      print(paste("cluster ", cl, " is not found in ", ROI, sep = ""))
      next
    } else {
      cluster_xy$Location_Center_X = round(cluster_xy$Location_Center_X)
      cluster_xy$Location_Center_Y = round(cluster_xy$Location_Center_Y)
      names(cluster_xy) = c("ObjectNumber","x","y","cluster")
      colours_in_mask <- inner_join(ROI1, cluster_xy[,c("x","y","cluster")])
      min = min(unique(colours_in_mask$value))
      print(min(unique(colours_in_mask$value)))
      colours_in_mask = colours_in_mask[-(which(colours_in_mask$value == min)),]
      ROI2[which(ROI1$value %in% colours_in_mask$value), "value"] = n
      ROI2$name[ROI2$value == n] = first(cluster_xy$cluster)
      n = n+1
    }
  }
  print("now we come to the bit that is suspicious")
  # ROI2$name = ifelse(ROI2$name == "1", "T/DC1",
  #                    ifelse(ROI2$name == "2", "T/M1",
  #                           ifelse(ROI2$name == "3", "T/reg",
  #                                  ifelse(ROI2$name == "4", "T/M2",
  #                                         ifelse(ROI2$name == "5", "T/NA","NA" )))))

  
  
  # #"T/DC1", "T/M1", "T/reg", "T/M2", "T/NA")
  # Recode ROI2$name using dplyr::recode
  ROI2$name <- dplyr::recode(ROI2$name,
                             "1" = "T/DC1",
                             "2" = "T/M1",
                             "3" = "T/reg",
                             "4" = "T/M2",
                             "5" = "T/NA",
                             .default = "NA")
  
  background = ROI2[which(ROI2$value < 1),"unique_px_ID"]
  ROI2[which(Outline$value == 1 ),"value"] = 1
  ROI2[which(ROI2$unique_px_ID %in% background),"value"] = 0
  ROI2$name[ROI2$value == 0] = 101
  ROI2$name[ROI2$value == 1] = 102
  
  
  ## Subset labels based on communities present in that image
  labels = labels[which(labels$cluster %in% unique(ROI2$name)),]
  
  # Change 101 & 102 values to blank
  labels$cluster[labels$cluster == 101 | labels$cluster == 102] <- " "
  
  p = ggplot(ROI2,aes(x=x,y=-y, fill=as.factor(value))) +
    geom_raster() +
    theme_void() +
    theme(legend.title=element_blank(),
          legend.text = element_text(colour = "black")) +
    scale_fill_manual(values = alpha(labels$colour,  1), labels = labels$cluster)
  
  filename = paste(ROI, ".pdf", sep = "")
  ggsave(plot = p, device = "pdf", width=5.6, height=5, dpi=300, path = path, filename = filename)
  print("plot saved")
  return(p)
}

# Apply cell_map function to each ROI in images1
for (ROI in images1) {
  print(ROI)
  p <- cell_map(ROI, "filename", Rphenograph20_240626, top5, out_dir, labels)
}



