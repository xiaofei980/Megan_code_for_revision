library(ggplot2) 
library(Rtsne) 
library(clustree) 
library(ggdendro) 
library(tiff) 
library(reshape2) 
library(gplots)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(gridExtra)
library(ggridges)
library(ggpubr)
library(scales)
library(clustree)
library(cowplot)
library(umap)
library(lme4)
library(Rphenograph)
library(ggpubr)
library(rlang)

input_path = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Data/Figures_input/"
output_path="/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Results/Figures_output/"

# Name of input and output directories:
Rphenograph20_240627 = read.csv(paste0(input_path, "/Rphenograph_Megan_26June_output_20clusters_k250_13ct_fractions.csv"),check.names=F)
# celltypepercTable include the frequency information of each source cell to other cell type. 
celltypepercTable<-read.csv(paste0(input_path, "/frequencies_25px_20240626.csv"),check.names=F)
##add community information into celltypepercTable
names(celltypepercTable)[names(celltypepercTable)=="source_ID"]<-"cell_ID"
names(celltypepercTable)
Rphenograph20_community <- Rphenograph20_240627[, c("cell_ID", "community","ExpGroup")]
celltypepercTable_community <- merge(x = celltypepercTable, y = Rphenograph20_community, by = "cell_ID")
print(head(celltypepercTable_community))

#colnames(Rphenograph21["cluster"])<-"community"
#colnames(Rphenograph21[18])<-community
Rphenograph20_240627$community = Rphenograph20_240627$cluster

Rphenograph20_240627$cell_ID = Rphenograph20_240627$source_ID
data1[!(data1$cell_ID %in% Rphenograph20_240627$cell_ID),]

Rphenograph20_240627_2 = inner_join(x=data1, y=Rphenograph20_240627, by= "cell_ID")
Rphenograph20_240627 = Rphenograph20_240627_2
rm(Rphenograph20_240627_2)


colours1 =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")


# Calculate the percentage of cell types 
## The next code can be used all the time ,add marker or delete maker threshold 
## code is long ,so check and choose best name for table for reusing ..
##record the cutoff in a slide and check with the raw image.  try to reuse /rewrite a code that map to FIJI 
Rphenograph20_240627$PDL1exp <- ifelse(Rphenograph20_240627$MI_PDL1 >= 0.5 , "PDL1+", "PDL1-")
Rphenograph20_240627$Casp3exp <- ifelse(Rphenograph20_240627$MI_Casp3 >= 0.8 , "Casp3+", "Casp3-") 
Rphenograph20_240627$Ki67exp <- ifelse(Rphenograph20_240627$MI_Ki67 >= 0.5 , "Ki67+", "Ki67-") 
Rphenograph20_240627$PD1exp <- ifelse(Rphenograph20_240627$MI_PD1 >= 0.15 , "PD1+", "PD1-") 
Rphenograph20_240627$LAG3exp <- ifelse(Rphenograph20_240627$MI_LAG3 >= 0.01 , "LAG3+", "LAG3-") 
Rphenograph20_240627$CD44exp <- ifelse(Rphenograph20_240627$MI_CD44 >= 0.4 , "CD44+", "CD44-") 
Rphenograph20_240627$MHCIIexp <- ifelse(Rphenograph20_240627$MI_MHCII >= 0.25 , "MHCII+", "MHCII-") 
Rphenograph20_240627$CD86exp <- ifelse(Rphenograph20_240627$MI_CD86 >= 0.25 , "CD86+", "CD86-") 

str(Rphenograph20_240627)
Rphenograph20_240627$ExpGroup <- as.factor(Rphenograph20_240627$ExpGroup)
Rphenograph20_240627$ROI_ID.x <- as.factor(Rphenograph20_240627$ROI_ID.x)
str(Rphenograph20_240627)
# Remove the unnamed column
Rphenograph20_240627 <- Rphenograph20_240627[, !colnames(Rphenograph20_240627) == ""]
# Verify the column removal
colnames(Rphenograph20_240627)

# limits = c("MRTX+PD1", "MRTX+PD1+CTLA-4")
# stat_summary_means_20_community=Rphenograph20_240627[Rphenograph20_240627$ExpGroup %in% limits,]%>%
#   group_by(community,ExpGroup,ROI_ID.x) %>%

# here just copy from celltype_group script .
stat_summary_means_20_community = Rphenograph20_240627 %>%
  group_by(ExpGroup,ROI_ID.x,community) %>%
  # stat_summary_means_20_group = Rphenograph20_240627 %>%
  # group_by(ROI_ID.x,ExpGroup) %>%
  dplyr::summarise(
    # marker mean intensity: PD1 ,LAG3,TIM3, PDL1, CD86,MHCII, CXCL9, PVR...... 
    MI_PD1 = mean(MI_PD1),
    MI_MHCII= mean(MI_MHCII),
    MI_PDL1 = mean(MI_PDL1),
    MI_LAG3 = mean(MI_LAG3),
    MI_CD86 = mean(MI_CD86),
    
    #each cell type count 
    totaln = n(),
    Tumourn = sum(annotation=="Tumour"),
    Tcelln = sum(annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    CD4n = sum(annotation=="T cells CD4"),
    CD8n = sum(annotation=="T cells CD8"),
    Tregn = sum(annotation=="T reg cells"),
    Bcelln = sum(annotation=="B cells"),
    Neutn = sum(annotation=="Neutrophils"),
    Maccelln = sum(annotation %in% c("Macrophages type 1","Macrophages type 2")),
    Mac1n = sum(annotation=="Macrophages type 1"),
    Mac2n = sum(annotation=="Macrophages type 2"),
    DCcelln = sum(annotation %in% c("Dendritic cells","Dendritic cells CD103")),#DCcelln means all DCs 
    DCn = sum(annotation=="Dendritic cells"), #DCn means DC others, here is a little bit confused, maybe change later, now to fit previous figure, 
    cDC1n = sum(annotation=="Dendritic cells CD103"),
    Fibn = sum(annotation=="Fibroblasts"),
    immunecellsn = sum(annotation %in% c("B cells","T cells CD4","T cells CD8","Dendritic cells","T reg cells","Neutrophils","Macrophages type 1","Macrophages type 2","Dendritic cells CD103")),
    #immunce cells: B T Mac DC Neutrophils
    
    #PD1, TIM3, LAG3 expression on different T cells,here 
    PD1Tcelln = sum(PD1exp=="PD1+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    PD1CD4Tcelln = sum(PD1exp=="PD1+" & annotation=="T cells CD4"),
    PD1CD8Tcelln = sum(PD1exp=="PD1+" & annotation=="T cells CD8"),
    PD1Tregn = sum(PD1exp=="PD1+" & annotation=="T reg cells"),
    
    LAG3Tcelln = sum(LAG3exp=="LAG3+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    LAG3CD8Tcelln = sum(LAG3exp=="LAG3+" & annotation=="T cells CD8"),
    LAG3CD4Tcelln = sum(LAG3exp=="LAG3+" & annotation=="T cells CD4"),
    LAG3Tregn = sum(LAG3exp=="LAG3+" & annotation=="T reg cells"),
    
    # I didnt include TIM3 in this analysis 
    # TIM3Tcelln = sum(TIM3exp=="TIM3+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    # TIM3CD8Tcelln = sum(TIM3exp=="TIM3+" & annotation=="T cells CD8"),
    # TIM3CD4Tcelln = sum(TIM3exp=="TIM3+" & annotation=="T cells CD4"),
    # TIM3Tregn = sum(TIM3exp=="LAG3+" & annotation=="T reg cells"),
    
    #CD44, MHCII, CD86, Sirpa (IL1b)  on myeloid cells
    CD44Maccelln = sum(CD44exp=="CD44+" & annotation %in% c("Macrophages type 1","Macrophages type 2")),
    CD44Mac1n = sum(CD44exp=="CD44+" & annotation=="Macrophages type 1"),
    CD44Mac2n = sum(CD44exp=="CD44+" & annotation=="Macrophages type 2"),
    CD44DCcelln = sum(CD44exp=="CD44+" & annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    CD44DCn = sum(CD44exp=="CD44+" & annotation=="Dendritic cells"),
    CD44cDC1n = sum(CD44exp=="CD44+" & annotation=="Dendritic cells CD103"),
    
    CD86Maccelln = sum(CD86exp=="CD86+" & annotation %in% c("Macrophages type 1","Macrophages type 2")),
    CD86Mac1n = sum(CD86exp=="CD86+" & annotation=="Macrophages type 1"),
    CD86Mac2n = sum(CD86exp=="CD86+" & annotation=="Macrophages type 2"),
    CD86DCcelln = sum(CD86exp=="CD86+" & annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    CD86DCn = sum(CD86exp=="CD86+" & annotation=="Dendritic cells"),
    CD86cDC1n = sum(CD86exp=="CD86+" & annotation=="Dendritic cells CD103"),
    
    MHCIIMaccelln = sum(MHCIIexp=="MHCII+" & annotation %in% c("Macrophages type 1","Macrophages type 2")),
    MHCIIMac1n = sum(MHCIIexp=="MHCII+" & annotation=="Macrophages type 1"),
    MHCIIMac2n = sum(MHCIIexp=="MHCII+" & annotation=="Macrophages type 2"),
    MHCIIDCcelln = sum(MHCIIexp=="MHCII+" & annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    MHCIIDCn = sum(MHCIIexp=="MHCII+" & annotation=="Dendritic cells"),
    MHCIIcDC1n = sum(MHCIIexp=="MHCII+" & annotation=="Dendritic cells CD103"),
    
    PDL1Maccelln = sum(PDL1exp=="PDL1+" & annotation %in% c("Macrophages type 1","Macrophages type 2")),
    PDL1Mac1n = sum(PDL1exp=="PDL1+" & annotation=="Macrophages type 1"),
    PDL1Mac2n = sum(PDL1exp=="PDL1+" & annotation=="Macrophages type 2"),
    PDL1DCcelln = sum(PDL1exp=="PDL1+" & annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    PDL1DCn = sum(PDL1exp=="PDL1+" & annotation=="Dendritic cells"),
    PDL1cDC1n = sum(PDL1exp=="PDL1+" & annotation=="Dendritic cells CD103")
    
    #casp3 ,Ki67,CD44 on tumor
    
    # SirpaMaccelln = sum(Sirpaexp=="Sirpa+" & annotation %in% c("Macrophages type 1","Macrophages type 2")),
    # SirpaMac1n = sum(Sirpaexp=="Sirpa+" & annotation=="Macrophages type 1"),
    # SirpaMac2n = sum(Sirpaexp=="Sirpa+" & annotation=="Macrophages type 2"),
    # SirpaDCcelln = sum(Sirpaexp=="Sirpa+" & annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    # SirpaDCn = sum(Sirpaexp=="Sirpa+" & annotation=="Dendritic cells"),
    # SirpacDC1n = sum(Sirpaexp=="Sirpa+" & annotation=="Dendritic cells CD103")
  ) %>%
  dplyr::mutate(
    # Each cell type in all cells
    Bcelln_totaln = Bcelln/totaln,
    Neutn_totaln = Neutn/totaln,
    Maccelln_totaln= Maccelln/totaln,
    Tcelln_totaln= Tcelln/totaln,
    DCcelln_totaln= DCcelln/totaln,
    Fibn_totaln = Fibn/totaln,
    
    
    # Each immune cell type % in  all immune cells 
    Bcelln_immunecells = Bcelln/immunecellsn,
    Neutn_immunecells = Neutn/immunecellsn,
    Maccelln_immunecells= Maccelln/immunecellsn,
    Tcelln_immunecells= Tcelln/immunecellsn,
    DCcelln_immunecells= DCcelln/immunecellsn,
    Fibn_immunecells = Fibn/immunecellsn,
    
    # Each sub cell type %  in each cell type 
    CD4n_Tcelln = CD4n/Tcelln,
    CD8n_Tcelln = CD8n/Tcelln,
    Mac2n_Maccelln = Mac2n/Maccelln,
    Mac1n_Maccelln = Mac1n/Maccelln,
    DCn_DCcelln = DCn/DCcelln,
    cDC1n_DCcelln = cDC1n/DCcelln,
    
    # Marker expression relative to cell type
    PD1Tcelln_Tcelln = PD1Tcelln/Tcelln,
    PD1CD8Tcelln_CD8n =  PD1CD8Tcelln/CD8n,
    PD1CD4Tcelln_CD4n = PD1CD4Tcelln/CD4n,
    PD1Tregn_Tregn = PD1Tregn/Tregn,
    LAG3Tcelln_Tcelln = PD1Tcelln/Tcelln,
    LAG3CD8Tcelln_CD8n = LAG3CD8Tcelln/CD8n,
    LAG3CD4Tcelln_CD4n = LAG3CD4Tcelln/CD4n, 
    LAG3Tregn_Tregn = LAG3Tregn/Tregn,
    # TIM3Tcelln_Tcelln = TIM3Tcelln/Tcelln,
    # TIM3CD8Tcelln_CD8n = TIM3CD8Tcelln/CD8n,
    # TIM3CD4Tcelln_CD4n = TIM3CD4Tcelln/CD4n, 
    # TIM3Tregn_Tregn = TIM3Tregn/Tregn,
    
    CD86DCn_DCn = CD86DCn/DCn,   
    CD86DCcelln_DCcelln =  CD86DCcelln/DCcelln,
    CD86cDC1n_cDC1n =  CD86cDC1n/cDC1n,
    MHCIIDCn_DCn = MHCIIDCn/DCn,   
    MHCIIDCcelln_DCcelln =  MHCIIDCcelln/DCcelln,
    MHCIIcDC1n_cDC1n =  MHCIIcDC1n/cDC1n,
    PDL1DCn_DCn = PDL1DCn/DCn,   
    PDL1DCcelln_DCcelln =  PDL1DCcelln/DCcelln,
    PDL1cDC1n_cDC1n =  PDL1cDC1n/cDC1n,
    
    MHCIIMaccelln_Maccelln =MHCIIMaccelln/Maccelln,            
    MHCIIMac1n_Mac1n = MHCIIMac1n/Mac1n,              
    MHCIIMac2n_Mac2n = MHCIIMac2n/Mac2n, 
    CD86Maccelln_Maccelln =CD86Maccelln/Maccelln,            
    CD86Mac1n_Mac1n = CD86Mac1n/Mac1n,              
    CD86Mac2n_Mac2n = CD86Mac2n/Mac2n,   
    PDL1Maccelln_Maccelln =PDL1Maccelln/Maccelln,            
    PDL1Mac1n_Mac1n = PDL1Mac1n/Mac1n,              
    PDL1Mac2n_Mac2n = PDL1Mac2n/Mac2n,   
    CD44Maccelln_Maccelln =CD44Maccelln/Maccelln,            
    CD44Mac1n_Mac1n = CD44Mac1n/Mac1n,              
    CD44Mac2n_Mac2n = CD44Mac2n/Mac2n
  )

filenm = paste(output_path, "stat_summary_means_20_community",  ".csv", sep = "")
write.csv(stat_summary_means_20_community, filenm)




##Lag3 quantity within communities, keep this code because of level of factor ,seperate adjust the x discrete. 
# limits = c("MRTX+PD1", "MRTX+PD1+CTLA-4")
# LAG3_summary_community=Rphenograph20_240627[Rphenograph20_240627$ExpGroup %in% limits,]%>%
#   group_by(community,ExpGroup,ROI_ID.x)%>%
#   dplyr::summarise(LAG3cells = sum(LAG3exp=="LAG3+"),
#                    Tcelln = sum(annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
#                    LAG3Tcelln = sum(LAG3exp=="LAG3+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
#                    immunecellsn = sum(annotation %in% c("B cells","T cells CD4","T cells CD8","Dendritic cells","T reg cells","Neutrophils","Macrophages type 1","Macrophages type 2","Dendritic cells CD103"))
#   )%>%
#   dplyr::mutate(LAG3Tcelln_Tcelln=LAG3Tcelln/Tcelln)
# filenm = paste(output_path, "LAG3_summary_community",  ".csv", sep = "")
# write.csv(LAG3_summary_community, filenm)
# # Assuming the community variable should be a factor with specific levels
# LAG3_summary_community$community <- factor(LAG3_summary_community$community, 
#                                            levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))



##240627 try function MarkerExpWithinCommunity
output_folder <- file.path(output_path, "/marker_within_community/")
MarkerExpWithinCommunity <- function(data, title, y_axis, y_label,filename) {
  p <- ggplot(data, aes(x = community, y = y_axis, fill = ExpGroup)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Hide outliers since we use points
    geom_point(aes(color = ExpGroup), position = position_dodge(width = 0.75), size = 1) +
    scale_x_discrete(limits = as.character(1:21)) +
    scale_y_continuous(labels = percent) +
    theme_classic() + 
    labs(
      title = title,
      x = "Community",
      y = y_label,
      fill = "Experimental Group",
      color = "Experimental Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  print(p)
  output_file <- file.path(output_folder, paste0(filename, ".png"))
  # Save the plot
  ggsave(
    plot = p,
    filename = output_file,
    device = "png",
    width =10 ,  # Adjust as needed
    height = 4,  # Adjust as needed
    dpi = 300
  )
}
stat_summary_means_20_community[is.na(stat_summary_means_20_community)] <- 0
str(stat_summary_means_20_community)
stat_summary_means_20_community$community <- as.factor(stat_summary_means_20_community$community)
stat_summary_means_20_community$ExpGroup <- as.factor(stat_summary_means_20_community$ExpGroup)
stat_summary_means_20_community$ROI_ID.x <- as.factor(stat_summary_means_20_community$ROI_ID.x)
# Example usage:
MarkerExpWithinCommunity( 
  data = PD1_summary_community,
  y_axis= PD1_summary_community$PD1Tcelln_Tcelln,
  title = "PD1+ T cell in T cell",
  y_label = "PD1+ T cell % in T cells",
  filename = "PD1+ T cell percentage within community"
)

# check MHCII,PDL1,CD86 expression within communities
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MHCIIcDC1n_cDC1n,
  title = "MHCII+cDC1 out of cDC1 cells within communities",
  y_label = "MHCII+cDC1% in cDC1 cells",
  filename = "MHCII+cDC1 out of cDC1 cells within communities"
)

MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$PD1CD4Tcelln_CD4n,
  title = "PD1+ CD4 out of CD4 T cells within communities",
  y_label = "PD1CD4Tcelln_CD4n",
  filename = "PD1+ CD4 out of CD4 T cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$PD1CD8Tcelln_CD8n,
  title = "PD1+ CD8 out of CD8 T cells within communities",
  y_label = "PD1CD8Tcelln_CD8n",
  filename = "PD1+ CD8 out of CD8 T cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$PD1CD8Tcelln_CD8n,
  title = "PD1+ CD8 out of CD8 T cells within communities",
  y_label = "PD1CD8Tcelln_CD8n",
  filename = "PD1+ CD8 out of CD8 T cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$PDL1cDC1n_cDC1n,
  title = "PDL1+cDC1 out of cDC1 cells within communities",
  y_label = "PDL1+cDC1% in cDC1 cells",
  filename = "PDL1+cDC1 out of cDC1 cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MHCIIcDC1n_cDC1n,
  title = "MHCIIcDC1 out of cDC1 cells within communities",
  y_label = "MHCII+cDC1% in cDC1 cells",
  filename = "MHCII+cDC1 out of cDC1 cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$CD86cDC1n_cDC1n,
  title = "CD86cDC1 out of cDC1 cells within communities",
  y_label = "CD86+cDC1% in cDC1 cells",
  filename = "CD86+cDC1 out of cDC1 cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$CD86DCn_DCn,
  title = "CD86DC out of DC cells within communities",
  y_label = "CD86+ DC % in DC cells",
  filename = "CD86+ DC out of DC cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MHCIIDCn_DCn,
  title = "MHCIIDC out of DC cells within communities",
  y_label = "MHCII+ DC % in DC cells",
  filename = "MHCII+ DC out of DC cells within communities"
)
## Mean intensity not percentage  . so set a new function 
MarkerExp_MI_WithinCommunity <- function(data, title, y_axis, y_label,filename) {
  p <- ggplot(data, aes(x = community, y = y_axis, fill = ExpGroup)) +
    geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Hide outliers since we use points
    geom_point(aes(color = ExpGroup), position = position_dodge(width = 0.75), size = 1) +
    scale_x_discrete(limits = as.character(1:21)) +
    scale_y_continuous() +
    theme_classic() + 
    labs(
      title = title,
      x = "Community",
      y = y_label,
      fill = "Experimental Group",
      color = "Experimental Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  print(p)
  output_file <- file.path(output_folder, paste0(filename, ".png"))
  # Save the plot
  ggsave(
    plot = p,
    filename = output_file,
    device = "png",
    width =10 ,  # Adjust as needed
    height = 4,  # Adjust as needed
    dpi = 300
  )
}
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MI_PD1,
  title = "MI_PD1 within community",
  y_label = "MI_PD1",
  filename = "MI_PD1 within community"
)
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MI_MHCII,
  title = "MI_MHCII within community",
  y_label = "MI_MHCII",
  filename = "MI_MHCII within community"
)
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MI_CD86,
  title = "MI_CD86 within community",
  y_label = "MI_CD86",
  filename = "MI_CD86 within community"
)
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_20_community,
  y_axis= stat_summary_means_20_community$MI_LAG3,
  title = "MI_LAG3 within community",
  y_label = "MI_LAG3",
  filename = "MI_LAG3 within community"
)

## Heatmap showing specific marker expression level 
### need to be optimized , showing beautiful images, also gray means over the limitation , also need to think 
heatmap_CD8andTreg <- stat_summary_means_20_community %>%
  pivot_longer(cols = starts_with("CD4n_totaln") | starts_with("CD8n_totaln") | 
                 starts_with("DCn_totaln") | starts_with("cDC1n_totaln") | 
                 starts_with("Mac2n_totaln") | starts_with("Mac1n_totaln") | 
                 starts_with("Tregn_totaln") | starts_with("Bcelln_totaln") | 
                 starts_with("Neutn_totaln") | starts_with("Maccelln_totaln") | 
                 starts_with("Tcelln_totaln") | starts_with("DCcelln_totaln"),
               names_to = "CellType", values_to = "Proportion")

p <- ggplot(heatmap_CD8andTreg, aes(x = community, y = CellType, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.1, limit = c(0, 0.2), space = "Lab", name = "Proportion") +
 # scale_x_discrete(scale_x_discrete(limits = as.character(1:21))) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Proportion of Each Cell Type per Community", x = "Community", y = "Cell Type")
print(p)
ggsave(filename = paste(output_path, "heatmap_CD8andTreg.pdf", sep = ""), plot = p, width = 10, height = 8)
ggsave(filename = paste(output_path, "heatmap_CD8andTreg.pdf", sep = ""), plot = p, width = 10, height = 8)


# barplot for CN4 and CN7 percentage  (all cells)

colours =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
             "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
             "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
             "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")
create.resultsdir <- function(path, name) {
  full_path <- file.path(path, name)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
  }
  return(full_path)
}





###heatmap per community
ROI_list=c("01_MRTX+PD1","02_MRTX+PD1","03_MRTX+PD1+CTLA-4" ,"04_MRTX+PD1+CTLA-4",
           "05_MRTX+PD1+CTLA-4" ,"07_MRTX+PD1+CTLA-4" ,"08_MRTX+PD1" )
heatmap_plot <- ggplot(Rphenograph20_240627, aes(x = factor(community), y = factor(ROI_ID.x), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Community", y = "Treatment Group", fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## CD8 T cell count 
# CD8 T cell in each community
CD8Tcell_count <- stat_summary_means_20_community %>%
  group_by(community, ExpGroup) %>%
  dplyr::summarise(CD8_community = sum(CD8n))

# CD8 t cell in each Group
total_CD8Tcell_count <- stat_summary_means_20_community %>%
  group_by(ExpGroup) %>%
  dplyr::summarise(total_CD8_eachgroup= sum(CD8n))
#calculate percentage
CD8Tcell_percentage <- CD8Tcell_count %>%
  left_join(total_CD8Tcell_count, by = "ExpGroup") %>%
  dplyr::mutate(percentage = (CD8_community / total_CD8_eachgroup))

filtered_data1<- CD8Tcell_percentage[CD8Tcell_percentage$community %in% c("1", "12","6") & CD8Tcell_percentage$ExpGroup == "MRTX+PD1", ]
total_percentage <- sum(filtered_data1$percentage, na.rm = TRUE)
print(total_percentage)

filtered_data2<- CD8Tcell_percentage[CD8Tcell_percentage$community %in% c("1", "2","9") & CD8Tcell_percentage$ExpGroup == "MRTX+PD1+CTLA-4", ]
total_percentage <- sum(filtered_data2$percentage, na.rm = TRUE)
print(total_percentage)


print(CD8Tcell_percentage)
ggplot(CD8Tcell_count,aes(x=community,y=CD8_community,colour=ExpGroup))+geom_boxplot()







