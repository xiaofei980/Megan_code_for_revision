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
# stat_summary_means_18=Rphenograph20_240627[Rphenograph20_240627$ExpGroup %in% limits,]%>%
#   group_by(community,ExpGroup,ROI_ID.x) %>%
  stat_summary_means_20 = Rphenograph20_240627 %>%
  group_by(ExpGroup,ROI_ID.x) %>%
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

filenm = paste(output_path, "stat_summary_means_20",  ".csv", sep = "")
write.csv(stat_summary_means_20, filenm)

# plot.boxplot function : wilcox test, percentage

plot.boxplot <- function(dataframe, colname, comparisons, title = "", subdir = "",
                         limits = c("MRTX+PD1", "MRTX+PD1+CTLA-4"),
                         autosave = T, percentage = T, binwidth = NULL) {
  symcolumn = sym(colname)
  directory = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Results/Figures_output/Organized_Figures/Boxplot/"
  if (subdir != "") {
    directory = file.path(directory, subdir)
    if (!dir.exists(directory)) dir.create(directory)
  }
  if (title == "") title = colname
  
  if (is.null(binwidth)) {
    # Calculate binwidth as 1/30 of the range of the data if not provided
    data_range <- range(dataframe[[colname]], na.rm = TRUE)
    binwidth <- (data_range[2] - data_range[1]) / 30
  }
  
  p = ggplot(dataframe[dataframe$ExpGroup %in% limits,],
             aes(x = ExpGroup, y = !!symcolumn, fill = ExpGroup)) +
    stat_compare_means(comparisons = comparisons, label = "p.signif",method = "wilcox.test") +
    geom_boxplot() +
    geom_dotplot(method = "histodot", binaxis='y', fill = "white",
                 stackdir='center', dotsize = 1.5) + scale_x_discrete(limits=limits) +
    { if (percentage) {
      list(scale_y_continuous(labels = percent, expand = expansion(mult = c(0.03, 0.1))))}
    } + scale_fill_manual(values = brewer.pal(n = 4, name = "Pastel1")) +
    theme_classic() + ggtitle(title) +
    theme(plot.title = element_text(size = 11, hjust = 0.5, face = "bold", vjust = 0),
          plot.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 15))
  
  if (autosave) {
    ggsave(filename = paste0("Boxplot_", colname, ".jpeg"), path = directory,
           bg = "white", width = 1300, height = 1500, units = "px")
  }
  
  return(p)
}
## Making plot test
plot.boxplot(dataframe = stat_summary_means_20, colname = "CD4n_Tcelln", comparisons = list(my_comparisons),
             title = "CD4 T cells % in T cells",limits = my_comparisons,autosave = T)
## making plot according to requirement 
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "CD8n_immunecells", comparisons = list(my_comparisons),
             title = "CD8 T cells % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Tregn_immunecells", comparisons = list(my_comparisons),
             title = "T reg cells % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "cDC1n_immunecells", comparisons = list(my_comparisons),
             title = "cDC1 % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Mac2n_immunecells", comparisons = list(my_comparisons),
             title = "Macrophages type 2 % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Mac1n_immunecells", comparisons = list(my_comparisons),
             title = "Macrophages type 1 % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Bcelln_immunecells", comparisons = list(my_comparisons),
             title = "B cells % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "DCn_immunecells", comparisons = list(my_comparisons),
             title = "DC others % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Neutn_immunecells", comparisons = list(my_comparisons),
             title = "Neutrophils % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Tcelln_immunecells", comparisons = list(my_comparisons),
             title = "T cells % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Maccelln_immunecells", comparisons = list(my_comparisons),
             title = "Macrophages % in immune cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "DCcelln_immunecells", comparisons = list(my_comparisons),
             title = "Dendritic cells % in immune cells", limits = my_comparisons,autosave = T)


plot.boxplot(dataframe = stat_summary_means_18_group, colname = "CD4n_Tcelln", comparisons = list(my_comparisons),
             title = "CD4 percentage in T cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "CD8n_Tcelln", comparisons = list(my_comparisons),
             title = "CD8 percentage in T cells", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Mac2n_Maccelln", comparisons = list(my_comparisons),
             title = "Mac2 percentage in Macrophages", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "Mac1n_Maccelln", comparisons = list(my_comparisons),
             title = "Mac1 percentage in Macrophages", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "DCn_DCcelln", comparisons = list(my_comparisons),
             title = "DC others percentage in DCs", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "cDC1n_DCcelln", comparisons = list(my_comparisons),
             title = "cDC1 percentage in DCs", limits = my_comparisons,autosave = T)

plot.boxplot(dataframe = stat_summary_means_18_group, colname = "CD86DCn_DCn", comparisons = list(my_comparisons),
             title = "CD86 DC others % in DC others",limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "CD86cDC1n_cDC1n", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "MHCIIDCn_DCn", comparisons = list(my_comparisons),
             title = "MHCII DC others % in DC others", limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "MHCIIcDC1n_cDC1n", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "PDL1DCn_DCn", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "PDL1cDC1n_cDC1n", comparisons = list(my_comparisons),
             title = "PDL1 DC others % in DC others", limits = my_comparisons,autosave = T)

plot.boxplot(dataframe = stat_summary_means_18_group, colname = "PD1Tcelln_Tcelln", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "LAG3Tcelln_Tcelln", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "PD1CD8Tcelln_CD8n", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "LAG3CD8Tcelln_CD8n", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "PD1CD4Tcelln_CD4n", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)
plot.boxplot(dataframe = stat_summary_means_18_group, colname = "LAG3CD4Tcelln_CD4n", comparisons = list(my_comparisons),
             limits = my_comparisons,autosave = T)


# plot.boxplot(dataframe = stat_summary_means, colname = "CD4n_immunecells", comparisons = list(my_comparisons),
#              title = "CD4 T cells percentage in immune cells", limits = my_comparisons)ggplot(stat_summary_means, aes(x=ExpGroup, y= CD4n_immunecells, colour = ROI_ID.x)) + 
#   geom_point() +
#   scale_y_continuous(labels = scales::percent) +
#  scale_x_discrete(limits=c("MRTX+PD1", "MRTX+PD1+CTLA-4"))


##
date <- Sys.Date()
output_file = paste0(output_path, date, "_Rphenograph20_MarkerExp.csv")
output_file
write.csv(stat_summary_means_20, file = output_file, row.names = F) 
