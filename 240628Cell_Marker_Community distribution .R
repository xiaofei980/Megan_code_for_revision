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
path = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Data/Figures_input/"
output_path="/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Results/Figures_output/"

# Name of input and output directories:
data1 = read.csv(paste0(path, "/240528celldata.csv"))
neighb1 = read.csv(paste0(path, "/20240528_allROIs_25px_Iris_neighbours_distanceCalculation_allCells.csv"))
Rphenograph21 = read.csv(paste0(path, "Rphenograph_Megan_28may_output_21clusters_k250_13ct_fractions.csv"),check.names=F)
#colnames(Rphenograph21["cluster"])<-"community"
#colnames(Rphenograph21[18])<-community
Rphenograph21$community = Rphenograph21$cluster

data1[which(data1$annotation == "Macrophages type1"), "annotation"] = "Macrophages type 1"
data1[which(data1$annotation == "Macrophages type2"), "annotation"] = "Macrophages type 2"
Rphenograph21[which(Rphenograph21$source_cluster == "Macrophages type1"), "source_cluster"] = "Macrophages type 1"
Rphenograph21[which(Rphenograph21$source_cluster == "Macrophages type2"), "source_cluster"] = "Macrophages type 2"
Rphenograph21$cell_ID = Rphenograph21$source_ID
data1[!(data1$cell_ID %in% Rphenograph21$cell_ID),]

Rphenograph21_2 = inner_join(x=data1, y=Rphenograph21, by= "cell_ID")
Rphenograph21 = Rphenograph21_2
rm(Rphenograph21_2)


colours1 =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")

p = ggplot(Rphenograph21, aes(
  x = reorder(factor(community), source_cluster == "T cells CD8", sum),
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
 geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution") +
 scale_fill_manual(values = colours1) 
p


#### per treatment group
p = ggplot(Rphenograph21[which(Rphenograph21$ExpGroup == "MRTX+PD1+CTLA-4"),], aes(
  x = reorder(factor(community), source_cluster == "T cells CD8", sum),
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
  #geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution") +
  geom_bar(stat = "count") + ylab("Cell count") + 
  scale_fill_manual(values = colours1) +
  ggtitle("MRTX+PD1+CTLA-4")
p

filename = "Community per cell type_MRTX+PD1+CTLA-4"
ggsave(
  plot = p,
  device = "png",
  width = 6.5,  
  height = 7.5,
  dpi = 300,
  path = output_path,
  filename = paste0(filename, ".png")
)

p = ggplot(Rphenograph21[which(Rphenograph21$ExpGroup == "MRTX+PD1"),], aes(
  x = reorder(factor(community), source_cluster == "T cells CD8", sum),
  fill = as.factor(source_cluster)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
  #geom_bar(position = "fill") + scale_y_continuous(labels = scales::percent) + ylab("Percentage distribution") +
  geom_bar(stat = "count") + ylab("Cell count") + 
  scale_fill_manual(values = colours1) +
  ggtitle("MRTX+PD1")
p

filename = "Community per cell type_MRTX+PD1"
ggsave(
  plot = p,
  device = "png",
  width = 6.5,  
  height = 7.5,
  dpi = 300,
  path = output_path,
  filename = paste0(filename, ".png")
)


# #barplot per ccount
# p = ggplot(Rphenograph21, aes(
#   x = reorder(factor(community), source_cluster == "T cells CD8", sum),
#   fill = as.factor(source_cluster)
# )) +
#   theme_classic() +
#   coord_flip() +
#   theme(
#     axis.title.x = element_text(size = 20),
#     axis.title.y = element_text(size = 20),
#     axis.text = element_text(size = 16)
#   ) +
#   theme(legend.text = element_text(size = 16)) +
#   scale_x_discrete(limits = rev) +
#   xlab("Community") +
#   labs(fill = "") +
#   geom_bar(stat = "count") + ylab("Cell count") + 
#   scale_fill_manual(values = colours1)
# p
# filename="barplot_per_count"
# 
# ggsave(
#   plot = p,
#   device = "png",
#   width = 6.5,
#   height = 7.5,
#   dpi = 300,
#   units="in",
#   path = output_path,
#   filename = paste0(filename, ".png")
# )

p = ggplot(Rphenograph21, aes(
  x = reorder(factor(community), source_cluster == "T cells CD8", sum),
  fill = as.factor(ExpGroup)
)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  theme(legend.text = element_text(size = 16)) +
  scale_x_discrete(limits = rev) +
  xlab("Community") +
  labs(fill = "") +
  
p
filename = "Community per treatment_"
ggsave(
  plot = p,
  device = "png",
  width = width,
  height = height,
  dpi = 300,
  path = output_path,
  filename = paste0(filename, ".png")
)



#compare exhaustion markers within community 7 between treatment groups - redo with counts above threshold
ggplot(Rphenograph21[which(Rphenograph21$community == "7" & Rphenograph21$annotation.x == "T cells CD8"),], aes(x= MI_CD8a, y=MI_PD1, colour=ExpGroup)) +
  geom_point()


# Check location of the communities

p = ggplot(Rphenograph21, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(size = 0.2, color = "black") +
  geom_point(data = Rphenograph21[which(Rphenograph21$community %in% c("7")),], size = 0.1, color = "red") +
  geom_point(data = Rphenograph21[which(Rphenograph21$community %in% c("5")),], size = 0.1, color = "blue") +
  geom_point(data = Rphenograph21[which(Rphenograph21$community %in% c("8")),], size = 0.1, color = "yellow") +
  geom_point(data = Rphenograph21[which(Rphenograph21$community %in% c("4")),], size = 0.1, color = "Green") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank()) +
  facet_wrap(~ ROI_ID.x, nrow = 2, ncol = 6) +
  ggtitle("CN7Red CN5Blue CN8Yellow CN4Green")
p

filename = "CN7_Red CN5_Blue CN8_Yellow CN4_Green"
ggsave(plot = p, device = "png", width=12, height=5, dpi=300, path = output_path, filename = filename)

#  
# p
p = ggplot(Rphenograph21, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(size = 0.2, color = "black") +
  geom_point(data = Rphenograph21[which(Rphenograph21$annotation.x %in% c("Tumour")),], size = 0.2, color = "red") +
  geom_point(data = Rphenograph21[which(Rphenograph21$annotation.x %in% c("Macrophages type 2")),], size = 0.2, color = "blue") +
  # geom_point(data = Rphenograph21[which(Rphenograph21$community %in% c("8")),], size = 0.2, color = "yellow") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(size = 0.5),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank()) +
  facet_wrap(~ ROI_ID.x, nrow = 2, ncol = 6) +
  ggtitle("Tumour_Red M2_Blue")
p
filename = "Tumour_Red M2_Blue"
ggsave(plot = p, device = "png", width=12, height=5, dpi=300, path = output_path, filename = filename)



# unnamed_cols <- which(names(Rphenograph18) == "")
# 
# # Removing the unnamed columns
# Rphenograph18 <- Rphenograph18[, -unnamed_cols]


# `summarise()` regrouping output by 'MouseID', 'ROI_name3', 'treatment', 'domain2' 
#Rphenograph21 = subset(Rphenograph21, select = -c(62))
#2040531 xiaofei  calculate the percentage of cell types 
##set PD1 ,LAG3,PDL1 cutoff to define positive and negative  PD1   0.15 PDL1 0.5 LAG3 0.01 Casp3  0.8 Ki67 0.5 MHCII
Rphenograph18$PDL1exp <- ifelse(Rphenograph18$MI_PDL1 >= 0.5 , "PDL1+", "PDL1-")
Rphenograph18$Casp3exp <- ifelse(Rphenograph18$MI_Casp3 >= 0.8 , "Casp3+", "Casp3-") 
Rphenograph18$Ki67exp <- ifelse(Rphenograph18$MI_Ki67 >= 0.5 , "Ki67+", "Ki67-") 
Rphenograph18$PD1exp <- ifelse(Rphenograph18$MI_PD1 >= 0.15 , "PD1+", "PD1-") 
Rphenograph18$LAG3exp <- ifelse(Rphenograph18$MI_LAG3 >= 0.01 , "LAG3+", "LAG3-") 
Rphenograph18$CD44exp <- ifelse(Rphenograph18$MI_CD44 >= 0.4 , "CD44+", "CD44-") 
Rphenograph18$MHCIIexp <- ifelse(Rphenograph18$MI_MHCII >= 0.25 , "MHCII+", "MHCII-") 
Rphenograph18$CD86exp <- ifelse(Rphenograph18$MI_CD86 >= 0.25 , "CD86+", "CD86-") 


##Lag3 quantity within communities
limits = c("MRTX+PD1", "MRTX+PD1+CTLA-4")
LAG3_summary_community=Rphenograph18[Rphenograph18$ExpGroup %in% limits,]%>%
  group_by(community,ExpGroup,ROI_ID.x)%>%
  dplyr::summarise(LAG3cells = sum(LAG3exp=="LAG3+"),
            Tcelln = sum(annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
            LAG3Tcelln = sum(LAG3exp=="LAG3+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
            immunecellsn = sum(annotation %in% c("B cells","T cells CD4","T cells CD8","Dendritic cells","T reg cells","Neutrophils","Macrophages type 1","Macrophages type 2","Dendritic cells CD103"))
            )%>%
  dplyr::mutate(LAG3Tcelln_Tcelln=LAG3Tcelln/Tcelln)
filenm = paste(output_path, "LAG3_summary_community",  ".csv", sep = "")
write.csv(LAG3_summary_community, filenm)
# Assuming the community variable should be a factor with specific levels
LAG3_summary_community$community <- factor(LAG3_summary_community$community, 
                                           levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))


# Ensure community is a factor with the correct levels
LAG3_summary_community$community <- factor(LAG3_summary_community$community, 
                                           levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))

# Plotting LAG3 cells in each community
ggplot(LAG3_summary_community, aes(x = community, y = LAG3cells, fill = ExpGroup)) +
  geom_boxplot() + 
  geom_dotplot(aes(fill = ExpGroup), binaxis = 'y', stackdir = 'center', dotsize = 0.4) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
 # scale_y_continuous(labels = percent) +  # Correct usage of the percent function
  theme_classic() + 
  ggtitle("LAG3 cells in each community") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
# Plotting LAG3 T cells in each community
ggplot(LAG3_summary_community, aes(x = community, y = LAG3cells, fill = ExpGroup)) +
  geom_boxplot() + 
  geom_dotplot(aes(fill = ExpGroup), binaxis = 'y', stackdir = 'center', dotsize = 0.4) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
  # scale_y_continuous(labels = percent) +  # Correct usage of the percent function
  theme_classic() + 
  ggtitle("LAG3 T cells in each community") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    #legend.justification = c(1,1),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
# Plotting LAG3Tcelln_Tcelln  in each community
ggplot(LAG3_summary_community, aes(x = community, y = LAG3Tcelln_Tcelln, fill = ExpGroup)) +
  geom_boxplot() + 
  geom_dotplot(aes(fill = ExpGroup), binaxis = 'y', stackdir = 'center', dotsize = 0.4) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
  scale_y_continuous(labels = percent) +  # Correct usage of the percent function
  theme_classic() + 
  ggtitle("LAG3 T cells % in T cells in each community") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

# check PD1 expression within communities
limits = c("MRTX+PD1", "MRTX+PD1+CTLA-4")
PD1_summary_community=Rphenograph18[Rphenograph18$ExpGroup %in% limits,]%>%
  group_by(community,ExpGroup,ROI_ID.x)%>%
  dplyr::summarise(PD1cells = sum(PD1exp=="PD1+"),
                   Tcelln = sum(annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
                   PD1Tcelln = sum(PD1exp=="PD1+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
                   immunecellsn = sum(annotation %in% c("B cells","T cells CD4","T cells CD8","Dendritic cells","T reg cells","Neutrophils","Macrophages type 1","Macrophages type 2","Dendritic cells CD103"))
  )%>%
  dplyr::mutate(PD1Tcelln_Tcelln=PD1Tcelln/Tcelln)
filenm = paste(output_path, "PD1_summary_community",  ".csv", sep = "")
write.csv(PD1_summary_community, filenm)
# # Assuming the community variable should be a factor with specific levels
# PD1_summary_community$community <- factor(PD1_summary_community$community, 
#                                            levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))
# 
# 
# # Ensure community is a factor with the correct levels
# PD1_summary_community$community <- factor(PD1_summary_community$community, 
#                                            levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))
# 
# # Plotting PD1 cells in each community
# ggplot(PD1_summary_community, aes(x = community, y = PD1cells, fill = ExpGroup)) +
#   geom_boxplot() + 
#   geom_dotplot(aes(fill = ExpGroup), binaxis = 'y', stackdir = 'center', dotsize = 0.4) +
#   scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
#   # scale_y_continuous(labels = percent) +  # Correct usage of the percent function
#   theme_classic() + 
#   ggtitle("PD1 cells in each community") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
#     plot.background = element_blank(),
#     axis.title = element_blank(),
#     axis.text = element_text(size = 14, face = "bold"),
#     legend.position = "right",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 10)
#   )
# # Plotting PD1 T cells in each community
# ggplot(PD1_summary_community, aes(x = community, y = PD1cells, fill = ExpGroup)) +
#   geom_boxplot() + 
#   geom_dotplot(aes(fill = ExpGroup), binaxis = 'y', stackdir = 'center', dotsize = 0.4) +
#   scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
#   # scale_y_continuous(labels = percent) +  # Correct usage of the percent function
#   theme_classic() + 
#   ggtitle("PD1 T cells in each community") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
#     plot.background = element_blank(),
#     axis.title = element_blank(),
#     axis.text = element_text(size = 14, face = "bold"),
#     legend.position = "right",
#     #legend.justification = c(1,1),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 10)
#   )
# # Plotting PD1Tcelln_Tcelln  in each community
# ggplot(PD1_summary_community, aes(x = community, y = PD1Tcelln_Tcelln, fill = ExpGroup)) +
#   geom_boxplot() + 
#   geom_dotplot(aes(fill = ExpGroup), binaxis = 'y', stackdir = 'center', dotsize = 0.4) +
#   scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
#   scale_y_continuous(labels = percent) +  # Correct usage of the percent function
#   theme_classic() + 
#   ggtitle("PD1 T cells % in T cells in each community") +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
#     plot.background = element_blank(),
#     axis.title = element_blank(),
#     axis.text = element_text(size = 14, face = "bold"),
#     legend.position = "right",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 10)
#   )


#######
# Assuming PD1_summary_community is already created and loaded

# Plotting PD1 cells in each community
ggplot(PD1_summary_community, aes(x = community, y = PD1cells, fill = ExpGroup)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width=0.7),size=1)+
  geom_jitter(aes(color = ExpGroup), width = 0.2, size = 1.5) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
  theme_classic() + 
  ggtitle("PD1 cells in each community") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )

####
p=ggplot(PD1_summary_community, aes(x = community, y = PD1cells, fill = ExpGroup)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Hide outliers since we use points
  geom_point(aes(color = ExpGroup), position = position_dodge(width = 0.75),size = 1) +
  scale_x_discrete(limits = as.character(1:21)) +
  theme_classic() + 
  labs(
    title = "PD1 Cells in Each Community",
    x = "Community",
    y = "Number of PD1 Cells",
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
filename = "PD1 Cells in Each Community_"


# Create the folder if it doesn't exist
output_folder <- file.path(output_path, "/marker_within_community/")
if (!file.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
# Full path including the new folder
output_file <- file.path(output_folder, paste0(filename, ".png"))

# Assuming 'p' is your ggplot object
ggsave(
  plot = p,
  filename = output_file,
  device = "png",
  width = 12,
  height = 4,
  dpi = 300
)





# Plotting PD1 T cells in each community
p=ggplot(PD1_summary_community, aes(x = community, y = PD1Tcelln, fill = ExpGroup)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Hide outliers since we use points
  geom_point(aes(color = ExpGroup), position = position_dodge(width = 0.75),size = 1) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
  theme_classic() + 
  labs(
    title = "PD1 T cells in each community",
    x = "Community",
    y = "Number of PD1 T Cells",
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
p
filename = "PD1 T cells in each community_"
output_folder <- file.path(output_path, "/marker_within_community/")
# Full path including the new folder
output_file <- file.path(output_folder, paste0(filename, ".png"))

# Assuming 'p' is your ggplot object
ggsave(
  plot = p,
  filename = output_file,
  device = "png",
  width = 12,
  height = 4,
  dpi = 300
)

# Plotting PD1Tcelln_Tcelln in each community
p=ggplot(PD1_summary_community, aes(x = community, y = PD1Tcelln_Tcelln, fill = ExpGroup)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +  # Hide outliers since we use points
  geom_point(aes(color = ExpGroup), position = position_dodge(width = 0.75),size = 1) +
  scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
  scale_y_continuous(labels = percent) +
  theme_classic() + 
  labs(
    title = "PD1 T cells % in T cells in each community",
    x = "Community",
    y = "Number of PD1 T Cells",
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
filename = "PD1 T cells % in T cells in each community_"
output_folder <- file.path(output_path, "/marker_within_community/")
if (!file.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}
# Full path including the new folder
output_file <- file.path(output_folder, paste0(filename, ".png"))

# Assuming 'p' is your ggplot object
ggsave(
  plot = p,
  filename = output_file,
  device = "png",
  width = 12,
  height = 4,
  dpi = 300
)

##xiaofei try function 
output_folder <- file.path(output_path, "/marker_within_community/")
MarkerExpWithinCommunity <- function(data, title, y_axis, y_label,filename) {
  
  # Plotting PD1 T cells in each community
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
stat_summary_means_18[is.na(stat_summary_means_18)] <- 0
str(stat_summary_means_18)
stat_summary_means_18$community <- as.factor(stat_summary_means_18$community)
stat_summary_means_18$ExpGroup <- as.factor(stat_summary_means_18$ExpGroup)
stat_summary_means_18$ROI_ID.x <- as.factor(stat_summary_means_18$ROI_ID.x)
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
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MHCIIcDC1n_cDC1n,
  title = "MHCII+cDC1 out of cDC1 cells within communities",
  y_label = "MHCII+cDC1% in cDC1 cells",
  filename = "MHCII+cDC1 out of cDC1 cells within communities"
)

MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$PD1CD4Tcelln_CD4n,
  title = "PD1+ CD4 out of CD4 T cells within communities",
  y_label = "PD1CD4Tcelln_CD4n",
  filename = "PD1+ CD4 out of CD4 T cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$PD1CD8Tcelln_CD8n,
  title = "PD1+ CD8 out of CD8 T cells within communities",
  y_label = "PD1CD8Tcelln_CD8n",
  filename = "PD1+ CD8 out of CD8 T cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$PD1CD8Tcelln_CD8n,
  title = "PD1+ CD8 out of CD8 T cells within communities",
  y_label = "PD1CD8Tcelln_CD8n",
  filename = "PD1+ CD8 out of CD8 T cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$PDL1cDC1n_cDC1n,
  title = "PDL1+cDC1 out of cDC1 cells within communities",
  y_label = "PDL1+cDC1% in cDC1 cells",
  filename = "PDL1+cDC1 out of cDC1 cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MHCIIcDC1n_cDC1n,
  title = "MHCIIcDC1 out of cDC1 cells within communities",
  y_label = "MHCII+cDC1% in cDC1 cells",
  filename = "MHCII+cDC1 out of cDC1 cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$CD86cDC1n_cDC1n,
  title = "CD86cDC1 out of cDC1 cells within communities",
  y_label = "CD86+cDC1% in cDC1 cells",
  filename = "CD86+cDC1 out of cDC1 cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$CD86DCn_DCn,
  title = "CD86DC out of DC cells within communities",
  y_label = "CD86+ DC % in DC cells",
  filename = "CD86+ DC out of DC cells within communities"
)
MarkerExpWithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MHCIIDCn_DCn,
  title = "MHCIIDC out of DC cells within communities",
  y_label = "MHCII+ DC % in DC cells",
  filename = "MHCII+ DC out of DC cells within communities"
)

MarkerExp_MI_WithinCommunity <- function(data, title, y_axis, y_label,filename) {
  
  # Plotting PD1 T cells in each community
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
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MI_PD1,
  title = "MI_PD1 within community",
  y_label = "MI_PD1",
  filename = "MI_PD1 within community"
)
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MI_MHCII,
  title = "MI_MHCII within community",
  y_label = "MI_MHCII",
  filename = "MI_MHCII within community"
)
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MI_CD86,
  title = "MI_CD86 within community",
  y_label = "MI_CD86",
  filename = "MI_CD86 within community"
)
MarkerExp_MI_WithinCommunity( 
  data = stat_summary_means_18,
  y_axis= stat_summary_means_18$MI_LAG3,
  title = "MI_LAG3 within community",
  y_label = "MI_LAG3",
  filename = "MI_LAG3 within community"
)



##plot dots on all original image

ggplot(neighb1, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(data=neighb1,size = .1, color = "black") +
  geom_point(data = neighb1[which(neighb1$MI_PD1 >= 0.5),],
             #           | celldata$MI_PECAM >=0.375)),],      
             size = 0.005, colour = "yellow") +
  #geom_point(data = celldata[which(celldata$MI_CD11c>0.2&celldata$kinv33 =="23"),], size = 0.01, colour = "red")+
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank())+
  facet_wrap(~ROI_ID,nrow=3)

##plot dots on one specific original image

ROI13<-subset(neighb1, ROI_ID == "01_MRTX+PD1") # "04_MRTX+PD1+CTLA-4"
ggplot(ROI13, aes(x= Location_Center_X, y = -Location_Center_Y)) +
  geom_point(size = 0.5, color = "black") +
  geom_point(data = ROI13[which(ROI13$MI_CD86 >= 0.3),], size = 0.5, colour = "red")+
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_rect(),
        panel.spacing = unit(0, "mm"),
        plot.background = element_blank())+
  facet_wrap(~ROI_ID,nrow=1)











# cd44  mac dc fibroblast  MHCII  myeloid      pD1 lag3 on T cell 
limits = c("MRTX+PD1", "MRTX+PD1+CTLA-4")
stat_summary_means_18=Rphenograph18[Rphenograph18$ExpGroup %in% limits,]%>%
# stat_summary_means_18 = Rphenograph18 %>%
  group_by(community,ExpGroup,ROI_ID.x) %>%
  # stat_summary_means_18_group = Rphenograph18 %>%
  # group_by(ROI_ID.x,ExpGroup) %>%
 
  dplyr::summarise(
    MI_PD1 = mean(MI_PD1),
    MI_MHCII= mean(MI_MHCII),
    MI_PDL1 = mean(MI_PDL1),
    MI_LAG3 = mean(MI_LAG3),
    MI_CD86 = mean(MI_CD86),
    totaln = n(),
    PD1Tcelln = sum(PD1exp=="PD1+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    LAG3Tcelln = sum(LAG3exp=="LAG3+" & annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    PD1CD8Tcelln = sum(PD1exp=="PD1+" & annotation=="T cells CD8"),
    LAG3CD8Tcelln = sum(LAG3exp=="LAG3+" & annotation=="T cells CD8"),
    PD1CD4Tcelln = sum(PD1exp=="PD1+" & annotation=="T cells CD4"),
    LAG3CD4Tcelln = sum(LAG3exp=="LAG3+" & annotation=="T cells CD4"),
    PD1Bcelln = sum(PD1exp=="PD1+" & annotation=="B cells"),
    LAG3Bcelln = sum(LAG3exp=="LAG3+" & annotation=="B cells"),

    CD44Maccelln = sum(CD44exp=="CD44+" & annotation %in% c("Macrophages type 1","Macrophages type 2")),
    CD44Mac1n = sum(CD44exp=="CD44+" & annotation=="Macrophages type 1"),
    CD44Mac2n = sum(CD44exp=="CD44+" & annotation=="Macrophages type 2"),
    CD44DCcelln = sum(CD44exp=="CD44+" & annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    CD44DCn = sum(CD44exp=="CD44+" & annotation=="Dendritic cells"),
    CD44cDC1n = sum(CD44exp=="CD44+" & annotation=="Dendritic cells CD103"),
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
    PDL1cDC1n = sum(PDL1exp=="PDL1+" & annotation=="Dendritic cells CD103"),
    
    
    CD4n = sum(annotation=="T cells CD4"),
    CD8n = sum(annotation=="T cells CD8"),
    DCn = sum(annotation=="Dendritic cells"),
    cDC1n = sum(annotation=="Dendritic cells CD103"),
    Mac2n = sum(annotation=="Macrophages type 2"),
    Mac1n = sum(annotation=="Macrophages type 1"),
    Tregn = sum(annotation=="T reg cells"),
    Bcelln = sum(annotation=="B cells"),
    Neutn = sum(annotation=="Neutrophils"),
    Fibn = sum(annotation=="Fibroblasts"),
    Tcelln = sum(annotation %in% c("T cells CD4","T cells CD8","T reg cells")),
    DCcelln = sum(annotation %in% c("Dendritic cells","Dendritic cells CD103")),
    Maccelln = sum(annotation %in% c("Macrophages type 1","Macrophages type 2")),
    immunecellsn = sum(annotation %in% c("B cells","T cells CD4","T cells CD8","Dendritic cells","T reg cells","Neutrophils","Macrophages type 1","Macrophages type 2","Dendritic cells CD103"))
   ) %>%
  dplyr::mutate(CD4n_immunecells = CD4n/immunecellsn,
         CD8n_immunecells = CD8n/immunecellsn,
         DCn_immunecells = DCn/immunecellsn,
         cDC1n_immunecells = cDC1n/immunecellsn,
         Mac2n_immunecells = Mac2n/immunecellsn,
         Mac1n_immunecells = Mac1n/immunecellsn,
         Tregn_immunecells = Tregn/immunecellsn,
         Bcelln_immunecells = Bcelln/immunecellsn,
         Neutn_immunecells = Neutn/immunecellsn,
         Maccelln_immunecells= Maccelln/immunecellsn,
         Tcelln_immunecells= Tcelln/immunecellsn,
         DCcelln_immunecells= DCcelln/immunecellsn,
         Fibn_immunecells = Fibn/immunecellsn,
         CD44Maccelln_immunecells =  CD44Maccelln/immunecellsn,
         CD44Mac1n_immunecells = CD44Mac1n/immunecellsn, 
         CD44Mac2n_immunecells =CD44Mac2n/immunecellsn,
         CD44DCcelln_immunecells = CD44DCcelln/immunecellsn,
         CD44DCn_immunecells = CD44DCn/immunecellsn,
         CD44cDC1n_immunecells = CD44cDC1n/immunecellsn,
         MHCIIMaccelln_immunecells =MHCIIMaccelln/immunecellsn,            
         MHCIIMac1n_immunecells = MHCIIMac1n/immunecellsn,              
         MHCIIMac2n_immunecells = MHCIIMac2n/immunecellsn,         
         MHCIIDCcelln_immunecells = MHCIIDCcelln/immunecellsn,             
         MHCIIDCn_immunecells = MHCIIDCn/immunecellsn,           
         MHCIIcDC1n_immunecells =  MHCIIcDC1n/immunecellsn,  
         PDL1Maccelln_immunecells =PDL1Maccelln/immunecellsn,            
         PDL1Mac1n_immunecells = PDL1Mac1n/immunecellsn,              
         PDL1Mac2n_immunecells = PDL1Mac2n/immunecellsn,         
         PDL1DCcelln_immunecells = PDL1DCcelln/immunecellsn,             
         PDL1DCn_immunecells = PDL1DCn/immunecellsn, 
         PDL1cDC1n_immunecells =  PDL1cDC1n/immunecellsn,  
         PD1Tcelln_immunecells = PD1Tcelln/immunecellsn,
         LAG3Tcelln_immunecells = PD1Tcelln/immunecellsn,
         PD1CD8Tcelln_immunecells =  PD1CD8Tcelln/immunecellsn,
         LAG3CD8Tcelln_immunecells = LAG3CD8Tcelln/immunecellsn,
         PD1CD4Tcelln_immunecells = PD1CD4Tcelln/immunecellsn,
         LAG3CD4Tcelln_immunecells = LAG3CD4Tcelln/immunecellsn,
         
         CD4n_Tcelln = CD4n/Tcelln,
         CD8n_Tcelln = CD8n/Tcelln,
         Mac2n_Maccelln = Mac2n/Maccelln,
         Mac1n_Maccelln = Mac1n/Maccelln,
         DCn_DCcelln = DCn/DCcelln,
         cDC1n_DCcelln = cDC1n/DCcelln,
         PD1Tcelln_Tcelln = PD1Tcelln/Tcelln,
         LAG3Tcelln_Tcelln = PD1Tcelln/Tcelln,
         PD1CD8Tcelln_Tcelln =  PD1CD8Tcelln/Tcelln,
         LAG3CD8Tcelln_Tcelln = LAG3CD8Tcelln/Tcelln,
         PD1CD4Tcelln_Tcelln = PD1CD4Tcelln/Tcelln,
         LAG3CD4Tcelln_Tcelln = LAG3CD4Tcelln/Tcelln,
         MHCIIMaccelln_Maccelln =MHCIIMaccelln/Maccelln,            
         MHCIIMac1n_Maccelln = MHCIIMac1n/Maccelln,              
         MHCIIMac2n_Maccelln = MHCIIMac2n/Maccelln,         
         MHCIIDCcelln_DCcelln = MHCIIDCcelln/DCcelln,             
         MHCIIDCn_DCcelln = MHCIIDCn/DCcelln,           
         MHCIIcDC1n_DCcelln =  MHCIIcDC1n/DCcelln,
         CD86DCn_DCn = CD86DCn/DCn,             
         CD86cDC1n_cDC1n = CD86cDC1n/cDC1n,  
         
         CD86cDC1n_DCcelln =  CD86cDC1n/DCcelln,
         CD86DCn_DCn = CD86DCn/DCn,           
         CD86cDC1n_cDC1n =  CD86cDC1n/cDC1n,
         MHCIIDCn_DCn = MHCIIDCn/DCn,           
         MHCIIcDC1n_cDC1n =  MHCIIcDC1n/cDC1n,
         PDL1DCn_DCn = PDL1DCn/DCn,           
         PDL1cDC1n_cDC1n =  PDL1cDC1n/cDC1n,
         PD1Tcelln_Tcelln = PD1Tcelln/Tcelln,
         LAG3Tcelln_Tcelln = PD1Tcelln/Tcelln,
         PD1CD8Tcelln_CD8n =  PD1CD8Tcelln/CD8n,
         LAG3CD8Tcelln_CD8n = LAG3CD8Tcelln/CD8n,
         PD1CD4Tcelln_CD4n = PD1CD4Tcelln/CD4n,
         LAG3CD4Tcelln_CD4n = LAG3CD4Tcelln/CD4n,
        
        
   

         CD4n_totaln = CD4n/totaln,
         CD8n_totaln = CD8n/totaln,
         DCn_totaln = DCn/totaln,
         cDC1n_totaln = cDC1n/totaln,
         Mac2n_totaln= Mac2n/totaln,
         Mac1n_totaln = Mac1n/totaln,
         Tregn_totaln = Tregn/totaln,
         Bcelln_totaln = Bcelln/totaln,
         Neutn_totaln = Neutn/totaln,
         Maccelln_totaln= Maccelln/totaln,
         Tcelln_totaln= Tcelln/totaln,
         DCcelln_totaln= DCcelln/totaln
         )
filenm = paste(output_path, "stat_summary_means_18",  ".csv", sep = "")
write.csv(stat_summary_means_18, filenm)

# plot.boxplot function

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

plot.boxplot(dataframe = stat_summary_means_18_group, colname = "CD4n_immunecells", comparisons = list(my_comparisons),
                  title = "CD4 T cells % in immune cells", limits = my_comparisons,autosave = T)
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




# plot.boxplot(dataframe = stat_summary_means, colname = "Fibn_immunecells", comparisons = list(my_comparisons),
#              title = "Fibroblasts % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD44Maccelln_immunecells", comparisons = list(my_comparisons),
#              title = "CD44+Macrophages % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD44Mac1n_immunecells", comparisons = list(my_comparisons),
#              title = "CD44+ Mac1 % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD44Mac2n_immunecells", comparisons = list(my_comparisons),
#              title = "CD44+ Mac2 % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD44DCcelln_immunecells", comparisons = list(my_comparisons),
#              title = "CD44+ DC % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD44DCn_immunecells", comparisons = list(my_comparisons),
#              title = "CD44+ DC others % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD44cDC1n_immunecells", comparisons = list(my_comparisons),
#              title = "CD44+ cDC1 % of immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "MHCIIMaccelln_immunecells", comparisons = list(my_comparisons),
#              title = "MHCII+ Macrophages % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "MHCIIMac2n_immunecells", comparisons = list(my_comparisons),
#              title = "MHCII+ type2 Macrophages % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "MHCIIMac1n_immunecells", comparisons = list(my_comparisons),
#              title = "MHCII+ type1 Macrophages % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "MHCIIDCcelln_immunecells", comparisons = list(my_comparisons),
#              title = "MHCII+ DC % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "MHCIIcDC1n_immunecells", comparisons = list(my_comparisons),
#              title = "MHCII+ cDC1 % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "MHCIIDCn_immunecells", comparisons = list(my_comparisons),
#              title = "MHCII+ DC other % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "PDL1Maccelln_immunecells", comparisons = list(my_comparisons),
#              title = "PDL1+ Macrophages % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "PDL1Mac2n_immunecells", comparisons = list(my_comparisons),
#              title = "PDL1+ type2 Macrophages % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "PDL1Mac1n_immunecells", comparisons = list(my_comparisons),
#              title = "PDL1+ type1 Macrophages % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "PDL1DCcelln_immunecells", comparisons = list(my_comparisons),
#              title = "PDL1+ DC % in immune cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "PDL1cDC1n_immunecells", comparisons = list(my_comparisons),
#              title = "PDL1+ cDC1 % in immune cells", limits = my_comparisons,autosave = T)


# plot.boxplot(dataframe = stat_summary_means, colname = "Fibn_immunecells", comparisons = list(my_comparisons),
#              title = "CD4 T cells percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "CD8n_totaln", comparisons = list(my_comparisons),
#              title = "CD8 T cells percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Tregn_totaln", comparisons = list(my_comparisons),
#              title = "T reg cells percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "cDC1n_totaln", comparisons = list(my_comparisons),
#              title = "cDC1 percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Mac2n_totaln", comparisons = list(my_comparisons),
#              title = "Macrophages type 2 percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Mac1n_totaln", comparisons = list(my_comparisons),
#              title = "Macrophages type 1 percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Bcelln_totaln", comparisons = list(my_comparisons),
#              title = "B cells percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "DCn_totaln", comparisons = list(my_comparisons),
#              title = "DC others percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Neutn_totaln", comparisons = list(my_comparisons),
#              title = "Neutrophils percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Tcelln_totaln", comparisons = list(my_comparisons),
#              title = "T cells percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "Maccelln_totaln", comparisons = list(my_comparisons),
#              title = "Macrophages percentage in total cells", limits = my_comparisons,autosave = T)
# plot.boxplot(dataframe = stat_summary_means, colname = "DCcelln_totaln", comparisons = list(my_comparisons),
#              title = "Dendritic cells percentage in total cells", limits = my_comparisons,autosave = T)


# plot.boxplot(dataframe = stat_summary_means, colname = "CD4n_immunecells", comparisons = list(my_comparisons),
#              title = "CD4 T cells percentage in immune cells", limits = my_comparisons)ggplot(stat_summary_means, aes(x=ExpGroup, y= CD4n_immunecells, colour = ROI_ID.x)) + 
#   geom_point() +
#   scale_y_continuous(labels = scales::percent) +
#  scale_x_discrete(limits=c("MRTX+PD1", "MRTX+PD1+CTLA-4"))



##
date <- Sys.Date()
output_file = paste0(output_path, date, "_Rphenograph18_new.csv")
output_file
write.csv(Rphenograph21, file = output_file, row.names = F) 

### stacked barplot for CN4 and CN7 percentage  (all cells)
colours =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
             "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
             "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
             "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")
plot.celltypeperc <- function(celldata, title, colours) {
  date <- str_sub(gsub("-", "", Sys.Date()), -6, -1)
  
  #Proportion of cell types in tumour domain
  p <- ggplot(celldata, 
              aes(x = factor(community, levels = unique(community)),
                  fill = factor(annotation, levels = unique(annotation)))) +
    geom_bar(position = "fill", colour = "black") +  # You can change the bar border color here
    scale_fill_manual(values = colours) +  # Assigning specific colors to annotations
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, size = 14, hjust = 1), 
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 6), 
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.title = element_blank())
  
  filename <- paste0(date, "_", title, ".jpeg")
  ggsave(plot = p, device = "jpeg", width = 2000, height = 2000, units = "px", dpi = 300,
         path = file.path(output_path, "Figures/"), filename = filename)
  return(p)
}

# Assuming neighb1 is defined somewhere before this code
community_4and7 <- filter(neighb1, community == 4 | community == 7)

# Call the function with the specified colors
plot.celltypeperc(celldata = community_4and7, title = "Stacked_community4/7 % of total cells", colours = colours)

plot.celltypeperc(celldata = community_4and7[!community_4and7$annotation %in% c("Tumour"),],
                  title = "Stacked community4/7_without_tumour", colours = colours)
plot.celltypeperc(celldata = community_4and7[!community_4and7$annotation %in% c("Tumour", "Epithelium", "Fibroblasts", "Endothelium"),],
                  title = "Stacked_commmunity4/7 immune_cells", colours = colours)
plot.celltypeperc(celldata = neighb1, title = "Stacked all community of total cells", colours = colours)
plot.celltypeperc(celldata = neighb1[!neighb1$annotation %in% c("Tumour", "Epithelium", "Fibroblasts", "Endothelium"),],
                  title = "Stacked all community of immune cells", colours = colours)

##  show composition of  T  or T reg proportion per community  
###1  calculate T cell percentage
###2  reshape to long data  melt()
###3  geom plot(geom_tile )
stat_summary_CD8TandTreg = Rphenograph21 %>%
  group_by(community) %>%
  summarise(
    totaln = n(),
    CD4n = sum(annotation.x=="T cells CD4"),
    CD8n = sum(annotation.x=="T cells CD8"),
    DCn = sum(annotation.x=="Dendritic cells"),
    cDC1n = sum(annotation.x=="Dendritic cells CD103"),
    Mac2n = sum(annotation.x=="Macrophages type 2"),
    Mac1n = sum(annotation.x=="Macrophages type 1"),
    Tregn = sum(annotation.x=="T reg cells"),
    Bcelln = sum(annotation.x=="B cells"),
    Neutn = sum(annotation.x=="Neutrophils"),
    Fibn = sum(annotation.x=="Fibroblasts"),
    Tcelln = sum(annotation.x %in% c("T cells CD4","T cells CD8","T reg cells")),
    DCcelln = sum(annotation.x %in% c("Dendritic cells","Dendritic cells CD103")),
    Maccelln = sum(annotation.x %in% c("Macrophages type 1","Macrophages type 2"))
  ) %>%
  mutate(
         CD4n_totaln = CD4n/totaln,
         CD8n_totaln = CD8n/totaln,
         DCn_totaln = DCn/totaln,
         cDC1n_totaln = cDC1n/totaln,
         Mac2n_totaln= Mac2n/totaln,
         Mac1n_totaln = Mac1n/totaln,
         Tregn_totaln = Tregn/totaln,
         Bcelln_totaln = Bcelln/totaln,
         Neutn_totaln = Neutn/totaln,
         Maccelln_totaln= Maccelln/totaln,
         Tcelln_totaln= Tcelln/totaln,
         DCcelln_totaln= DCcelln/totaln
  )
filenm = paste(output_path, "stat_summary_CD8TandTreg",  ".csv", sep = "")
write.csv(stat_summary_CD8TandTreg, filenm,row.names = FALSE)

heatmap_CD8andTreg <- stat_summary_CD8TandTreg %>%
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
scale_x_discrete(limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "Proportion of Each Cell Type per Community", x = "Community", y = "Cell Type")

print(p)

ggsave(filename = paste(output_path, "heatmap_CD8andTreg.pdf", sep = ""), plot = p, width = 10, height = 8)


ggsave(filename = paste(output_path, "heatmap_CD8andTreg.pdf", sep = ""), plot = p, width = 10, height = 8)



###heatmap per community
ROI_list=c("01_MRTX+PD1","02_MRTX+PD1","03_MRTX+PD1+CTLA-4" ,"04_MRTX+PD1+CTLA-4",
           "05_MRTX+PD1+CTLA-4" ,"06_ETP","07_MRTX+PD1+CTLA-4" ,"08_MRTX+PD1" )
heatmap_plot <- ggplot(Rphenograph21, aes(x = factor(community), y = factor(ROI_ID.x), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Community", y = "Treatment Group", fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## CD8 T cell count 
# CD8 T cell in each community
CD8Tcell_count <- stat_summary_means_18 %>%
  group_by(community, ExpGroup) %>%
  dplyr::summarise(CD8_community = sum(CD8n))

# CD8 t cell in each Group
total_CD8Tcell_count <- stat_summary_means_18 %>%
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

# 查看结果
print(CD8Tcell_percentage)
ggplot(CD8Tcell_count,aes(x=community,y=CD8_community,colour=ExpGroup))+geom_boxplot()




