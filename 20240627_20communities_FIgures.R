# compare  K200,250,300
library(ggplot2) 
library(dplyr)
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
base = "/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/"

# Name of input and output directories:
input_path = file.path(base, "Projectfile/Data/Figures_input/")
output_path="/mnt/Data1/groupfebe/runs/Xiaofei/Megan_IMC/Projectfile/Results/Figures_output/"


# Name of input and output directories:
data1 = read.csv(paste0(input_path, "/2024_06_26_celldata.csv"))
Rphenograph20_240626 = read.csv(paste0(input_path, "Rphenograph_Megan_26June_output_20clusters_k250_13ct_fractions.csv"),check.names=F)
unique(names(Rphenograph20_240626))
unique(data1$annotation)
unique(names(data1))

names(Rphenograph20_240626)[names(Rphenograph20_240626) == "cluster"] <- "community"
# data1[which(Rphenograph19$annotation == "Macrophages type1"), "annotation"] = "Macrophages type 1"

Rphenograph20_240626$cell_ID = Rphenograph20_240626$source_ID

data1[!(data1$cell_ID %in% Rphenograph20_240626$cell_ID),]
Rphenograph19_2 = inner_join(x=data1, y=Rphenograph20_240626, by= "cell_ID")
Rphenograph20_240626 = Rphenograph19_2
rm(Rphenograph19_2)
unique(names(Rphenograph20_240626))

write.csv(Rphenograph20_240626,file= paste0(output_path,"Rphenograph20_240626.csv"))
# Rphenograph19_new<-read.csv(file = paste0(output_path,"Rphenograhph19_new.csv"),check.names = F)
str(Rphenograph20_240626)
length(unique(Rphenograph20_240626$community))
colours1 =  c("B cells" = "#945931",  "Dendritic cells" = "#FF9933", "Dendritic cells CD103" = "#FFCC66", "Endothelium" = "#FFEA42",
              "Epithelium" = "#FFF4A4FF", "Fibroblasts" = "#ABDDA4FF", "Macrophages other" = "#618F75FF", "Macrophages type 1" = "#336666",
              "Macrophages type 2" = "#4E79A7FF", "Neutrophils" = "#A0CBE8FF", "NK cells" = "#938ABBFF", "T cells CD4" = "#B07AA1FF",
              "T cells CD8" = "#FF9D9AFF", "T reg cells" = "#CC6666", "Tumour" = "#704850FF", "Unclassified" = "#8C8C8C")
####### per treatment group--percentage 

p = ggplot(Rphenograph19_new, aes(
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
filename = "CD8_distribution_over_communities_Rphenograph19_new"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 6.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized Figures/",filename, ".png")
)
##MRTX+PD1+CTLA4 group percentage 
p = ggplot(Rphenograph19_new[which(Rphenograph19_new$ExpGroup == "MRTX+PD1+CTLA-4"),], aes(
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
filename = "CD8_distribution_over_communities_MRTX+PD1+CTLA-4_Rphenograph19_new"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 6.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized_Figures/",filename, ".png")
)
##MRTX+PD1 group percentage 
p = ggplot(Rphenograph19_new[which(Rphenograph19_new$ExpGroup == "MRTX+PD1"),], aes(
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
filename = "CD8_distribution_over_communities_MRTX+PD1_Rphenograph19_new"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 6.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized_Figures/",filename, ".png")
)


####Show how T cells are distributed among the communities


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
com_colours = getPalette(2)
com_colours = c("1" = "#E41A1C", "2" = "#9E425A", "3" = "#596A98", "4" = "#3B87A2", "5"= "#449B75", "6"= "#4DAF4A", "7"=  "#6B886D",
                "8" = "#896191", "9"= "#AC5782", "10" = "#D56B41", "11"= "#FF7F00", "12"= "#FFB214", "13"= "#FFE528", "14"= "#EDDD30",
                "15" = "#C9992C", "16"= "#A65628", "17"= "#C66764", "18" = "#E678A0", "19" = "#E485B7", "20"= "#BE8FA8", "21"= "#999999")
p = ggplot(Rphenograph18[which(Rphenograph31$annotation == "T cells CD8" & Rphenograph18$ExpGroup %in% c("MRTX+PD1")),], aes(
  x = ExpGroup,
  fill = reorder(factor(community), source_cluster == "T cells CD8", sum))) +
  theme_classic() +
  #coord_flip() +
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
  scale_fill_manual(values = com_colours)
p
filename = "CD8_distribution_over_communities_MRTX+PD1_18community"
ggsave(
  plot = p,
  device = "png",
  width = 4.5,  
  height = 6.5,
  dpi = 300,
  path = output_path,
  filename = paste0(filename, ".png")
)

#### per treatment group--cell count
p = ggplot(Rphenograph20_240626[which(Rphenograph20_240626$ExpGroup == "MRTX+PD1+CTLA-4"),], aes(
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

filename = "Community per cell type_MRTX+PD1+CTLA-4 Rphenograph20_240626"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 7.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized_Figures/",filename, ".png")
)

p = ggplot(Rphenograph20_240626[which(Rphenograph20_240626$ExpGroup == "MRTX+PD1"),], aes(
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

filename = "Community per cell type_MRTX+PD1 Rphenograph20_240626"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 7.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized_Figures/",filename, ".png")
)

#### making figures for top5 community in each group 
##MRTX+PD1 

##MRTX+PD1+CTLA4 group percentage 
Rphenograph20_240626_CTLA4ExpGroup<-Rphenograph20_240626[which((Rphenograph20_240626$ExpGroup == "MRTX+PD1+CTLA-4") & (Rphenograph20_240626$community %in% c("9","3","8","7","1"))) ,]
str(Rphenograph20_240626_CTLA4ExpGroup)
p = ggplot(Rphenograph20_240626_CTLA4ExpGroup[which(Rphenograph20_240626_CTLA4ExpGroup$ExpGroup == "MRTX+PD1+CTLA-4"),], aes(
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
filename = "top5_CD8_distribution_over_communities_MRTX+PD1+CTLA-4_Rphenograph20_240626_CTLA4ExpGroup"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 6.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized_Figures/",filename, ".png")
)
##MRTX+PD1 group percentage 
Rphenograph20_240626_PD1ExpGroup<-Rphenograph20_240626[which((Rphenograph20_240626$ExpGroup == "MRTX+PD1") & (Rphenograph20_240626$community %in% c("7","4","5","2","9"))) ,]
str(Rphenograph20_240626_PD1ExpGroup)
p = ggplot(Rphenograph20_240626_PD1ExpGroup[which(Rphenograph20_240626_PD1ExpGroup$ExpGroup == "MRTX+PD1"),], aes(
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
filename = "top5_CD8_distribution_over_communities_MRTX+PD1_Rphenograph20_240626_PD1ExpGroup"
ggsave(
  plot = p,
  device = "png",
  width = 8,  
  height = 6.5,
  dpi = 300,
  path = output_path,
  filename = paste0("/Organized_Figures/",filename, ".png")
)

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





