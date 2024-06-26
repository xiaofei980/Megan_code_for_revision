library(exiftoolr)
library(stringr)

# Script used to extract filename, height and width of images
base = "/mnt/Data1/imcyto/runs/Xiaofei/Megan_IMC/Projectfile/Data/"
image_path = file.path(base,"imagemetadata_input/full_stack.ome")
output_path = file.path(base, "Normalization_input/imagemetadata.csv")


image_metadata = exif_read(path = file.path(image_path, list.files(image_path)))
image_metadata = image_metadata[,c("FileName", "ImageWidth", "ImageHeight")]
colnames(image_metadata) = c("filename", "width", "height")
image_metadata$filename = str_remove(string = image_metadata$filename,
                                     pattern = "_full_stack.ome.tiff")
View(image_metadata)

image_metadata$filename <- gsub("CONTROL", "MRTX+PD1", image_metadata$filename )
image_metadata$filename <- gsub("TREATMENT", "MRTX+PD1+CTLA-4", image_metadata$filename )
write.csv(image_metadata, file = output_path, row.names = F)






# Script to remove unused images, and paste together split image
imagedata = image_metadata
index1 = grep(pattern = "MOC2_CCR2KO_2R_4_ROI1_x1", x = imagedata$filename)
index2 = grep(pattern = "MOC2_CCR2KO_2R_4_ROI1_x2", x = imagedata$filename)
width = imagedata$width[index1]
height = imagedata$height[index1] + imagedata$height[index2]
imagedata = rbind(imagedata, c("MOC2_CCR2KO_2R_4_ROI1", width, height))
imagedata = imagedata[order(imagedata$filename),]

remove.row = function(filenamelist, imagedata) {
  for (filename in filenamelist) {
    exclude = grepl(pattern = filename, x = imagedata[,"filename"])
    imagedata = imagedata[!exclude,]
  }l
  return(imagedata)
}

filenamelist = c("MOC1_MOCAF_1B_2_ROI1_1", "MOC2_CCR2KO_2R_4_ROI1_1",
                 "MOC2_WT_1R_4_ROI1_1", "MOC2_CCR2KO_2R_4_ROI1_x1", 
                 "MOC2_CCR2KO_2R_4_ROI1_x2" )

imagedata = remove.row(filenamelist = filenamelist, imagedata = imagedata)
output_path = "/home/imcyto/analysis/Xiaofei/Sahai_IMC/imagemetadataNEWADJ.csv"
write.csv(x = imagedata, file = output_path, row.names = F)
