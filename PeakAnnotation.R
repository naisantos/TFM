### Este script me permite obtener de los picos identificados por condicion los picos unicos y los comunes outliers, annotados con respecto al genoma

library("ChIPseeker")
TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene"
library(TxDb,character.only = TRUE)
library(TxDb,character.only = TRUE)
library(clusterProfiler)
library(org.Mm.eg.db)
mm <- org.Mm.eg.db
library(scuttle)
library(ggplot2)

library("writexl")
library(hrbrthemes)
library(dplyr)


outFolder <- "/media/sequentia/synology_office/Naiara/H3K4me3/UniquePeaks/"

#==========================================================================
# MALES

peakFile <- "/media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/peaks/blacklist/IGVM/Peaks.bed"

peakFile <- readPeakFile(peakFile)

#con esta funcion buscamos genes que estan a x distancia maxima de los picos,
#donde normalmente buscamos genes activadores
peakAnnotation <- annotatePeak(peakFile,tssRegion=c(-3000, 2000),TxDb=get(TxDb))

#Para anotar la ubicación de un pico dado en términos de características genómicas, 
#annotatePeakasigne picos a la anotación genómica en la columna de "anotación" de la 
#salida, que incluye si un pico está en TSS, Exon, 5' UTR, 3' UTR, Intronic o Intergénico.
plotAnnoBar(peakAnnotation,title = "Males")
# ggsave("/home/sradio/Intergenerational_project/ChIPseq_H3K4me3/IntersectionPeaks/v3/Male/Condition/PeakAnnotationPlot.png",width = 785,height = 284,units="mm")


peakAnnotation <- as.data.frame(peakAnnotation@anno)

entrezids <- unique(peakAnnotation$geneId) #tomar los que no estén repetidos


#Biological Id TRanslator -> bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
#(included in clusterProfiler)
gene.df <- bitr(entrezids, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL","GENENAME"),
                OrgDb = mm)
#ENSEMBL: ensembl ID of the nearest gene
#SYMBOL: gene symbol
#GENENAME: full gene name

#merge(x, y, # Data frames u objectos a ser transformados
# by.x = by, by.y = by, # Columnas usadas para unir)
#permite fusionar o unir dos data frames por columnas comunes o por nombres de fila
peakAnnotation <- merge(peakAnnotation,gene.df,by.x = "geneId", by.y = "ENTREZID")

#cambiar los nombres de algunas variables y borrar la que no nos interesa
peakAnnotation <- peakAnnotation %>%
  rename(PeakName = V4)

peakAnnotation$V5 <- NULL

peakAnnotation <- peakAnnotation %>%
  rename(scoreControl = V6)

peakAnnotation <- peakAnnotation %>%
  rename(scoreTreatment = V7)

peakAnnotation <- peakAnnotation %>%
  rename(scoreRate = V8)

peakAnnotation <- peakAnnotation %>%
  rename(FE_control = V9)

peakAnnotation <- peakAnnotation %>%
  rename(FE_treatment = V10)

peakAnnotation <- peakAnnotation %>%
  rename(FERate = V11)

peakAnnotation$V12 <- NULL


peakAnnotation <- peakAnnotation %>%
  rename(replicates_Control = V13)

peakAnnotation <- peakAnnotation %>%
  rename(replicates_Treated = V14)

peakAnnotation <- peakAnnotation %>%
  rename(Classification = V15)


peakAnnotation$geneId <- NULL
peakAnnotation$geneChr <- NULL
peakAnnotation$geneStart <- NULL
peakAnnotation$geneEnd <- NULL
peakAnnotation$geneLength <- NULL
peakAnnotation$geneStrand <- NULL
peakAnnotation$transcriptId <- NULL
peakAnnotation$strand <- NULL

UniquePeaks <- subset(peakAnnotation,Classification!="Common")
#hacer un subconjunto con los picos unicos


write.table(x = UniquePeaks, file = paste0(outFolder,"Male_PeaksAnnotation.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


#### Find outliers -> En el trabajo los dejamos a un lado

commonPeaks <- subset(peakAnnotation,Classification=="Common")
#Hacer un subconjutno con los picos que aparecen ambos tratamientos

write.table(x = commonPeaks, file = "/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/Male_PeaksAnnotation.tsv",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

#Convenience function to determine which values in a numeric vector are outliers 
#based on the median absolute deviation
outlier_values <- isOutlier(commonPeaks$FERate, type = "both")
commonPeaks$is_outlier <- as.logical(outlier_values) #para ver si son TRUE o FALSE

#Get or set specific attributes of an object.
#attr(x, which, exact = FALSE)
lower.value <- attr(outlier_values,"thresholds")[1] 
higher.value <- attr(outlier_values,"thresholds")[2]

# lower    higher 
#  0.694176  1.405824  



commonPeaks$Type <- ifelse(commonPeaks$FERate > 1, "Treated","Control")

# Outliers_RatioFE
p <- ggplot(commonPeaks, aes(x=log(FERate))) + geom_histogram(alpha=.9,binwidth = 0.05,bins = 500,color="black", fill="lightblue") + geom_segment(aes(x=log(lower.value),xend=log(lower.value),y=0,yend=100),show.legend = FALSE,colour="BLUE", linetype = "dotted")+ geom_segment(aes(x=log(higher.value),xend=log(higher.value),y=0,yend=100),show.legend = FALSE,colour="BLUE", linetype = "dotted")   + theme_ipsum(grid = FALSE) + theme(legend.position = "none") + ylab("# Peaks")

ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_OutlierHistogram.png",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_OutlierHistogram.pdf",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_OutlierHistogram.tiff",p,width = 1200,height = 1200,units = "px")


p <- ggplot(commonPeaks, aes(x=abs(log(FERate)), fill=Type)) + geom_density(alpha=.5) + theme_ipsum(grid = FALSE) +xlim(c(-1,3))

ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_OutlierDensity.png",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_OutlierDensity.pdf",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_OutlierDensity.tiff",p,width = 1200,height = 1200,units = "px")

commonPeaks <- subset(commonPeaks,is_outlier==TRUE) #tomar solo los outliers

commonPeaks$is_outlier <- NULL
 
write.table(x = commonPeaks, file = "/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Male_PeaksAnnotation.tsv",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


#========================================================================================
# FEMALE


outFolder <- "/media/sequentia/synology_office/Naiara/H3K4me3/UniquePeaks/"

peakFile <- "/media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/peaks/blacklist/IGVF/Peaks.bed"
peakFile <- readPeakFile(peakFile)

peakAnnotation <- annotatePeak(peakFile,tssRegion=c(-3000, 2000),TxDb=get(TxDb))

plotAnnoBar(peakAnnotation,title = "Females")
# ggsave("/home/sradio/Intergenerational_project/ChIPseq_H3K4me3/IntersectionPeaks/v3/Male/Condition/PeakAnnotationPlot.png",width = 785,height = 284,units="mm")

peakAnnotation <- as.data.frame(peakAnnotation@anno)

entrezids <- unique(peakAnnotation$geneId)
gene.df <- bitr(entrezids, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL","GENENAME"),
                OrgDb = mm)

peakAnnotation <- merge(peakAnnotation,gene.df,by.x = "geneId", by.y = "ENTREZID")

peakAnnotation <- peakAnnotation %>%
  rename(PeakName = V4)

peakAnnotation$V5 <- NULL

peakAnnotation <- peakAnnotation %>%
  rename(scoreControl = V6)

peakAnnotation <- peakAnnotation %>%
  rename(scoreTreatment = V7)

peakAnnotation <- peakAnnotation %>%
  rename(scoreRate = V8)

peakAnnotation <- peakAnnotation %>%
  rename(FE_control = V9)

peakAnnotation <- peakAnnotation %>%
  rename(FE_treatment = V10)

peakAnnotation <- peakAnnotation %>%
  rename(FERate = V11)

peakAnnotation$V12 <- NULL


peakAnnotation <- peakAnnotation %>%
  rename(replicates_Control = V13)

peakAnnotation <- peakAnnotation %>%
  rename(replicates_Treated = V14)

peakAnnotation <- peakAnnotation %>%
  rename(Classification = V15)


peakAnnotation$geneId <- NULL
peakAnnotation$geneChr <- NULL
peakAnnotation$geneStart <- NULL
peakAnnotation$geneEnd <- NULL
peakAnnotation$geneLength <- NULL
peakAnnotation$geneStrand <- NULL
peakAnnotation$transcriptId <- NULL
peakAnnotation$strand <- NULL

UniquePeaks <- subset(peakAnnotation,Classification!="Common") 
#Hacer un subconjunto con los picos unicos


write.table(x = UniquePeaks, file = paste0(outFolder,"Female_PeaksAnnotation.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


#### Find outliers 

commonPeaks <- subset(peakAnnotation,Classification=="Common")

write.table(x = commonPeaks, file = "/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/Female_PeaksAnnotation.tsv",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


outlier_values <- isOutlier(commonPeaks$FERate, type = "both")
commonPeaks$is_outlier <- as.logical(outlier_values)


lower.value <- attr(outlier_values,"thresholds")[1]
higher.value <- attr(outlier_values,"thresholds")[2]

# lower    higher 
#  0.556937 1.313063 



commonPeaks$Type <- ifelse(commonPeaks$FERate > 1, "Treated","Control")

# Outliers_RatioFE
p <- ggplot(commonPeaks, aes(x=log(FERate))) + geom_histogram(alpha=.9,binwidth = 0.05,bins = 500,color="black", fill="lightblue") + geom_segment(aes(x=log(lower.value),xend=log(lower.value),y=0,yend=500),show.legend = FALSE,colour="BLUE", linetype = "dotted")+ geom_segment(aes(x=log(higher.value),xend=log(higher.value),y=0,yend=500),show.legend = FALSE,colour="BLUE", linetype = "dotted")   + theme_ipsum(grid = FALSE) + theme(legend.position = "none") + ylab("# Peaks")

ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_OutlierHistogram.png",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_OutlierHistogram.pdf",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_OutlierHistogram.tiff",p,width = 1200,height = 1200,units = "px")


p <- ggplot(commonPeaks, aes(x=abs(log(FERate)), fill=Type)) + geom_density(alpha=.5) + theme_ipsum(grid = FALSE) +xlim(c(-1,3))

ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_OutlierDensity.png",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_OutlierDensity.pdf",p,width = 1200,height = 1200,units = "px")
ggsave("/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_OutlierDensity.tiff",p,width = 1200,height = 1200,units = "px")

commonPeaks <- subset(commonPeaks,is_outlier==TRUE)

commonPeaks$is_outlier <- NULL

write.table(x = commonPeaks, file = "/media/sequentia/synology_office/Naiara/H3K4me3/CommonEnrichedPeaks/Female_PeaksAnnotation.tsv",quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



