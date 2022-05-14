##Filtrar los picos que van unidos a regiones reguladoras (promotor, 5' utr, 1 intron y 1 exon)
#para saber cuantos son

library("ChIPseeker")
library(clusterProfiler)
library(scuttle)
library(reshape)
library(hrbrthemes)
library(dplyr)
library(devtools)

outFolder <- "/media/sequentia/synology_office/Naiara/H3K4me3/FilterPeaks/"


UniqueMale <- readPeakFile("/media/sequentia/synology_office/Naiara/H3K4me3/UniquePeaks/Male_PeaksAnnotation.tsv")
UniqueMale <- as.data.frame(UniqueMale)
df1 <- subset(UniqueMale,substr(UniqueMale$annotation, 1, 8) == "Promoter" |  substr(UniqueMale$annotation, 1, 1) == "5" )
df2 <- subset(UniqueMale, grepl("intron 1 of", UniqueMale$annotation) )
df3 <- subset(UniqueMale, grepl("exon 1 of", UniqueMale$annotation))
Funimale <- merge_recurse(list(df1, df2, df3))


UniqueFemale <- readPeakFile("/media/sequentia/synology_office/Naiara/H3K4me3/UniquePeaks/Female_PeaksAnnotation.tsv")
UniqueFemale <- as.data.frame(UniqueFemale)
df1 <- subset(UniqueFemale,substr(UniqueFemale$annotation, 1, 1) == "P" |  substr(UniqueFemale$annotation, 1, 1) == "5" )
df2 <- subset(UniqueFemale, grepl("intron 1 of", UniqueFemale$annotation) )
df3 <- subset(UniqueFemale, grepl("exon 1 of", UniqueFemale$annotation))
Funifemale <- merge_recurse(list(df1, df2, df3))


CommonMale <- readPeakFile("/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/Male_PeaksAnnotation.tsv")
CommonMale <- as.data.frame(CommonMale)
df1 <- subset(CommonMale,substr(CommonMale$annotation, 1, 1) == "P" |  substr(CommonMale$annotation, 1, 1) == "5" )
df2 <- subset(CommonMale, grepl("intron 1 of", CommonMale$annotation) )
df3 <- subset(CommonMale, grepl("exon 1 of", CommonMale$annotation))
Fcommale <- merge_recurse(list(df1, df2, df3))
write.table(x = Fcommale, file = paste0(outFolder,"Male_Common.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(Fcommale, SYMBOL, GENENAME, )
write.table(subset1, file = paste0(outFolder,"Male_Common_SYMBOLGENENAME.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(Fcommale, seqnames, start,end,width,scoreRate)
write.table(subset1, file = paste0(outFolder,"Male_Common_IGV.bed",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



CommonFemale <- readPeakFile("/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/Female_PeaksAnnotation.tsv")
CommonFemale <- as.data.frame(CommonFemale)
df1 <- subset(CommonFemale,substr(CommonFemale$annotation, 1, 1) == "P" |  substr(CommonFemale$annotation, 1, 1) == "5" )
df2 <- subset(CommonFemale, grepl("intron 1 of", CommonFemale$annotation) )
df3 <- subset(CommonFemale, grepl("exon 1 of", CommonFemale$annotation))
Fcomfemale <- merge_recurse(list(df1, df2, df3))
write.table(x = Fcomfemale, file = paste0(outFolder,"Female_Common.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(Fcomfemale, SYMBOL, GENENAME, )
write.table(subset1, file = paste0(outFolder,"Female_Common_SYMBOLGENENAME.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(Fcomfemale, seqnames, start,end,width,scoreRate)
write.table(subset1, file = paste0(outFolder,"Female_Common_IGV.bed",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



#venn diagram
library(VennDiagram)

#MALES

Genescontrol <- subset(Funimale, Funimale$Classification == "Control")
Genestreated <- subset(Funimale, Funimale$Classification == "Treated")

venn.diagram(
  x = list(Genescontrol$ENSEMBL, Fcommale$ENSEMBL, Genestreated$ENSEMBL),
  category.names = c("Genescontrol" , "Genesencomun " , "Genestreated"), filename = 'venn_diagrammMALES.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c("#D8B70A", "#FDD262", "#DD8D29"))

##Tomar solo los genes que son unicos, es decir, que solo esten en Treated o Control
x <- c()

for (i in 1:nrow(Genescontrol)) {
  if (sum(Fcommale$ENSEMBL== Genescontrol$ENSEMBL[i], na.rm = TRUE) == 0 & sum(Genestreated$ENSEMBL== Genescontrol$ENSEMBL[i], na.rm = TRUE) == 0) {
    x <- append(x, i)
  }
}

MaControl <- Genescontrol[x, ]
write.table(x = MaControl$SYMBOL, file = paste0(outFolder,"Male_ControlSYMBOL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(x = MaControl, file = paste0(outFolder,"Male_Control.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

subset1 <- select(MaControl, Classification, ENSEMBL)
write.table(subset1, file = paste0(outFolder,"Male_Control_ENSEMBL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(MaControl, ENSEMBL, GENENAME)
write.table(subset1, file = paste0(outFolder,"Male_Control_GENENAME.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(MaControl, seqnames, start,end,width,scoreRate)
write.table(subset1, file = paste0(outFolder,"Male_Control_IGV.bed",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



x <- c()
for (i in 1:nrow(Genestreated)) {
  if (sum(Fcommale$ENSEMBL== Genestreated$ENSEMBL[i], na.rm = TRUE) == 0 & sum(Genescontrol$ENSEMBL== Genestreated$ENSEMBL[i], na.rm = TRUE) == 0 ) {
    x <- append(x, i)
  }
}

MaTreated <- Genestreated[x, ]
write.table(x = MaTreated$SYMBOL, file = paste0(outFolder,"Male_TreatedSYMBOL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(x = MaTreated, file = paste0(outFolder,"Male_Treated.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

subset1 <- select(MaTreated, Classification, ENSEMBL)
write.table(subset1, file = paste0(outFolder,"Male_Treated_ENSEMBL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(MaTreated, ENSEMBL, GENENAME)
write.table(subset1, file = paste0(outFolder,"Male_Treated_GENENAME.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(MaTreated, seqnames, start,end,width,scoreRate)
write.table(subset1, file = paste0(outFolder,"Male_Treated_IGV.bed",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)




#FEMALES

Genescontrol <- subset(Funifemale, Funifemale$Classification == "Control")
Genestreated <- subset(Funifemale, Funifemale$Classification == "Treated")

venn.diagram(
  x = list(Genescontrol$ENSEMBL, Fcomfemale$ENSEMBL, Genestreated$ENSEMBL),
  category.names = c("Genescontrol" , "Genesencomun " , "Genestreated"), filename = 'venn_diagrammFEMALES.png',
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c("#D8B70A", "#FDD262", "#DD8D29"))


##Tomar solo los genes que son unicos, es decir, que solo esten en Treated o Control
x <- c()

for (i in 1:nrow(Genescontrol)) {
  if (sum(Fcomfemale$ENSEMBL== Genescontrol$ENSEMBL[i], na.rm = TRUE) == 0 & sum(Genestreated$ENSEMBL== Genescontrol$ENSEMBL[i], na.rm = TRUE) == 0) {
    x <- append(x, i)
  }
}

FemControl <- Genescontrol[x, ]
write.table(x = FemControl$SYMBOL, file = paste0(outFolder,"Female_ControlSYMBOL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(x = FemControl, file = paste0(outFolder,"Female_Control.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

subset1 <- select(FemControl, Classification, ENSEMBL)
write.table(subset1, file = paste0(outFolder,"Female_Control_ENSEMBL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(FemControl, ENSEMBL, GENENAME)
write.table(subset1, file = paste0(outFolder,"Female_Control_GENENAME.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(FemControl, seqnames, start,end,width,scoreRate)
write.table(subset1, file = paste0(outFolder,"Female_Control_IGV.bed",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



x <- c()
for (i in 1:nrow(Genestreated)) {
  if (sum(Fcomfemale$ENSEMBL== Genestreated$ENSEMBL[i], na.rm = TRUE) == 0 & sum(Genescontrol$ENSEMBL== Genestreated$ENSEMBL[i], na.rm = TRUE) == 0) {
    x <- append(x, i)
  }
}

FemTreated <- Genestreated[x, ]
write.table(FemTreated$SYMBOL, file = paste0(outFolder,"Female_TreatedSYMBOL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(x = FemTreated, file = paste0(outFolder,"Female_Treated.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


subset1 <- select(FemTreated, Classification, ENSEMBL)
write.table(subset1, file = paste0(outFolder,"Female_Treated_ENSEMBL.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(FemTreated, ENSEMBL, GENENAME)
write.table(subset1, file = paste0(outFolder,"Female_Treated_GENENAME.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)
subset1 <- select(FemTreated, seqnames, start,end,width,scoreRate)
write.table(subset1, file = paste0(outFolder,"Female_Treated_IGV.bed",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

#==========================================================================================
## Hay que hacer un ggplot del Fold Enrichment de los picos unicos de males y females


p <- ggplot(MaControl, aes(x= FE_control)) + geom_histogram(alpha=.9,binwidth = 0.05,bins = 500,color="lightblue", fill="black") 
p

p <- ggplot(MaTreated, aes(x= FE_control)) + geom_histogram(alpha=.9,binwidth = 0.05,bins = 500,color="lightblue", fill="black") 
p

p <- ggplot(FemControl, aes(x= FE_control)) + geom_histogram(alpha=.9,binwidth = 0.05,bins = 500,color="lightblue", fill="black") 
p

p <- ggplot(FemTreated, aes(x= FE_control)) + geom_histogram(alpha=.9,binwidth = 0.05,bins = 500,color="lightblue", fill="black") 
p







