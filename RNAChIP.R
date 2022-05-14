
##Incorporar los datos ChIP con lo de RNAseq

##FEMALE

library(readxl)
DEGs_Female <- read_excel('/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/DEGs_Female.xlsx')
Female_Common <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/female_Regulated_Peaks_Genes.tsv", header=TRUE)

names(Female_Common)[4] = "Gene ID" #cambiar el nombre de la columna para juntar las dos bases de datos

FemRNAChIP <- merge(x = DEGs_Female, y = Female_Common, by = c("Gene ID" ))


##MALE

DEGs_Male <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/DEGs_Male.xlsx")
Male_Common <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/male_Regulated_Peaks_Genes.tsv", header=TRUE)

names(Male_Common)[4] = "Gene ID" 

MaleRNAChIP <- merge(x = DEGs_Male, y = Male_Common, by = c("Gene ID" ))

##============================================================================================
# Hay que hacer una tabla que nos indique cuantos se co-regulan, down con down, up con up y mixtos

##Nos fijamos en el signo de la columna logFC de cada uno

##FEMALE
##el primero es que aparece UP en el RNA

FemUPUP <- FemRNAChIP[FemRNAChIP$logFC.x > 0 & FemRNAChIP$logFC.y > 0,] #9
FemUPDOWN <- FemRNAChIP[FemRNAChIP$logFC.x > 0 & FemRNAChIP$logFC.y < 0,] #21
FemDOWNUP <- FemRNAChIP[FemRNAChIP$logFC.x < 0 & FemRNAChIP$logFC.y > 0,] #1
FemDOWNDOWN <- FemRNAChIP[FemRNAChIP$logFC.x < 0 & FemRNAChIP$logFC.y < 0,] #63


##MALE
MaleUPUP <- MaleRNAChIP[MaleRNAChIP$logFC.x > 0 & MaleRNAChIP$logFC.y > 0,] #62
MaleUPDOWN <- MaleRNAChIP[MaleRNAChIP$logFC.x > 0 & MaleRNAChIP$logFC.y < 0,] #6
MaleDOWNUP <- MaleRNAChIP[MaleRNAChIP$logFC.x < 0 & MaleRNAChIP$logFC.y > 0,] #4
MaleDOWNDOWN <- MaleRNAChIP[MaleRNAChIP$logFC.x < 0 & MaleRNAChIP$logFC.y < 0,] #49

outFolder <- '/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/'

write.table(x = MaleUPUP$`Gene name`, file = paste0(outFolder,"Male_UPUPsymbol.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


##Ahora mismo solo estamos interesados en MALES y en UPUP y DOWNDOWN, y es por ello
##que los juntamos los dos en una misma base de datos

Males <- merge(x = MaleUPUP, y = MaleDOWNDOWN, all = TRUE)
write.table(x = Males, file = paste0(outFolder,"MalesRNAChIP.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



#crear la nueva columna entrezID


Males$entrez = mapIds(org.Mm.eg.db,
                     keys=Males$"Gene ID", #Column containing Ensembl gene ids
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


Males <- Males %>% distinct()

library(dplyr)

Males <- data.frame(Males)

MalesEntrezid <- dplyr::select(Males, Gene.ID, entrez)
MalesFC <- dplyr::select(Males,  Gene.ID, logFC.x)


write_tsv(MalesEntrezid,file = '/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/Male_entrez.tsv',col_names = FALSE)
write_tsv(MalesFC,file = '/media/sequentia/synology_office/Naiara/H3K4me3/RNAseq/Male_FC.tsv',col_names = FALSE)













