##Poder poner los simbolos de los genes en los archivos de GO y KEGG
#El procedimiento es el mismo para 8 archivos,  UPmalesunique, DOWNmalesunique, UPmalescommon, DOWNmalescommon y sus respectivos en hembras

outFolder <- "/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/"

library("org.Mm.eg.db")
library(dplyr)

## COMMON FEMALES 

GOEA_UPfemale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_UPfemale.filt.txt", header=FALSE)

l <- as.list(scan(text=GOEA_UPfemale.filt[,11], what="",  sep=';')) 
#creamos una lista con los identificadores de los genes

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
       keys=as.character(df$V1), #Column containing Ensembl gene ids
       column="SYMBOL",
       keytype="ENSEMBL",
       multiVals="first") #a partir de los identificadores de ENSEBL, creamos una base de datos que contiene el simbolo de cada uno

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3] #tomamos solo las columnas de ENSEMBL y SYMBOL

write.table(df2, file = paste0(outFolder,"GOFemaleUp.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


GOEA_DOWNfemale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_DOWNfemale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_DOWNfemale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOFemaleDOWN.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


## COMMON MALES

GOEA_UPmale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_UPmale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_UPmale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOmaleUp.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


GOEA_DOWNmale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_DOWNmale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_DOWNmale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOmaleDOWN.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)



## UNIQUE FEMALES

GOEA_Controlfemale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_Controlfemale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_Controlfemale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOFemaleControl.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


GOEA_Treatedfemale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_Treatedfemale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_Treatedfemale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOFemaleTreated.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)




 
## UNIQUE MALES


GOEA_Controlmale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_Controlmale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_Controlmale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOmaleControl.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


GOEA_Treatedmale.filt <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/GOEA_Treatedmale.filt.txt", header=FALSE)
l <- as.list(scan(text=GOEA_Treatedmale.filt[,11], what="",  sep=';'))

df <- as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))


df1 <- mapIds(org.Mm.eg.db,
              keys=as.character(df$V1), #Column containing Ensembl gene ids
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")

df1 <- as.data.frame(df1)

df2 <- merge(df,df1, by = "row.names")

df2 <- df2[!duplicated(df2$V1), ] ##los enseblen deben aparecer solo 1 vez
df2 <- df2[,2:3]

write.table(df2, file = paste0(outFolder,"GOmaleTreated.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

























