##Poder poner los simbolos de los genes en los archivos de GO y KEGG
#Es el mismo procedimiento que con los GO

outFolder <- "/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/"

library("org.Mm.eg.db")
library(dplyr)
library(readxl)


## COMMON FEMALES 

KEGG_UPfemale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_UPfemale.xlsx")
KEGG_UPfemale <- as.data.frame(KEGG_UPfemale)

l <- as.list(scan(text=as.character(KEGG_UPfemale[,11]), what="",  sep=';'))


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

write.table(df2, file = paste0(outFolder,"UPfemaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

KEGG_DOWNfemale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_DOWNfemale.xlsx")
KEGG_DOWNfemale <- as.data.frame(KEGG_DOWNfemale)
l <- as.list(scan(text=as.character(KEGG_DOWNfemale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"DOWNfemaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)


## COMMON MALES

KEGG_UPmale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_UPmale.xlsx")
KEGG_UPmale <- as.data.frame(KEGG_UPmale)
l <- as.list(scan(text=as.character(KEGG_UPmale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"UPmaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)





KEGG_DOWNmale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_DOWNmale.xlsx")
KEGG_DOWNmale <- as.data.frame(KEGG_DOWNmale)
l <- as.list(scan(text=as.character(KEGG_DOWNmale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"DOWNmaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)




## UNIQUE FEMALE

KEGG_Controlfemale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_Controlfemale.xlsx")
KEGG_Controlfemale <- as.data.frame(KEGG_Controlfemale)
l <- as.list(scan(text=as.character(KEGG_Controlfemale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"ControlfemaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)





KEGG_Treatedfemale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_Treatedfemale.xlsx")
KEGG_Treatedfemale <- as.data.frame(KEGG_Treatedfemale)
l <- as.list(scan(text=as.character(KEGG_Treatedfemale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"TreatedfemaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)




## UNIQUE MALE

KEGG_Controlmale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_Controlmale.xlsx")
KEGG_Controlmale <- as.data.frame(KEGG_Controlmale)
l <- as.list(scan(text=as.character(KEGG_Controlmale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"ControlmaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)





KEGG_Treatedmale <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/GOKEGGnames/KEGG_Treatedmale.xlsx")
KEGG_Treatedmale <- as.data.frame(KEGG_Treatedmale)
l <- as.list(scan(text=as.character(KEGG_Treatedmale[,11]), what="",  sep=';'))

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

write.table(df2, file = paste0(outFolder,"TreatedmaleKEGG.tsv",sep=""),quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)











