library("AnnotationDbi")
library("org.Mm.eg.db")
library(dplyr)
library(tidyr)
library(readr)

###Males

male <- read.delim2("/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/male_GroupedGenes.tsv", header=FALSE)

#Creamos una nueva columna llamada ENTREZID, que contiene numeros enteros como idenfitifcadores de genes
male$entrez = mapIds(org.Mm.eg.db,
                     keys=male$V1, #Column containing Ensembl gene ids
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


male$V2 <- NULL #eliminamos la columna de FC

male <- male %>% distinct() #Tomamos solamente los genes una única vez

male <- drop_na(male)

write_tsv(male,file = "/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/Male_ens2gb.tsv",col_names = FALSE)

###Females


female <- read.delim2("/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/female_GroupedGenes.tsv", header=FALSE)

female$entrez = mapIds(org.Mm.eg.db,
                     keys=female$V1, #Column containing Ensembl gene ids
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


female$V2 <- NULL

female <- female %>% distinct()

female <- drop_na(female)

write_tsv(female,file = "/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/Female_ens2gb.tsv",col_names = FALSE)
