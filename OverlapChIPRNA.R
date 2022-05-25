#El overlap de adultos y neonatales de RNA y ChIP
library(readxl)

NeoRNAChIP <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP_RNA/MalesRNAChIP.xlsx")

AdultosRNAChIP <- read_excel("/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP_RNA/MalesTreated_vs_Control_UP.xlsx")
AdultosRNAChIP <- subset(AdultosRNAChIP, AdultosRNAChIP$seqnames != '-' &  AdultosRNAChIP$Fold != '-')
AdultosRNAChIP <- subset(AdultosRNAChIP, AdultosRNAChIP$RNASeq_p.value < 0.05)



names(AdultosRNAChIP)[1] = "Gene ID"

NeoAdultChIPRNA <- merge(x = NeoRNAChIP, y = AdultosRNAChIP, by = c("Gene ID" ))

MaleUPUP <- NeoAdultChIPRNA[NeoAdultChIPRNA$logFC.x > 0 & NeoAdultChIPRNA$RNASeq_logFC > 0 & NeoAdultChIPRNA$logFC.y > 0 & NeoAdultChIPRNA$Fold > 0 ,] #8
MaleUPDOWN <- NeoAdultChIPRNA[NeoAdultChIPRNA$logFC.x > 0 & NeoAdultChIPRNA$RNASeq_logFC < 0 & NeoAdultChIPRNA$logFC.y > 0 & NeoAdultChIPRNA$Fold < 0,] #2


outFolder <- '/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP_RNA/'


write.table(x = NeoAdultChIPRNA, file = paste0(outFolder,"NeoAdultChIPRNA.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

write.table(x = MaleUPUP$`Gene name`, file = paste0(outFolder,"Male_UPUPsymbol.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

library("writexl")

write.table(x = MaleUPUP, file = paste0(outFolder,"NeoAdultChIPRNAUPUP.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

write_xlsx(MaleUPUP, "/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP_RNA/NeoAdultChIPRNAUPUP.xlsx")
write_xlsx(NeoAdultChIPRNA, "/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP_RNA/NeoAdultChIPRNA.xlsx")


library(VennDiagram)
library(wesanderson)
library(gridExtra)


draw.pairwise.venn(area1 = 651, area2 = 111	, cross.area = 8,alpha = c(0.5,0.5), category = c("Adults ChIP and RNA-seq", "Neonates ChIP and RNA-seq"), 
                   fill = c("#D8B70A", "#FDD262") , cex=1, cat.cex=1, cat.fontfamily = rep("serif", 2), margin = 0.1)






############################################################################################

neoUp <- NeoRNAChIP[NeoRNAChIP$logFC.x > 0 & NeoRNAChIP$logFC.y > 0, ]
AduUP <- AdultosRNAChIP[AdultosRNAChIP$RNASeq_logFC> 0 & AdultosRNAChIP$Fold > 0, ]



UP <- merge(x = neoUp, y = AduUP, by = c("Gene ID" )) #8


AduDOWN <- AdultosRNAChIP[AdultosRNAChIP$RNASeq_logFC< 0 & AdultosRNAChIP$Fold < 0, ]
write.table(x = AduDOWN$Symbol, file = paste0(outFolder,"down.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)














