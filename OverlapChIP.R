#El overlap de adultos y neonatales de ChIP
library(readxl)

NeoChIP <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP/ChIPNeoMales.tsv", header = TRUE)
AdultosChIP <- read.delim("/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP/ChIPAdultMales.tsv", header = TRUE)


NeoAdultChIP <- merge(x = NeoChIP, y = AdultosChIP, by = c("GeneID" ))

MaleUPUP <- NeoAdultChIP[NeoAdultChIP$logFC.x > 0 & NeoAdultChIP$logFC.y > 0 ,] #410

outFolder <- '/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP/'

write.table(x = MaleUPUP$Symbol.x, file = paste0(outFolder,"Male_UPUPsymbol.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

library("writexl")

write.table(x = MaleUPUP, file = paste0(outFolder,"NeoAdultChIP.tsv",sep="")
            ,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE)

write_xlsx(MaleUPUP, "/media/sequentia/synology_office/Naiara/H3K4me3/Overlap_ChIP/NeoAdultChIP.xlsx")























