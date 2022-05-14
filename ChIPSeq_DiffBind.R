library(R.utils)
library(DiffBind)
library("writexl")
library("ChIPseeker")
library("org.Mm.eg.db")
library("clusterProfiler")
library(parallel)


TxDb <- "TxDb.Mmusculus.UCSC.mm10.knownGene" #El archivo que contiene el genoma de ratón

library(TxDb,character.only = TRUE)


#read config file
samples <- dba(sampleSheet="/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/malesdata_DiffBind.csv")

#let diffBind count number of reads on peaks
#Hacer conteo de los datos mirando en regiones iguales en todos los samples(C y T)
samples <- dba.count(samples,summits = FALSE)


# Normalize 
samples_norm <- dba.normalize(samples,method = DBA_EDGER,normalize = "TMM")
#los diferentes metodos de normalizacion: 
#DBA_EDGER, DBA_DESEQ2, DBA_ALL_METHODS

#let diffBind know you want to perform the analysis on conditions
samples_norm <- dba.contrast(samples_norm, contrast=c("Condition","P","C"))

#perform the differential analysis
#Performs default generation of a consensus peakset, read counting, 
#normalization, and setting up of contrasts if they have not been specified.
samples_norm <- dba.analyze(samples_norm,method=DBA_DESEQ2, bBlacklist=FALSE, bGreylist=FALSE)
#Para ver mas o menos marcas de metilalcion(H3K4me3) en ambos tratamientos ->
#para saber si tiene mas probabildad de que se regulen mas

dba.show(samples_norm, bContrasts=TRUE)

# Plot 
plot(samples_norm)
dev.off()

#get Differentially Binding Sites (DBS)
#dbs <- dba.report(samples)
dbs <- dba.report(samples_norm, th=.01, bUsePval=TRUE, fold=0)
#th -> significance threshold
#bUsePval	 -> logical indicating whether to use FDR (FALSE) or p-value (TRUE) for thresholding.
dba.plotVolcano(samples_norm,th = .01,bUsePval = TRUE, fold=0)



#peak annotation regarding gene context
peakAnnoDBS <- annotatePeak(
  dbs,
  tssRegion=c(-3000, 3000),
  TxDb=get(TxDb),
  annoDb="org.Mm.eg.db"
)

#En un intervalo que nosotros concretamos, en nuestro caso (-3000,3000), busca genes que estén
#annoDb	 -> annotation package

genesDBS <- as.data.frame(peakAnnoDBS)
write.table(genesDBS, file="~/DiffBindGenes_Treated_vs_Control_male.tsv",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)


genesDBS_down <- subset(genesDBS,Fold<(0))$geneId
genesDBS_up <- subset(genesDBS,Fold>0)$geneId

# length(genesDBS_down) -> 265


# length(genesDBS_up) -> 612


#compareCluster -> Given a list of gene set, this function will compute 
#profiles of each gene cluster.

compGO <- compareCluster(
  geneCluster=list("up"=unique(genesDBS_up), "down"=unique(genesDBS_down)),
  fun="enrichGO",
  pvalueCutoff=0.05,
  OrgDb=org.Mm.eg.db,
  ont="BP"
)

dotplot(compGO, showCategory = 10, title = "GO BP",font.size=8,label_format = 50)

#==================================================================================
##### FEMALES

#read config file

samples <- dba(sampleSheet="/media/sequentia/synology_office/Naiara/H3K4me3/CommonPeaks/femalesdata_DiffBind.csv")

#let diffBind count number of reads on peaks
Fsamples <- dba.count(samples,summits = FALSE)

# Normalize 
FNsamples <- dba.normalize(Fsamples,method = DBA_EDGER,normalize = "TMM")



#let diffBind know you want to perform the analysis on conditions

FNsamples <- dba.contrast(FNsamples, contrast=c("Condition","P","C"))
#minMembers	 -> when automatically generating contrasts, minimum number of unique 
#samples in a group. 

#perform the differential analysis
FNsamples <- dba.analyze(FNsamples,method=DBA_DESEQ2, bBlacklist=FALSE, bGreylist=FALSE)

dba.show(FNsamples, bContrasts=TRUE,)

plot(FNsamples, contrast=1,)
dev.off()


#get Differentially Binding Sites (DBS)
#dbs <- dba.report(samples)
dbs <- dba.report(FNsamples, th=.01, bUsePval=TRUE, fold=0)

dba.plotVolcano(FNsamples,th = .01,bUsePval = TRUE, fold=0)
dev.off()



#peak annotation regarding gene context
peakAnnoDBS <- annotatePeak(
  dbs,
  tssRegion=c(-3000, 3000),
  TxDb=get(TxDb),
  annoDb="org.Mm.eg.db"
)


genesDBS <- as.data.frame(peakAnnoDBS)
write.table(genesDBS, file="/home/sradio/DiffBindGenes_Treated_vs_Control_female.tsv",quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

genesDBS_down <- subset(genesDBS,Fold<(0))$geneId
genesDBS_up <- subset(genesDBS,Fold>0)$geneId


compGO <- compareCluster(
  geneCluster=list("up"=unique(genesDBS_up), "down"=unique(genesDBS_down)),
  fun="enrichGO",
  pvalueCutoff=0.05,
  OrgDb=org.Mm.eg.db,
  ont="BP"
)

dotplot(compGO, showCategory = 10, title = "GO BP",font.size=8,label_format = 50)





