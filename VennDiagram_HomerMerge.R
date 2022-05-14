## This script is to build the replicates Venn Diagram Merge

library(VennDiagram)
library(wesanderson)
library(gridExtra)
# CM
draw.triple.venn(area1 = 26046, area2 = 24193	, area3 = 16666	, n12 = 22339	, n23 = 16552		, n13 = 16567	, n123 = 16508,alpha = c(0.5,0.5,0.5), category = c("CM1", "CM2", "CM3")
                 , fill = c("#D8B70A", "#FDD262", "#DD8D29") , cex=1.5, cat.cex=1.5, cat.fontfamily = rep("serif", 3))

# CF 
draw.triple.venn(area1 = 24378	,area2 =24956	, area3 = 22222	, n12 = 22098	, n23 = 21013	, n13 = 21179	, n123 = 20494, alpha = c(0.5,0.5,0.5), 
                 category =  c("CF1", "CF2", "CF3"),fill =c("#D8B70A", "#FDD262", "#DD8D29") , cex=1.5, cat.cex=1.5, cat.fontfamily = rep("serif", 3))

#PM

draw.triple.venn(area1 = 22546	, area2 = 15925	, area3 = 22255	, n12 = 15699	, n23 = 15695	, n13 = 20485	, n123 = 15615,alpha = c(0.5,0.5,0.5),
                 category = c("PM1", "PM2", "PM3"), fill = c("#D8B70A", "#FDD262", "#DD8D29") , cex=1.5, cat.cex=1.5, cat.fontfamily = rep("serif", 3))

# PF

draw.triple.venn(area1 = 21305	, area2 = 21938	, area3 = 23503	, n12 = 19875	, n23 = 20695	, n13 = 20379	, n123 = 19579,alpha = c(0.5,0.5,0.5), category = c("PF1", "PF2", "PF3")
                 , fill = c("#D8B70A", "#FDD262", "#DD8D29") , cex=1.5, cat.cex=1.5, cat.fontfamily = rep("serif", 3))

###Conditions



### CM vs PM (partial 2 or more replicates)

draw.pairwise.venn(area1 = 27422	,area2 =23782	, cross.area = 22443,alpha = c(0.5,0.5), category = c("CM", "PM"), fill = c("#D8B70A", "#02401B") , cex=1.5, cat.cex=1.5, cat.fontfamily = rep("serif", 2),
                          lty =2, filename = NULL)

### CF vs PF (partial 2 or more replicates)

draw.pairwise.venn(area1 = 26714	,area2 =25066	, cross.area = 23380,alpha = c(0.5,0.5), category = c("CF", "PF"),  fill = c("#02401B","#D8B70A"), cex=1.5, cat.cex=1.5, cat.fontfamily = rep("serif", 2),
                   lty =2, filename = NULL)
                   
