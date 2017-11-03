source("http://bioconductor.org/biocLite.R")
biocLite("GOstats")
library("arabidopsis.db0")
biocLite("org.At.tair.db")
library("org.At.tair.db")
library("GO.db")
library("AnnotationDbi")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")
library("Biobase")
biocLite("Rgraphviz")
library(DBI)



dirname <- "~/Documents/Research/300/Twisst/"
#analysisname <- "DIPvsTETexpression"
analysisname <- "PutSelEvents"
datafile <-  "~/Documents/Research/300/AccessoryInfo/Gene_Universe.txt"
geneUniverse <- read.table( datafile, header=FALSE)
subsetname='twisst3c.topo3'
#subsetname='overTET'
datafile2 <-  "~/Documents/Research/300/Twisst/Twisst3c_topo3_0.99outliers_ThalIDs.txt"

geneSelected1 <- read.table( "~/Documents/Research/300/Twisst/Twisst3c_topo3_0.99outliers_ThalIDs.txt", header=FALSE)
geneSelected2 <- read.table( "~/Documents/Research/300/Twisst/Twisst3c_topo1_0.99outliers_ThalIDs.txt", header=FALSE)
geneSelected3 <- read.table( "~/Documents/Research/300/Twisst/Twisst4c_topo3_0.99outliers_ThalIDs.txt", header=FALSE)
geneSelected4 <- read.table( "~/Documents/Research/300/Twisst/Twisst4c_topo1_0.99outliers_ThalIDs.txt", header=FALSE)

geneSelected1$V1=as.character(geneSelected1$V1)
geneSelected2$V1=as.character(geneSelected2$V1)
geneSelected3$V1=as.character(geneSelected3$V1)
geneSelected4$V1=as.character(geneSelected4$V1)

geneSelected5 = geneSelected1[geneSelected1$V1 %in% geneSelected3$V1,]
geneSelected6 = geneSelected2[geneSelected2$V1 %in% geneSelected4$V1,]

geneSelectedList1 <- as.vector(geneSelected1$V1)
geneSelectedList2 <- as.vector(geneSelected2$V1)
geneSelectedList3 <- as.vector(geneSelected3$V1)
geneSelectedList4 <- as.vector(geneSelected4$V1)
geneSelectedList5 <- as.vector(geneSelected5)
geneSelectedList6 <- as.vector(geneSelected6)



geneUniverseList <- as.vector(geneUniverse$V1)



hgCutoff <-  0.05
params <-  new("GOHyperGParams", geneIds=geneSelectedList1,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c_topo3_0.99outliers_GOreport.html")

params <-  new("GOHyperGParams", geneIds=geneSelectedList2,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c_topo1_0.99outliers_GOreport.html")

params <-  new("GOHyperGParams", geneIds=geneSelectedList3,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst4c_topo3_0.99outliers_GOreport.html")

params <-  new("GOHyperGParams", geneIds=geneSelectedList4,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst4c_topo1_0.99outliers_GOreport.html")

params <-  new("GOHyperGParams", geneIds=geneSelectedList5,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c4c_topo3_0.99outliers_GOreport.html")

params <-  new("GOHyperGParams", geneIds=geneSelectedList6,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = hgCutoff, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c4c_topo1_0.99outliers_GOreport.html")


## Bonferonni corrected ##

params <-  new("GOHyperGParams", geneIds=geneSelectedList1,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = 0.05/973, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c_topo3_0.99outliers_GOreport_bonf.txt")

params <-  new("GOHyperGParams", geneIds=geneSelectedList2,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = 0.05/1143, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c_topo1_0.99outliers_GOreport_bonf.txt")

params <-  new("GOHyperGParams", geneIds=geneSelectedList3,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = 0.05/775, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst4c_topo3_0.99outliers_GOreport_bonf.txt")

params <-  new("GOHyperGParams", geneIds=geneSelectedList4,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = 0.05/1019, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst4c_topo1_0.99outliers_GOreport_bonf.txt")

params <-  new("GOHyperGParams", geneIds=geneSelectedList5,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = 0.05/294, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c4c_topo3_0.99outliers_GOreport_bonf.txt")

params <-  new("GOHyperGParams", geneIds=geneSelectedList6,
 universeGeneIds=geneUniverseList ,annotation = "org.At.tair",
 ontology = "BP", pvalueCutoff = 0.05/849, conditional = FALSE,
 testDirection = "over")

conditional(params) <-  TRUE
hgCondOver <-  hyperGTest(params)
hgCondOver 
htmlReport(hgCondOver, file="~/Documents/Research/300/Twisst/Twisst3c4c_topo1_0.99outliers_GOreport_bonf.txt")

# # Get genes involved with specified GO
GOname <- "0050826"
res1 <- select(org.At.tair.db, keys= paste("GO:", GOname, sep = "") , columns=c("TAIR","GO"), 
               keytype="GO")
WGtairids <- as.vector(res1$TAIR)
LOI <- as.vector(geneSelectedList[geneSelectedList %in% WGtairids ])
LOI

write.table(LOI[1:length(LOI)], paste(dirname, "GO",GOname,"_genelist_", subsetname,"_",analysisname, ".txt", sep = ""), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 

# # # PLUG IN THERE THE NUMBER OF GO WITH p < 0.05
vectampon <- pvalues(hgCondOver)[1:95 ]
vecnames <- as.vector(names(vectampon))
write.table(paste(vecnames, sep = '\t',vectampon) , paste(dirname, "GOofInterest_pvals_", subsetname,"_",analysisname, ".txt", sep = ""), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 

