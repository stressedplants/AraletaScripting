library(ggplot2)
library(gridExtra)
library(pheatmap)
library(Seurat)
library(RColorBrewer)

##############Note:
#Please download GSM5591756_Seurat_dedoublet
# from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184511
# and place in data/scRNAseq/

###############Universal variables and helper functions##############
#Note that the time points are really wacky, if you look at data table S1 in Hickman et al., 2017
ZTs=c("00:30",
      "00:45",
      "01:10",
      "02:00",
      "02:50",
      "03:30",
      "04:50",
      "06:10",
      "07:30",
      "08:50",
      "10:00",
      "11:10",
      "13:30",
      "15:50",
      "21:20")
clustLab=c("epidermal_0",
           "mesophyll_1",
           "stress_mesophyll_2", 
           "mesophyll_3",
           "mesophyll_4",
           "bundlesheath_5",
           "phloem_6",
           "vascular_7",
           "stress_epidermal_8",
           "companion_9",
           "unknown_10",
           "unknown_11",
           "guard_12",
           "sieve_13",
           "hydathode_14",
           "myrosinase_15",
           "unknown_16"
)

clustLabPretty=c("epidermal",
                 "mesophyll group 1",
                 "stressed mesophyll", 
                 "mesophyll group 2",
                 "mesophyll group 3",
                 "bundle sheath",
                 "phloem",
                 "vascular",
                 "stressed epidermal",
                 "companion",
                 "unknown group 1",
                 "unknown group 2",
                 "guard",
                 "sieve",
                 "hydathode",
                 "myrosinase",
                 "unknown group 3"
)
names(clustLab)=paste("cl", c(0:16), sep="")
names(clustLabPretty)=paste("cl", c(0:16), sep="")

zscore <- function(table){
  apply(table, 1, function(i){
    if(max(i)==0){
      i
    }else{
      (i-mean(i))/sd(i)
    }
    
  })
}


#####################################
# Step 1:  Accuracy of CIBERSORTx on Procko et al. 2022 data set (single cell type)
file="outputs/CIBERSORToutputs/outputs_mixturePreds/simulatedLeafParts_eckerChory/CIBERSORTx_Job72_Results_leafPartsEC_EckerChory.txt"
mat=read.table(file, header=T, row.names = 1)
pheatmap(mat[,1:17], cluster_rows=F, cluster_cols=F)

#Ugh, re-order for making it pretty:
reorder_col=c(4:7, 2:3, 1, 10, 12, 8:9, 11, 15:16, 14, 13, 17)
reorder_row=c(1:4, 19:34, 5:18)
pheatmap(mat[reorder_row,reorder_col], cluster_rows=F, cluster_cols=F)

#Re-name rows and columns 
mat_renamed=mat[reorder_row,reorder_col]
colnames(mat_renamed)=clustLabPretty

rownames(mat_renamed)=unlist(lapply(clustLabPretty, function(i){ paste(i, c("rep 1", "rep 2"))}))


pheatmap(mat_renamed, cluster_rows=F, cluster_cols=F, filename="figures/TissueTimer_FigureS1_pureTissues.pdf", width=6, height=6)
dev.off()


########
# Step 2:  Accuracy of CIBERSORTx on Procko et al. 2022 data set (simulated mixed cell types)
file="outputs/CIBERSORToutputs/outputs_mixturePreds/simulatedLeaves_eckerChory/CIBERSORTx_Job69_Results_simLeavesEC_EckerChory.txt"
mat=read.table(file, header=T, row.names = 1)

#What are the real ratios?
file="outputs/CIBERSORToutputs/outputs_mixturePreds/simulatedLeaves_eckerChory/leafratios_real.txt"
realratios=read.table(file, header=F)
rownames(realratios)=paste("cl", c(0:16), sep="")

sapply(rownames(realratios), function(i){
  plot(as.numeric(realratios[i,]), as.numeric(mat[,i]), main=clustLabPretty[i], xlab="true tissue proportion", ylab="predicted tissue proportion")
  abline(c(0,1))
})

cors=sapply(rownames(realratios), function(i){
  cor(as.numeric(realratios[i,]), as.numeric(mat[,i]))
})

barplot(cors, xlab="cluster", ylab="Pearson's correlation", ylim=c(0,1))

slopes=sapply(rownames(realratios), function(i){
  a=as.numeric(realratios[i,])
  b=as.numeric(mat[,i])
  lm(b ~ a)$coefficients[2]
})
names(slopes)=names(cors)
barplot(slopes, xlab="cluster", ylab="Slope of line of best fit")

#Now make it pretty
pdf("figures/TissueTimer_FigureS2_mixTissuesScatter.pdf", height=15, width=12)
par(mfcol = c(5,4), oma = c(0.1, 0.1, 0.1, 0.1))
sapply(rownames(realratios), function(i){
  plot(as.numeric(realratios[i,]), as.numeric(mat[,i]), main=clustLabPretty[i], xlab="true tissue proportion", ylab="predicted tissue proportion", pch=20)
  abline(c(0,1))
})
dev.off()

names(cors)=clustLabPretty
names(slopes)=clustLabPretty
pdf("figures/TissueTimer_Figure1B_mixTissuesSummary.pdf", height=10, width=6)
par(mfcol = c(2,1), oma = c(1, 0.5, 0.5, 0.5), omi=c(1, 0.5, 0.5, 0.5))

barplot(cors[order(cors)], las=3, ylab="Pearson's correlation", ylim=c(0,1))
abline(h=0.8, lty=2)
barplot(slopes[order(cors)], las=3, ylab="Slope of line of best fit")
abline(h=1, lty=2)


dev.off()


###################
# Step 3:  Accuracy of CIBERSORTx on Ecker/Chory signatures on microdissection data

### Micro-dissection data:
file="outputs/CIBERSORToutputs/outputs_mixturePreds/simulatedLeafParts_dissection/CIBERSORTx_Job99_eckerChory_allLeaves.txt"
mat=read.table(file, header=T, row.names = 1)
pheatmap(mat[,1:17], cluster_rows=F, cluster_cols=F)

types=substring(row.names(mat), 1, 5)

lapply(unique(types), function(i){
  pheatmap(mat[which(types==i),1:17], cluster_rows=F, cluster_cols=T)
})


colnames(mat)=clustLabPretty[colnames(mat)]
titles=c("bulkR", "Guard", "Lower Epidermal", "Palisade Mesophyll", 'Spongy Mesophyll', "Upper Epidermal", "Vascular")
names(titles)=unique(types)
hm=lapply(unique(types)[-1][c(2, 5, 3,4,6,1)], function(i){
  pheatmap(mat[which(types==i),1:17], cluster_rows=T, cluster_cols=T, show_rownames = F, main=titles[i])[[4]]
})
pdf("figures/TissueTimer_Figure1C_S3_labelMicrodissectedSamples.pdf", height=15, width=10)

grid.arrange(grobs=hm, col=3)

dev.off()

########
# Step 4: Gene expression in a developmental time course

devExp=read.table("outputs/CIBERSORToutputs/outputs_mixturePreds/development/CIBERSORTx_Job62_Results_devEckerChory.txt", header=T, row.names=1)

colnames(devExp)[1:17]=clustLabPretty[colnames(devExp)[1:17]]
pheatmap(devExp[,1:17], cluster_rows=F, cluster_cols=T)

pheatmap(zscore(t(devExp[,1:17])), cluster_rows=F, cluster_cols=T)

annot=data.frame(rep=c(rep("1", 14), rep("2", 14)), 
                 stage=c(rep("Growth", 5), rep("Mature", 4), rep("Senescent", 5),
                         rep("Growth", 5), rep("Mature", 4), rep("Senescent", 5)))
rownames(annot)=rownames(devExp)
pheatmap(zscore(t(devExp[,1:17])), cluster_rows=F, cluster_cols=T, annotation_row = annot, gaps_row=c(5, 9, 14, 19, 23), labels_row=c(seq(4, 30, 2), seq(4, 30, 2)), file="figures/TissueTimer_Figure2A_developmental.pdf", width=6, height=6)
dev.off()

pheatmap(devExp[,1:17], cluster_rows=F, cluster_cols=T, annotation_row = annot, gaps_row=c(5, 9, 14, 19, 23), labels_row=c(seq(4, 30, 2), seq(4, 30, 2)), file="figures/TissueTimer_FigureS5A_developmental_unscaled.pdf", width=6, height=6)
dev.off()

########
# Step 5: Tissue ratios over a diurnal timecourse and with MeJA treatment

file="outputs/CIBERSORToutputs/outputs_mixturePreds/diurnal/CIBERSORTx_Job63_Results_JAEckerChory.txt"
mat=read.table(file, header=T, row.names = 1)

colnames(mat)[1:17]=clustLabPretty[colnames(mat)[1:17]]
pheatmap(mat[,1:17], cluster_rows=F, cluster_cols=T)

pheatmap(zscore(t(mat[,1:17])), cluster_rows=F, cluster_cols=T)




#average across reps

###could try having reps as columns
groups=sapply(strsplit(rownames(mat), "_"), function(i){
  i[length(i)]
})
splitByRep=lapply(unique(groups), function(i){
  mat[which(groups==i),]
})
mergeReps=splitByRep[[1]]+splitByRep[[2]]+splitByRep[[3]]+splitByRep[[4]]  #do.call(cbind, splitByRep)

rownames(mergeReps)[1:15]=ZTs

pheatmap(zscore(t(mergeReps[1:15,1:17])), cluster_rows=F, cluster_cols=T)

###Make pretty

annot=data.frame(light=c(rep("day", 11), rep("night", 4)), row.names = rownames(mergeReps)[1:15])

annot_col=list(light=c("lightyellow", "black"))
names(annot_col[[1]])=c("day", "night")

pheatmap(zscore(t(mergeReps[1:15,1:17])), cluster_rows=F, cluster_cols=T, annotation_row = annot, annotation_colors = annot_col, gaps_row=c(6, 11))

#for supplements
pheatmap(mergeReps[1:15,1:17]/4, cluster_rows=F, cluster_cols=T, annotation_row = annot, annotation_colors = annot_col, gaps_row=c(6, 11), file="figures/TissueTimer_FigureS5B_diurnal_unscaled.pdf", width=6, height=6)
dev.off()

pheatmap(zscore(t(mergeReps[1:15,1:17])), cluster_rows=F, cluster_cols=T, annotation_row = annot, annotation_colors = annot_col, gaps_row=c(6, 11), file="figures/TissueTimer_Figure2B_diurnal.pdf", width=6, height=6)
dev.off()


#compare JA and Mock
pheatmap(-log(1+mergeReps[2:15,1:17])+log(1+mergeReps[16:29,1:17]), cluster_rows=F, annotation_row = annot, annotation_colors = annot_col, gaps_row=c(5, 10), file="figures/TissueTimer_FigureS6_MeJA.pdf", width=6, height=6)


########
# Step 6: Analysis of light sensitivity of each cell type

lightGeneList=read.table("data/otherData/lightInductionGeneList.txt", header=T, row.names=1, stringsAsFactors = T)
###Downloadable from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184511
a=readRDS("data/scRNAseq/GSM5591756_Seurat_dedoublet")
###use cluster designations from paper
asub=a 

b=UMAPPlot(asub)
clusters=b$data$ident

cells=which(clusters %in% c(0:4))
genes=rownames(lightGeneList)[which(rownames(lightGeneList) %in% rownames(asub))]

a=sapply(paste(c(0:16)), function(i){
  temp=rowSums(asub[["SCT"]]@data[genes,which(paste(clusters)==i)])/length(which(paste(clusters)==i))
  
})
a=a[which(rowSums(a)>1),]
pheatmap(zscore(a))
lightGeneList$induced=as.factor(lightGeneList$induced)
anncol=list(time = c("red", "lightyellow", "grey"), induced=c("cornflowerblue", "orange"))
names(anncol$time)=c("night", "day", "none")
names(anncol$induced)=c("FALSE", "TRUE")
colnames(a)=clustLabPretty
singleCellExpression=a
pheatmap(zscore(a), annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
ab=zscore(a)

breaksList=seq(-3, 4, 0.1)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="day" & lightGeneList[colnames(ab),1]=="TRUE")], filename = "figures/TissueTimer_Figure3_DT.pdf", width=6, height=4, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="night" & lightGeneList[colnames(ab),1]=="TRUE")],filename = "figures/TissueTimer_Figure3_NT.pdf", width=6, height=4, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="day" & lightGeneList[colnames(ab),1]=="FALSE")],filename = "figures/TissueTimer_Figure3_DF.pdf", width=6, height=4,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="night" & lightGeneList[colnames(ab),1]=="FALSE")],filename = "figures/TissueTimer_Figure3_NF.pdf", width=6, height=4, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)

pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="day" & lightGeneList[colnames(ab),1]=="TRUE")], filename = "figures/TissueTimer_Figure3_DT.pdf", width=6, height=4, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="night" & lightGeneList[colnames(ab),1]=="TRUE")],filename = "figures/TissueTimer_Figure3_NT.pdf", width=6, height=4, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="day" & lightGeneList[colnames(ab),1]=="FALSE")],filename = "figures/TissueTimer_Figure3_DF.pdf", width=6, height=4,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)
pheatmap(zscore(a)[,which(lightGeneList[colnames(ab),2]=="night" & lightGeneList[colnames(ab),1]=="FALSE")],filename = "figures/TissueTimer_Figure3_NF.pdf", width=6, height=4, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)#, annotation_col = data.frame(lightGeneList), annotation_colors=anncol)

########
# Step 7: cell type predictions for plants during bolting

bolted=read.table("outputs/CIBERSORToutputs/outputs_mixturePreds/bolted/CIBERSORTx_Job101_Results_bolted.txt", header=T, row.names=1)
unbolted=read.table("outputs/CIBERSORToutputs/outputs_mixturePreds/bolted/CIBERSORTx_Job102_Results_notbolted.txt", header=T, row.names=1)
boltingstat=rbind(bolted, unbolted)

boltPhenotypes=read.table("data/otherData/phenotype_data_flowering_variability.csv", header=T, row.names = 1, sep=",")

#pheatmap ordered by 

annot=data.frame(bolted=c(rep("yes", length(rownames(bolted))), rep("no", length(rownames(unbolted)))), row.names = rownames(boltingstat))

colnames(boltingstat)[1:17]=clustLabPretty[colnames(boltingstat)[1:17]]
pheatmap(boltingstat[,1:17], cluster_rows=T, cluster_cols=T, annotation_row = annot) ####note: try ordering by pseudotime
tempPseudoTime=as.character(read.table("data/bulkRNAseq/bolting/TPM_pseudotime.csv", sep=",", stringsAsFactors = F)[1,])[-1]
orderLeafSize=rownames(boltPhenotypes)[order(boltPhenotypes[,"leaf_avg"])]
orderLeafSize=orderLeafSize[which(orderLeafSize %in% tempPseudoTime)]
orderBiomass=rownames(boltPhenotypes)[order(boltPhenotypes[,"biomass"])]
orderBiomass=orderBiomass[which(orderBiomass %in% tempPseudoTime)]


pheatmap(boltingstat[tempPseudoTime,1:17], cluster_rows=F, cluster_cols=T, annotation_row = annot)
pheatmap(boltingstat[tempPseudoTime,1:17], cluster_rows=F, cluster_cols=T, annotation_row = annot, file="figures/TissueTimer_FigureS7_boltingheatmapunscaled.pdf", width=4, height=10)
dev.off()
pheatmap(boltingstat[orderLeafSize,1:17], cluster_rows=F, cluster_cols=T, annotation_row = annot)
pheatmap(boltingstat[orderBiomass,1:17], cluster_rows=F, cluster_cols=T, annotation_row = annot)


#remove empty columns
boltingFiltered=boltingstat[tempPseudoTime,1:17][,which(colSums(boltingstat[tempPseudoTime,1:17])>1.5)]

pheatmap(boltingFiltered, cluster_rows=F, cluster_cols=T, annotation_row = annot, scale="column") #, file="TissueTimer_Figure3_changeInPseudotime.pdf", width=6, height=6) ####noute: try ordering by pseudotime
pheatmap(boltingFiltered, cluster_rows=F, cluster_cols=T, annotation_row = annot, scale="column", file="figures/TissueTimer_Figure4A_boltingheatmapscaled.pdf", width=4, height=10) #, file="TissueTimer_Figure3_changeInPseudotime.pdf", width=6, height=6) ####noute: try ordering by pseudotime
dev.off()
pheatmap(boltingFiltered[orderLeafSize,], cluster_rows=F, cluster_cols=T, annotation_row = annot, scale="column")
pheatmap(boltingFiltered[orderBiomass,], cluster_rows=F, cluster_cols=T, annotation_row = annot, scale="column")

#Let's look at Spearman correlations, as well as a few line plots
spearmanPseudoTime=apply(boltingFiltered, 2, function(i){
  cor(c(1:length(i)), i, method="spearman")
})
spearmanLeafSize=apply(boltingFiltered[orderLeafSize,], 2, function(i){
  cor(c(1:length(i)), i, method="spearman")
})
spearmanBiomass=apply(boltingFiltered[orderBiomass,], 2, function(i){
  cor(c(1:length(i)), i, method="spearman")
})
#tidy
spearmanPseudoTime=data.frame("Spearman"=spearmanPseudoTime, "Cell"=names(spearmanPseudoTime), "Type"=rep("pseudotime", length(spearmanPseudoTime)))
spearmanLeafSize=data.frame("Spearman"=spearmanLeafSize, "Cell"=names(spearmanLeafSize), "Type"=rep("leaf size", length(spearmanLeafSize)))
spearmanBiomass=data.frame("Spearman"=spearmanBiomass, "Cell"=names(spearmanBiomass), "Type"=rep("biomass", length(spearmanBiomass)))

spearmanTidy=rbind(spearmanPseudoTime, spearmanLeafSize, spearmanBiomass)

spearmanTidy$Cell=factor(spearmanTidy$Cell, levels=c("vascular", "bundle sheath", "sieve", "companion", "mesophyll group 1", "mesophyll group 2", "mesophyll group 3", "epidermal", "myrosinase"))
ggplot(spearmanTidy, aes(fill=Type, y=Spearman, x=Cell)) + 
  geom_bar(position="dodge", stat="identity") + theme_minimal() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c('cornflowerblue', 'lightgreen', 'orange')) +xlab("Cell Type") + ylab("Spearman's Correlation") +ggsave("figures/TissueTimer_Figure4C_spearmanCorr.pdf", width=5)

#show some examples as line plots:
sapply(c(1:9), function(i){
  plot(c(), xlim=c(min(boltingFiltered[,i]),max(boltingFiltered[,i])), ylim=c(0,1), xlab="Plants in order", ylab="Proportion of cell type", main=colnames(boltingFiltered)[i])
  points(boltingFiltered[orderBiomass,i], boltPhenotypes[orderBiomass,"biomass"]/max(boltPhenotypes[orderBiomass,"biomass"]), col="cornflowerblue", pch=20)
  points(boltingFiltered[orderBiomass,i], boltPhenotypes[orderBiomass,"leaf_avg"]/max(boltPhenotypes[orderBiomass,"leaf_avg"]), col="lightgreen", pch=20)
  points(boltingFiltered[,i], c(1:65)/65, col="orange", pch=20)
  
})


##########
# Step 8:  Which transcription factors are in which cell types?

#get DAP-seq network
setwd("data/otherData/dap_data_v4/")
in.list <- list.files(recursive = T, full.names = T)
makeEntireDAPseqNetwork=lapply(in.list, function(i){
  read.table(i, header=T)
})
setwd("../../..")
makeEntireDAPseqNetwork=do.call(rbind, makeEntireDAPseqNetwork)

dapseqGenes=unique(c(makeEntireDAPseqNetwork[,1], makeEntireDAPseqNetwork[,2]))

##############make less filtered singleCellExpressionTable

a=readRDS("data/scRNAseq/GSM5591756_Seurat_dedoublet")

###use cluster designations from paper
asub=a #[genes,]
b=UMAPPlot(asub)
clusters=b$data$ident
a=sapply(paste(c(0:16)), function(i){
  temp=rowSums(asub[["SCT"]]@data[,which(paste(clusters)==i)])/length(which(paste(clusters)==i))
})
colnames(a)=clustLabPretty

#for the TFs, what are their cellular localisations?
tfsInDAPseq=unique(makeEntireDAPseqNetwork[,1])
tfsInDAPseqAndSC=tfsInDAPseq[which(tfsInDAPseq %in% rownames(a))]

aTFsDAPseq=a[tfsInDAPseqAndSC,]
pheatmap(aTFsDAPseq[which(rowSums(aTFsDAPseq)>1),])

#tissue specific gene network, as determined by cellular concentrations
thresh=0.25
subNetSameCell=apply(makeEntireDAPseqNetwork, 1, function(i){
  if(i[1] %in% rownames(a) & i[2] %in% rownames(a)){
    a[i[1], ]>thresh & a[i[2],]>thresh
  }else{rep(0,length(a[1,]))}
})
ids=which(colSums(subNetSameCell)>0)
subNetSameCell=cbind(t(subNetSameCell)[ids,], makeEntireDAPseqNetwork[ids,])
colSums(subNetSameCell[,1:17])
plot(table(table(subNetSameCell[,19])))

sapply(c(1:17), function(i){
  write.table(subNetSameCell[which(subNetSameCell[,i]==1),c(18,19)], file=paste("outputs/filteredNetworks/CellTypeNetworkThresh_0.25_cellType_", clustLabPretty[i], ".tsv", sep=""), quote=F, row.names = F)
})

#save it as list:
cellSpecificNetworks=lapply(c(1:17), function(i){
  subNetSameCell[which(subNetSameCell[,i]==1),c(18,19)]
})

####################################################
# Step 9: CIBERSORTx can impute realistic gene expression values

#generate expected gene expression by cell type
#also compare to imputed values to see if that's any closer?
imputedBolted=read.table("outputs/CIBERSORToutputs/outputs_imputed/CIBERSORTx_Job112_boltingImputed/CIBERSORTxGEP_Job112_GEPs.txt", header=T, row.names = 1)
imputedUnbolted=read.table("outputs/CIBERSORToutputs/outputs_imputed/CIBERSORTx_Job113_unboltedImputed/CIBERSORTxGEP_Job113_GEPs.txt", header = T, row.names = 1)
imputedBoth=read.table("outputs/CIBERSORToutputs/outputs_imputed/CIBERSORTx_Job114_allboltedImputed/CIBERSORTxGEP_Job114_GEPs.txt", header = T, row.names = 1)

colnames(imputedBolted)=clustLabPretty[colnames(imputedBolted)]
colnames(imputedUnbolted)=clustLabPretty[colnames(imputedUnbolted)]
colnames(imputedBoth)=clustLabPretty[colnames(imputedBoth)]
genesInBoth=rownames(imputedBolted)[which(rownames(imputedBolted) %in% rownames(a))]

pdf("figures/TissueTimer_FigureS9_imputeVSscRNA-seq.pdf", height=15, width=12)
par(mfcol = c(3,3), oma = c(0.1, 0.1, 0.1, 0.1))
sapply(colnames(imputedBolted)[2:10], function(i){
  plot(log(a[genesInBoth,i]), log(imputedBoth[genesInBoth,i]), main=i, xlab="log(expression from scRNA-seq)", ylab="log(imputed expression)")
})
dev.off()


#compare predicted values for imputed bolted and imputed unbolted
sapply(colnames(imputedBolted)[2:10], function(i){
  plot(log(imputedUnbolted[genesInBoth,i]), log(imputedBolted[genesInBoth,i]), col=rgb(0.1, 0.1, 0.1, 0.3), main=i, xlab="log(predicted unbolted expression)", ylab="log(predicted bolted expression)")
})

################################
# Step 10: For genes that are differentially expressed in bolting plants: what tissues are they expressed in?

boltingExpression=read.table("data/bulkRNAseq/bolting/TPM_pseudotime.csv", sep=",", stringsAsFactors = F, header=T, row.names = 1)

imputedHighResUp=lapply(paste("cl", c(0,1,3,4,5,7,9,13,15), sep=""), function(i){
  t=read.table(paste("outputs/CIBERSORToutputs/outputs_imputed/CIBERSORTx_Job116_highRes_UP_allBolted/CIBERSORTxHiRes_Job116_",i,"_Window35.txt", sep=""), header=T, row.names = 1, sep="\t")
  t[which(!apply(t, 1, function(j){is.na(j[1]) | sd(j)==0})),]
})

names(imputedHighResUp)=clustLabPretty[paste("cl", c(0,1,3,4,5,7,9,13,15), sep="")]

uniqueDiffExp=unique(unlist(sapply(imputedHighResUp, function(i){rownames(i)})))
write.table(uniqueDiffExp, "outputs/geneLists/uniqueDiffExpCellTypeSpecificBolting.txt", quote=F, row.names = F)

#what is the expression profile of these genes in scRNA-seq experiment?
pheatmap(a[uniqueDiffExp[which(uniqueDiffExp %in% rownames(a))],names(imputedHighResUp)], scale='row')

#what is the expression profile of these genes in scRNA-seq experiment?
heatmapBoltingExp=pheatmap(boltingExpression[uniqueDiffExp[which(uniqueDiffExp %in% rownames(boltingExpression))],tempPseudoTime], scale='row', cluster_cols = F)
pheatmap(boltingExpression[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                           tempPseudoTime], scale='row',cluster_cols=F, cluster_rows=F)
pheatmap(t(zscore(boltingExpression[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                                    tempPseudoTime])), cluster_cols=F, cluster_rows=F)
pheatmap(a[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
           names(imputedHighResUp)], scale='row',cluster_cols=F, cluster_rows=F)

#now how about the predicted cell types
predictedCellTypes=sapply(heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))], function(j){
  sapply(imputedHighResUp, function(i){
    if(j %in% rownames(i)){
      if(sum(i[j, c(1:3)])>sum(i[j, c(67:69)])){-3}else{3}
      
    }else{0}
  })
})

pheatmap(t(predictedCellTypes), cluster_rows = F)


#combine these figures:
sub1=t(zscore(boltingExpression[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                                tempPseudoTime]))
sub2=t(predictedCellTypes)
sub3=t(zscore(a[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                names(imputedHighResUp)]))

breaks=c(max(which(sub2[,1]!=0)),
         max(which(sub2[,2]!=0)),
         max(which(sub2[,3]!=0)),
         max(which(sub2[,5]!=0)),
         max(which(sub2[,6]!=0)))

pheatmap(sub1, cluster_rows=F, cluster_cols=F, show_rownames =F, show_colnames=F, gaps_row=breaks, filename = "figures/TissueTimer_Figure5A_imputedUp_expressionInBoltRNAseq.pdf", width=5)

pheatmap(abs(sub2), cluster_rows=F, cluster_cols=F, show_rownames =F, gaps_row=breaks, filename = "figures/TissueTimer_Figure5BimputedUp_predictedTissue.pdf", width=2.5)
write.table(sub2, file="outputs/geneLists/TissueTimer_TableSX_predictedCellTypeUp.csv", sep=",")
pheatmap(sub3, cluster_rows=F, cluster_cols=F, show_rownames =F, gaps_row=breaks, filename = "figures/TissueTimer_Figure5CimputedUp_scRNAseq.pdf", width=2.5)
dev.off()
#####Let's repeat this excercise with genes that decrease their expression after bolting

imputedHighResDown=lapply(paste("cl", c(0,1,3,4,5,7,9,13,15), sep=""), function(i){
  t=read.table(paste("outputs/CIBERSORToutputs/outputs_imputed/CIBERSORTx_Job117_output_hiresDOWN_bolting/CIBERSORTxHiRes_Job117_",i,"_Window35.txt", sep=""), header=T, row.names = 1, sep="\t")
  t[which(!apply(t, 1, function(j){is.na(j[1]) | sd(j)==0})),]
})

names(imputedHighResDown)=clustLabPretty[paste("cl", c(0,1,3,4,5,7,9,13,15), sep="")]

uniqueDiffExp=unique(unlist(sapply(imputedHighResDown, function(i){rownames(i)})))

#what is the expression profile of these genes in scRNA-seq experiment?
pheatmap(a[uniqueDiffExp[which(uniqueDiffExp %in% rownames(a))],names(imputedHighResDown)], scale='row')

#what is the expression profile of these genes in scRNA-seq experiment?
heatmapBoltingExp=pheatmap(boltingExpression[uniqueDiffExp[which(uniqueDiffExp %in% rownames(boltingExpression))],tempPseudoTime], scale='row', cluster_cols = F)
pheatmap(boltingExpression[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                           tempPseudoTime], scale='row',cluster_cols=F, cluster_rows=F)
pheatmap(t(zscore(boltingExpression[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                                    tempPseudoTime])), cluster_cols=F, cluster_rows=F)
pheatmap(a[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
           names(imputedHighResDown)], scale='row',cluster_cols=F, cluster_rows=F)

#now how about the predicted cell types
predictedCellTypes=sapply(heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))], function(j){
  sapply(imputedHighResDown, function(i){
    if(j %in% rownames(i)){
      if(sum(i[j, c(1:3)])>sum(i[j, c(67:69)])){-3}else{3}
      
    }else{0}
  })
})

pheatmap(t(predictedCellTypes), cluster_rows = F)


#combine these figures:
sub1=t(zscore(boltingExpression[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                                tempPseudoTime]))
sub2=t(predictedCellTypes)
sub3=t(zscore(a[heatmapBoltingExp$tree_row$labels[which(heatmapBoltingExp$tree_row$labels %in% rownames(a))],
                names(imputedHighResDown)]))

breaks=c(max(which(sub2[,1]!=0)),
         max(which(sub2[,2]!=0)),
         max(which(sub2[,3]!=0)),
         max(which(sub2[,4]!=0)))

pheatmap(sub1, cluster_rows=F, cluster_cols=F, show_rownames =F, show_colnames=F, gaps_row=breaks, filename = "figures/TissueTimer_FigureS8A_imputedDown_expressionInBoltRNAseq.pdf", width=5)

pheatmap(abs(sub2), cluster_rows=F, cluster_cols=F, show_rownames =F, gaps_row=breaks, filename = "figures/TissueTimer_FigureS8BimputedDown_predictedTissue.pdf", width=2.5)
write.table(sub2, file="outputs/geneLists/TissueTimer_Table_predictedCellTypeDown.csv", sep=",")
pheatmap(sub3, cluster_rows=F, cluster_cols=F, show_rownames =F, gaps_row=breaks, filename = "figures/TissueTimer_FigureS8CimputedDown_scRNAseq.pdf", width=2.5)


