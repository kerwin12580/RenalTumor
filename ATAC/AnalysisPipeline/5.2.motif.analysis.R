#' @description: identify the tumor-specific TFs

library(Signac)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(101)
library(ggpubr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(matrixStats)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G

source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")
setwd("/data/active_data/lzl/RenalTumor-20200713/DataAnalysis-20210803/scATAC")

scATAC.data <- readRDS("scATAC.data.pro.rds")
motif.info <- data.frame(originName = names(scATAC.data@assays$Peaks@motifs@motif.names), TF = unlist(scATAC.data@assays$Peaks@motifs@motif.names))
rownames(motif.info) <- NULL
motif.info$originName <- gsub("_", "-", motif.info$originName)

#####基于chromVAR数据分析
cellType.motifs.chromVAR <- readRDS("5.Motif/cellType.motifs.chromVAR.human_pwms_v2.rds")

###################1.plot motifs deviation score heatmap and number ratioplot
top5 <- sapply(cellType.motifs.chromVAR, function(x){
    return(x$gene[1:5])
})
top5 <- unique(as.character(top5))
####Plot --- cellType specific motif heatmap
sig.motifs <- cellType.motifs.chromVAR %>% 
                                bind_rows %>%
                                filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
                                select(gene) %>%
                                distinct()
idx <- match(sig.motifs$gene, motif.info$TF)
sig.motifNames <- motif.info$originName[idx]

motifs.avgExp <- AverageExpression(scATAC.data, features = sig.motifNames, assays = "chromvar")
motifs.avgExp <- motifs.avgExp$chromvar
rownames(motifs.avgExp) <- sig.motifs$gene
zScore <- function(x){(x - mean(x)) /sd(x)}
motifs.avgExp.scale <- apply(motifs.avgExp, 1, zScore) %>% t() # row: TF; column: cell type
mark.idx <- match(top5, rownames(motifs.avgExp.scale))
pdf("5.Motif/Analysis/chromVAR.cellType.sig.motifs.heatmap.pdf")
Heatmap(motifs.avgExp.scale, name = "Deviation score", 
        width = unit(10, "cm"), height = unit(12, "cm"),
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        show_row_dend = F, show_column_dend = F, show_column_names = T, show_row_names = F, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 8), by_row = T)) + 
rowAnnotation(link = anno_mark(at = mark.idx, labels = top5, link_width = unit(2, "mm"), labels_gp = gpar(fontsize = 5), padding = unit(1, "mm")))
dev.off()

####Plot --- top5 cell-type specific TFs
cellType.sig.motifs <- lapply(names(cellType.motifs.chromVAR), function(x){
    sig <- cellType.motifs.chromVAR[[x]] %>% filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% mutate(celltype = x)
    return(sig)
})
names(cellType.sig.motifs) <- names(cellType.motifs.chromVAR)
sig.motifs.num <- cellType.sig.motifs %>% bind_rows()
sig.motifs.num <- as.data.frame(table(sig.motifs.num$celltype))
pdf("5.Motif/Analysis/chromVAR.cellType.sig.motifs.barplot.pdf")
ggbarplot(sig.motifs.num, x="Var1", y="Freq", fill = "Var1", color = "Var1",
          sort.by.groups=FALSE, sort.val = "desc", #不按组排序
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") + rotate_x_text(30)
dev.off()

###################2.Screen cell-specific enriched motifs
source(file = "/home/longzhilin/Analysis_Code/code/analysis.diff.survival.TCGA.R")
DESeq2.normalized_counts <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.normalized_counts.rds")
DESeq2.normalized_counts <- log2(DESeq2.normalized_counts+1)
DESeq2.result <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/DESeq2.result.rds")
clin.data <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/TCGA/KIRC/Result/clin.data.rds")

# extreact the avg_log2FC and FDR value
order.TFs <- cellType.motifs.chromVAR[[1]]$gene
motifs.FC <- sapply(cellType.motifs.chromVAR, function(x){
  index <- match(order.TFs, x$gene)
  return(round(x$avg_log2FC[index], 2))
})
rownames(motifs.FC) <- order.TFs

motifs.fdr <- sapply(cellType.motifs.chromVAR, function(x){
  index <- match(order.TFs, x$gene)
  return(x$p_val_adj[index])
})
rownames(motifs.fdr) <- order.TFs
a <- motifs.fdr
a[which(a==0)] <- 2
motifs.fdr[which(motifs.fdr==0)] <- min(a)*0.001 ###4.940656e-324, Multiply the minimum value by 0.01, instead of 0
motifs.fdr <- round(-log10(motifs.fdr), 2)

#3.deviation score
all.motifs.avgExp <- AverageExpression(scATAC.data, features = motif.info$originName, assays = "chromvar")
all.motifs.avgExp <- all.motifs.avgExp$chromvar
rownames(all.motifs.avgExp) <- motif.info$TF

################################## screening Strategy
#1.deviation score: sd > median(sd) + 4*mad
#2.Tumor cells: fdr> 0.0001 & avg_log2FC>=4
source(file = "/home/longzhilin/Analysis_Code/DataScience/MAD.R")
avg.sd <- rowSds(all.motifs.avgExp)
names(avg.sd) <- rownames(all.motifs.avgExp)
sd.mad <- DoubleMAD(avg.sd)
mad.threshold <- median(avg.sd) + 4*sd.mad[2] #Take the right

#Only in cell type: FDR < 0.0001 & log2FC>=4; log2FC<1 in other cell types
sig.label <- sapply(order.TFs, function(x){
     if(avg.sd[x] > mad.threshold){
        fdr <- motifs.fdr[x,]
        FC <- motifs.FC[x,]
        sig.fdr <- which(fdr > -log10(0.0001))
        sig.fc <- which(FC >= 4)  
        others <- which(FC >= 1)

        sig.idx <- intersect(sig.fdr, sig.fc)
        common.len <- length(sig.idx)

        if(common.len==1 & length(others)==1){
        #if(common.len==1){
            return("specific")
        }else{
            return("common")
        }
    }else{
        return("low variation")
    }
})
cellType.high.specific.TF <- names(which(sig.label=="specific")) # 49
heatmapEM.fdr <- motifs.fdr[cellType.high.specific.TF,]
heatmapEM.FC <- motifs.FC[cellType.high.specific.TF,]
cellType.sd <- avg.sd[cellType.high.specific.TF]
write.xlsx(list(FDR = heatmapEM.fdr, FC = heatmapEM.FC), file = "5.Motif/Analysis/cellType.specific.TFs.xlsx", sheetName = c("FDR", "Log2FC"), rowNames = T)

#### Tumor specific TFs
idx <- which(colnames(heatmapEM.FC) == "Tumor")
index <- which(heatmapEM.FC[, idx] >= 4)
tumor.specific.TFs <- data.frame(Name = rownames(heatmapEM.fdr)[index], FDR = heatmapEM.fdr[index, idx], avg_log2FC = heatmapEM.FC[index, idx], sd = cellType.sd[index])
heatmapEM.fdr.tumor <- heatmapEM.fdr[index,]
heatmapEM.FC.tumor <- heatmapEM.FC[index,]
saveRDS(tumor.specific.TFs, file = "5.Motif/Analysis/tumor.specific.TFs.rds")

# survival in TCGA-KIRC data
write.xlsx(tumor.specific.TFs, file = "5.Motif/Analysis/tumor.specific.TFs.xlsx", sheetName = c("Tumor specific TFs"), rowNames = T)
pdf("5.Motif/Analysis/tumor.specific.TFs.survival.pdf")
tumor.TFs.signature.res <- analysis.diff.survival.TCGA(interest.gene = tumor.specific.TFs$Name, diff.gene.pro = DESeq2.result, exp.data.process = DESeq2.normalized_counts, clin.data = clin.data, EnhancedVolcano.plot = F, Box.plot = F, main = "tumor.specific.TFs", meta.signature = T, single.signature = T)
dev.off()

# surcivial in ICB data
library(ggpubr)
source(file = "/home/longzhilin/Analysis_Code/SurvivalAnalysis/Cox.function.R")
source(file = "/home/longzhilin/Analysis_Code/code/RCC.ICB.analysis.R")
normalized_expression <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/normalized_expression.rds")
patient.info.RNA <- readRDS("/data/active_data/lzl/RenalTumor-20200713/Data/ICB.therapy/patient.info.RNA.rds")
#sex：F - 0; M - 1
patient.info.RNA$Sex <- gsub("^F|FEMALE|Female$", 0, patient.info.RNA$Sex)
patient.info.RNA$Sex <- gsub("^M|Male|MALE$", 1, patient.info.RNA$Sex)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("PRIMARY", 0, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
patient.info.RNA$Tumor_Sample_Primary_or_Metastasis <- gsub("METASTASIS", 1, patient.info.RNA$Tumor_Sample_Primary_or_Metastasis)
single <- as.list(tumor.specific.TFs$Name)
names(single) <- tumor.specific.TFs$Name
TF.list <- c(list(All = tumor.specific.TFs$Name), single)
pdf("5.Motif/Analysis/tumor.specific.TFs.ICB.pdf")
cox.res <- RCC.icb.analysis(signature.list = TF.list , expresionMatrix = normalized_expression, clincal.info = patient.info.RNA)
dev.off()


####footprint plot
DefaultAssay(scATAC.data) <- "Peaks"
scATAC.data <- Footprint(
  object = scATAC.data,
  motif.name = c("OTP", "ISL1", "VENTX", "HOXC5"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)
pdf("5.Motif/Analysis/Footprint.pdf")
PlotFootprint(scATAC.data, features = "OTP")
PlotFootprint(scATAC.data, features = "ISL1")
PlotFootprint(scATAC.data, features = "VENTX")
PlotFootprint(scATAC.data, features = "HOXC5")
dev.off()

library(ggplot2)
library(ggseqlogo)
PWMatrixToProbMatrix <- function(x){
        if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
        m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
        m <- t(t(m)/colSums(m))
        m
}
library(chromVARmotifs)
data("human_pwms_v2")
index <- grep("OTP|ISL1|VENTX|HOXC5", names(human_pwms_v2))
PPM.list <- lapply(index, function(y){
    PPM <- PWMatrixToProbMatrix(human_pwms_v2[[y]])
})
pdf("5.Motif/Analysis/tumor.TF.logo.pdf")
names(PPM.list) <- as.character(unlist(scATAC.data@assays$Peaks@motifs@motif.names[index]))
p <- ggseqlogo(PPM.list) + theme(plot.title = element_text(hjust = 0.5))
print(p)
dev.off()