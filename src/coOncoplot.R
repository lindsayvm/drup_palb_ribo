#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("maftools") 
library(maftools)
library(data.table)
library(RColorBrewer)
library(randomcoloR)

setwd("/home/l.leek/drup_palb_ribo/")

#Load mutation data
mut.df = fread("data/mut_wgs_final_06062022.txt", data.table = F)
#Load metadata
clinical.df = fread("data/Clinical_13052022.txt", data.table = F)
clinical.df$tumorType[clinical.df$tumorType %in% c("Neuroendorine carcinoma","pNET","Atypical carcinoid")] = "Neuroendocrine carcinoma"
clinical.df$tumorType[clinical.df$tumorType %in% c("Bladder cancer")] = "Urothelial cell carcinoma"
clinical.df$tumorType[clinical.df$tumorType %in% c("RCC")] = "Renal cell carcinoma"
clinical.df$tumorType = gsub("cancer","carcinoma",clinical.df$tumorType)

#Generate MAF object
data.maf = read.maf(maf = mut.df, 
                    clinicalData = clinical.df,
                    #  cnTable = custom.cn.data,
                    verbose = T)

#Load data 
mut_ngs.df = fread("data/mut_ngs_final.txt",data.table = F)

#Generate MAF object
ngs.maf = read.maf(maf = mut_ngs.df, clinicalData = clinical.df)


#Select color for variants
x = 'Set2' #https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
varClass = unique(data.maf@data$Variant_Classification)
col = RColorBrewer::brewer.pal(n = length(varClass), name = x)
names(col) = varClass

#color for clinical features
unique(data.maf@clinical.data$clinicalResponse)
clinicalResponse.v = c("SD", "PD", "NE", "Non CR/PD")
clinicalResponse = c("#7CFC00","#FFA500", "#FF0000", "#808080")
names(clinicalResponse) = clinicalResponse.v 

#color for tumors
tumorType.v = unique(data.maf@clinical.data$tumorType)
tumorType = distinctColorPalette(k = length(tumorType.v), altCol = FALSE, runTsne = FALSE)
names(tumorType) = tumorType.v

#Final colored list. Names of the list elements should match those in clinicalFeatures arguments 
annotation_col = list(clinicalResponse = clinicalResponse, 
                      tumorType = tumorType)





#NGS + WGS ~ Comparing two cohorts (MAFs). 
wgsNgs = mafCompare(m1 = data.maf, m2 = ngs.maf, 
                    m1Name = 'WGS', m2Name = 'NGS', minMut = 1)

# top 25 in both
top25_WGS = wgsNgs$results$Hugo_Symbol[order(wgsNgs$results$WGS, decreasing = T)][1:25]
top25_NGS = wgsNgs$results$Hugo_Symbol[order(wgsNgs$results$NGS, decreasing = T)][1:25]
topgenes = unique(c(top25_WGS, top25_NGS))

top25_WGSNGS = wgsNgs$results[wgsNgs$results$Hugo_Symbol %in% topgenes, ]
top25_WGSNGS = top25_WGSNGS[order(top25_WGSNGS$WGS, decreasing = T), c(1:3)]
#write.table(x = top25_WGSNGS, file = "Desktop/top25_WGSNGS.tsv", sep = "\t", row.names = F, col.names = T)
#25 meest voorkomende WGS en 25 meest vorkomende NGS

topgenes =  c("CDKN2A","CDK4","CDK6","CCND1","CCND2","CCND3", 
              topgenes)
topgenes = topgenes[!duplicated(topgenes)]
if(length(topgenes) >= 25){
  topgenes = topgenes[1:25]
}


png(file="results/wgs_ngs_coOncoplot.png",width=4600,height=4500, res=300)

coOncoplot(m1 = data.maf, m2 = ngs.maf, 
           m1Name = 'WGS', m2Name = 'NGS', 
           genes = topgenes, 
           removeNonMutated = TRUE,
           bgCol = "white", borderCol = "grey", colors = col,
           clinicalFeatures1 = c('clinicalResponse', "tumorType"),
           clinicalFeatures2 = c('clinicalResponse', "tumorType"),
           sortByAnnotation1 = TRUE, sortByAnnotation2 = TRUE,
           annotationColor1 = annotation_col,
           annotationColor2 = annotation_col,
           showSampleNames = F)
dev.off()

