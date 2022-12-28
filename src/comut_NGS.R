#https://www.bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
#https://bioconductor.statistik.tu-dortmund.de/packages/3.5/bioc/vignettes/maftools/inst/doc/maftools.html

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("maftools") 
library(maftools)
library(data.table)
library(RColorBrewer)
library(randomcoloR)

setwd("/home/l.leek/drup_palb_ribo/")

#Load mutation data
mut.df = fread("data/mut_ngs_final.txt", data.table = F)
mut.df = mut.df[!mut.df$Tumor_Sample_Barcode %in% c("01-0294","02-0067","07-0168"), ]
#Load metadata
clinical.df = fread("data/Clinical_13052022.txt", data.table = F)

#Load CN status
#cn.df = fread("data/Voorbeeld_CN.tsv", data.table = F)

#Generate MAF object
data.maf = read.maf(maf = mut.df, 
                    clinicalData = clinical.df,
                    #  cnTable = custom.cn.data,
                    verbose = T)

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

#top 50
freq.df = as.data.frame(sort(table(mut.df$Hugo_Symbol), decreasing = T))
genes = as.character(freq.df$Var1)
genes =  c("CDKN2A","CDK4","CDK6","CCND1","CCND2","CCND3", 
           genes)
genes = genes[!duplicated(genes)]
if(length(genes) >= 50){
  genes = genes[1:50]
}

#Make plot
png(file="results/myOncoplot_NGS.png",width=4600,height=4500, res=300)

oncoplot(maf = data.maf,   
         bgCol = "white", borderCol = "grey", colors = col,
         showTumorSampleBarcodes = F,
         sortByMutation = FALSE,
         sortByAnnotation = TRUE,
         genes = genes, 
         altered = T,
         clinicalFeatures = c('clinicalResponse', 'tumorType'),
         annotationColor = annotation_col,
         keepGeneOrder = T) 

while (!is.null(dev.list()))  dev.off()


