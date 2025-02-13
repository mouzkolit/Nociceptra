
#load all the packages
library(DEGreport)
library(dplyr)
library(tibble)
library(xml2)
library(clusterProfiler)
library(DOSE)
library(data.table)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(org.Mm.eg.db)
library("biomaRt")
library("DESeq2")
library("RColorBrewer")
library("gplots")
dev.off()
directory = "../files/raw_counts_mrna" # please set here the directory where files are located/ or mirnas
setwd(directory) # setting the directory


outputPrefix <- "mRNA_seq_analysis" #change for renaming
sampleFiles<- c("19101-0001.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0002.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0003.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0004.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0005.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0006.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0007.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0008.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0009.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0010.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0011.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0012.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0013.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0014.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0015.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0016.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0017.fastqAligned.sortedByCoord.out.bam.counts",
                "19101-0018.fastqAligned.sortedByCoord.out.bam.counts"
) # please add here the filenames --> either all samples or cell line specific samples

sampleNames <- c("AD3_1.1",
                 "AD3_1.2",
                 "AD3_1.3",
                 "AD3_5.1",
                 "AD3_5.2",
                 "AD3_5.3",
                 "AD3_9.1",
                 "AD3_9.2",
                 "AD3_9.3",
                 "AD3_16.1",
                 "AD3_16.2",
                 "AD3_16.3",
                 "AD3_26.1",
                 "AD3_26.2",
                 "AD3_26.3",
                 "AD3_36.1",
                 "AD3_36.2",
                 "AD3_36.3"
                 ) # please add here the name of the samples

#sampleCondition = rep(c("AD2", "840", "AD3"), each = 18)
time = rep(c("DAY0", "DAY5", "DAY9", "DAY16", "DAY26","DAY36"), each = 3) # add here your time covariates either for one cell line as represented here or for change for all cell lines
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = time)

treatments = c("DAY0","DAY5","DAY9","DAY16","DAY26","DAY36") # same


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)


# run the deseq model over the time
dds <- DESeq(ddsHTSeq, test = "LRT", reduced=~1)
res <-results(dds)
res <- lfcShrink(dds=dds, coef=2, res=res) # for differential expression D36 AD3/AD2 we used ashr shrinkage

mat_rld  = assay(dds)
des = colData(ddsHTSeq)
# Subset the LRT results to return genes with padj < 0.05
ressig = subset(res, padj < 0.05) # for miRNAs set to 0.01
# save data results and normalized reads to csv
# insert gene_name_chromosome to the table 
#make the table for all and for signficant genes
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata1 <- merge(as.data.frame(ressig), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata1)[1] <- 'ensembl_gene_id_version' # set the name


# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "_normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "-test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plotDispEsts(dds)
dev.copy(svg, paste0(outputPrefix, "-outlier_detection_cooks.svg"))
dev.off()


#gene dispersion plot 
plotDispEsts(dds)
dev.copy(svg, paste0(outputPrefix, "-Gene_dispersion_plots.svg"))
dev.off()

# check and replace outliers
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
  
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "-replaceoutliers.csv"))
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]


#write the output
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "-replaceoutliers-results.csv"))

# MA plot of RNAseq data for entire dataset
plotMA(res, ylim=c(-8,8),main = "mRNAseq experiment", alpha = 0.05)
dev.copy(svg, paste0(outputPrefix, "-MAplot_initial_analysis.svg"))
dev.off()

# transform raw counts into normalized values with variance stablization or regularized log
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# put them into the dataframe and save them for concatenation of later used genes  
rld_counts = as.data.frame(assay(rld))
vst_counts = as.data.frame(assay(vsd))
write.csv(rld_counts, paste0(outputPrefix, "-rlog-transformed-counts.csv"))
write.csv(vst_counts,file = paste0(outputPrefix, "-vst-transformed-counts.csv"))


# save normalized values as input for the clustering ect.
write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "-rlog-transformed-counts.csv")))
write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "-vst-transformed-counts.csv")))

# plot to show effect of transformation
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(svg,paste0(outputPrefix, "-variance_stabilizing.svg"))
dev.off()

# clustering analysis
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, sampleNames, sep=" : "))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
dev.copy(svg, paste0(outputPrefix, "-clustering.svg"))
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(10,10))
dev.off()

# save the PCA
dev.copy(svg, paste0(outputPrefix, "-pca.svg"))
degPCA(log2(counts(dds)+0.5), colData(dds), # here one could also use vst counts or rlog counts
       condition="condition", name="condition", shape="condition")
dev.off()



# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data

# 100 top expressed miRNAs with heatmap.2
select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:100]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.3, labRow=F,
          main="100 Top Expressed mRNAs")
dev.copy(svg, paste0(outputPrefix, "-HEATMAP.svg"))
dev.off()


# save both files in the corresponding csv tables the normalized ouput for significant genes
write.csv(resdata1, file = paste0(outputPrefix, "-results-with-normalized_significant.csv"))

# the overall expression of genes
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized.csv"))