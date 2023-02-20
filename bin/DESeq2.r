#!/usr/bin/env Rscript
if (!require("data.table")){
  install.packages("data.table", dependencies=TRUE)
  library("data.table")
}
if (!require("reshape2")) {
  install.packages("reshape2", dependencies=TRUE)
  library("reshape2")
}
if (!require("writexl")) {
  install.packages("writexl", dependencies=TRUE)
  library("writexl")
}
if (!require("readxl")) {
  install.packages("readxl", dependencies=TRUE)
  library("readxl")
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE, repos='http://cloud.r-project.org/')
}

suppressMessages(library(DESeq2))

suppressMessages(library(dplyr))

suppressMessages(library(gplots))

suppressMessages(library(ggplot2))

suppressMessages(library(ggrepel))

suppressMessages(library(RColorBrewer))

suppressMessages(library(EnhancedVolcano))

suppressMessages(library(reshape2))

suppressMessages(library(writexl))

suppressMessages(library(readxl))
args <- commandArgs(trailingOnly = TRUE)
#arg 1 - input file name
input <- as.character(args[3:(length(args))])

#arg 2 - category file 
category_file= args[1]
paired_samples=args[2]
print(paired_samples)
filelist<-list()
filelist<-input
header<-names(filelist)

# Read the categories file
categories=read.csv(category_file,header=FALSE)

# Read the categories and format the filenames
categories$V2=as.factor(categories$V2)
number_levels=nlevels(as.factor(categories$V2))
categories$V1 =gsub(".fq.gz","",categories$V1)
categories$V1 =gsub(".fastq.gz","",categories$V1)
  
# Prepare the combined data frame with gene ID as rownames and sample ID as colname
data<-do.call("cbind", lapply(filelist, fread, header=FALSE, select=c(3)))
unmapped<-do.call("cbind", lapply(filelist, fread, header=FALSE, select=c(4)))
data<-as.data.frame(data)
unmapped<-as.data.frame(unmapped)
temp <- fread(filelist[1],header=FALSE, select=c(1))
rownames(data)<-temp$V1
rownames(unmapped)<-temp$V1
filenames=data.frame(filelist)
for(category in categories$V1){
  match=grep(category,filenames[,1])
  filenames[match,]=category
}
colnames(data)<-filenames[,1]
colnames(unmapped)<-filenames[,1]

data<-data[rownames(data)!="*",,drop=FALSE]
unmapped<-unmapped[rownames(unmapped)=="*",,drop=FALSE]

# Write the summary table of unmapped reads
write.table(unmapped,file="mature_unmapped_read_counts.txt",sep='\t',quote=FALSE)

# Remove genes with 0 reads in all samples
row_sub = apply(data, 1, function(row) all(row ==0 ))
# Only subset if at least one sample is remaining
nr_keep <- table(row_sub)
nr_keep <- as.numeric(nr_keep[names(nr_keep) == FALSE])
if (!is.null(nr_keep) & length(nr_keep)>0 & nr_keep > 0){
  data<-data[!row_sub,]
}

write.csv(data,file="normalized_counts.csv")
RNAseq_data=data


# Check if two categories are present because DESeq2 requires two categories
if(number_levels!=2)
{
  stop("Two condtions needed in the category file")
}

# Take the first category in the category file as base level and second as the level to compare
categories$V2=  factor(categories$V2, levels=as.character(unique(categories$V2)))
base_level=levels(categories$V2)[1]
level_to_compare=levels(categories$V2)[2]

# Format RNA seq count data to fit the DESeqDataSetFromMatrix countData format
RNAseq_data <- RNAseq_data[-1, ] 

# Creat col_data to match DESeqDataSetFromMatrix colData 
col_data=data.frame(categories$V2)
colnames(col_data)=c("condition")
rownames(col_data)=categories$V1

# Order the col_data to fullfil rownames(colData) == colnames(countData) condition in DESeqDataSetFromMatrix
colOrder =(rownames(col_data))
RNAseq_data = RNAseq_data[,colOrder]

#  Format data for DESeqDataSetFromMatrix
data= as.data.frame(apply(RNAseq_data, 2, as.numeric))
rownames(data)=rownames(RNAseq_data)
data=as.matrix(data)

#paired analysis
if(paired_samples){
  if(ncol(categories)<3){
    stop("Please specify the sample names, conditions and subject details for paired analysis")
  }
  # DESeq2 analysis for paired samples, takes in un-normalised counts
  col_data$subject=categories$V3
  
  dds <- DESeqDataSetFromMatrix(countData=data, col_data, formula(~ subject+condition))
}else{
  # DESeq2 analysis for unpaired samples, takes in un-normalised counts
  dds <- DESeqDataSetFromMatrix(countData=data, col_data, formula(~ condition))
}

dds <- DESeq(dds)

# Nomralised counts 
normalised_counts=counts(dds,normalized=TRUE)
un_nomralised_counts=counts(dds,normalized=FALSE)

res <- results(dds, contrast = c("condition",level_to_compare,base_level), cooksCutoff=FALSE)

results.DF <- data.frame(res)

res <- res[order(res$log2FoldChange, decreasing=TRUE),]
results.DF <- results.DF %>% arrange(pvalue)
resultsDF=merge(results.DF,normalised_counts,by="row.names")
results.subset.DF <- results.DF %>% filter(pvalue < 0.05)

up <- (resultsDF[!is.na(resultsDF$pvalue) & resultsDF$pvalue <= 0.05 &
                   resultsDF$log2FoldChange >= 0, ])
down <- (resultsDF[!is.na(resultsDF$pvalue) & resultsDF$pvalue <= 0.05 &    
                     resultsDF$log2FoldChange <= 0, ]) 
write_to_file=cbind(resultsDF[,1],resultsDF[ , c("log2FoldChange", "pvalue", "padj")] , resultsDF[,8:ncol(resultsDF)])
names(write_to_file)[names(write_to_file) == 'rownames(results.DF)'] <- 'ID'
file_name=paste(base_level,"_vs_",level_to_compare,".csv",sep="")
write.csv(write_to_file,file=file_name)




sample_numbers=categories %>% group_by(V2) %>% summarise(num=length(V2))
#density plot
RPM <- t(t(data)/colSums(data))*1e6
inlog <- log(RPM)
colLabel <- c(rep("#E41A1C",sample_numbers$num[1]), rep("#377EB8", sample_numbers$num[2]))
colTy <- c(rep(1, sample_numbers$num[1]), rep(2,sample_numbers$num[2]))
### Create a density plot
max.density <- 0
for(i in 1:ncol(RPM)){
  my.density <- density(inlog[,i])
  if(max(my.density$y) > max.density){
    max.density <- max(my.density$y)
  }
}
inlog.density.ymax <- max.density*1.2 # Get y max and add 20% 
width <- 6
height <- 6
density_plot_name=paste0(base_level,"_vs_",level_to_compare,"_density_plot.pdf", sep="")
pdf(density_plot_name, height=height, width=width)
density_plot=plot(density(inlog[,1]),
                  #ylim=c(0,0.2), 
                  ylim=c(0,inlog.density.ymax),
                  main="Density plot of counts per gene",
                  lty=colTy[1],
                  xlab="Log of RPM per gene",
                  ylab="Density",
                  col=colLabel[1])
for(i in 2:ncol(RPM)){
  lines(density(inlog[,i]), lty=colTy[i], col=colLabel[i])
}
legend("topright", legend=colnames(RPM), lty=colTy, col=colLabel)
dev.off()


rld <- varianceStabilizingTransformation(dds)
#rld <- rlogTransformation(dds)  # Alternative

### Convert the rld transformation into a matrix
rld.matrx <- assay(rld)   
rld.df <- data.frame(rld.matrx)
### Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
sampleDists <- data.frame(sampleDists)
#par(mar=c(6,4,4,5)+0.1) 
width <- 6
height <- 6
distance_matrix_plot_name=paste0(base_level,"_vs_",level_to_compare,"_distance_matrix_plot.pdf", sep="")


pdf(distance_matrix_plot_name, height=height, width=width)
heatmap_seq=heatmap.2(as.matrix(sampleDists), key=F, trace="none",
                      col=colorpanel(100, "black", "white"),
                      #ColSideColors=mycols[file.names], 
                      #RowSideColors=mycols[file.names],
                      #cexRow = 0.8,
                      #cexCol = 0.8,
                      margins=c(12,10),
                      srtCol=45,
                      main="Sample Distance Matrix")
dev.off()



### Principal component analysis
#Tpm PCA plot (not DESeq2)
d <- dist(t(RPM))
fit=cmdscale(d, eig=TRUE, k=2)
x=fit$points[,1]
y=fit$points[,2]

### ggplot2 scatterplot
gg.df <- data.frame(cbind("PC1" = x, "PC2" = y))
gg.df$Group <- c(rep(sample_numbers$V2[1], sample_numbers$num[1]), rep(sample_numbers$V2[2], sample_numbers$num[2]))
pca_plot=ggplot(gg.df, aes(PC1, PC2, color = Group, shape = Group)) + 
  geom_point() + 
  geom_text_repel(label=rownames(gg.df), show.legend = F) +
  ggtitle("Principal Component Analysis")


width <- 6
height <- 6
pca_plot_name=paste0(base_level,"_vs_",level_to_compare,"_PCA_plot.pdf", sep="")
ggsave(pca_plot_name, pca_plot, width = width, height = height)

