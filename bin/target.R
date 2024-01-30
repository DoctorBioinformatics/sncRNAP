#!/usr/bin/env Rscript

# add miRNA targets from mirDB to results file 

targte_file= sub("\\.zip$", "", args[1])
#arg[1] should be miRDB_v6.0_human_target_prediction_result.tsv
targets=read.csv(unz(args[1], paste0(targte_file,".tsv")),sep="\t",header=F) #"args[1]=miRDB_v6.0_human_target_prediction_result"
colnames(targets)=c("miRNA","target","score")

#arg[2] should be the csv result file from DEseq2
read_counts=read.csv(arg[2])#C_vs_T.csv"
read_counts$miRNA <- sub("_.*$", "", read_counts$resultsDF...1.)
read_counts=merge(read_counts,targets,by.x = "miRNA",by.y = "miRNA", all.x = T)
read_counts$miRNA=NULL
write.csv(filename=arg[2], read_counts)