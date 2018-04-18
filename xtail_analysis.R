#source("http://bioconductor.org/biocLite.R")
#biocLite("devtools")    # only if devtools not yet installed

library('xtail')

args <- commandArgs(trailingOnly = TRUE)
RPFs<-strsplit(args[1], ",")[[1]]
RNAs<-strsplit(args[2], ",")[[1]]
condition<-strsplit(args[3], ",")[[1]]
base_dir<-args[4]


datalist = lapply(RPFs, function(x){read.table(file=x,header=FALSE,col.names=c("Gene", sub(".trim.uncontam.cnt", "", sub(getwd(),"",x))))})
rpfm <- Reduce(function(...) merge(..., by=1), datalist)
rownames(rpfm) <- rpfm[,1]
rpfm <- rpfm[,-1]
write.table(rpfm, file.path(base_dir,"counttable_rpfm.txt"), quote=F, sep="\t", col.names=NA )

datalist = lapply(RNAs, function(x){read.table(file=x,header=FALSE,col.names=c("Gene", sub(".trim.uncontam.cnt", "", sub(getwd(),"",x))))})
rnam <- Reduce(function(...) merge(..., by=1), datalist)
rownames(rnam) <- rnam[,1]
rnam <- rnam[,-1]
write.table(rnam, file.path(base_dir,"counttable_rnam.txt"), quote=F, sep="\t", col.names=NA )

test.results<-xtail(rnam,rpfm,condition)
summary(test.results)
write.xtail(test.results, file.path(base_dir,"xtail_analysis_table.txt"), quote=F, sep="\t", col.names=NA )

test.tab <- resultsTable(test.results, log2FCs=T, log2Rs=T)
write.table(test.tab, file.path(base_dir,"xtail_analysis_TE.txt"), quote=F, sep="\t", col.names=NA )
