#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

args = commandArgs(trailingOnly=TRUE)

packages <- c("GenomicRanges","ggplot2")

zzz<-lapply(packages, function(xxx) suppressMessages(library(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

source(file.path(script.dirname,"functionMul.R"))

# Rscript $HOME/Aimin/IS-Seq-python3/R/ScriptGetIScluster.R ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/multihit_allSites.rds ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE_multihitData/Results.RData ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE/Results.RData ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80 20

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_multi_hit=args[1]
  input_multi_hit_MLE=args[2]
  input_uniq_hit_MLE=args[3]
  output.dir=args[4]
  num.is.clusters=args[5]
}

GetJC <- function(df) {
  d <- sapply(names(df), function(x) sapply(names(df), function(y) {
    n <- length(intersect(df[[x]],df[[y]]))
    m <- length(union(df[[x]],df[[y]]))
    d <- n/m
    d
  }))
  d
}

GetReaId4MultihitIs <- function(input_multi_hit,input_multi_hit_MLE,input_uniq_hit_MLE,output.dir,num.is.clusters) {
  
  #input_multi_hit <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/multihit_allSites.rds'
  
  input_multi_hit <- '/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/multihit_allSites.rds'
  
  multihit.allSites <- readRDS(input_multi_hit)
  multihit.allSites.reduced <- flank(multihit.allSites, -1, start=TRUE)
  
  
  GetMerged <- function(multihit.allSites.reduced) {
    multihit.allSites.reduced.1 <- GenomicRanges::reduce(multihit.allSites.reduced, min.gapwidth = 7, with.revmap=TRUE)
    multihit.allSites.reduced.1$counts <- sapply(multihit.allSites.reduced.1$revmap, length)
    multihit.allSites.reduced.1
  }
  
  
  GetEachIs <- function(multihit.allSites.reduced) {
    multihit.allSites.reduced.1 <- GenomicRanges::reduce(multihit.allSites.reduced, min.gapwidth = 0, with.revmap=TRUE)
    multihit.allSites.reduced.1$counts <- sapply(multihit.allSites.reduced.1$revmap, length)
    multihit.allSites.reduced.1
  }
  
  IS.merged <- GetMerged(multihit.allSites.reduced) 
  
  IS.each <- GetEachIs(multihit.allSites.reduced)
  
  input_uniq_hit <- '/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/allSites.rds'
  
  uniq.hit.allSites <- readRDS(input_uniq_hit)
  uniq.hit.allSites.reduced <- flank(uniq.hit.allSites, -1, start=TRUE)
  
  IS.merged.uniq <- GetMerged(uniq.hit.allSites.reduced) 
  
  IS.each.uniq <- GetEachIs(uniq.hit.allSites.reduced)
  
  
  multihit.allSites.reduced$IS.index <- paste0(seqnames(multihit.allSites.reduced),'_',start(multihit.allSites.reduced),'_',strand(multihit.allSites.reduced))
  
  IS.index <- unique(multihit.allSites.reduced$IS.index)
  
  read.ID <- lapply(IS.index, function(u){
    ID <- multihit.allSites.reduced[multihit.allSites.reduced$IS.index==u,]$ID
    ID
  })
  
  names(read.ID) <- IS.index
  
  d <- 1-GetJC(read.ID)
  
  #IS.name <- colnames(d)
  
  #colnames(d) <- seq(1,782,1)
  #rownames(d) <- seq(1,782,1)
  
  hc.rd <- hclust(as.dist(d))
  plot(hc.rd, labels = colnames(d), xlab="",ylab = "Share # reads", main = "Cluster Dendrogram between ISs")
  
  group.index<- cutree(hc.rd,k=num.is.clusters)
  
  N <- as.data.frame(table(group.index))
  
  readName <- lapply(N$group.index, function(u){
    
    z <- names(group.index)[which(group.index==u)]
    
  })
  
  #input_multi_hit_MLE <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE_multihitData/Results.RData'
  
  input_multi_hit_MLE <- "/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData"
  
  multi.hit.MLE <- GetIsRel(input_multi_hit_MLE)
  
  #IS.name[which(group.index==1)]
  
  multi.hit.MLE$is.index <- paste0(seqnames(multi.hit.MLE),'_',start(multi.hit.MLE),'_',strand(multi.hit.MLE)) 
  
  group.index.id <- unique(group.index)
  
  clusters <- lapply(1:length(group.index.id), function(u){
    
    print(group.index.id[[u]])
    
    reads.key <- names(group.index)[which(group.index==group.index.id[[u]])]
    
    reads.key.1 <- as.data.frame(do.call(rbind,strsplit(reads.key,split = "_")))
    
    colnames(reads.key.1) <- c('seqnames','start','strand')
    
    reads.key.1$end <- reads.key.1$start
    
    reads.key.1.gr <- makeGRangesFromDataFrame(reads.key.1)
    
    ov <- findOverlaps(reads.key.1.gr,multi.hit.MLE,maxgap=7)
    

    cluster <- multi.hit.MLE[subjectHits(ov)]
    
    cluster$Frag.MLE.cluster <- sum(cluster$Frag.MLE)
    
    cluster$cluster.ID <- paste0('Cluster',u)
    
    cluster
  })
  
  clusters.all <- do.call(c,clusters)
  
  #input_uniq_hit_MLE <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE/Results.RData'
  
  #input_multi_hit_MLE <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE_multihitData/Results.RData'
  
  uniq.hit.MLE <- GetIsRel(input_uniq_hit_MLE)
  
  uniq.hit.MLE$is.index <- paste0(seqnames(uniq.hit.MLE),'_',start(uniq.hit.MLE),'_',strand(uniq.hit.MLE)) 
  
  uniq.hit.MLE$Frag.MLE.cluster <- uniq.hit.MLE$Frag.MLE
    
  uniq.hit.MLE$cluster.ID <- uniq.hit.MLE$is.index
  
  all.IS <- c(uniq.hit.MLE,clusters.all)
  
  write.table(as.data.frame(all.IS),file.path(output.dir,'CombinedUniqCluster.txt'),quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)
  
  #res <- data.frame(IS=all.IS$cluster.ID, abundance=all.IS$Frag.MLE.cluster)
  
  #res.1 <- unique(res)
  
  #res.2 <- res.1[order(res.1$abundance,decreasing = T),]
  
  #d1 <- d[upper.tri(d, diag = FALSE)] 
  
  #library(rstatix)
  
  #d2 <- d %>% pull_upper_triangle()
  
  #CDF <- ecdf(d1)
  #plot(CDF)
  
}

GetReaId4MultihitIs(input_multi_hit,input_multi_hit_MLE,input_uniq_hit_MLE,output.dir,num.is.clusters)

