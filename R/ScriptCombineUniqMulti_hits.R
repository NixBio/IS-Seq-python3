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

# Rscript $HOME/Aimin/IS-Seq-python3/R/ScriptCombineUniqMulti_hits.R /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0

# Rscript $HOME/Aimin/IS-Seq-python3/R/ScriptCombineUniqMulti_hits.R /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/FragMLE/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/FragMLE_multihitData/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0

# Rscript $HOME/Aimin/IS-Seq-python3/R/ScriptCombineUniqMulti_hits.R /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE_multihitData/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80

# Rscript $HOME/Aimin/IS-Seq-python3/R/ScriptCombineUniqMulti_hits.R /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/FragMLE/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/FragMLE_multihitData/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95

#/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0

# Rscript $HOME/Aimin/IS-Seq-python3/R/ScriptCombineUniqMulti_hits.R ~/Aimin/ISseqOutput/vcn/Uniq/MOI30CLC9_IS0/Results.RData ~/Aimin/ISseqOutput/vcn/Uniq/MOI30CLC9_IS0/Results.RData ~/Aimin/ISseqOutput/vcn/UniqAndMulti/MOI30CLC9_IS0

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_uniq_hit=args[1]
  input_multi_hit=args[2]
  output.dir=args[3]
}

library(ggplot2)

#multihit.95 <- readRDS('/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/multihitData.rds')

#saveRDS(multihit.95$unclusteredMultihits,'/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/multihit_allSites.rds')

#input.95.multihit <- "/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData"

#load("/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData")

GetIsRel <- function(input) {
  load(input)
  IS.grouped$REL <- IS.grouped[,4]/sum(IS.grouped[,4])
  colnames(IS.grouped)[4] <- 'Frag.MLE'
  IS.grouped.sorted <- IS.grouped[order(IS.grouped$Frag.MLE,decreasing = T),]
  IS.grouped.sorted
  IS.grouped.gr <- makeGRangesFromDataFrame(data.frame(chr=IS.grouped.sorted$seqnames,start=IS.grouped.sorted$start,end=IS.grouped.sorted$start,strand=IS.grouped.sorted$strand,Frag.MLE=IS.grouped.sorted$Frag.MLE,REL=IS.grouped.sorted$REL),keep.extra.column=T)
  IS.grouped.gr
}

CombineUniqAndMulti1 <- function(IS.all, IS.multihit) {
  
  ov <- suppressWarnings(findOverlaps(IS.all, IS.multihit,maxgap=7))
  
  #IS.multihit[subjectHits(ov)]
  #IS.multihit[-subjectHits(ov)]
  
  IS.common <- IS.all[queryHits(ov)]
  #IS.common
  IS.common$Frag.MLE <- IS.all[queryHits(ov)]$Frag.MLE+IS.multihit[subjectHits(ov)]$Frag.MLE
  
  #IS.common
  IS.single.multihit <- suppressWarnings(c(IS.all[-queryHits(ov)],IS.common,IS.multihit[-subjectHits(ov)]))
  
  IS.single.multihit
  IS.single.multihit$REL <- IS.single.multihit$Frag.MLE/sum(IS.single.multihit$Frag.MLE)
  
  IS.single.multihit[order(IS.single.multihit$REL,decreasing = T)]
  
  IS.all.with.multihit <- c(IS.all[-queryHits(ov)],IS.common)
  IS.all.with.multihit
  
  IS.all.with.multihit$REL <- IS.all.with.multihit$Frag.MLE/sum(IS.all.with.multihit$Frag.MLE)
  
  IS.all.with.multihit[order(IS.all.with.multihit$REL,decreasing = T)]
  IS.all.with.multihit
}

CombineUniqAndMultiUseSpecificIS <- function(IS.all, IS.multihit) {
  
 
  #IS.multihit[subjectHits(ov)]
  #IS.multihit[-subjectHits(ov)]
  
  #IS.all <- IS.uniq.hit
  #IS.multihit <-  IS.multihit 
  
  ov <- suppressWarnings(findOverlaps(IS.all, IS.multihit,maxgap=7))
  
  IS.all$is.index <- paste0(seqnames(IS.all),'_',start(IS.all),'_',strand(IS.all))
  
  one.uniq.IS <-  IS.all[which(IS.all$is.index %in% c('chr7_75504993_-'))]
  
  IS.common <- IS.multihit[subjectHits(ov)]
  IS.common$is.index <- paste0(seqnames(IS.common),'_',start(IS.common),'_',strand(IS.common))
  one.IS.common <- IS.common[which(IS.common$is.index %in% c('chr7_75504993_-'))]
  
  IS.all[which(IS.all$is.index %in% c('chr7_75504993_-'))]$Frag.MLE <- one.uniq.IS$Frag.MLE+one.IS.common$Frag.MLE
  
  IS.all$REL <- IS.all$Frag.MLE/sum(IS.all$Frag.MLE)
  
  IS.all
  
}

MakeDataSet <- function(IS.all.with.multihit) {
  
  n <- length(IS.all.with.multihit)
  
  if(n>10){
    top.10.IS.df <- as.data.frame(IS.all.with.multihit[order(IS.all.with.multihit$REL,decreasing = T)][c(1:10)])
    
    other.IS.df <- as.data.frame(IS.all.with.multihit[order(IS.all.with.multihit$REL,decreasing = T)][-c(1:10)])
    
    top.10.IS.df.1 <- data.frame(IS=paste0(top.10.IS.df$seqnames,'_',top.10.IS.df$start,'_',top.10.IS.df$strand),Frag.MLE=top.10.IS.df$Frag.MLE,REL=top.10.IS.df$REL)
    
    other.IS.df.1 <- data.frame(IS='Other.IS',Frag.MLE=sum(other.IS.df$Frag.MLE),REL=sum(other.IS.df$REL))
    
    all.IS <- rbind(top.10.IS.df.1,other.IS.df.1)
    
  }else
  {
    top.10.IS.df <- as.data.frame(IS.all.with.multihit[order(IS.all.with.multihit$REL,decreasing = T)][c(1:n)])
    top.10.IS.df.1 <- data.frame(IS=paste0(top.10.IS.df$seqnames,'_',top.10.IS.df$start,'_',top.10.IS.df$strand),Frag.MLE=top.10.IS.df$Frag.MLE,REL=top.10.IS.df$REL)
    all.IS <- top.10.IS.df.1
  }
  
  all.IS
  
}

GenerateBarPlot <- function(all.IS,outout.dir,title) {
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  library(tidyverse)
  all.IS$REL <- round(100*all.IS$REL,2)
  
  sp <- ggplot(all.IS, aes(x=fct_inorder(IS), y = REL)) +
    geom_bar(stat="identity",position = "dodge")  + geom_text(aes(label=REL),vjust = -0.5, colour = "black")+
    theme(plot.title=element_text(hjust=0.5)) + ylim(0,100) +
    ylab("Relative Abundance")  +  theme(legend.position="none")+xlab("IS")
  
  sp <- sp + theme_bw() + ggtitle(title)+theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7,color="black"),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(output.dir,paste0(title,".png")),plot = sp)
  
  sp
  
}

#input_multi_hit <- '/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/FragMLE_multihitData/Results.RData' 

IS.multihit <- GetIsRel(input_multi_hit)
IS.multihit.top.10 <- MakeDataSet(IS.multihit)
GenerateBarPlot(IS.multihit.top.10,output.dir,'Multihit')

#input.0.CL6 <- "/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE/Results.RData"

#input_uniq_hit <- '/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/FragMLE/Results.RData'

IS.uniq.hit <- GetIsRel(input_uniq_hit)
IS.unique.top.10 <- MakeDataSet(IS.uniq.hit)
GenerateBarPlot(IS.unique.top.10,output.dir,'Uniqhit')

IS.uniq.multihit <- CombineUniqAndMulti1(IS.uniq.hit,IS.multihit)
IS.uniq.multihit.top.10 <- MakeDataSet(IS.uniq.multihit)
GenerateBarPlot(IS.uniq.multihit.top.10,output.dir,"UniqAndMultihit")

IS.uniq.multihit.4.one.IS <- CombineUniqAndMultiUseSpecificIS(IS.uniq.hit,IS.multihit)
IS.uniq.multihit.4.one.IS.top.10 <- MakeDataSet(IS.uniq.multihit.4.one.IS )
GenerateBarPlot(IS.uniq.multihit.4.one.IS.top.10,output.dir,"UniqAndMultihit4OneIS")



