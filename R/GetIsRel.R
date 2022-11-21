output.dir <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80'

input <- "~/Library/CloudStorage/OneDrive-SanaBiotechnology/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE/Results.RData"

GetIsRel <- function(input) {
  load(input)
  IS.grouped$REL <- IS.grouped[,4]/sum(IS.grouped[,4])
  colnames(IS.grouped)[4] <- 'Frag.MLE'
  IS.grouped.sorted <- IS.grouped[order(IS.grouped$Frag.MLE,decreasing = T),]
  IS.grouped.sorted
  IS.grouped.gr <- makeGRangesFromDataFrame(data.frame(chr=IS.grouped.sorted$seqnames,start=IS.grouped.sorted$start,end=IS.grouped.sorted$start,strand=IS.grouped.sorted$strand,Frag.MLE=IS.grouped.sorted$Frag.MLE,REL=IS.grouped.sorted$REL),keep.extra.column=T)
  IS.grouped.gr
}

IS.all <- GetIsRel(input)

input.95 <- "~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/FragMLE/Results.RData"

IS.all.95 <- GetIsRel(input.95)

output.dir.95 <- "~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95"

IS.unique.mapping.95 <- MakeDataSet(IS.all.95)
GenerateBarPlot(IS.unique.mapping.95,output.dir.95,'Uniqhit')

multihit.95 <- readRDS('~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/multihitData.rds')
saveRDS(multihit.95$unclusteredMultihits,'~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/multihit_allSites.rds')

input.95.multihit <- "~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/FragMLE_multihitData/Results.RData"
IS.multihit.mapping.95 <- GetIsRel(input.95.multihit)

IS.multihit.mapping.95.1 <- MakeDataSet(IS.multihit.mapping.95)

GenerateBarPlot(IS.multihit.mapping.95.1,output.dir.95,'Multihit')




input.multihit <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/FragMLE_multihitData/Results.RData'

IS.multihit <- GetIsRel(input.multihit)

CombineUniqAndMulti <- function(IS.all, IS.multihit) {
  
  ov <- findOverlaps(IS.all,IS.multihit,maxgap=7)
  
  #IS.multihit[subjectHits(ov)]
  #IS.multihit[-subjectHits(ov)]
  
  IS.common <- IS.all[queryHits(ov)]
  
  IS.common$Frag.MLE <- IS.all[queryHits(ov)]$Frag.MLE+IS.multihit[subjectHits(ov)]$Frag.MLE
  
  IS.single.multihit <- c(IS.all[-queryHits(ov)],IS.common,IS.multihit[-subjectHits(ov)])
  
  IS.single.multihit$REL <- IS.single.multihit$Frag.MLE/sum(IS.single.multihit$Frag.MLE)
  
  IS.single.multihit[order(IS.single.multihit$REL,decreasing = T)]
  
  IS.all.with.multihit <- c(IS.all[-queryHits(ov)],IS.common)
  
  
  IS.all.with.multihit$REL <- IS.all.with.multihit$Frag.MLE/sum(IS.all.with.multihit$Frag.MLE)
  
  IS.all.with.multihit[order(IS.all.with.multihit$REL,decreasing = T)]
  IS.all.with.multihi
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

IS.all.with.multihit.95 <- CombineUniqAndMulti1(IS.all.95,IS.multihit.mapping.95)

IS.all.with.multihit.95.top.10 <- MakeDataSet(IS.all.with.multihit.95)
GenerateBarPlot(IS.all.with.multihit.95.top.10,output.dir,"UniqAndMultihit")


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

library(ggplot2)

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


all.IS <- MakeDataSet(IS.all.with.multihit)
GenerateBarPlot(all.IS,output.dir,"UniqAndMultihit")

IS.unique.mapping <- MakeDataSet(IS.all)
GenerateBarPlot(IS.unique.mapping,output.dir,'Uniqhit')

IS.multihit.mapping <- MakeDataSet(IS.multihit)
GenerateBarPlot(IS.multihit.mapping,output.dir,'Multihit')


allSites.input <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/allSites.rds'

GetIsReadId <- function(allSites.input,index) {
  allSites <- readRDS(allSites.input)
  sites.reduced <- flank(allSites, -1, start=TRUE)
  ov1 <- findOverlaps(IS.common[index],sites.reduced)
  reads.ID <- sites.reduced[subjectHits(ov1)]$ID
  reads.ID
}


IS.common

multihit.allSites.input <- '~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/multihit_allSites.rds'

num.common.reads <- lapply(1:11, function(u){
  
  read.from.single <- GetIsReadId(allSites.input,u)
  read.from.multihit <- GetIsReadId(multihit.allSites.input,u)
  
  common.reads <- intersect(read.from.multihit,read.from.single)
  
  n <- length(common.reads)
  
  n
})

allSites <- readRDS(allSites.input)

GetReaId4Is <- function() {
  
  multihit.allSites <- readRDS('~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/multihit_allSites.rds')
  
  multihit.allSites.reduced <- flank(multihit.allSites, -1, start=TRUE)
  
  multihit.allSites.reduced$IS.index <- paste0(seqnames(multihit.allSites.reduced),'_',start(multihit.allSites.reduced),'_',strand(multihit.allSites.reduced))
         
  IS.index <- unique(multihit.allSites.reduced$IS.index)
  
read.ID <- lapply(IS.index, function(u){
    ID <- multihit.allSites.reduced[multihit.allSites.reduced$IS.index==u,]$ID
    ID
  })
 
names(read.ID) <- IS.index

GetJC <- function(df) {
  d <- sapply(names(df), function(x) sapply(names(df), function(y) {
      n <- length(intersect(df[[x]],df[[y]]))
      m <- length(union(df[[x]],df[[y]]))
      d <- n/m
      d
  }))
  d
}

d <- 1-GetJC(read.ID)

IS.name <- colnames(d)

colnames(d) <- seq(1,782,1)
rownames(d) <- seq(1,782,1)

hc.rd <- hclust(as.dist(d))
plot(hc.rd, labels = colnames(d), xlab="",ylab = "Share # reads", main = "Cluster Dendrogram between ISs")

group.index<- cutree(hc.rd,k=20)

N <- lapply(1:20, function(u){
  
  n <- length(IS.name[which(group.index==u)])
  
})


IS.name[which(group.index==1)]

d1 <- d[upper.tri(d, diag = FALSE)] 

library(rstatix)

d2 <- d %>% pull_upper_triangle()

CDF <- ecdf(d1)
plot(CDF)

}


