#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

packages <- c("optparse","stringr","dplyr","sonicLength","GenomicRanges","dplyr")

zzz<-lapply(packages, function(xxx) suppressMessages(library(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  input.multi.hit=args[1]
  input.uniq.hit=args[2]
  sampleName=args[3]
  SST=args[4]
  output.dir=args[5]
}

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R /local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa/MOI30CLB7/IS0/multihitData.rds /local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa/MOI30CLB7/IS0

#$HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE/Results.RData $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/allSites.rds MOI30CLB7 SST80 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS80/Uniq_cluster

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/allSites.rds MOI30CLB7 SST95 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS95/Uniq_cluster



# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS95/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS95/allSites.rds MOI50CLH6 SST95 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS95/Uniq_cluster

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS80/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS80/allSites.rds MOI50CLH6 SST80 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS80/Uniq_cluster

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS0/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS0/allSites.rds MOI50CLH6 SST0 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI50CLH6/IS0/Uniq_cluster

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/allSites.rds CL6 SST0 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/Uniq_cluster

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS80/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS80/allSites.rds CL6 SST80 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS80/Uniq_cluster

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReads.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS95/multihitData.rds $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS95/allSites.rds CL6 SST95 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS95/Uniq_cluster

#MOI50CLH6


#$HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0


#print(input.multi.hit)

#input.multi.hit <- '$HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/multihitData.rds'

#input.uniq.hit <- '$HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/allSites.rds'

#sampleName1 <- basename(dirname(dirname(input.multi.hit)))
#sampleName2 <- basename(dirname(input.multi.hit))

cmd <- paste0('Rscript ', file.path(script.dirname,'GetFragMLE.R '),input.uniq.hit,' ',paste0(sampleName,'_',SST),' ',output.dir)

cat(cmd,"\n")

system(cmd)

res.uniq <- readRDS(file.path(output.dir,'IS.grouped.rds'))

colnames(res.uniq) <- c('seqnames','start','strand','MLE')


multihit <- readRDS(input.multi.hit)
#library(GenomicRanges)

num.is <- length(multihit$unclusteredMultihits)

if(num.is!=0){

multihit.allSites.reduced <- flank(multihit$unclusteredMultihits, -1, start=TRUE)

n <- length(multihit$clusteredMultihitPositions)

# cluster <- lapply(1:n, function(u){
#   
#   x <- multihit$clusteredMultihitPositions[[u]]
#   
#   ov <- findOverlaps(x,multihit.allSites.reduced,maxgap=0)
#   
#   #multihit$clusteredMultihitPositions$`1`
#   
#   #multihit$clusteredMultihitPositions$`1`[queryHits(ov)]
#   
#   #multihit.allSites.reduced[subjectHits(ov)]
#   
#   res <- as.data.frame(table(width(multihit$unclusteredMultihits[match(names(multihit.allSites.reduced[subjectHits(ov)]),names(multihit$unclusteredMultihits))])))
#   
#   res
#   
# })


n <- length(multihit$clusteredMultihitPositions)
# library(dplyr)

# cluster1 <- lapply(1:n, function(u){
#   
#   x <- multihit$clusteredMultihitPositions[[u]]
#   
#   ov <- findOverlaps(x,multihit.allSites.reduced,maxgap=0)
#   
#   #multihit$clusteredMultihitPositions$`1`
#   
#   #multihit$clusteredMultihitPositions$`1`[queryHits(ov)]
#   
#   #multihit.allSites.reduced[subjectHits(ov)]
#   
#   res <- multihit$unclusteredMultihits[subjectHits(ov)]
#   
#   res1 <- flank(res, -1, start=TRUE)
#   
#   res1$Frag <- width(res)
#   
#   res2 <- data.frame(IS=paste0(seqnames(res1),'_',start(res1),'_',strand(res1)),Frag=res1$Frag)
#   
#   
#   #res2 %>% group_by(IS) %>% length(unique(Frag))
#   
#   res3 <- as.data.frame(res2 %>% group_by(IS) %>% summarise(n()))
#   
#   #res4 <- mean(res3$`n()`)
#   
#   #ov <- findOverlaps(res1,x,maxgap=0)
#   
#   # 
#   # res2 <- as.data.frame(table(width(multihit$unclusteredMultihits[match(names(res1[subjectHits(ov)]),names(res))])))
#   # 
#   # res2
#   
# })

cluster2 <- lapply(1:n, function(u){
  
  x <- multihit$clusteredMultihitPositions[[u]]
  
  ov <- findOverlaps(x,multihit.allSites.reduced,maxgap=0)
  
  #multihit$clusteredMultihitPositions$`1`
  
  #multihit$clusteredMultihitPositions$`1`[queryHits(ov)]
  
  #multihit.allSites.reduced[subjectHits(ov)]
  
  res <- multihit$unclusteredMultihits[subjectHits(ov)]
  
  res1 <- flank(res, -1, start=TRUE)
  
  res1$Frag <- width(res)
  
  #res2 <- as.data.frame(res1)
  
  m <- length(res1)
  
  res3 <- data.frame(id=paste(rep(paste0('chr_',paste0('cluster',u)),m),rep(paste0('start_',paste0('cluster',u)),m),rep(paste0('strand_',paste0('cluster',u)),m)),Frag=res1$Frag)
             
  res3
  #res2 <- data.frame(IS=paste0(seqnames(res1),'_',start(res1),'_',strand(res1)),Frag=res1$Frag)
  
  
  #res2 %>% group_by(IS) %>% length(unique(Frag))
  
  #res3 <- as.data.frame(res2 %>% group_by(IS) %>% summarise(n()))
  
  #res4 <- mean(res3$`n()`)
  
  #ov <- findOverlaps(res1,x,maxgap=0)
  
  # 
  # res2 <- as.data.frame(table(width(multihit$unclusteredMultihits[match(names(res1[subjectHits(ov)]),names(res))])))
  # 
  # res2
  
})

#names(cluster2) <- NULL

estAbund.df <- do.call(rbind,cluster2)
  
estAbund.df <- unique(estAbund.df)

id.estAbund.df <- estAbund.df$id

fit.real.INSP <- estAbund(id.estAbund.df, estAbund.df$Frag,min.length=min(estAbund.df$Frag))

res <- fit.real.INSP$theta[unique(id.estAbund.df)]
res <- as.data.frame(res)
colnames(res) <- 'MLE'
#res

res.multi <- data.frame(do.call(rbind,strsplit(row.names(res),split=" ")),MLE=as.integer(res$MLE))
colnames(res.multi) <- c('seqnames','start','strand','MLE')

res.uniq.multi <- rbind(res.uniq,res.multi)

}else{
  res.uniq.multi <- res.uniq
}

#temp1 <- multihit.allSites.reduced[subjectHits(ov)][which(start(multihit.allSites.reduced[subjectHits(ov)])==48277454)]

#temp2 <- multihit.allSites.reduced[subjectHits(ov)][which(start(multihit.allSites.reduced[subjectHits(ov)])==59587069)]

#saveRDS(multihit$unclusteredMultihits,file.path(output.dir,paste0(sampleName1,"_",sampleName2,"_","multihit_allSites.rds")))




#              )
#Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetFragMLE.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/allSites.rds MOI30CLB7_IS0 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/uniq


# uniqhit <- readRDS(input.uniq.hit)
# 
# uniqhit.allSites.reduced <- flank(uniqhit, -1, start=TRUE)
# 
# estAbund.df <- data.frame(Chromosome=seqnames(uniqhit.allSites.reduced),Position=start(uniqhit.allSites.reduced),Ort=strand(uniqhit.allSites.reduced),Frag=width(uniqhit))
# 
# estAbund.df <- unique(estAbund.df)
# id.estAbund.df <- with(estAbund.df, paste(Chromosome, Position, Ort))
# 
# fit.real.INSP <- estAbund(id.estAbund.df, estAbund.df$Frag,min.length=min(estAbund.df$Frag))
# res.uniq <- fit.real.INSP$theta[unique(id.estAbund.df)]
# res.uniq <- as.data.frame(res.uniq)
# colnames(res.uniq) <- 'MLE'
# res.uniq


# sampleName1 <- basename(dirname(dirname(input.multi.hit)))
# sampleName2 <- basename(dirname(input.multi.hit))
# 
# multihit <- readRDS(input.multi.hit)
# library(GenomicRanges)
# 
# multihit.allSites.reduced <- flank(multihit$unclusteredMultihits, -1, start=TRUE)



#output.directory.ISSeq.paper <- '$HOME/IS-Seq_output/IS-Seq/submission/final_from_luca/ReSubmission'

#output.file <- file.path(output.directory.ISSeq.paper,paste0('MOI30CLB7_uniq_cluster_Frag_MLE_SST0','.xls'))

#output.file <- file.path(output.directory.ISSeq.paper,paste0('MOI30CLB7_uniq_cluster_Frag_MLE_SST0','.xls'))

output.file <- file.path(output.dir,paste0(sampleName,'_uniq_cluster_Frag_MLE_',SST,'.xls'))

write.table(res.uniq.multi,file = output.file, append = FALSE, quote=FALSE,row.names=FALSE,col.names = TRUE,sep = "\t")

# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetFragMLE.R $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/allSites.rds MOI30CLB7_IS0 $HOME/IS-Seq_output/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/MOI30CLB7/IS0/uniq


