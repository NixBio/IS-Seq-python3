#!/usr/bin/env Rscript

packages <- c("optparse","stringr","curl","rtracklayer")

zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="input file", metavar="character"),
  make_option(c("-g", "--gene_annotation"), type="character", default=NULL,
          help="input gene annotation", metavar="character"),
  make_option(c("-r", "--rmsk"), type="character", default=NULL,
              help="input rmsk", metavar="character"),
  make_option(c("-m", "--chrom"), type="character", default=NULL,
              help="input chromInfo", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default=NULL,
              help="output file", metavar="character")
);

example.use <- "Example: Rscript $HOME/IS-Seq-python3/R/makeREFIndex1.R -i https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz -g https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz -r https://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/database/rmsk.txt.gz -m https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz -o /home/ayan/user/ispipe/utilsRefData/hg38/GRCh38.primary_assembly.genome.fa\n"

opt_parser = OptionParser(option_list=option_list,epilogue=example.use);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input count file)", call.=FALSE)
}

input.file <- opt$input_file
input.gene.annotation <- opt$gene_annotation
input.rmsk <- opt$rmsk
input.chrom <- opt$chrom

output.file <- opt$out_file

out.dir.name <- dirname(output.file)
if (!dir.exists(out.dir.name)){dir.create(out.dir.name, recursive = TRUE)}

t.input <- file.info(input.file)$mtime
t.output <- file.info(output.file)$mtime

x <- dirname(output.file)
y <- basename(x)
z <- file.path(x,paste0(y,"ChrOnly.fa"))

if(is.na(t.output)){
  
  cat("Download fasta file\n")
  cmd0 <- paste0('wget ',input.file,' -P ',x)
  system(cmd0)
  
  cat("gunzip file\n")
  cmd1 <- paste0('gunzip ',paste0(output.file,'.gz'))
  system(cmd1)
  
  cat("samtools faidx fasta file\n")
  cmd2 <- paste0('samtools faidx ',output.file)
  system(cmd2)
  
  
  cat("get chr1-22 and X, Y,M fasta file only\n")
  chr<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
  
  
  z.output <- file.info(z)$mtime
  
  if(is.na(z.output)){
     
  null <- lapply(chr, function(u){
    
    cmd3 <- paste0('samtools faidx ',output.file,' ',u,' >> ',z)
    system(cmd3)
    
  })
    
  }
  
  cat("make bwa index\n")
  cmd4= paste0('bwa-mem2 index ',z)
  system(cmd4)
  
}

w <- file.path(x,paste0('repeatMasker',y,'BED'))

t.w <- file.info(w)$mtime

if(is.na(t.w)){
  
cat("get repeatMaskerBED\n")

  url = input.rmsk
  tmp <- tempfile()
  curl_download(url, tmp)
  
  zz=gzfile(tmp,'rt')
  
  dat=read.table(zz,header=F,sep = "\t")
  dat=dat[,c(6,7,8,11,2,10)]
  
  colnames(dat) <- c("genoName","genoStart","genoEnd","repName","swScore","strand")
  
  dat=dat[order(dat$genoName,dat$genoStart),]
  
  write.table(dat,file = w, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

}

library(readr)

zz <- "/local_scratch/ISseqOutput1/vector_data_for_IS-Seq/PLAS000066_LVV_clean_bed.csv"

dir.name <- '/local_scratch/ISseqOutput1/vector_data_for_IS-Seq/vector'

if(!dir.exists(dir.name)){
  
  dir.create(dir.name)
}

y <- 'MND-GFP'
gene.file <- file.path(dir.name,paste0(y,'_genesKNOWN_sorted.bed'))


gene.file.time <- file.info(gene.file)$mtime

if(is.na(gene.file.time)){
  
  cat("get genesKNOWN_sorted.bed\n")
  
  url = input.gene.annotation
  tmp <- tempfile()
  curl_download(url, tmp)
  
  zz=gzfile(tmp,'rt')
  
  dat=read.table(zz,header=F,sep = ",")
  
  dat.gene <- dat[dat$V3=='gene',]
  
  gene.name <- str_sub(unlist(lapply(base::strsplit(as.character(dat.gene$V9),';'),function(u)u[3])),start = 12L)
  
  #MND-GFP
  
  gene.hg38.table <- data.frame(chr='MND-GFP',start=dat$V2,end=dat$V3,strand=dat$V4,geneName=gsub(" ","",dat$V1),type=rep('KNOWN',length(dat$V1)))
  
  write.table(gene.hg38.table,file = gene.file, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

}

outputFile1 <- file.path(x,"hg38.genome.sorted.txt")
outputFile2 <- file.path(x,"hg38.genome.num.txt")

hg38.genome.num.time <- file.info(outputFile2)$mtime

if(is.na(hg38.genome.num.time)){
  
cat("get hg38.genome.num file\n")

con <- gzcon(url(input.chrom))
txt <- readLines(con)
hg38.genome <- read.csv(
  textConnection(txt),
  sep = "\t",
  header = FALSE,
  row.names = 1,
)

colnames(hg38.genome) <- c("seqlengths", "bit_file")
hg38.genome.size <- data.frame(chr=row.names(hg38.genome),length=hg38.genome$seqlengths)

makeChrInt <- function(hg38.genome.size) {
  
  chr2 <- lapply(hg38.genome.size$chr,function(u){
    simple.names = gsub("^chr", "",u)
    rs <- data.frame(chrInt=simple.names,chr=u)
  })
  
  chr3 <- do.call(rbind,chr2)
  
  chr3[chr3$chrInt=="X",]$chrInt <- '23'
  
  chr3[chr3$chrInt=="Y",]$chrInt <- '24'
  
  chr3[chr3$chrInt=="M",]$chrInt <- '25'
  
  chr3.p1 <- chr3[which(!is.na(as.integer(chr3$chrInt))),]
  
  chr3.p1$chrInt <- as.integer(chr3.p1$chrInt)
  
  chr3.p2 <- chr3[which(is.na(as.integer(chr3$chrInt))),]
  
  chr3.p2$chrInt <- seq(max(chr3.p1$chrInt)+1,max(chr3.p1$chrInt)+dim(chr3.p2)[1],1)

  chr4 <- rbind(chr3.p1,chr3.p2)
  
  chr5  <- chr4[order(chr4$chrInt),]
  
  chr5
  
}

oldw <- getOption("warn")

options(warn = -1)

hg38.genome.num <- makeChrInt(hg38.genome.size)

options(warn = oldw)

hg38.genome.size.sorted <- hg38.genome.size[match(hg38.genome.num$chr,hg38.genome.size$chr),]

write.table(hg38.genome.size.sorted,file = outputFile1, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

write.table(hg38.genome.num,file = outputFile2, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)

}

docker run --rm -v /local_scratch/ISseqOutput/utilsRefData/vector:/in aiminy/isseq:2.6 bwa-mem2 index /in/pCDY.MND.GFP-KAN-BB-Aldevron-5LTR-3SIN-LTR.fa


total_vector_host_sam <- readRDS("/local_scratch/ISseqOutput1/total_vector_host_sam.rds")

total_vector_host_sam_test7 <- readRDS("/local_scratch/ISseqOutput1/test7_fastq/total_vector_host_sam.rds")


test6_total_vector_host_sam 
total_vector_host_sam_test7


boxplot(data.frame(Test7_VectorMapping=total_vector_host_sam_test7[,7],
           Test7_HostMapping=total_vector_host_sam_test7[,8],
           Test7_None=total_vector_host_sam_test7[,9]))
        
        
boxplot(data.frame(Test6_VectorMapping=test6_total_vector_host_sam[,7],
        Test6_HostMapping=test6_total_vector_host_sam[,8],
        Test6_None=test6_total_vector_host_sam[,9]))


aws s3 sync /local_scratch/ISseqOutput1/test7Read/ s3://sana-hsc-gene-therapy/ISA/ISseqOutput/Luca/Aimin/test7Read/
  
aws s3 sync /local_scratch/ISseqOutput1/test6Read/ s3://sana-hsc-gene-therapy/ISA/ISseqOutput/Luca/Aimin/test6Read/


  