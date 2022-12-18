#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

packages <- c("optparse","stringr","GenomicRanges","ShortRead","tidysq")

zzz<-lapply(packages, function(xxx) suppressMessages(library(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  host.fa=args[1]
  num.of.IS=args[2]
  read.name=args[3]
  read_out_prefix=args[4]
  output.dir=args[5]
  output.fasta=args[6]
  number.of.read=args[7]
  #output.sam=args[8]
}

# Example:

# Rscript /home/aimin.yan/Aimin/IS-Seq-python3/R/GetFragmentByLocationScript.R /home/ubuntu/DEMO/IS-Seq/utilsRefData/IsSeq/hg38/hg38ChrOnly.fa 1000 IS1000 Simulation1000IS /local_scratch/ISseqOutput/Simulation fragment.fasta 200000000000

#host.fa <-'/local_scratch/ISseqOutput/utilsRefData/hg38/hg38ChrOnly.fa'

rfa.host <- readFasta(host.fa)

hg38.df <- cbind(as.data.frame(rfa.host@id),as.data.frame(width(rfa.host@sread)))

colnames(hg38.df) <- c('chr','length')

cat(num.of.IS,"\n")

hg38.df.is <- data.frame(hg38.df,num_is=round((as.integer(num.of.IS)/sum(hg38.df$length))*hg38.df$length,0))

set.seed(1)
is <- lapply(1:dim(hg38.df.is)[1], function(u){
  
  z <- sample(seq(1,hg38.df.is[u,]$length,1),hg38.df.is[u,]$num_is)
  std <- sample(c('+','-'),length(z),replace = T)
  z1 <- data.frame(chr=rep(hg38.df.is[u,]$chr,length(z)),as.data.frame(z),std=as.data.frame(std))    
  z1
})
names(is) <- hg38.df.is$chr
is1 <- do.call(rbind,is)

positive.IS <- data.frame(chr='chr19',z=49461738,std='-')   

is1 <- rbind(is1,positive.IS)

#is1[is1$chr=='chr19'&is1$z==24958819,]$z <- 49461738
#is1[is1$chr=='chr19'&is1$z==49461738,]$std <- '-'

chr.index <- data.frame(as.data.frame(unique(is1$chr)),index=rownames(as.data.frame(unique(is1$chr))))

colnames(chr.index) <- c('chr','index')

is2 <- merge(is1,chr.index,by='chr',sort=FALSE)

colnames(is2) <- c('chr','position','std','chr.indx')

getFrag <- function(rfa.host,chr.index=19,fragment.strand='negative',fragment.start=49461738,fragment.length=3000,fragment.name='chr19_49461738_-',output.dir,output.fasta) {
  
  if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}
  
  rfa.h <- import_sq(sread(rfa.host))
  
  rfa.host.chr19 <- rfa.h[chr.index,]$sq
  
  if(fragment.strand=='-'){
    
    cat(fragment.start,'\n')
    cat(fragment.length,'\n')
    
    fragment.end = as.integer(fragment.start)-as.integer(fragment.length)
    
    IS.up.3000.fragment <- bite(rfa.host.chr19,  fragment.end:fragment.start)
    IS.up.3000.fragment.DNA <- export_sq(IS.up.3000.fragment, "Biostrings::DNAStringSet")
    IS.up.3000.fragment.DNA.sq <- import_sq(IS.up.3000.fragment.DNA)
    IS.up.3000.fragment.DNA.sq.reverse.complement <- complement(reverse(IS.up.3000.fragment.DNA.sq$sq))
    temp.rfa.h <- rfa.h
    temp.rfa.h[chr.index,]$sq <- IS.up.3000.fragment.DNA.sq.reverse.complement
    
  }
  
  if(fragment.strand=='+'){
    
    fragment.end = as.integer(fragment.start)+as.integer(fragment.length)
    
    IS.up.3000.fragment <- bite(rfa.host.chr19,  fragment.start:fragment.end)
    print(IS.up.3000.fragment)
    IS.up.3000.fragment.DNA <- export_sq(IS.up.3000.fragment, "Biostrings::DNAStringSet")
    print(IS.up.3000.fragment.DNA)
    IS.up.3000.fragment.DNA.sq <- import_sq(IS.up.3000.fragment.DNA)
    print(IS.up.3000.fragment.DNA.sq)
    temp.rfa.h <- rfa.h
    print(temp.rfa.h)
    temp.rfa.h[chr.index,]$sq <-  IS.up.3000.fragment.DNA.sq$sq
    print(temp.rfa.h[chr.index,]$sq)
    
  }
  
  rfa.h.DNA <- export_sq(temp.rfa.h$sq, "Biostrings::DNAStringSet")
  
  rfa.host@sread <- rfa.h.DNA
  rfa.host@id[chr.index] <- fragment.name
  
  Fragment.out <- DNAStringSet(sread(rfa.host)[chr.index])
  names(Fragment.out) <- fragment.name
  
  #output.file=file.path(out.dir.name,output.fasta)
  #writeFasta(rfa.host[19], file=output.file, mode="w")
  
  writeXStringSet(Fragment.out, file.path(output.dir,output.fasta),append=TRUE)
  
}

#output.dir <- '~/SHARE/Aimin/Simulation100IS'
#output.fasta <- 'fragment.fasta'

if(!file.exists(file.path(output.dir,output.fasta))){
  
  null <- lapply(1:dim(is2)[1],function(u){
    
    #u <- 1
    
    chr.index <- as.integer(is2[u,]$chr.indx)
    fragment.strand <- is2[u,]$std
    fragment.start <- is2[u,]$position
    fragment.name <- paste0(is2[u,]$chr,'_',fragment.start,'_',fragment.strand)
    
    getFrag(rfa.host,chr.index=chr.index,fragment.strand=fragment.strand,fragment.start=fragment.start,fragment.length=3000,fragment.name=fragment.name,output.dir,output.fasta)
    
  })
  
}

#/out/Simulation1000IS

#read_out_prefix <- 'Simulation100IS'

out <- file.path(output.dir,read_out_prefix)

temp <- paste0(out,'1.fq')

if(!file.exists(temp)){
  
  output.file <- file.path(output.dir,output.fasta)
  
  number.of.read <- as.integer(number.of.read)
  
  cat(number.of.read,'\n')
  
  cmd= paste0('art_illumina -ss MSv3 -p -i ',output.file,' -l 250 -c ',number.of.read,' -m 1000 -s 300 -d ',read.name,' -o ',out,' -na')

  cat(cmd,'\n')
  system(cmd)
  
}

#if(fragment.type!="random"){

output.sam <- file.path(output.dir,output.sam)

  if(!file.exists(temp)){
  
  cmd1=paste0('bwa-mem2 mem -t 8 ',host.fa,' ',paste0(out,'1.fq'),' ',paste0(out,'2.fq'),' > ',output.sam)

  cat(cmd1,'\n')

  system(cmd1)
  #/home/ubuntu/DEMO/IS-Seq/utilsRefData/IsSeq/hg38/hg38ChrOnly.fa

  cmd2=paste0('Rscript ',file.path(script.dirname,"sam2filterNo.R"),' ',output.sam,' ',file.path(dirname(host.fa),'repeatMaskerhg38BED'),' ',0,' ',output.dir,' POOL-ISA-AVRO-6-Preclin')

  cat(cmd2,'\n')
  system(cmd2)

}

#Rscript ~/ispipe/R/sam2filterNo.R /home/ubuntu/SHARE/ISseqOutput/Dec282021/CutAdapt/align/R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/repeatMaskerhg38BED 0 /home/ubuntu/SHARE/Aimin/TestSimulation/CL6 POOL-ISA-AVRO-6-Preclin



#/home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/hg38ChrOnly.fa /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/simulationUp_3000_49461738_Frag1.fq /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/simulationUp_3000_49461738_Frag2.fq > /home/ubuntu/SHARE/Aimin/TestSimulation/UpFrag49461738R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam


#bwa mem -t 8 /home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/hg38ChrOnly.fa /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/simulationUp_3000_49461738_Frag1.fq /home/ubuntu/SHARE/Aimin/INSPIIRED_test_output/simulationUp_3000_49461738_Frag2.fq > /home/ubuntu/SHARE/Aimin/TestSimulation/UpFrag49461738R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam



#read_out_prefix <- 'Simulation100IS200kRead'
#out <- file.path(output.dir,read_out_prefix)

#number.of.read <- 200000
#read.name <- 'IS100'

#cmd= paste0('art_illumina -ss MSv3 -p -i ',output.file,' -l 250 -c ',as.integer(number.of.read),' -m 1000 -s 300 -d ',read.name,' -o ',out)

#cat(cmd,'\n')

