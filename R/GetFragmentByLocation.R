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
  fragment.type=args[2]  	
  fragment.length=args[3]
  fragment.strand=args[4]
  fragment.start=args[5]
  fragment.name=args[6]
  read.name=args[7]
  read_out_prefix=args[8]
  output.dir=args[9]
  output.fasta=args[10]
  number.of.read=args[11]
  output.sam=args[12]
  #output.file=args[7]
}

# Example:

# Rscript ~/ispipe/R/GetRandomFragment.R /home/ubuntu/DEMO/IS-Seq/utilsRefData/IsSeq/hg38/hg38ChrOnly.fa "target" 3000 "negative" 49461738 "chr19_49461738" "target" "chr19_49461738_" /home/ubuntu/SHARE/Aimin/Test_3000_negative_fragment hg38Chr19OnlyUpIS3000Fragment.fa 2000 UpFrag49461738R1_R2_Barcode_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_aligned_mem.sam


if (!dir.exists(output.dir)){dir.create(output.dir, recursive = TRUE)}

if(fragment.type=="random"){
  
  #Fragment.random <- paste(DNA_BASES[sample(1:4,3000,replace=T)], collapse="")
  
  Fragment.random <- paste(DNA_BASES[sample(1:4,fragment.length,replace=T)], collapse="")
  
  Fragment.random.1 <- DNAStringSet(Fragment.random)
  
  #output.dir <- '/home/ubuntu/SHARE/Aimin/RandomFragment'
  
  
  #names(Fragment.random.1) <- 'ChrRandom'
  
  names(Fragment.random.1) <- fragment.name
  
  Fragment.out <- Fragment.random.1
  
  #writeXStringSet(Fragment.random.1, file.path(output.dir,'RandomFragment_1.fa'))
  
  
  
  
  #out <- file.path(output.dir,'Random_1')
  
  #output.file <- file.path(output.dir,'RandomFragment_1.fa')
  
  #cmd= paste0('art_illumina -ss MSv3 -p -i ',output.file,' -l 250 -c 20000000 -m 1000 -s 300 -d "Random" -o ',out)
  
  #cat(cmd)
  
}else{

  #host.fa <- '/home/ubuntu/SHARE/D32_Platform_Development/MANUSCRIPTS_DATA/ISAtest/MiSeqTest/utilsRefData/hg38/hg38ChrOnly.fa'
  
  host.fa <-'/local_scratch/ISseqOutput/utilsRefData/hg38/hg38ChrOnly.fa'
  
  rfa.host <- readFasta(host.fa)
  
  hg38.df <- cbind(as.data.frame(rfa.host@id),as.data.frame(width(rfa.host@sread)))
  
  colnames(hg38.df) <- c('chr','length')
  
  num.of.IS <- 1000
  hg38.df.is <- data.frame(hg38.df,num_is=round((num.of.IS/sum(hg38.df$length))*hg38.df$length,0))
  
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
  
  
  rfa.h <- import_sq(sread(rfa.host))
  
  rfa.host.chr19 <- rfa.h[19,]$sq
  
  if(fragment.strand=='negative'){
    
    
    cat(fragment.start,'\n')
    cat(fragment.length,'\n')
    
    fragment.end = as.integer(fragment.start)-as.integer(fragment.length)
    
    IS.up.3000.fragment <- bite(rfa.host.chr19,  fragment.end:fragment.start)
    IS.up.3000.fragment.DNA <- export_sq(IS.up.3000.fragment, "Biostrings::DNAStringSet")
    IS.up.3000.fragment.DNA.sq <- import_sq(IS.up.3000.fragment.DNA)
    IS.up.3000.fragment.DNA.sq.reverse.complement <- complement(reverse(IS.up.3000.fragment.DNA.sq$sq))
    temp.rfa.h <- rfa.h
    temp.rfa.h[19,]$sq <- IS.up.3000.fragment.DNA.sq.reverse.complement
    
  }
  
  if(fragment.strand=='positive'){
    
    fragment.end = as.integer(fragment.start)+as.integer(fragment.length)
    
    IS.up.3000.fragment <- bite(rfa.host.chr19,  fragment.start:fragment.end)
    print(IS.up.3000.fragment)
    IS.up.3000.fragment.DNA <- export_sq(IS.up.3000.fragment, "Biostrings::DNAStringSet")
    print(IS.up.3000.fragment.DNA)
    IS.up.3000.fragment.DNA.sq <- import_sq(IS.up.3000.fragment.DNA)
    print(IS.up.3000.fragment.DNA.sq)
    temp.rfa.h <- rfa.h
    print(temp.rfa.h)
    temp.rfa.h[19,]$sq <-  IS.up.3000.fragment.DNA.sq$sq
    print(temp.rfa.h[19,]$sq)
    
  }
  
  rfa.h.DNA <- export_sq(temp.rfa.h$sq, "Biostrings::DNAStringSet")
  
  rfa.host@sread <- rfa.h.DNA
  rfa.host@id[19] <- fragment.name
  
  Fragment.out <- DNAStringSet(sread(rfa.host)[19])
  names(Fragment.out) <- fragment.name
  
  #output.file=file.path(out.dir.name,output.fasta)
  #writeFasta(rfa.host[19], file=output.file, mode="w")
  
}

writeXStringSet(Fragment.out, file.path(output.dir,output.fasta))

#out <- paste0(output.dir,'Random')

out <- file.path(output.dir,read_out_prefix)

#output.file <- file.path(output.dir,'RandomFragment.fa')

output.file <- file.path(output.dir,output.fasta)

#cmd= paste0('art_illumina -ss MSv3 -p -i ',output.file,' -l 250 -c 17274461 -m 1000 -s 300 -d "Random" -o ',out)

number.of.read <- as.integer(number.of.read)

cmd= paste0('art_illumina -ss MSv3 -p -i ',output.file,' -l 250 -c ',number.of.read,' -m 1000 -s 300 -d ',read.name,' -o ',out)

#read.name

cat(cmd,'\n')

system(cmd)

if(fragment.type!="random"){
  
output.sam <- file.path(output.dir,output.sam)

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







