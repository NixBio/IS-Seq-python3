#!/usr/bin/env Rscript

packages <- c("optparse","stringr","GenomicRanges")

zzz<-lapply(packages, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_sam=args[1]  	
  input_repeatMasker=args[2]
  input_MAPQ=args[3]
  output.dir=args[4]
  ISA_run=args[5]
}

print(input_sam)
print(input_repeatMasker)
print(input_MAPQ)
print(output.dir)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

output.bam <- file.path(output.dir,'BAM')

if(!dir.exists(output.bam)){dir.create(output.bam,recursive = TRUE)}

sampleParser <- function(input_sam) {
  x <- basename(input_sam)  
  s0 <- unlist(lapply(str_locate_all(x,"aligned"),"[",c(1)))-2
  sampleName <- str_sub(x,1,s0)
  sampleName
}

sampleName <- sampleParser(input_sam)

output.file.0 <- file.path(output.bam,paste0(sampleName,'_aligned_mem.bam'))
  
t.input <- file.info(input_sam)$mtime
t.output.0 <- file.info(output.file.0)$mtime

if(is.na(t.output.0)|(t.input>t.output.0)){
cmd= paste0('samtools view -bS ',input_sam,' > ',output.file.0)  
cat(cmd,'\n')
system(cmd)
}

output.file.1<- file.path(output.bam,paste0(sampleName,'_aligned_mem_mapped.bam'))

t.output.1 <- file.info(output.file.1)$mtime

if(is.na(t.output.1)|(t.output.0>t.output.1)){
  cmd= paste0('samtools view -b -F 4 ',output.file.0,' > ',output.file.1)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.2<- file.path(output.bam,paste0(sampleName,'_aligned_mem_mapped_primary.bam'))

t.output.2 <- file.info(output.file.2)$mtime

if(is.na(t.output.2)|(t.output.1>t.output.2)){
  cmd= paste0('samtools view -b -F 256 ',output.file.1,' > ',output.file.2)  
  cat(cmd,'\n')
  system(cmd)
}

output.bam.1 <- file.path(output.dir,'BAMSorted')

if(!dir.exists(output.bam.1)){dir.create(output.bam.1,recursive = TRUE)}

output.file.3<- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_mapped_primary_sort.bam'))

t.output.3 <- file.info(output.file.3)$mtime

if(is.na(t.output.3)|(t.output.2>t.output.3)){
  cmd= paste0('samtools sort ',output.file.2,' -o ',output.file.3)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.4<- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_mapped_primary_sort.bam.index'))

t.output.4 <- file.info(output.file.4)$mtime

if(is.na(t.output.4)|(t.output.3>t.output.4)){
  cmd= paste0('samtools index ',output.file.3)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.5 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_sort_inMask.bam'))

t.output.5 <- file.info(output.file.5)$mtime

if(is.na(t.output.5)|(t.output.3>t.output.5)){
  cmd= paste0('bedtools intersect -abam ',output.file.3,' -b ',input_repeatMasker,' > ',output.file.5)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.6 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_sort_nonMask.bam'))

t.output.6 <- file.info(output.file.6)$mtime

if(is.na(t.output.6)|(t.output.3>t.output.6)){
  cmd= paste0('bedtools intersect -v -abam ',output.file.3,' -b ',input_repeatMasker,' > ',output.file.6)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.7 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_sort_inMask_qual.bam'))

t.output.7 <- file.info(output.file.7)$mtime

if(is.na(t.output.7)|(t.output.3>t.output.7)){
  cmd= paste0('samtools view -bq ',input_MAPQ,' ',output.file.5,' > ',output.file.7)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.8 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter.bam'))

t.output.8 <- file.info(output.file.8)$mtime

if(is.na(t.output.8)|(t.output.7>t.output.8)){
  cmd= paste0('samtools merge ',output.file.8,' ',output.file.6,' ',output.file.7)  
  cat(cmd,'\n')
  system(cmd)
}

output.file.9 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead.bam'))

t.output.9 <- file.info(output.file.9)$mtime

if(is.na(t.output.9)|(t.output.8>t.output.9)){
  cmd= paste0('PicardCommandLine AddOrReplaceReadGroups I=',output.file.8,' O=',output.file.9,' SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU=shoot RGSM=DePristo')
  cat(cmd,'\n')
  system(cmd)
}

output.file.10 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead.bam.index'))

t.output.10 <- file.info(output.file.10)$mtime

if(is.na(t.output.10)|(t.output.9>t.output.10)){
  cmd= paste0('samtools index ',output.file.9)
  cat(cmd,'\n')
  system(cmd)
}

output.file.11 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt.bam'))

t.output.11 <- file.info(output.file.11)$mtime

if(is.na(t.output.11)|(t.output.9>t.output.11)){
  cmd= paste0('python /usr/src/IS-Seq-python3/utils/pysam_parse.py ',output.file.9)
  cat(cmd,'\n')
  system(cmd)
}

output.file.12 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_unmapped.bam'))

t.output.12 <- file.info(output.file.12)$mtime

if(is.na(t.output.12)|(t.output.11>t.output.12)){
  cmd= paste0('samtools view -b -f 4 ',output.file.11,' > ',output.file.12)
  cat(cmd,'\n')
  system(cmd)
}

output.file.13 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_supplementary.bam'))

t.output.13 <- file.info(output.file.13)$mtime

if(is.na(t.output.13)|(t.output.11>t.output.13)){
  cmd= paste0('samtools view -b -f 2048 ',output.file.11,' > ',output.file.13)
  cat(cmd,'\n')
  system(cmd)
}

output.file.14 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam'))

t.output.14 <- file.info(output.file.14)$mtime

if(is.na(t.output.14)|(t.output.11>t.output.14)){
  cmd= paste0('samtools view -b -F 2048 ',output.file.11,' > ',output.file.14)
  cat(cmd,'\n')
  system(cmd)
}

output.file.15 <- file.path(output.bam.1,paste0(sampleName,'_aligned_mem_allFilter_rehead_exact3nt_nonSupplementary.bam.index'))

t.output.15 <- file.info(output.file.15)$mtime

if(is.na(t.output.15)|(t.output.14>t.output.15)){
  cmd= paste0('samtools index ',output.file.14)
  cat(cmd,'\n')
  system(cmd)
}

output.file.16 <- file.path(output.dir,paste0(ISA_run,'_FB-P5-Rd1-LTR.9_FB-P7-Rd2-LC.9_final_parse_filterNo.txt'))

t.output.16 <- file.info(output.file.16)$mtime

if(is.na(t.output.16)|(t.output.14>t.output.16)){
  cmd= paste0('python /usr/src/IS-Seq-python3/utils/try_pysam.py ',output.file.14,' ',ISA_run)
  cat(cmd,'\n')
  system(cmd)
}
