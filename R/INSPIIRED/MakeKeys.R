#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input_R1=args[1]
  input_R2=args[2]
  sampleName=args[3]
  output.dir=args[4]
}

print(input_R1)
print(input_R2)
print(output.dir)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

library(data.table)

cmd1 <- paste0("cat < ",input_R1," | grep '@M'")
cmd2 <- paste0("cat < ",input_R2," | grep '@M'")

R1.readName <- data.table::fread(cmd1,header=F)
R2.readName <- data.table::fread(cmd2,header=F)

R1Name <- sapply(strsplit(as.character(R1.readName$V1), "@"), "[[", 2)
R1 <- data.frame(R1=seq(1:length(R1Name)),names=R1Name)

R2Name <- sapply(strsplit(as.character(R2.readName$V1), "@"), "[[", 2)
R2 <- data.frame(R2=seq(1:length(R2Name)),names=R2Name)

cat('R1','\n')
print(head(R1))

cat('R2','\n')
print(head(R2))

keys <- merge(R2,R1,by='names')

keys$readPairKey <- paste0(keys$R2,'_',keys$R1)

keys <- keys[,c(2,3,1,4)]

keys$names <- paste0(sampleName,'%',keys$names)

cat('R1R2','\n')
print(head(keys))

saveRDS(keys,file.path(output.dir,'keys.rds'))
