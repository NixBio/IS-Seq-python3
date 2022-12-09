#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

packages <- c("optparse","stringr")

zzz<-lapply(packages, function(xxx) suppressMessages(library(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input.dir=args[1]
  #output.dir=args[2]
  #IS.type=args[3]
}

# Example: Rscript /home/aimin.yan/Aimin/IS-Seq-python3/R/CombineUniqMultiBash.R /home/aimin.yan/Aimin/ISseqOutput

#input.dir <- "~/Aimin/ISseqOutput"

input.files <- list.files(input.dir,pattern = "Results.RData",recursive = TRUE,full.names=TRUE)

df1 <- data.frame(sampleName=basename(dirname(input.files[grep('Uniq',input.files)])),UniqfileName=input.files[grep('Uniq',input.files)])

df2 <- data.frame(sampleName=basename(dirname(input.files[-grep('Uniq',input.files)])),MultiName=input.files[-grep('Uniq',input.files)])

input.data <- merge(df1,df2,by="sampleName")

input.files[-grep('Uniq',input.files)]

#if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

null <- lapply(1:dim(input.data)[1], function(u){
  
  sampleName <- input.data[u,1]
  UniqfileName <- input.data[u,2]
  MultiName <- input.data[u,3]
  
  cat(sampleName,"\t",UniqfileName,"\t",MultiName,"\n")
  
  cmd <- paste0(paste0("nohup docker run --rm -v ",paste0(UniqfileName,":/in1")," --rm -v ",paste0(MultiName,":/in2")," --rm -v ",paste0(input.dir,"/vcn/UniqAndMulti/",sampleName,"/:/out"),paste0(" aiminy/isseq:2.1 Rscript /usr/src/IS-Seq-python3/R/ScriptCombineUniqMulti_hits.R /in1 /in2 /out"),paste0(" > ~/Aimin/log/logCombineUniqWithMulti.txt 2>&1 &")))
  
  cat(cmd,"\n")
  
  system(cmd)
  
})

#temp <- strsplit("MOI30CLB6_IS0_multihit_allSites.rds",split="_")




