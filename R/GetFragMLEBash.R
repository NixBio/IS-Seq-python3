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
  output.dir=args[2]
  IS.type=args[3]
}

# Example: Rscript /home/aimin.yan/Aimin/IS-Seq-python3/R/GetFragMLEBash.R /local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa ~/Aimin/ISseqOutput/vcn/Uniq Uniq

input.files <- list.files(input.dir,pattern = "*allSites.rds",recursive = TRUE,full.names=TRUE)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

null <- lapply(input.files, function(u,IS.type){
  
  
  
  #u <- '/home/aimin.yan/Aimin/ISseqOutput/MOI50CLH6_IS95_multihit_allSites.rds'
  
  #input.dir <- '/home/aimin.yan/Aimin/ISseqOutput'
  #output.dir <- '/home/aimin.yan/Aimin/ISseqOutput'
  
   cat(u,"\n")
  
  if(IS.type=="Uniq"){
    #u <- '/local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa/MOI50CLH6/IS95/allSites.rds'
    x <- basename(u)
    ss1 <- basename(dirname(u))
    ss2 <- basename(dirname(dirname(u)))
    
    #ss1 <- unique(as.character(sapply(strsplit(x,"_"), `[`, 1)))
    #ss2 <- unique(as.character(sapply(strsplit(x,"_"), `[`, 2)))
    y <- paste0(ss2,'_',ss1)
    
    cat(dirname(u),'\n')
    cmd <- paste0('nohup docker run --rm -v ',paste0(dirname(u),':/in')," --rm -v ",paste0(output.dir,":/out "),"aiminy/isseq:2.1 Rscript /usr/src/IS-Seq-python3/R/GetFragMLE.R ", paste0("/in/",x),' ',y,' ',paste0("/out/",y)," > ~/Aimin/log/logGetMultiHitsFragMLE.txt 2>&1 &")
    cat(cmd,"\n")
    system(cmd)
    
  } else{
    
    x <- basename(u)
    ss1 <- unique(as.character(sapply(strsplit(x,"_"), `[`, 1)))
    ss2 <- unique(as.character(sapply(strsplit(x,"_"), `[`, 2)))
    y <- paste0(ss1,'_',ss2)
    
    cmd <- paste0('nohup docker run --rm -v ',paste0(input.dir,':/in')," --rm -v ",paste0(output.dir,":/out "),"aiminy/isseq:2.1 Rscript /usr/src/IS-Seq-python3/R/GetFragMLE.R ", paste0("/in/",x),' ',y,' ',paste0("/out/",y)," > ~/Aimin/log/logGetMultiHitsFragMLE.txt 2>&1 &")
    
    cat(cmd,"\n")
    system(cmd)
  }
  
  
  
  
  
},IS.type)

#temp <- strsplit("MOI30CLB6_IS0_multihit_allSites.rds",split="_")




