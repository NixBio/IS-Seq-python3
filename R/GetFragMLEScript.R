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
}

cmd <- paste0('nohup docker run --rm -v ',paste0(input.dir,':/in')," --rm -v ",paste0(output.dir,":/out "),"aiminy/isseq:2.1 Rscript /usr/src/IS-Seq-python3/R/GetFragMLE.R", paste0("/in/","CL6_IS0_multihit_allSites.rds"),CL6_IS0,paste0("/out/",CL6_IS0)," > ~/Aimin/log/logGetMultiHitsFragMLE.txt 2>&1 &")
                                                                                                                                               cat(cmd,"\n")
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                                                               
                                                                                                        
                                                                                                                                               
                                                                                                                                               
                                                                                                                                                                                     
  
  