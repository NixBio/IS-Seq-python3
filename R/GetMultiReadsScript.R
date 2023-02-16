# Rscript $HOME/IS-Seq/IS-Seq-python3/R/GetMultiReadsScript.R /local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

packages <- c("optparse","stringr")

zzz<-lapply(packages, function(xxx) suppressMessages(library(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

args = commandArgs(trailingOnly=TRUE)

#print(args[1])

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  #input_uniq_hit=args[1]
  input.dir=args[1]
  output.dir=args[2]
}

if(!dir.exists(output.dir)){
  dir.create(output.dir,recursive = TRUE)
}

input.files <- list.files(input.dir,pattern = "multihitData.rds",recursive = TRUE,full.names=TRUE)

null <- lapply(input.files, function(u){
  cmd=paste0("Rscript ",file.path(script.dirname,"GetMultiReads.R ") ,u,' ',output.dir)
  cat(cmd,'\n')
  system(cmd)
})
