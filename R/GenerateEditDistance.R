# Examples:

# Rscript /home/home/aimin.yan/Aimin/IS-Seq-python3/R/GenerateEditDistance.R /home/local_scratch/ISseqOutput/230221M086400005000000000KWN5YRerun

initial.options <- commandArgs(trailingOnly = FALSE)
print(initial.options)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
print(script.dirname)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
  input.dir=args[1]
}

input.files <- list.files(input.dir,pattern = "exact3nt_nonSupplementary.bam$",recursive = TRUE,full.names=TRUE)

null <- lapply(input.files, function(u){
  cmd <- paste0('python /home/home/aimin.yan/Aimin/IS-Seq-python3/CompareEditDistanceDefinitions.py ',u)
  cat(cmd,'\n')
  system(cmd)
})


input.ed.files <- list.files('/local_scratch/ISseqOutput/230221M086400005000000000KWN5YRerun',pattern = "*_edit_distance.csv$",recursive = TRUE,full.names=TRUE)

res <- lapply(input.ed.files, function(u){
  
  x <- read.csv(u,header = T)
  x
})
  

input.ed.files <- list.files('/local_scratch/ISseqOutput/230221M086400005000000000KWN5YRerun',pattern = "*_edit_distance_1.csv$",recursive = TRUE,full.names=TRUE)

res <- lapply(input.ed.files, function(u){
  
  x <- read.csv(u,header = T)
  x <- x[,-1]
  x
})

names(res) <-  basename(input.ed.files)

