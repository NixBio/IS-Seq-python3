GetData <- function(input.dir) {
  
  Is.out <- list.files(path=input.dir, pattern="*CollisionClean_CollisionTable.txt", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
  
  Is.df <- lapply(Is.out,function(u){
    z <- read.table(u,header = T)
    z <- z[order(z$chr,z$pos),]
    row.names(z) <- seq(1,dim(z)[1],1)
    z
  })
  
  Is.df
}

input.dir <- "/local_scratch/ISseqOutput/Simulation100allChr"

Is.df <- GetData(input.dir)

Is.df.old <- GetData('/local_scratch/simulation100IS')

identical(Is.df[[1]][,c(1:4,6)],Is.df.old[[1]][,c(1:4,6)])

