#!/usr/bin/env Rscript

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
  input.file=args[1]
  input.dir=args[2]
  output.dir=args[3]
  ref.genome=args[4]
}

print(input.file)
print(input.dir)
print(output.dir)
print(ref.genome)

if(!dir.exists(output.dir)){dir.create(output.dir,recursive = TRUE)}

input <- read.table(input.file,sep=',',header=T)


input.R1 = file.path(input.dir,'R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLcDEMULTIPLEXINGTofq',paste0('R1_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_',input$Fusion.Primer.LTR..ID,'.fq_trimwithCutAdapt'))


input.R1R2 = file.path(input.dir, paste0('CutAdapt/R1_R2_Barcode_',input$Fusion.Primer.LTR..ID,'_',input$Fusion.Primer.LC..ID,'_trimmedID'))

input.R2 = file.path(input.dir,'RandomBarcodRemovalOutPut4R2CuReRun',paste0('R2_fastq_trim12nt_qcTrimmed_MatchBlastLtrLc_Barcode_',input$Fusion.Primer.LC..ID,'.fq_trimwithCutAdapt'))


if(TRUE){

cat('Start conversion\n')

input.all <- data.frame(R1=input.R1,R1R2=input.R1R2,R2=input.R2,SampleName=gsub('-','',input$Sample.Type))

null <- lapply(1:dim(input.all)[1], function(u){
 
  
  sample.dir <- file.path(output.dir,input.all[u,]$SampleName)

  if(!dir.exists(sample.dir)){dir.create(sample.dir,recursive = TRUE)}

  cmd.R1 <- paste0('Rscript ',file.path(script.dirname,'FqToFa.R '),input.all[u,]$R1,' ',input.all[u,]$R1R2,' ',input.all[u,]$SampleName,' ', sample.dir,' ',ref.genome)

  print(cmd.R1)
  system(cmd.R1)
  
  cmd.R2 <- paste0('Rscript ',file.path(script.dirname,'FqToFa.R '),input.all[u,]$R2,' ',input.all[u,]$R1R2,' ',input.all[u,]$SampleName,' ', sample.dir,' ',ref.genome)

  print(cmd.R2)
  system(cmd.R2)
  
})

}

