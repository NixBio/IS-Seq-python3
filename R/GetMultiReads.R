# Rscript $HOME/Aimin/IS-Seq-python3/R/GetMultiReads.R /local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa/MOI30CLB7/IS0/multihitData.rds /local_scratch/ISseqOutput/vcn/IsaByINSPIIRED/fa/MOI30CLB7/IS0

#/Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE/Results.RData /Users/c-aimin.yan/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0/FragMLE_multihitData/Results.RData ~/OneDrive/Aimin/Vcn_INSPIIRED/share/ISseqOutput/Dec282021/IsaByINSPIIREDTimeTest/fa/CL6/IS0


args = commandArgs(trailingOnly=TRUE)

#print(args)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  input.multi.hit=args[1]
  output.dir=args[2]
}

#print(input.multi.hit)

multihit <- readRDS(input.multi.hit)
saveRDS(multihit$unclusteredMultihits,file.path(output.dir,"multihit_allSites.rds"))


