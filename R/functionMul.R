GetIsRel <- function(input) {
  load(input)
  IS.grouped$REL <- IS.grouped[,4]/sum(IS.grouped[,4])
  colnames(IS.grouped)[4] <- 'Frag.MLE'
  IS.grouped.sorted <- IS.grouped[order(IS.grouped$Frag.MLE,decreasing = T),]
  IS.grouped.sorted
  IS.grouped.gr <- makeGRangesFromDataFrame(data.frame(chr=IS.grouped.sorted$seqnames,start=IS.grouped.sorted$start,end=IS.grouped.sorted$start,strand=IS.grouped.sorted$strand,Frag.MLE=IS.grouped.sorted$Frag.MLE,REL=IS.grouped.sorted$REL),keep.extra.column=T)
  IS.grouped.gr
}