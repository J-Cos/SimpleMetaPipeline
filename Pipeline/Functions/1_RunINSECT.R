#RunINSECT(  SeqDataTable=LuluOutput2,
#            classifier="INSECT_COIclassifier.rds")
            


RunINSECT<-function(SeqDataTable, classifier, cores=1) {

    seqs<-SeqDataTable$Sequence
    tree <-readRDS(file.path(path, "Data", classifier))
    longdf<-classify(seqs,tree,cores=cores)

    SeqDataTable<-cbind(SeqDataTable, longdf[c(3:4, 6:11)])
    return(SeqDataTable)
}

#write.csv( SeqDataTable,"SeqDataTable.csv")
