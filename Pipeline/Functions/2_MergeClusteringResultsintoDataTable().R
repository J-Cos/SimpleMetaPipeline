MergeClusteringResultsintoDataTable<-function(TableToMergeTo, OTULabeledSeqs) {
    # combine OTU labelled seqs with ESVtable
    # merge by sequence  this is all that is required when swarm performed on ESVs (dada only)
    SeqDataTable<-merge(TableToMergeTo, OTULabeledSeqs, by="Sequence", all=TRUE)
        
    # but when performed on cESVs (dada then lulu1) (i.e. when some ESVs have been removed before swarm)
    # we need to then back fill the correct otu assignments for each ESV depending on which cESV ESVs have been grouped into. 
    removedESVs<-which(is.na(SeqDataTable["OTU"]))
    if (length(removedESVs)!=0) {
        for (i in 1: length(removedESVs) ) {
            parentESVrow<-which(SeqDataTable$ESV == SeqDataTable[removedESVs[i],]$curatedESV) 
            SeqDataTable[removedESVs[i],]["OTU"]<- SeqDataTable [ parentESVrow ,]$OTU
            SeqDataTable[removedESVs[i],]["OTUrepresentativeSequence"]<- FALSE
        }
    }

    return(SeqDataTable)

}