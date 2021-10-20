RunLULU<-function(TableToMergeTo) {

        #run lulu
            curated_result<-lulu(OtuTableForLulu, MatchListForLulu)

        #combine output into seq datatable
            curated_OTU_map<-curated_result$otu_map
            curated_OTU_map$OTU<-rownames(curated_OTU_map)
            curated_OTU_map<-curated_OTU_map[,c(3,4,6)]
            names(curated_OTU_map)<-c("curatedOTU", "CuratedOTURepresentativeSequence", "OTU")
            curated_OTU_map$CuratedOTURepresentativeSequence<-curated_OTU_map$CuratedOTURepresentativeSequence=="parent"


            SeqDataTable<-merge(TableToMergeTo, curated_OTU_map, by= "OTU", all=TRUE)
            SeqDataTable$CuratedOTURepresentativeSequence<-SeqDataTable$CuratedOTURepresentativeSequence==TRUE & SeqDataTable$OTUrepresentativeSequence==TRUE

        return(SeqDataTable)


}