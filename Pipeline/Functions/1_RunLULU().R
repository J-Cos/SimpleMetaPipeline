RunLULU<-function(TableToMergeTo, MatchRate, MinRelativeCo, RatioType) {

        #run lulu
            curated_result<-lulu(OtuTableForLulu, MatchListForLulu, minimum_ratio_type =RatioType, minimum_match =MatchRate , minimum_relative_cooccurence = MinRelativeCo)

        #combine output into seq datatable
            curated_OTU_map<-curated_result$otu_map
            curated_OTU_map$OTU<-rownames(curated_OTU_map)
            curated_OTU_map<-curated_OTU_map[,c(3,4,6)]
            names(curated_OTU_map)<-c("curatedOTU", "CuratedOTURepresentativeSequence", "OTU")
            curated_OTU_map$CuratedOTURepresentativeSequence<-curated_OTU_map$CuratedOTURepresentativeSequence=="parent"

           #merge with seq data table
            SeqDataTable<-merge(TableToMergeTo, curated_OTU_map, by= "OTU", all=TRUE)
            SeqDataTable$CuratedOTURepresentativeSequence<-SeqDataTable$CuratedOTURepresentativeSequence==TRUE & SeqDataTable$OTUrepresentativeSequence==TRUE

        return(SeqDataTable)


}