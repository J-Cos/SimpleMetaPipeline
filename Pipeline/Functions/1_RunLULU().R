RunLULU<-function(TableToMergeTo, MatchRate, MinRelativeCo, RatioType, clustering) {

    # run lulu
    curated_result<-lulu(OtuTableForLulu, MatchListForLulu, minimum_ratio_type =RatioType, minimum_match =MatchRate , minimum_relative_cooccurence = MinRelativeCo)

    # combine output into seq datatable
    curated_OTU_map<-curated_result$otu_map
    curated_OTU_map$OTU<-rownames(curated_OTU_map)
    curated_OTU_map<-curated_OTU_map[,c(3,4,6)]
    names(curated_OTU_map)<-c(paste0("curated", clustering), paste0("Curated",clustering,"RepresentativeSequence"), clustering)
    curated_OTU_map[[ paste0("Curated",clustering,"RepresentativeSequence") ]]<-curated_OTU_map[[ paste0("Curated",clustering,"RepresentativeSequence") ]]=="parent"

    # merge with seq data table
    SeqDataTable<-merge(TableToMergeTo, curated_OTU_map, by= clustering, all=TRUE)

    # label esv which are not OTU representatives as also not being cOTUrepresentatives
    if (clustering=="OTU") { 
        SeqDataTable$CuratedOTURepresentativeSequence<-SeqDataTable$CuratedOTURepresentativeSequence==TRUE & SeqDataTable$OTUrepresentativeSequence==TRUE
    }

    return(SeqDataTable)

}