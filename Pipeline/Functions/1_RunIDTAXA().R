# you can specify what clustering level you want to assign to

#expects classifier trained with seq names in following format, number of ranks can vary: "Root; Arthropoda; Arachnida; Araneae; Uloboridae; Zosis; geniculata" 


RunIdtaxa<-function(Type, trainingSet, TableToMergeTo, SeqsToAssign=SeqsToAssign, threshold) {
    if (Type == "Assign") {
        #SeqsToAssign should modify which seqs are selected here
            if (SeqsToAssign=="ESVs"){
                dna <- Biostrings::DNAStringSet(TableToMergeTo[,2])
                names(dna)<-TableToMergeTo$ESV
            } else if (SeqsToAssign == "OTUs"){
                dna <- Biostrings::DNAStringSet(TableToMergeTo[TableToMergeTo$OTUrepresentativeSequence==TRUE,2])
                names(dna)<-TableToMergeTo$OTU

            } else if (SeqsToAssign== "cOTUs") {
                dna <- Biostrings::DNAStringSet(TableToMergeTo[TableToMergeTo$CuratedOTURepresentativeSequence==TRUE,2])
                names(dna)<-TableToMergeTo$curatedOTU

            } else {print ("RunIdTaxa does not recognise SeqsToAssign argument") }

        #load the specified trainingset
        load(file.path(path, "Data", "Classifiers", trainingSet))

        #classify
            ids <- IdTaxa(dna,
                        trainingSet,
                        type="extended",
                        strand="both",
                        threshold=threshold,
                        processors=4)

           ExtractFromIds<-function(list_item, category){
                ListItemLength<-length(list_item[[category]])
                if (ListItemLength < Nranks) {
                   list_item[[category]][(ListItemLength+1):Nranks] <-list_item[[category]][ListItemLength]
                }
                return(list_item[[category]])
            }

            Nranks<-max( unlist ( lapply ( lapply( ids, '[[', 1), length)))
            TaxonVecList<-lapply(ids, ExtractFromIds, category="taxon")
            TaxonDf<-plyr::ldply(TaxonVecList, rbind)

            ConfVecList<-lapply(ids, ExtractFromIds, category="confidence")
            ConfDf<-plyr::ldply(ConfVecList, rbind)

            IdtaxaDf<-merge(TaxonDf, ConfDf, ".id")
            names(IdtaxaDf)<-c("ESV", paste0("Rank_", 1:Nranks), paste0("Rank_", 1:Nranks, "_Confidence"))

        #merge with seq data table
            SeqDataTable<-merge(TableToMergeTo,IdtaxaDf, by= "ESV", all=TRUE)
            
            return(list(SeqDataTable=SeqDataTable, PlotData=ids))
    } else {
        print("No assignment performed")
        return(NULL)
    }
}
