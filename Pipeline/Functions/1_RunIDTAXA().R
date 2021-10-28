# you can specify what clustering level you want to assign to


RunIdtaxa<-function(Type, trainingSet, TableToMergeTo=LuluOutput, SeqsToAssign=SeqsToAssign, threshold=0) {
    if (Type == "Load" | Type == "Create") {
            #SeqsToAssign should modify which seqs are selected here
            if (SeqsToAssign=="ESVs"){
                dna <- Biostrings::DNAStringSet(LuluOutput[,2])
            } else if (SeqsToAssign == "OTUs"){
                dna <- Biostrings::DNAStringSet(LuluOutput[LuluOutput$OTUrepresentativeSequence==TRUE,2])
            } else if (SeqsToAssign== "cOTUs") {
                dna <- Biostrings::DNAStringSet(LuluOutput[LuluOutput$CuratedOTURepresentativeSequence==TRUE,2])
            } else {print ("RunIdTaxa does not recognise SeqsToAssign argument") }


            ids <- IdTaxa(dna,
                        trainingSet,
                        type="extended",
                        strand="top",
                        threshold=threshold,
                        processors=4)


            #convert output to neat dataframe dependent on whether trainingset includes rank information (derived from a taxid file)
            # or not, if it doesn't then currently rank assignments are assumed (not ideal).
            
            desiredRanks<-c("domain", "kingdom", "phylum", "class", "order", "family", "genus")

            if (is.null(trainingSet$rank)) {

                IDtaxa_df<-bind_rows(lapply(ids,unlist), .id = "ESV")
                IDtaxa_df<-cbind(IDtaxa_df$ESV, select(IDtaxa_df, starts_with( "taxon")), select(IDtaxa_df, starts_with( "confidence")))
                names(IDtaxa_df)[1]<-"ESV"
                names(IDtaxa_df)[2:  (length(desiredRanks)+1)]<-desiredRanks
                names(IDtaxa_df)[ (length(desiredRanks)+2) : (length(desiredRanks)*2+1) ]<-paste0("confidence ", desiredRanks)
                IDtaxa_df<-IDtaxa_df[ 1: (length(desiredRanks)*2+1) ]

            } else {
                
                IDtaxa_df<-data.frame(matrix(NA, ncol=2*length(desiredRanks), nrow=0))
                names(IDtaxa_df)<-c(desiredRanks, paste0(desiredRanks, "_confidence"))
                
                for (i in 1: length(ids)) {
                    selection<-ids[[i]]$rank %in% desiredRanks
                    insertions<-desiredRanks %in% ids[[i]]$rank
                    IDtaxa_df[i,1:length(desiredRanks)][insertions]<-ids[[i]]$taxon[selection]
                    IDtaxa_df[i, (length(desiredRanks)+1) : (length(desiredRanks)*2) ] [insertions] <-ids[[i]]$confidence[selection]
                }
                IDtaxa_df$ESV<-rownames(IDtaxa_df)

            }

            #merge with seq data table
                TableToMergeTo$ESV<-rownames(TableToMergeTo)
                SeqDataTable<-merge(TableToMergeTo, IDtaxa_df, by= "ESV", all=TRUE)
                SeqDataTable$ESV<-paste0("ESV_", SeqDataTable$ESV)


            return(list(SeqDataTable=SeqDataTable, IDTAXAplotdata=(ids), trainingSet=trainingSet))
    } else {
        print("No assignment performed")
        return(NULL)
    }
}
