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
                        strand="top",
                        threshold=threshold,
                        processors=4)

            #convert output to neat dataframe
                IDtaxa_df<-bind_rows(lapply(ids,unlist), .id = "ESV")
                names(IDtaxa_df)[2:8]<-c("rootrank", "domain", "phylum", "class", "order", "family", "genus")
                names(IDtaxa_df)[9:15]<-paste0("confidence ", c("rootrank", "domain", "phylum", "class", "order", "family", "genus"))
                IDtaxa_df<-IDtaxa_df[-c(16:22)]

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
