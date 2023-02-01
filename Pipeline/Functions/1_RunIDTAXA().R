# you can specify what clustering level you want to assign to

#expects classifier trained with seq names in following format, number of ranks can vary: "Root; Arthropoda; Arachnida; Araneae; Uloboridae; Zosis; geniculata" 

RunIdtaxa<-function(Type, trainingSet, TableToMergeTo, SeqsToAssign=SeqsToAssign, threshold, QuerySequenceChunkSize=500, parallel=FALSE, multithread=multithread) {
    
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

        # splitquery sequences into chunks
            dna_list<-split(dna, ceiling(seq_along(dna)/QuerySequenceChunkSize))
            print(paste0("Query sequences split into ", length(dna_list), " chunks of ", QuerySequenceChunkSize, " sequences"))
        
        if (multithread==TRUE){
            mc.cores<-250
        } else {
            mc.cores<-multithread
        }


        if (parallel){
            # parallel
                print(paste0("Starting parallel classification with ", as.integer(mc.cores), " cores"))

                ids<-parallel::mclapply(dna_list, mc.cores=mc.cores, function(dnachunk){
                                                        IdTaxa(dnachunk,
                                                                trainingSet,
                                                                type="extended",
                                                                strand="both",
                                                                threshold=threshold,
                                                                processors=1)
                                                        }
                                    )
        } else {
            #linear
                #first check if any chunks completed in a previous exited run, if so load these and append  additional chunks to this object
                    if(file.exists(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_ids.RDS")))) {
                        ids<-readRDS(file=file.path(path, "IntermediateOutputs", paste0(dataname,"_ids.RDS") ))
                    } else {
                        ids<-list()
                    }
                
                print(paste0("Starting linear classification at chunk ", length(ids)+1))

                #loop over incomplete chunks
                for (chunkNumber in (length(ids)+1):length(dna_list)) {
                    #classify
                    ids[[chunkNumber]] <- IdTaxa(dna_list[[chunkNumber]],
                                            trainingSet,
                                            type="extended",
                                            strand="both",
                                            threshold=threshold,
                                            processors=NULL)
                    print(paste0("Chunk ", chunkNumber, " complete"))
                    CacheOutput(ids)
                }
        }
                
        print(paste0("All chunks complete"))
        ids<-purrr::reduce(ids, c)
        print(paste0("All chunks merged"))

        #if rank in ids
            if ('rank' %in% colnames(as.data.frame(ids[1][[1]]))) {

                ReformatIds<-   function (id) {
                                    x<-c(id$taxon, id$confidence)
                                    names(x)<-c(id$rank, paste0(id$rank, "_confidence"))
                                    x<-x[c(desiredranks, paste0(desiredranks, "_confidence"))]
                                    names(x)<-c(desiredranks, paste0(desiredranks, "_confidence"))
                                    return(x)
                                }

                output_list <- lapply(ids, ReformatIds)
                
                IdtaxaDf<-plyr::ldply(output_list, rbind)
                names(IdtaxaDf)[1]<-"ESV"

        #if rank not in ids
            } else {

                ExtractFromIds<-function(list_item, category){
                        ListItemLength<-length(list_item[[category]])
                        if (ListItemLength < Nranks) {
                        list_item[[category]][(ListItemLength+1):Nranks] <-NA
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

            }


        #merge with seq data table
            SeqDataTable<-merge(TableToMergeTo,IdtaxaDf, by= "ESV", all=TRUE)
            
            return(list(SeqDataTable=SeqDataTable, PlotData=ids))
    } else {
        print("No assignment performed")
        return(NULL)
    }
}
