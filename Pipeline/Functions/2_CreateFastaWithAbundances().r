            CreateFastaWithAbundances<-function(SeqDataTable, clustering){
                
                abundances <-SeqDataTable %>%
                                group_by ( {{ clustering }} ) %>%
                                summarise_if(is.numeric,sum) %>%
                                select(starts_with("Sample")) %>% 
                                rowSums()
        
                sequences_abundances<-paste0(SeqDataTable$Sequence, "_",abundances)

                write.fasta(sequences=as.list(SeqDataTable$Sequence), names=sequences_abundances,
                        file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_ESVSequencesWithAbundances.fasta")))    
            }