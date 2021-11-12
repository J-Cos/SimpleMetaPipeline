            CreateFastaWithAbundances<-function(SeqDataTable, clustering){
                
                #extract representative sequences
                if (clustering=="ESV") {
                    seqs<-SeqDataTable

                    abundances <-SeqDataTable %>%
                        group_by ( ESV ) %>%
                        summarise_if(is.numeric,sum) %>%
                        select(starts_with("Sample")) %>% 
                        rowSums()

                } else if (clustering == "curatedESV"){
                    seqs<-SeqDataTable[which(SeqDataTable$CuratedESVRepresentativeSequence==TRUE),]

                    abundances <-SeqDataTable %>%
                        group_by ( curatedESV ) %>%
                        summarise_if(is.numeric,sum) %>%
                        select(starts_with("Sample")) %>% 
                        rowSums()

                }

                sequences_abundances<-paste0(seqs$Sequence, "_",abundances)

                write.fasta(sequences=as.list(seqs$Sequence), names=sequences_abundances,
                        file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_ESVSequencesWithAbundances.fasta")))    
            }