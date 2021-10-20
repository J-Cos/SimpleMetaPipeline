            CreateFastaWithAbundances<-function(ESVtable,ESVsequences){
                sequences_abundances<-c()
                abundances<-colSums(ESVtable)
                for (r in 1:length(colnames(ESVtable))) {
                    sequences_abundances[r]<-paste0(colnames(ESVtable)[r], "_",abundances[r])
                }

                write.fasta(sequences=as.list(ESVsequences), names=sequences_abundances,
                        file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_ESVsequences.fasta")))    
            }