
CreateMatchlistForLulu<-function(Input, MatchRate, clustering, HPC) {
            #get clustering labelled fasta for lulu
            
            if (clustering=="ESV") {
                seqs<-Input
            } else if (clustering == "OTU"){
                seqs<-Input[which(Input$OTUrepresentativeSequence==TRUE),]
            }


            write.fasta(sequences=as.list(seqs$Sequence), names=seqs[[clustering]],
                                file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_", clustering, "_sequences.fasta")))

            #create matchlist from fasta
            sequencesfile<-paste0(dataname, "_", clustering, "_sequences.fasta")

            if (HPC==FALSE) { #slighlty different vsearch behaviour on cluster - not sure why - maxhits 10 vs maxhits 0
                system(command= paste0("cd ", file.path(path, "IntermediateOutputs"), "&& vsearch --usearch_global ", 
                                        sequencesfile ," --db ", sequencesfile, 
                                        " --self --id .", MatchRate," --iddef 1 --userout ", dataname ,"_", clustering ,"_match_list.txt -userfields query+target+id --maxaccepts 0 --maxrejects 0 --query_cov .9 --maxhits 10"))
            } else if (HPC==TRUE) { #slighlty different vsearch behaviour on cluster - not sure why - maxhits 10 vs maxhits 0
                system(command= paste0("cd ", file.path(path, "IntermediateOutputs"), "&& vsearch --usearch_global ", 
                        sequencesfile ," --db ", sequencesfile, 
                        " --self --id .", MatchRate," --iddef 1 --userout ", dataname ,"_", clustering ,"_match_list.txt -userfields query+target+id --maxaccepts 0 --maxrejects 0 --query_cov .9 --maxhits 0"))
            }

            matchlist <- read.table(file.path(path, "IntermediateOutputs", paste0(dataname ,"_", clustering ,"_match_list.txt")), sep="\t", header=FALSE)

            return(matchlist)
}
