
CreateMatchlistForLulu<-function(Input=SwarmOutput, MatchlistRate) {
            #get otu labelled fasta for lulu
            OTUseqs<-Input[which(Input$OTUrepresentativeSequence==TRUE),]
            write.fasta(sequences=as.list(OTUseqs$sequence), names=OTUseqs$OTU,
                                file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_OTUsequences.fasta")))

            #create matchlist from fasta
            OTUsequencesfile<-paste0(dataname, "_OTUsequences.fasta")
            system(command= paste0("cd ", file.path(path, "IntermediateOutputs"), "&& vsearch --usearch_global ", 
                                    OTUsequencesfile ," --db ", OTUsequencesfile, 
                                    " --self --id .", MatchlistRate," --iddef 1 --userout ", dataname ,"_OTUmatch_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"))

            matchlist <- read.table(file.path(path, "IntermediateOutputs", paste0(dataname ,"_OTUmatch_list.txt")), sep="\t", header=FALSE)

            return(matchlist)
}
