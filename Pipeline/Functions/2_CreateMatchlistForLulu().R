CreateMatchlistForLulu<-function(Input, MatchRate, clustering, multithread=multithread) {
    # get clustering labelled fasta for lulu
    if (multithread==TRUE){
        mc.cores<-250
    } else {
        mc.cores<-multithread
    }

    if (clustering=="ESV") {
        seqs<-Input
    } else if (clustering == "OTU"){
        seqs<-Input[which(Input$OTUrepresentativeSequence==TRUE),]
    }

    write.fasta(
        sequences=as.list(seqs$Sequence), 
        names=seqs[[clustering]],
        file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_", clustering, "_sequences.fasta")))

    # create matchlist from fasta
    sequencesfile<-paste0(dataname, "_", clustering, "_sequences.fasta")

    system(command= paste0(
        "cd ", 
        file.path(path, "IntermediateOutputs"), 
        "&& vsearch --usearch_global ", 
        sequencesfile ,
        " --db ", 
        sequencesfile, 
        " --threads ", 
        mc.cores,
        " --self --id .", 
        MatchRate,
        " --iddef 1 --userout ", 
        dataname ,
        "_", 
        clustering ,
        "_match_list.txt -userfields query+target+id --maxaccepts 0 --maxrejects 0 --query_cov .9 --maxhits 10"))

    matchlist <- read.table(file.path(path, "IntermediateOutputs", paste0(dataname ,"_", clustering ,"_match_list.txt")), sep="\t", header=FALSE)

    return(matchlist)
}
