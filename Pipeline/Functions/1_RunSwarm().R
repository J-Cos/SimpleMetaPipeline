RunSwarm<-function(differences, threads, TableToMergeTo){         
         
    # enabling sudo on desktop and non-sudo on cluster
    system(command=paste0(
        "swarm", 
        " -d ", 
        differences,
        " -f ", 
        " -t ", 
        threads, 
        " ", 
        file.path(path, "IntermediateOutputs", paste0(dataname, "_ESVSequencesWithAbundances.fasta")),
        " -s ", 
        file.path(path, "IntermediateOutputs", paste0(dataname, "_SwarmStats.csv")),
        " -w ", 
        file.path(path, "IntermediateOutputs", paste0(dataname, "SwarmSequences.fasta")),
        " -o ", 
        file.path(path, "IntermediateOutputs", paste0(dataname, "_SwarmAllAmplicons.txt")), 
        " -i ", 
        file.path(path, "IntermediateOutputs", paste0(dataname, "_SwarmInternalStructures.txt"))))

    # read in data created by swarm
    Allamplicons <- readLines(file.path(path, "IntermediateOutputs", paste0(dataname, "_SwarmAllAmplicons.txt")))
    
    # create seq data table containing, i) sequence id ii) otu assignment, iii) whether a sequence is representative for that OTU
    SeqDataList<-list()
    for (i in 1:(length(Allamplicons))) {
        OTU<-data.frame(Sequence=unlist(strsplit(Allamplicons[i]," ")), OTU=paste0("OTU_", i))
        #OTU$abundance<-as.numeric(sapply(strsplit(OTU$sequence, "_"), `[`, 2))
        OTU$Sequence<-sapply(strsplit(as.character(OTU$Sequence), "_"), `[`, 1)
        OTU$OTUrepresentativeSequence<-FALSE
        OTU$OTUrepresentativeSequence[1]<-TRUE # for each of the esvs in a single otu the first is the representative seq (this is standard swarm output)
        SeqDataList[[i]]<-OTU
    }
    SeqDataTable<-as.data.frame(data.table::rbindlist(SeqDataList))

    SeqDataTable<-MergeClusteringResultsintoDataTable(TableToMergeTo=TableToMergeTo, OTULabeledSeqs=SeqDataTable)

    return(SeqDataTable)
}