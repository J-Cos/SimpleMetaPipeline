RunVSEARCH<-function(SimilarityThreshold=0.97, clustering="curatedESV", TableToMergeTo=LuluOutput1) { 

    #get representative seqs for clustering level
    if (clustering=="ESV") {
        Sequences<-as.list(TableToMergeTo$Sequence)
    } else if (clustering =="curatedESV") {
        Sequences<-as.list(TableToMergeTo$Sequence[TableToMergeTo$CuratedESVRepresentativeSequence])
    } else  { (print("unrecognised clustering request")) }

    write.fasta(sequences=Sequences, names=Sequences,
                        file.out=file.path(path, "IntermediateOutputs", paste0(dataname, "_", clustering, "_sequences_NamedAsSeqs.fasta")))

    #create matchlist from fasta
    sequencesfile<-paste0(dataname, "_", clustering, "_sequences_NamedAsSeqs.fasta")
    system(command= paste0("mkdir -p ", file.path(path,"IntermediateOutputs", paste0(dataname , "_vsearchOTUs")),
                            " && vsearch --cluster_size ",
                            file.path(path,"IntermediateOutputs",sequencesfile),
                            " --centroids ", file.path(path,"IntermediateOutputs","centroidsequences.fasta"),   
                            " --clusterout_id ",
                            " --clusters ", file.path(path,"IntermediateOutputs", paste0(dataname , "_vsearchOTUs"), "Cluster"),
                                " --uc ", file.path(path,"IntermediateOutputs", paste0(dataname , "_vsearchdata")),
                            "   --id ", SimilarityThreshold,
                            " --relabel OTU   ")) 

    vsearchclusters<-as.data.frame(read.table(file.path(path,"IntermediateOutputs", paste0(dataname , "_vsearchdata"))))
    vsearchclusters<-vsearchclusters[which(vsearchclusters$V1!="C"),c(1,2,9)]

    OTULabeledSeqs<-data.frame(OTU=paste0("OTU_",vsearchclusters$V2), Sequence=vsearchclusters$V9,
                                OTUrepresentativeSequence= (vsearchclusters$V1=="S") )

    SeqDataTable<-MergeClusteringResultsintoDataTable(TableToMergeTo=TableToMergeTo, 
                                                    OTULabeledSeqs=OTULabeledSeqs)

    return(SeqDataTable)
}