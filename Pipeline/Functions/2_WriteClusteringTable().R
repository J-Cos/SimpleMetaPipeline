WriteClusteringTable<-function(FinalOutput, pipeline){
    if (pipeline=="DSL"){
        ClusteringTable<-data.frame(  
            "Number of ESVs"= dim(FinalOutput)[1], 
            "Number of OTUs" = length(unique(FinalOutput$OTU)), 
            "Number of cOTUs" = length(unique(FinalOutput$curatedOTU)))
    } else if (pipeline=="DLSL"){
        ClusteringTable<-data.frame(  
            "Number of ESVs"= dim(FinalOutput)[1], 
            "Number of cESVs"= length(unique(FinalOutput$curatedESV)),  
            "Number of OTUs" = length(unique(FinalOutput$OTU)), 
            "Number of cOTUs" = length(unique(FinalOutput$curatedOTU)))
    }
    write.csv( ClusteringTable, file=file.path(path,"Results",paste0(dataname,"_ClusteringTable.csv")) )
}            
