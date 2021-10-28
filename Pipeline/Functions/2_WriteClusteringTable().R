            
WriteClusteringTable<-function(){
    ClusteringTable<-data.frame(  
    "Number of ESVs"= dim(LuluOutput)[1], 
    "Number of OTUs" = dim(unique(LuluOutput[1]))[1], 
    "Number of cOTUs" = dim(unique(LuluOutput[6]))[1])

    write.csv( ClusteringTable, file=file.path(path,"Results",paste0(dataname,"_ClusteringTable.csv")) )
}            
