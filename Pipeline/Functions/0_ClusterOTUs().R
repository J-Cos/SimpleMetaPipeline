ClusterOTUs<-function(linkage, #either complete (vsearch) or single (swarm)
                    differences,
                    threads,
                    TableToMergeTo,
                    HPC,
                    SimilarityThreshold,
                    clustering 
                        ) {

    if (linkage=="single") {

        CreateFastaWithAbundances(SeqDataTable=TableToMergeTo, clustering=clustering)

        ClusterOutput<-RunSwarm(differences=differences, 
                threads=threads, 
                TableToMergeTo=TableToMergeTo, 
                HPC=HPC)

    }else if (linkage=="complete"){

        ClusterOutput<-RunVSEARCH(SimilarityThreshold=SimilarityThreshold, 
                clustering=clustering, 
                TableToMergeTo=TableToMergeTo)

    } else { 
        print("unknown linkage specified. Must be either 'single' (swarm) or 'complete' (vsearch) ")
        return(NULL)
    }

    return(ClusterOutput)
}
