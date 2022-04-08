

CreateAssignmentAgreementHeatmaps<-function(d) {
    #  capitalise unclassifieds in idtaxa output to match blast unclassifieds
        rankcols<-c("Rank_1", "Rank_2", "Rank_3", "Rank_4", "Rank_5", "Rank_6", "Rank_7", "Rank_8")
        for (i in 1: length(rankcols)) {
            d[[rankcols[i]]]<-paste0(toupper(substr(d[[rankcols[i]]], 1, 1)), substr(d[[rankcols[i]]], 2, nchar(d[[rankcols[i]]])), sep="")
        }

    # ensure we are only looking at blast assignable sequences loses ca. 1500 seqs
        d<-d[!is.na(d$Blast_percentIdentical),]


    #create heatmaps for all ranks fro all values of blast parameters (query coverage and percent identical)
        coverages<-seq(from=50,to=100, by=1)
        similarities<-seq(from=50,to=100, by=0.5)
        ranks<-list(    c("Kingdom", "Rank_2"),
                        c("Phylum", "Rank_3"),
                        c("Class", "Rank_4"),
                        c("Order", "Rank_5"),
                        c("Family", "Rank_6"),
                        c("Genus", "Rank_7"),
                        c("Species", "Rank_8")
        )

        PercentSharedIDs_plotlist<-list()                            

        for (rank in 1:length(ranks)) {

            PercentSharedIDs<-matrix(   nrow =length(coverages),
                                    ncol=length(similarities),
                                    dimnames=list(coverages, similarities) )
            PercentSeqsAssigned<-PercentSharedIDs

            for (qcovIndex in 1:length(coverages)) {
                for (simIndex in 1:length(similarities) ) {

                dblast<-d[   d$Blast_query_coverage>=coverages[qcovIndex] & d$Blast_percentIdentical>=similarities[simIndex] , ]
                PercentSharedIDs[qcovIndex, simIndex]<-  sum(dblast[[ ranks[[rank]][1] ]] == dblast[[ ranks[[rank]][2] ]], na.rm=TRUE)/dim(dblast)[1]
                PercentSeqsAssigned[qcovIndex, simIndex]<-dim(dblast)[1]/dim(d)[1]

                PercentSharedIDs_plotlist[[rank]]<-melt(PercentSharedIDs) %>%
                                                    ggplot(aes(Var1, Var2)) +
                                                        geom_tile(aes(fill = value)) + 
                                                        #geom_text(aes(label = round(value, 3))) +
                                                        scale_fill_viridis(limits=c(0, 1) )+
                                                        xlab("Query Coverage") +
                                                        ylab("Percent Identical") +
                                                        ggtitle(ranks[[rank]])
                }
            }
        }

                PercentSeqsAssigned_plot<-melt(PercentSeqsAssigned) %>%
                                                    ggplot(aes(Var1, Var2)) +
                                                        geom_tile(aes(fill = value)) + 
                                                        #geom_text(aes(label = round(value, 3))) +
                                                        scale_fill_viridis( limits=c(0, 1) )+
                                                        xlab("Query Coverage") +
                                                        ylab("Percent Identical") +
                                                        ggtitle("Proportion of blast assignable seqs with a blast top hit at these parameter settings")

    return(PercentSharedIDs_plotlist)
}