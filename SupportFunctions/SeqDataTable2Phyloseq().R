# function to take output of bioinformatic pipeline into phyloseq
    #example
        #SeqDataTable2Phyloseq(SeqDataTablefile="test_SeqDataTable.RDS", clustering="ESV", Metadata=NULL)

SeqDataTable2Phyloseq<-function(SeqDataTablefile, clustering, Metadata=NULL){
    numcols<-17 # number of columns after the end of the samples

    SeqDataTable<-readRDS(file.path(path, "Results", SeqDataTablefile))

    if (clustering == "ESV") {

        #otumat
        otumat<-SeqDataTable %>% 
            group_by( ESV ) %>% 
            summarise_if(is.numeric,sum) %>% 
            column_to_rownames("ESV")%>% 
            as.matrix() 
        #taxmat
        taxmat<-SeqDataTable[,(dim(SeqDataTable)[2]-13):(dim(SeqDataTable)[2]-7)]
        taxmat$cluster<-SeqDataTable$ESV
        taxmat<-arrange(taxmat, cluster) %>%
            column_to_rownames("cluster") %>%
            as.matrix()
        #refseqs
        refseqs_df<-SeqDataTable[,c(3,1)] %>%
            arrange(, ESV) %>%
            column_to_rownames("ESV")
        refseqs <- DNAStringSet(refseqs_df$sequence)
        names(refseqs)<-rownames(refseqs_df)


    } else if (clustering == "OTU") {

        #otumat
        otumat<-SeqDataTable %>% 
            group_by( OTU ) %>% 
            summarise_if(is.numeric,sum) %>% 
            column_to_rownames("OTU")%>% 
            as.matrix() 
        #taxmat
        taxmat<-SeqDataTable[SeqDataTable$OTUrepresentativeSequence==TRUE,][,(dim(SeqDataTable)[2]-13):(dim(SeqDataTable)[2]-7)]
        taxmat$cluster<-SeqDataTable$OTU[SeqDataTable$OTUrepresentativeSequence==TRUE]
        taxmat<-arrange(taxmat, cluster) %>%
            column_to_rownames("cluster") %>%
            as.matrix()
        #refseqs
        refseqs_df<-SeqDataTable[SeqDataTable$OTUrepresentativeSequence==TRUE,c(3,2)] %>%
            arrange(, OTU) %>%
            column_to_rownames("OTU")
        refseqs <- DNAStringSet(refseqs_df$sequence)
        names(refseqs)<-rownames(refseqs_df)

    } else if (clustering =="cOTU") {

        #otumat
        otumat<-SeqDataTable %>% 
            group_by( curatedOTU ) %>% 
            summarise_if(is.numeric,sum) %>% 
            column_to_rownames("curatedOTU")%>% 
            as.matrix() 
        #taxmat
        taxmat<-SeqDataTable[SeqDataTable$CuratedOTURepresentativeSequence==TRUE,][,(dim(SeqDataTable)[2]-13):(dim(SeqDataTable)[2]-7)]
        taxmat$cluster<-SeqDataTable$curatedOTU[SeqDataTable$CuratedOTURepresentativeSequence==TRUE]
        taxmat<-arrange(taxmat, cluster) %>%
            column_to_rownames("cluster") %>%
            as.matrix()
        #refseqs
        refseqs_df<-SeqDataTable[SeqDataTable$CuratedOTURepresentativeSequence==TRUE,c(3,7)] %>%
            arrange(, curatedOTU) %>%
            column_to_rownames("curatedOTU")
        refseqs <- DNAStringSet(refseqs_df$sequence)
        names(refseqs)<-rownames(refseqs_df)

    } else { 
        print("unrecognised clustering") 
    }


    if (is.null(Metadata)){
        ps<-phyloseq(
            # sample_data(),
            otu_table(otumat, taxa_are_rows=TRUE),
            tax_table(taxmat),
            #phy_tree(),
            refseq(refseqs)
        )
    } else {
        ps<-phyloseq(
            sample_data(Metadata),
            otu_table(otumat, taxa_are_rows=TRUE),
            tax_table(taxmat),
            #phy_tree(),
            refseq(refseqs)
        ) 
    }


    return(ps)
}
