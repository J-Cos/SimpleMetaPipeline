#dependencies 
    #phyloseq

# function to take output of bioinformatic pipeline into phyloseq
    #example
        #example<-SeqDataTable2Phyloseq(SeqDataTablefile="23s_test_SeqDataTable.RDS", clustering="cOTU", Metadata=NULL)

SeqDataTable2Phyloseq<-function(SeqDataTablePath, clustering, Metadata=NULL){

    SeqDataTable<-readRDS(SeqDataTablePath)

    # very explicit column math to ensure 
        #1) that this function works with variable numbers of samples
        #2) that it is clear how 
        numberstandardcolumns<- 6
        numSamples<-which(names(SeqDataTable)=="CuratedOTURepresentativeSequence") -numberstandardcolumns

        firstTaxCol<-numberstandardcolumns + numSamples + 1
        lastConfCol<-length(names(SeqDataTable))
        numRanks<-(lastConfCol-firstTaxCol+1)/2
        lastTaxCol<-firstTaxCol+numRanks-1
        firstConfCol<-lastConfCol-numRanks

    #main function if statements creating ps components
        if (clustering == "ESV") {

            #otumat
                otumat<-SeqDataTable %>% 
                    group_by( ESV ) %>% 
                    summarise_if(is.numeric,sum) %>% 
                    column_to_rownames("ESV")%>% 
                    as.matrix() 

                reformatSampleNames<-function(list_item) {
                    string<-unlist(strsplit(list_item, "_"))
                    newstring<-string[c(-1, (-length(string)+1):-length(string))]
                    newstring<-paste(newstring, collapse="_")
                    return(newstring)
                }
                colnames(otumat)<-lapply(colnames(otumat), reformatSampleNames)


            if(firstTaxCol!=lastConfCol) {
                #taxmat
                taxmat<-SeqDataTable[,firstTaxCol:lastTaxCol]
                taxmat$cluster<-SeqDataTable$ESV
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
            }

            #refseqs
            refseqs_df<-SeqDataTable[,c("Sequence","ESV")] %>%
                arrange(, ESV) %>%
                column_to_rownames("ESV")
            refseqs <- DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)


        } else if (clustering == "OTU") {

            #otumat
            otumat<-SeqDataTable %>% 
                group_by( OTU ) %>% 
                summarise_if(is.numeric,sum) %>% 
                column_to_rownames("OTU")%>% 
                as.matrix() 

            if(firstTaxCol!=lastConfCol) {
                #taxmat
                taxmat<-SeqDataTable[SeqDataTable$OTUrepresentativeSequence==TRUE,][,firstTaxCol:lastTaxCol]
                taxmat$cluster<-SeqDataTable$OTU[SeqDataTable$OTUrepresentativeSequence==TRUE]
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
            }

            #refseqs
            refseqs_df<-SeqDataTable[SeqDataTable$OTUrepresentativeSequence==TRUE,c("sequence","OTU")] %>%
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
            
            if(firstTaxCol!=lastConfCol) {
                #taxmat
                taxmat<-SeqDataTable[SeqDataTable$CuratedOTURepresentativeSequence==TRUE,][,firstTaxCol:lastTaxCol]
                taxmat$cluster<-SeqDataTable$curatedOTU[SeqDataTable$CuratedOTURepresentativeSequence==TRUE]
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
            }

            #refseqs
            refseqs_df<-SeqDataTable[SeqDataTable$CuratedOTURepresentativeSequence==TRUE,c("sequence","curatedOTU")] %>%
                arrange(, curatedOTU) %>%
                column_to_rownames("curatedOTU")
            refseqs <- DNAStringSet(refseqs_df$sequence)
            names(refseqs)<-rownames(refseqs_df)

        } else { 
            print("unrecognised clustering") 
        }

    #build ps if statements
        if (!exists("taxmat")) {
            if (is.null(Metadata)){
                ps<-phyloseq(
                    # sample_data(),
                    otu_table(otumat, taxa_are_rows=TRUE),
                    #tax_table(taxmat),
                    #phy_tree(),
                    refseq(refseqs)
                )
            } else {
                ps<-phyloseq(
                    sample_data(Metadata),
                    otu_table(otumat, taxa_are_rows=TRUE),
                    #tax_table(taxmat),
                    #phy_tree(),
                    refseq(refseqs)
                ) 
            }
        } else {
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

        }

    return(ps)
}
