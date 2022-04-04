#dependencies 
    #library(phyloseq)
    #library(tidyverse)
    #library(DECIPHER)
               
# function to take output of bioinformatic pipeline into phyloseq
    #example
        #example<-SeqDataTable2Phyloseq(SeqDataTablefile="23s_test_SeqDataTable.RDS", clustering="curatedOTU", Metadata=NULL, assignment="BLAST")

#clustering can be ESV, curatedESV, OTU, curatedOTU
#assignment can be BLAST, Idtaxa or None

SeqDataTable2Phyloseq<-function(SeqDataTablePath, clustering, Metadata=NULL, assignment="None"){

    #check clustering valid, break if not
        validclusterings<-c("ESV", "curatedESV", "OTU", "curatedOTU")
        if (!clustering  %in%  validclusterings) { print(paste0(clustering, " clustering choice invalid please choose one of: ", validclusterings)) ; return(NULL)}
        validassignments<-c("BLAST", "Idtaxa", "None")
        if (!assignment  %in%  validassignments) { print(paste0( assignment, " assignment choice invalid please choose one of: ", validassignments)) ; return(NULL) }

    #preparatory steps

        SeqDataTable<-readRDS(SeqDataTablePath)

        #rename this column to match all others - should be fixed upstream in the pipeline later
        names(SeqDataTable)[names(SeqDataTable)=="OTUrepresentativeSequence"]<-"OTURepresentativeSequence"


        #get seq datatable and reorder columns
            clusteringcolumns<- c(  "ESV",
                                    "Sequence",
                                    "curatedESV",                      
                                    "CuratedESVRepresentativeSequence",
                                    "OTU",
                                    "OTURepresentativeSequence",       
                                    "curatedOTU",                       
                                    "CuratedOTURepresentativeSequence")
            
            #ensure standard order to clustering cols as above
            SeqDataTable<-SeqDataTable %>% relocate(clusteringcolumns)

        #recode clustering variables to work with tidyverse as arguments
            symclustering<-sym(clustering)
            RepresentativeSequenceColumnName<-paste0(toupper(substr(clustering, 1, 1)), substr(clustering, 2, nchar(clustering)), sep="", "RepresentativeSequence")

        ## get indices
            MainIndices<-which(names(SeqDataTable) %in% clusteringcolumns)
            SampleIndices<-grep("Sample_", colnames(SeqDataTable))
            if (assignment =="BLAST") {
                AssignmentIndices<-( grep("Blast_query_coverage", colnames(SeqDataTable)) +1 ) : dim(SeqDataTable)[2]
            } else if (assignment== "Idtaxa") {
                AssignmentIndices <- ( last(SampleIndices) + 1 ) : ( grep("_Confidence", colnames(SeqDataTable))[1] -1 )
            }


    #create input matrices
        #otumat - with sample names reformating to match standard in metadata
        otumat<-SeqDataTable %>% 
            group_by( !!symclustering ) %>% 
            summarise_if(is.numeric,sum) %>% 
            column_to_rownames(clustering)%>% 
            as.matrix()
            
        reformatSampleNames<-function(list_item) {
            string<-unlist(strsplit(list_item, "_"))
            newstring<-string[c(-1, (-length(string)+1):-length(string))]
            newstring<-paste(newstring, collapse="_")
            return(newstring)
        }
        colnames(otumat)<-lapply(colnames(otumat), reformatSampleNames)

        #taxmat
        if (clustering=="ESV") {
              taxmat<-SeqDataTable[,AssignmentIndices]
                taxmat$cluster<-SeqDataTable$ESV
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
        } else {
            taxmat<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,][,AssignmentIndices]
            taxmat$cluster<-SeqDataTable[[clustering]][SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE]
            taxmat<-arrange(taxmat, cluster) %>%
                column_to_rownames("cluster") %>%
                as.matrix()
        }


        #refseqs
        if (clustering=="ESV") {
              refseqs_df<-SeqDataTable[,c("Sequence",clustering)] %>%
            arrange(, !!symclustering ) %>%
            column_to_rownames(clustering)
            refseqs <- DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        } else {
            refseqs_df<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,c("Sequence",clustering)] %>%
            arrange(, !!symclustering ) %>%
            column_to_rownames(clustering)
            refseqs <- DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        }


    #create phyloseq object depending on which matrices present
        if (assignment=="None" & is.null(Metadata)) {
            ps<-phyloseq(
                        #sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        #tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        } else if (assignment!="None" & is.null(Metadata)) {
            ps<-phyloseq(
                        #sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        } else if (assignment=="None" & !is.null(Metadata)) {
            ps<-phyloseq(
                        sample_data(Metadata),
                        otu_table(otumat, taxa_are_rows=TRUE),
                        #tax_table(taxmat),
                        #phy_tree(),
                        refseq(refseqs) 
                    )
        } else if (assignment!="None" & !is.null(Metadata)) {
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
