#dependencies 
    #library(phyloseq)
    #library(tidyverse)
    #library(DECIPHER)
               
# function to take output of bioinformatic pipeline into phyloseq
    #example
        #example<-SeqDataTable2Phyloseq(SeqDataTablePath="23s_test_SeqDataTable.RDS", clustering="curatedOTU", Metadata=NULL, assignment="BLAST", ClusterAssignment="RepresentativeSequence")
        # SeqDataTable2Phyloseq(SeqDataTablePath="~/Dropbox/BioinformaticPipeline_Env/Results/16s_multirun_test_SeqDataTable.RDS", clustering="ESV", Metadata=NULL, assignment="Idtaxa", ClusterAssignment="RepresentativeSequence")

#Metadata argument deprecated, now overwritten by run data. Instead add metadata afterward.
#clustering can be ESV, curatedESV, OTU, curatedOTU
#assignment can be BLAST, Idtaxa or None
#ClusterAssignment can be "RepresentativeSequence" or between 0 and 1 to represent the proportion of reads that must share an assignment for OTU to receive that assignment
#ReformatSampleNames removes "001" from end and 'Sample_' from begining of all sample names (this is added by the sequencing/pipeline as standard hence this argument defalts to true)

SeqDataTable2Phyloseq<-function(SeqDataTablePath, clustering, Metadata=NULL, assignment="None", BLASTThreshold=NULL, ClusterAssignment="RepresentativeSequence", ReformatSampleNames=TRUE){


    require(phyloseq)
    require(tidyverse)
    require(DECIPHER)

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
                
                #filter out blast top hits below threshold
                if (!is.null(BLASTThreshold)) {
                    ESVsBelowBlastThreshold<-SeqDataTable$Blast_percentIdentical < BLASTThreshold | SeqDataTable$Blast_query_coverage < BLASTThreshold
                    ESVsBelowBlastThresholdOrNa<-ESVsBelowBlastThreshold
                    ESVsBelowBlastThresholdOrNa[is.na(ESVsBelowBlastThresholdOrNa)] <-TRUE
                    SeqDataTable[ESVsBelowBlastThresholdOrNa,AssignmentIndices]<-NA
                }

            } else if (assignment== "Idtaxa") {
                AssignmentIndices <- ( max(SampleIndices) + 1 ) : ( grep("_confidence", ignore.case=TRUE, colnames(SeqDataTable))[1] -1 )
            }


    #create input matrices
        #otumat - with sample names reformating to match standard in metadata
            #create
                if (clustering=="ESV"){
                    otumat<-SeqDataTable %>% 
                        select("ESV", SampleIndices) %>%
                        column_to_rownames("ESV")%>% 
                        as.matrix()

                } else {        
                    otumat<-SeqDataTable %>% 
                        group_by( !!symclustering ) %>% 
                        summarise_at(colnames(SeqDataTable) [SampleIndices],sum) %>% 
                        column_to_rownames(clustering)%>% 
                        as.matrix()
                }

            # if sample names have run appended 
                #1) pick deepest version of any duplicate samples
                StrippedSampleNames_original<-otumat %>% colnames %>%
                            strsplit(., "__") %>%
                            unlist  %>%
                            `[`(c(TRUE, FALSE))
                
                DuplicatedSamples<-which(table(StrippedSampleNames_original)>1) %>% names
                for (DuplicatedSample in DuplicatedSamples) {
                    DuplicateSampleRunNames<-colnames(otumat)[grep(DuplicatedSample, colnames(otumat),fixed=TRUE)]
                    SizeOfDuplicates<-otumat[,DuplicateSampleRunNames] %>% colSums
                    DuplicateToKeep<-sample(SizeOfDuplicates[SizeOfDuplicates==max(SizeOfDuplicates)], 1) %>% names #keep biggest duplicate sample in case there are two equally large)
                    DuplicatesToRemove<-DuplicateSampleRunNames[DuplicateSampleRunNames!=DuplicateToKeep] #
                    otumat<-otumat[,-which(colnames(otumat)==DuplicatesToRemove)] #remove unwanted duplicate colums from otumat
                }
                
                # 2) strip '__RunX' from end of sample name and record 'RunX' in metadata 
                SampsAndRuns<-otumat %>% colnames %>%
                            strsplit(., "__") %>%
                            unlist 
                StrippedSampleNames<-SampsAndRuns %>%
                            `[`(c(TRUE, FALSE))
                RunIdentifier<-SampsAndRuns %>%
                            `[`(!c(TRUE, FALSE))
                
                colnames(otumat)<-StrippedSampleNames


                #remove "Sample_" and added to front of sample names by pipeline and "001" added by sequencing
                    if(ReformatSampleNames){  
                        reformatSampleNames<-function(list_item) {
                            string<-unlist(strsplit(list_item, "_"))
                            newstring<-string[c(-1, (-length(string)+1):-length(string))]
                            newstring<-paste(newstring, collapse="_")
                            return(newstring)
                        }
                        colnames(otumat)<-lapply(colnames(otumat), reformatSampleNames)
                    }
                
                Metadata<-data.frame(row.names=colnames(otumat),RunIdentifier) %>% sample_data


        #taxmat
        if (clustering=="ESV") {
              taxmat<-SeqDataTable[,AssignmentIndices]
                taxmat$cluster<-SeqDataTable$ESV
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
        } else {

            if (ClusterAssignment=="RepresentativeSequence") { #take assignments from representative seqs
                taxmat<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,][,AssignmentIndices]
                taxmat$cluster<-SeqDataTable[[clustering]][SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE]
                taxmat<-arrange(taxmat, cluster) %>%
                    column_to_rownames("cluster") %>%
                    as.matrix()
            } else { # take assignments based on most common assignments to ESVs within cluster (with threshold set as function argument)
                taxmat<-c()
                for (cluster in rownames(otumat)) { #loop through clusters
                    ClusterDataTable<-SeqDataTable[SeqDataTable[[symclustering]]==cluster,]
                    ClusterAssignments<-cbind(  ClusterDataTable[,AssignmentIndices],
                            Proportions=rowSums(ClusterDataTable[,SampleIndices])/sum(ClusterDataTable[,SampleIndices]) )
                    
                    clusterRow<-cluster
                    Ranks<-colnames(ClusterAssignments)[-length(colnames(ClusterAssignments))]
                    for (col in Ranks ) { #loop through ranks for the cluster
                        if (length(unique(ClusterAssignments[[col]]))==1 ) {
                            clusterRow <- c( clusterRow,unique(ClusterAssignments[[col]]) )
                        } else {
                            AbundancePerAssignment<-ClusterAssignments %>%
                                                        group_by(!!sym(col)) %>%
                                                        summarise(Proportions=sum(Proportions))
                            if (max(AbundancePerAssignment$Proportions) >ClusterAssignment) {
                                #colAssignment<-as.character(AbundancePerAssignment[AbundancePerAssignment$Proportions==max(AbundancePerAssignment$Proportions),1])
                                colAssignment<-as.character(AbundancePerAssignment[order(-AbundancePerAssignment$Proportions),][1,1]) #avoids problem of two equally abundant options
                                clusterRow <- c( clusterRow,colAssignment )
                            } else {
                                clusterRow<- c( clusterRow,NA )
                            }
                        }
                    }
                names(clusterRow)<-c("clusterName",Ranks)
                taxmat<-as.data.frame(rbind(taxmat, clusterRow))
                }
            rownames(taxmat)<-NULL
            taxmat<-column_to_rownames(taxmat,"clusterName" )
            taxmat<-as.matrix(taxmat)
            }
        }

        #refseqs
        if (clustering=="ESV") {
            refseqs_df<-SeqDataTable[,c("Sequence","ESV")]
            refseqs_df<- arrange(refseqs_df, ESV ) %>%
            column_to_rownames("ESV")
            refseqs <- DNAStringSet(refseqs_df$Sequence)
            names(refseqs)<-rownames(refseqs_df)
        } else {
            refseqs_df<-SeqDataTable[SeqDataTable[[RepresentativeSequenceColumnName]]==TRUE,c("Sequence",clustering)]
            refseqs_df<- arrange(refseqs_df, clustering ) %>%
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
