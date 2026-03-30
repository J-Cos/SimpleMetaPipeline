###############################################################
# simple script to assess phyla level assignment to COI datasets
###############################################################

#0) load tidyverse / phyloseq / pipeline
library(tidyverse)
library(phyloseq)

#1) make a results dataframe that also contains the input parameters used to generate results
out_df<-data.frame(
        Paths=c(rep("~/Dropbox/BioinformaticPipeline_Env/Results/Margaux2019-21_COI/Margaux2019-21_COI_SeqDataTable.RDS",11),
                rep("~/Dropbox/BioinformaticPipeline_Env/Results/Rosie2021_COI/Rosie2021_COI_SeqDataTable.RDS",11),
                rep("~/Dropbox/CrossPacific_Paper/Data/CrossPacific_2015_SeqDataTable.RDS",11) ),
        Data=c(rep("Margaux", 11), rep("Rosie", 11), rep("EmmaSubset", 11)),
        Clustering=rep(c("ESV", "OTU", "OTU", "OTU", "OTU","OTU", "curatedOTU", "curatedOTU", "curatedOTU", "curatedOTU", "curatedOTU"),3),
        ClusterAssignment=rep(c("RepresentativeSequence", "RepresentativeSequence", 0.1, 0.5, 0.75, 0.9, "RepresentativeSequence", 0.1, 0.5, 0.75, 0.9), 3),
        ProportionReadsUnassignedAtPhyla=NA,
        ProportionESVsOrClustersUnassignedAtPhyla=NA)

#2) make our function
    GetProportionAssignedToPhyla<-function(ps) {

        #unclassified options
        UnclassifiedOptions<-c("unclassified_Eukaryota", "Unclassified")

        #convert NAs to "Unclassified"
        tax_df<-as.data.frame(tax_table(ps) ) 
        tax_df[is.na(tax_df)]<-UnclassifiedOptions[2]
        tax_table(ps)<-as.matrix(tax_df)

        #get prop ESVs unclassified at phyla
        propESVs<-sum(tax_table(ps)[,3]==UnclassifiedOptions[1] | tax_table(ps)[,3]==UnclassifiedOptions[2])/length(taxa_names(ps))

        #get prop reads unclassified at phyla
        numReads<-tapply(taxa_sums(ps), factor(tax_table(ps)[, "Phylum"]), sum) #not using tax_glom (potentially related to phyloseq issue 223) #Phylum or Rank_3
        df<-data.frame(numReads)
        df$propReads<-df$numReads/sum(df$numReads)

        return(c(   sum(df$propReads[rownames(df)==UnclassifiedOptions[1]], df$propReads[rownames(df)==UnclassifiedOptions[2]]) ,
                    propESVs))
    }

#3 run our function on each row in the results dataframe through a loop
    for (row in 1:dim(out_df)[1]) {

        #function
        out_df[row,5:6] <-  SeqDataTable2Phyloseq(  SeqDataTablePath=out_df$Paths[row], 
                                                    clustering=out_df$Clustering[row], 
                                                    Metadata=NULL, 
                                                    assignment="BLAST", #IDTAXA or BLAST
                                                    BLASTThreshold=85,
                                                    ClusterAssignment=out_df$ClusterAssignment[row]) %>%
                            GetProportionAssignedToPhyla()
        
        #progress updater
        print(paste0("Row ", row, " complete"))
    }


out_df$ClusteringDetails<-paste0(out_df$Clustering, "_", out_df$ClusterAssignment)
out_df %>%
    ggplot( aes(x=ProportionReadsUnassignedAtPhyla, y=ProportionESVsOrClustersUnassignedAtPhyla)) +
        geom_point(aes(color=ClusterAssignment, shape=Clustering), size=5) +
        stat_ellipse(aes(group=Data)) +
        geom_text(x=0.07, y=0.05, label="Emma Subset") +
        geom_text(x=0.12, y=0.14, label="Margaux") +
        geom_text(x=0.2, y=0.22, label="Rosie") +
        #ylim(c(0, 0.3))+
        #xlim(c(0, 0.35))+
        ggtitle("Assignment at phyla for COI datasets (Midori-Biocode database, BLAST algorithm (no thresholds))")


out_df2<-(out_df)#idtaxa
out_df3<-(out_df)#blast unthresholded

combined_out<-(rbind(rbind(out_df, out_df2), out_df3))
combined_out$Approach<-c(rep("BLASTthresholded",33), rep("IDTAXA", 33), rep("BLASTunthresholded", 33))

combined_out %>%
    ggplot( aes(x=ProportionReadsUnassignedAtPhyla, y=ProportionESVsOrClustersUnassignedAtPhyla)) +
        geom_point(aes(color=ClusterAssignment, shape=Clustering), size=2) +
        stat_ellipse(aes(group=Data)) +
        #geom_text(x=0.07, y=0.05, label="Emma Subset") +
        #geom_text(x=0.12, y=0.14, label="Margaux") +
        #geom_text(x=0.2, y=0.22, label="Rosie") +
        #ylim(c(0, 0.3))+
        #xlim(c(0, 0.35))+
        ggtitle("Assignment at phyla for COI datasets (Midori-Biocode database, BLAST algorithm (no thresholds))")+
        facet_wrap(~Approach)








write.csv(combined_out, file="combined_out.csv")

ps<-SeqDataTable2Phyloseq(  SeqDataTablePath=out_df$Paths[row], 
                                                    clustering=out_df$Clustering[row], 
                                                    Metadata=NULL, 
                                                    assignment="BLAST",
                                                    ClusterAssignment=out_df$ClusterAssignment[row]) 