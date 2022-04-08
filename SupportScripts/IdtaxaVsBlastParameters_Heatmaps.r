# a script to explore how blast parameters effect assignment agreement with idtaxa
# 
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(reshape2)
library(viridis)


#data load
    b<-SeqDataTable2Phyloseq(SeqDataTablePath="/home/j/Dropbox/BioinformaticPipeline_Env/Results/Margaux2019_COI_210322/Margaux2019_COI_SeqDataTable.RDS", clustering="curatedOTU", Metadata=NULL, assignment="BLAST")
    i<-SeqDataTable2Phyloseq(SeqDataTablePath="/home/j/Dropbox/BioinformaticPipeline_Env/Results/Margaux2019_COI_210322/Margaux2019_COI_SeqDataTable.RDS", clustering="curatedOTU", Metadata=NULL, assignment="Idtaxa")
    d<-readRDS("/home/j/Dropbox/BioinformaticPipeline_Env/Results/Margaux2019_COI_210322/Margaux2019_COI_SeqDataTable.RDS")
    dtest_new<-readRDS("/home/j/Dropbox/BioinformaticPipeline_Env/Results/COI_test_terrestrial subsampled_0confidence/COI_test_SeqDataTable.RDS")
    dtest_original<-readRDS("/home/j/Dropbox/BioinformaticPipeline_Env/Results/COI_test_originalMidoriBiocodeClassifier/COI_test_SeqDataTable.RDS")
    dtest_new40<-readRDS('/home/j/Dropbox/BioinformaticPipeline_Env/Results/COI_test_SeqDataTable.RDS')

#data explore
    names(d)
    dblast<-(d[d$Blast_query_coverage>99,])
    sum(d$Blast_percentIdentical>85, na.rm=TRUE)/ dim(d) [1]
    tail(sort(table(d$Phylum)),10)
    tail(sort(table(dblast$Phylum)),10)
    tail(sort(table(d$Rank_3)),12)





#create some heatmaps
    PercentSharedIDs_plotlist<-CreateAssignmentAgreementHeatmaps(d)
    grid.arrange(grobs = PercentSharedIDs_plotlist , ncol = 3, top=textGrob("BIOT COI (full), Original Midori-Biocode - Proportion of agreed assignments - out of all blast assignable sequences with a blast top hit at these parameter settings",gp=gpar(fontsize=20,font=3) ) ) ## display plot

    PercentSharedIDs_plotlist_test_original<-CreateAssignmentAgreementHeatmaps(dtest_original)
    grid.arrange(grobs = PercentSharedIDs_plotlist_test_original , ncol = 3, top=textGrob("BIOT COI (test), Original Midori-Biocode - Proportion of agreed assignments - out of all blast assignable sequences with a blast top hit at these parameter settings",gp=gpar(fontsize=20,font=3) ) ) ## display plot

    PercentSharedIDs_plotlist_test_new<-CreateAssignmentAgreementHeatmaps(dtest_new)
    grid.arrange(grobs = PercentSharedIDs_plotlist_test_new , ncol = 3, top=textGrob("BIOT COI (test), New terrestrial subsampled Midori-Biocode - 0 confidence - Proportion of agreed assignments - out of all blast assignable sequences with a blast top hit at these parameter settings",gp=gpar(fontsize=20,font=3) ) ) ## display plot

    PercentSharedIDs_plotlist_test_new<-CreateAssignmentAgreementHeatmaps(dtest_new40)


pdf(file = "/home/j/Desktop/BlastIdtaxaComparisonHeatmaps.pdf", width=20, height =8 )      

    grid.arrange(grobs = PercentSharedIDs_plotlist_test_new , ncol = 3, top=textGrob(
" Heatmaps showing proportion of sequences with agreed Idtaxa and BLAST assignments out of all BLAST assignable sequences with a BLAST top hit at these parameter settings \n
Idtaxa 40% confidence threshold used throughout (moderate confidence) ; BIOT COI (Subsampled data) ; New terrestrial order subsampled Midori-Biocode reference library ",
                                                                                        gp=gpar(fontsize=12,font=1) ) ) ## display plot

    grid.arrange(grobs = PercentSharedIDs_plotlist_test_original , ncol = 3, top=textGrob(
" Heatmaps showing proportion of sequences with agreed Idtaxa and BLAST assignments out of all BLAST assignable sequences with a BLAST top hit at these parameter settings \n
Idtaxa 40% confidence threshold used throughout (moderate confidence) ; BIOT COI (Subsampled data) ; Original Midori-Biocode reference library ",
                                                                                        gp=gpar(fontsize=12,font=1) ) ) ## display plot

dev.off()


PercentSharedIDs_plotlist_test_new[[4]]







# add a heatmap showing proportion of sequences with blast tophit at these parameters
rankcols<-c("Rank_1", "Rank_2", "Rank_3", "Rank_4", "Rank_5", "Rank_6", "Rank_7", "Rank_8")
for (i in 1: length(rankcols)) {
    d[[rankcols[i]]]<-paste0(toupper(substr(d[[rankcols[i]]], 1, 1)), substr(d[[rankcols[i]]], 2, nchar(d[[rankcols[i]]])), sep="")
}
d<-dtest_new40


dblastassignmentexists<-d[!is.na(d$Blast_percentIdentical),]

dblast<-dblastassignmentexists[   dblastassignmentexists$Blast_query_coverage>=95 & dblastassignmentexists$Blast_percentIdentical>=95 , ]

dblast$Rank_1==dblast$Root


tail(sort(table(dblast$Rank_8)),10)

tail(sort(table(dblast$Species)),10)

sum(dblast$Rank_8 == dblast$Species, na.rm=TRUE) /dim(dblast)[1]

dblast$Blast_percentIdentical

    head(cbind(dblast$Rank_7, dblast$Genus),10)


dblastassignmentexists$Blast_percentIdentical [dblastassignmentexists$Blast_percentIdentical>=99]