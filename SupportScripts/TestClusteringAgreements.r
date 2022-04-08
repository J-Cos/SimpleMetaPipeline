path<- "/home/j/Dropbox/BioinformaticPipeline_Env"

#2  #R packages
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( 
                    phyloseq)

    #Functions
        function_files<-list.files(file.path(path, "BioinformaticPipeline","SupportFunctions", ""))
        sapply(file.path(path, "BioinformaticPipeline","SupportFunctions",function_files),source)


df1<-read.csv(file.path(path, "Results", "PNGFullTest_23s_SeqDataTable.csv"))
df2<-readRDS(file.path(path, "Results", "PNGFullTest_23s", "PNGFullTest_23s_SeqDataTable.RDS"))
GetAgreementStats(df1)
GetAgreementStats(df)

GetAgreementStats<-function(df){
    cESVs<- data.frame(cESVs=unique(df$curatedESV),
                        SameAssignment=rep(NA, length(unique(df$curatedESV)))
                        )
    for (cESV in cESVs$cESVs) {
        cESVs$SameAssignment[cESVs$cESVs==cESV]<-length(unique(df[df$curatedESV==cESV,]$Rank_5))==1
        cESVs$NumSeqs[cESVs$cESVs==cESV]<-length(df[df$curatedESV==cESV,]$Rank_5)
        cESVs$NumAssignments[cESVs$cESVs==cESV]<-length(unique(df[df$curatedESV==cESV,]$Rank_5))
    }

    cESVs_morethan1<-cESVs[cESVs$NumSeqs>1,]

    plot(density(cESVs_morethan1$NumAssignments))

    SameAssignment_Percent_cESVs<-sum(cESVs_morethan1$SameAssignment)/dim(cESVs_morethan1)[1]
    SameAssignment_Percent_ESVs<-sum(as.numeric(cESVs_morethan1$SameAssignment)*cESVs_morethan1$NumSeqs)/sum(cESVs_morethan1$NumSeqs)

    return(list(SameAssignment_Percent_cESVs,SameAssignment_Percent_ESVs ))
}

 
#agreement between kmers and cesvs at different lulu matchrates

FilePaths<-file.path(path, "Results", paste0(90:99,"_23s_test_SeqDataTable.RDS"))

data<-lapply(FilePaths, readRDS)
Stats<-lapply(data, GetAgreementStats)
Stats<-matrix(unlist(Stats), nrow=2)
rownames(Stats)<-c("Percent_cESVs", "Percent_ESVS")
colnames(Stats)<-c(90:99)

data[1]




names(df1)










#agreement between cesv, otus and cotus


df<-read.csv('/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/latestCOIrun.csv')
length(unique(df$curatedESV))
sum(df$CuratedESVRepresentativeSequence)


    cESVs<-data.frame(cESVs=unique(df$curatedESV),
                    OTUagreement=NA,
                    cOTUagreement=NA)

    for (i in 1: dim(cESVs)[1]) {
        df_cESV<-df[df$curatedESV==cESVs[i,1],]
        cESVs$OTUagreement[i]<-max(table(df_cESV$OTU))  /length(df_cESV$OTU)
        cESVs$cOTUagreement[i]<-max(table(df_cESV$curatedOTU))  /length(df_cESV$curatedOTU)
    }

    summary(cESVs)
 ####
    cOTUs<-data.frame(cOTUs=unique(df$curatedOTU),
                    OTUagreement=NA,
                    cESVagreement=NA)

    for (i in 1: dim(cOTUs)[1]) {
        df_cOTU<-df[df$curatedOTU==cOTUs[i,1],]
        cOTUs$OTUagreement[i]<-max(table(df_cOTU$OTU))  /length(df_cOTU$OTU)
        cOTUs$cESVagreement[i]<-max(table(df_cOTU$curatedESV))  /length(df_cOTU$curatedESV)
    }

    summary(cOTUs)

 ####
    OTUs<-data.frame(OTUs=unique(df$OTU),
                    cOTUagreement=NA,
                    cESVagreement=NA)

    for (i in 1: dim(OTUs)[1]) {
        df_OTU<-df[df$OTU==OTUs[i,1],]
        OTUs$cOTUagreement[i]<-max(table(df_OTU$curatedOTU))  /length(df_OTU$curatedOTU)
        OTUs$cESVagreement[i]<-max(table(df_OTU$curatedESV))  /length(df_OTU$curatedESV)
    }

    summary(OTUs)
