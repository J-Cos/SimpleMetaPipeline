CreateOtuTableForLulu<-function(Input=SwarmOutput){
            
            nsamples<-dim(DadaOutput$ESVtable)[1]


            #extract from SeqDataTable
            OTUtable<-Input[,2:(2+nsamples)]%>%
                group_by(OTU) %>%
                    summarise_each(funs(sum))
            
            #convert to df and set OTU names to rownames
            OTUtable<-as.data.frame(OTUtable)
            rownames(OTUtable)<-OTUtable$OTU
            OTUtable<-OTUtable[,-1]

            return(OTUtable)
}