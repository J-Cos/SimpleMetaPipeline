CreateOtuTableForLulu<-function(Input=SwarmOutput){
            
            nsamples<-DadaOutput$SeqDataTable %>%
                        select(starts_with("Sample")) %>% 
                        dim() %>%
                        '['(2)


            #extract from SeqDataTable
            OTUtable<-Input%>%
                group_by(OTU) %>%
                select(starts_with("Sample")) %>% 
                    summarise_each(funs(sum))
            
            #convert to df and set OTU names to rownames
            OTUtable<-as.data.frame(OTUtable)
            rownames(OTUtable)<-OTUtable$OTU
            OTUtable<-OTUtable[,-1]

            return(OTUtable)
}