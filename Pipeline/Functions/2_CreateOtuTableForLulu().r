CreateOtuTableForLulu<-function(Input, clustering){
            
            nsamples<-DadaOutput$SeqDataTable %>%
                        select(starts_with("Sample")) %>% 
                        dim() %>%
                        '['(2)


            #extract from SeqDataTable
            if (clustering=="ESV") {
                OTUtable<-Input%>%
                    group_by( ESV ) %>%
                    select(starts_with("Sample")) %>% 
                        summarise_each(funs(sum))
            } else if (clustering=="OTU") {
                OTUtable<-Input%>%
                    group_by( OTU ) %>%
                    select(starts_with("Sample")) %>% 
                        summarise_each(funs(sum))
            }
            #convert to df and set OTU names to rownames
            OTUtable<-as.data.frame(OTUtable)
            rownames(OTUtable)<-OTUtable[[clustering]]
            OTUtable<-OTUtable[,-1]

            return(OTUtable)
}