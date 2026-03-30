#create taxa dataframe from biocode fasta
#E.G. x<-CreateTaxaDfFromBiocodeFasta(path="~/Dropbox/BioinformaticPipeline_Env/Data/Raw", BiocodeFasta="BIOCODE_JAN21_ER_edit.fasta")

    CreateTaxaDfFromBiocodeFasta<-function(path, BiocodeFasta) {
        
        #readfasta as strings
        BiocodeFasta <- readLines(file.path(path,BiocodeFasta))

        # if odd number of rows then the last is a fragement and should be removed
        if (length(BiocodeFasta)%%2==1) {
            BiocodeFasta<-BiocodeFasta[1:length(BiocodeFasta)-1]
        }


        BiocodeFastaNames<-BiocodeFasta[seq(from=1, to=length(BiocodeFasta), by=2)]
        
        Tax_df<-as.data.frame( matrix(NA, ncol =8 , nrow = length(BiocodeFastaNames)) )
        names(Tax_df)<-c(    "NumID", "ID", "Phylum", "Class", "Order", "Family", "Genus", "Species")
        
        for (i in 1:length(BiocodeFastaNames)) {

            string<-sub('.', '', BiocodeFastaNames[i])
            string<-strsplit(string, split= "\\|")
            string<-unlist(string)

            Tax_df$NumID[i]<-string[[3]]
            Tax_df$ID[i]<-string[1]

            x<-unlist(strsplit(string[2], split=","))
            x<-unlist(lapply(FUN=strsplit, x, ":", 2))
            ranks<-rep(NA, 6)
            for (j in 1: length(ranks)) {
                ranks[j]<-x[ (j*2) ]
            }

            Tax_df[i,3:8]<-ranks
        }

        return(Tax_df)
    }