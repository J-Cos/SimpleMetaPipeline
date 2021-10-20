SaveSequenceDataTableToDataBase<-function(Input){

    SeqDataTable<-Input


    # Create a connection to an on-disk SQLite 
        dbConn <- dbConnect(SQLite(), file.path(path, "Results", paste0(dataname, ".sqlite"))) # path to new database file

    #prepare sequences for saving as database   

        #preparing sequencees in Biostrigns format
            dna <- Biostrings::DNAStringSet(SeqDataTable$sequence)
            
            #sequence alignment, two options either direct if number of sequences small or via database if many sequences
            DNA<-AlignSeqs(dna)
                
            #set names to the rownames from seq table so they share unique identifiers 
            #(this works as the order is reserved by the previous steps)
            names(DNA)<-1:length(DNA)
                

            # Import sequences created in clustering script
            DECIPHER::Seqs2DB(DNA, type="XStringSet", dbFile=dbConn, identifier=paste0(dataname, "_ESVs"))

    #prepare table for saving as database
        # saving as database reuquires non numeric first column characters, so check this and modify if required
            nsamples<-dim(SeqDataTable)[2]-5
            
            if ( sum(substring(names(SeqDataTable),1,1) %in% as.character(0:9)) >0 ) { #check if any of the first characters of the columns are numbers
                names(SeqDataTable)[2+1:(2+nsamples)]<-paste0("Sample_",names(SeqDataTable)[2+1:(2+nsamples)])

            } else {print("Sample Names Non-numeric")}

        #add TO DB seqdatatable created in clustering script
            Add2DB(SeqDataTable, dbConn)
}