#test taxonomies succesfully merge - i.e there are not duplicate entries - at higher ranks in particular

#1 set seed for replication
    set.seed(0.1)

#2 dependencies    
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( tidyverse,
                    seqinr,
                    DECIPHER,
                    Biostrings
                    )

    combinefasta<-"ARMS_classifier_combinedSeqs.fasta"
    path<-"../../BioinformaticPipeline_Env"
    seqs <- readDNAStringSet(file.path(path, "Data", "Classifiers", combinefasta) )




    splits<-strsplit(names(seqs), ";")
    Tax_df<-plyr::ldply(splits, rbind )




#change the number to see different ranks
tail(sort(table(Tax_df[6])), 20)